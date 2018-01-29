# Packages needed for parallel processing need to be loaded "twice":s
using Bio.Seq: BioSequence, DNASequence, DNAAlphabet, DNAKmer, canonical, FASTASeqRecord, FASTAReader, each, neighbors
@everywhere using Bio.Seq: BioSequence, DNASequence, DNAAlphabet, DNAKmer, canonical, FASTASeqRecord, FASTAReader, each, neighbors
using DataStructures: DefaultDict
@everywhere using DataStructures: DefaultDict
using OpenGene: fastq_open, fastq_read
@everywhere using OpenGene: fastq_open, fastq_read

import GZip
import JLD: load, save
import FileIO: File, @format_str
import Blosc
using FastaIO

## Helper structs:
type AlleleCoverage
  allele::Int16
  votes::Int64
  depth::Int16
  covered_kmers::Int16
  uncovered_kmers::Int16
  gap_list::Vector{Tuple{Int16,Int16}}
end

type NovelAllele
  template_allele::Int
  sequence::String
  n_mutations::Int8
  mutations_list::Array{Any}
  depth::Int16
  uncorrected_gaps::Array{Tuple{Int32,Int32}}
end

function EmptyNovelAllele()
  return NovelAllele(-1,"",0,[],0,[])
end

type AlleleCall
  allele::String
  flag::String
  coverage::Float16
  depth::Int
  report_txt::String
  novel_allele::NovelAllele
  alleles_to_check::Vector{Tuple{String,String, String}} # tuples of allele label, sequence, report.
end

include("db_graph.jl")

function check_files(files)
  dont_exist = [file for file in files if !isfile(file)]
  if length(dont_exist) > 0
    Lumberjack.warn("The following input file(s) could not be found: $(join(dont_exist,',')), aborting ...")
    exit(-1)
  end
end

function complement_alleles(vector, m)
  comp_vector = Int16[]
  expected::Int16 = 1
  i::Int16 = 1
  l::Int16 = length(vector)
  while i <= l
    diff = vector[i] - expected
    if diff < 0
      i += 1
      continue
    end
    while diff > 0
      push!(comp_vector, expected)
      expected += 1
      diff -= 1
    end
    i += 1
    expected +=1
  end
  while expected <= m
    push!(comp_vector, expected)
    expected += 1
  end
  return comp_vector
end

function combine_loci_classification(k, results, loci)
  # Arrays to store the kmer classification votes:
  # all kmers in one big list:
  kmer_list = String[]

  # list with elements n1, l1, n2, l2, ..., where n1 is the number of kmers on
  # the kmer_list that corresponds to the the locus l1
  loci_list = Int32[]

  # weight of each kmer:
  weight_list = Int16[]

  # n1, a_1, ..., a_n1, n2, b_1, ..., b_n2, ...
  # n1 is the number of alleles corresponding to the next kmer, a_1,...,a_n1 are
  # the alleles getting the vote for the kmer;
  alleles_list = Int16[]

  # Store the allele ids for each locus; In the DB, we index the alleles from 1 to n,
  # so we have to store the original ID to translate back in the output:
  # same format as alleles list: a list size (# of alleles), then the list of alleles;
  allele_ids_per_locus = Int[]

  # locus to int.
  l2int = Dict{String,Int16}(locus => idx for (idx, locus) in enumerate(loci))
  weight::Int16 = 0

  for (locus,(kmer_class, allele_ids, kmer_weights)) in zip(loci,results)
    n_alleles = length(allele_ids)
    half = n_alleles/2
    n_kmers_in_class = 0
    for (kmer, allele_list) in kmer_class
      l_as = length(allele_list)
      if l_as == n_alleles
        continue
      elseif l_as > half
        allele_list = complement_alleles(allele_list, n_alleles)
        if isempty(allele_list)
          continue
        end
        weight = -get(kmer_weights, kmer, 1)
      else
        weight = get(kmer_weights, kmer, 1)
      end
      # kmer as string:
      push!(kmer_list, "$(convert(DNAKmer{k}, kmer))")
      push!(weight_list, weight)
      push!(alleles_list, length(allele_list))
      append!(alleles_list, allele_list)
      n_kmers_in_class += 1
    end
    push!(loci_list, n_kmers_in_class)
    push!(loci_list, l2int[locus])
    push!(allele_ids_per_locus, n_alleles)
    append!(allele_ids_per_locus, allele_ids)
  end
  return (loci_list, weight_list, alleles_list, kmer_list, allele_ids_per_locus)
end

function kmer_class_for_each_locus(k::Int8, files::Vector{String}, compress::Bool)
  loci = [splitext(basename(file))[1] for file in files]
  results = pmap(file->build_db_graph(DNAKmer{k}, file), files)
  return results, loci
end

function save_db(k, kmer_db, loci, filename, profile, args)

  loci_list, weight_list, alleles_list, kmer_list, allele_ids_per_locus = kmer_db
  d = Dict(
    "loci_list"=> Blosc.compress(loci_list),
    "weight_list" => Blosc.compress(weight_list),
    "allele_ids_per_locus" => Blosc.compress(allele_ids_per_locus),
    "kmer_list" => join(kmer_list,""),
    "loci"=>loci,
    "args"=>JSON.json(args)
  )
  # alleles list: potentially larger than 2G, so check limit:
  limit = floor(Int,2000000000/sizeof(alleles_list[1])) # Blosc limit is 2147483631, slightly less than 2G (2147483648). Using 2000000000 just to make it easier on the eye. :)
  if sizeof(alleles_list) <= limit
    d["alleles_list"] = Blosc.compress(alleles_list)
  else
    # we have to break this into smaller chunks:
    n_chunks = ceil(Int, length(alleles_list)/limit)
    l = length(alleles_list)
    for i = 0:n_chunks-1
      name = "alleles_list_$i"
      s = i*limit + 1
      e = min((i+1)*limit,l)
      d[name] = Blosc.compress(alleles_list[s:e])
    end
  end
  # mkdir:
  mkpath(dirname(filename))
  # database:
  JLD.save(File(format"JLD", "$filename"), d)
  # Profile:
  if profile != nothing
    cp(profile, "$filename.profile", remove_destination=true)
  end
end
function open_db(filename)
  # profile:
  profile = nothing
  if isfile("$filename.profile")
    types = Array{String}[]
    open("$filename.profile") do f
      header = split(readline(f),"\t")
      for l in eachline(f)
        values = split(strip(l),"\t")
        push!(types, values)
      end
      profile = Dict("header"=>header, "types"=>types)
    end
  end

  # Compressed database, open and decompress/decode in memory:
  d = JLD.load("$filename")
  # alleles_list might be split into smaller parts:
  alleles_list = []
  if haskey(d, "alleles_list")
    alleles_list = Blosc.decompress(Int16, d["alleles_list"])
  else
    idx = 0
    while haskey(d, "alleles_list_$idx")
      append!(alleles_list,Blosc.decompress(Int16, d["alleles_list_$idx"]))
      idx += 1
    end
  end
  build_args = JSON.parse(d["args"])
  k = build_args["k"]
  loci = d["loci"]
  loci_list = Blosc.decompress(Int32, d["loci_list"])
  # weight_list = Blosc.decompress(Int8, d["weight_list"])
  weight_list = Blosc.decompress(Int16, d["weight_list"])
  allele_ids_per_locus = Blosc.decompress(Int, d["allele_ids_per_locus"])
  kmer_str = d["kmer_list"]
  # build a dict to transform allele idx (1,2,...) to original allele ids:
  loci2alleles = Dict{Int16, Vector{Int16}}(idx => Int16[] for (idx,locus) in enumerate(loci))
  allele_ids_idx = 1
  locus_idx = 1
  while allele_ids_idx < length(allele_ids_per_locus)
    n_alleles = allele_ids_per_locus[allele_ids_idx]
    append!(loci2alleles[locus_idx], allele_ids_per_locus[allele_ids_idx+1:allele_ids_idx + n_alleles])
    allele_ids_idx += 1 + n_alleles
    locus_idx += 1
  end
  # build the kmer db in the usual format:
  # kmer_classification = DefaultDict{DNAKmer{k}, Vector{Tuple{Int16, Int16, Vector{Int16}}}}(() -> Vector{Tuple{Int16, Int16, Vector{Int16}}}())
  kmer_classification = Dict{DNAKmer{k}, Vector{Tuple{Int16, Int16, Vector{Int16}}}}()
# tuple is locus idx, weight, and list of alleles;
  loci_list_idx = 1
  allele_list_idx = 1
  kmer_idx = 1
  weight_idx = 1
  # for each loci, the number of kmers and the locus idx;
  while loci_list_idx < length(loci_list)
    n_kmers = loci_list[loci_list_idx]
    locus_idx = loci_list[loci_list_idx+1]
    loci_list_idx += 2
    for i in 1:n_kmers
      # get current kmer
      kmer = DNAKmer{k}(kmer_str[kmer_idx:kmer_idx+k-1])
      kmer_idx += k
      # get list of alleles for this kmer:
      n_alleles = alleles_list[allele_list_idx]
      current_allele_list = alleles_list[allele_list_idx+1:allele_list_idx+n_alleles]
      allele_list_idx += n_alleles + 1
      weight = weight_list[weight_idx]
      weight_idx += 1
      # save in db
      if !haskey(kmer_classification, kmer)
        kmer_classification[kmer] = Tuple{Int16, Int16, Vector{Int16}}[]
      end
      push!(kmer_classification[kmer], (locus_idx, weight, current_allele_list))
    end
  end
  return kmer_classification, loci, loci2alleles, k, profile, build_args
end

function _find_profile(alleles, profile)
  if profile == nothing
    return 0,0
  end
  l = length(alleles)
    alleles_str = map(string,alleles)
  for genotype in profile["types"]
    # first element is the type id, all the rest is the actual genotype:
    if alleles_str == genotype[2:1+l]
      # assuming Clonal Complex after genotype, when present.
      return genotype[1], length(genotype) >= l + 2 ? genotype[l+2] : ""
    end
  end
  return 0,0
end

function read_alleles(fastafile, ids)
  alleles = Dict{Int, String}()
  idx = 1
  for (name, seq) in FastaReader(fastafile)
    if in(idx, ids)
      alleles[idx] = seq
    end
    idx += 1
  end
  return alleles
end


function LocusNotPresentCall()
  return AlleleCall("0", "", 0.0, 0, "Not present, no kmers found.", EmptyNovelAllele(), [])
end

function LocusPresentCall(call, depth, txt)
  return AlleleCall(call, "", 1.0, depth, txt, EmptyNovelAllele(), [])
end

function MultiplePresentCall(call, depth, txt, alleles_to_check)
  return AlleleCall(call, "+", 1.0, depth, txt, EmptyNovelAllele(), alleles_to_check)
end

function UncoveredAlleleCall(report_txt, alleles_to_check)
  # return AlleleCall("0", "?", 0.0, 0, report_txt, nothing, alleles_to_check)
  return AlleleCall("0", "?", 0.0, 0, report_txt, EmptyNovelAllele(), alleles_to_check)
end

function NovelAlleleCall(report_txt, novel_allele, alleles_to_check)
  return AlleleCall("N", "", 1.0, novel_allele.depth, report_txt, novel_allele, alleles_to_check)
end

function PartiallyCoveredAllele(call, coverage, depth, report_txt, alleles_to_check)
  return AlleleCall(call, "-", coverage, depth, report_txt, EmptyNovelAllele(), alleles_to_check)
end


function call_alleles(k, kmer_count, votes, loci_votes, loci, loci2alleles, fasta_files, kmer_thr, max_mutations, output_votes)

  # Function to call for each allele
  function call_allele(votes_per_allele, locus_votes, locus, locus2alleles, fasta_file)
    if locus_votes == 0   # 1st - if no locus votes, no presence:
      return LocusNotPresentCall()
    end
    # Votes: Dict{Int16,Dict{Int16,Int64}}; locus_idx -> { allele_idx -> votes};
    sorted_voted_alleles = sort(collect(votes_per_allele), by=x->-x[2]) # sort by max votes for this locus;
    allele_set = Set([al for (al, votes) in sorted_voted_alleles[1:min(10,end)]])
    allele_seqs = read_alleles(fasta_file, allele_set)
    # TODO: instead of 10-20 from the sort, something a bit smarter? All within some margin from the best?
    allele_coverage = [AlleleCoverage(al, votes, sequence_coverage(DNAKmer{k}, allele_seqs[al], kmer_count, kmer_thr)...) for (al, votes) in sorted_voted_alleles[1:min(10,end)]]
    # filter to find fully covered alleles:
    covered = [x for x in allele_coverage if (x.depth >= kmer_thr)] # if I remove alleles with negative votes, I might the reconstruct the same allele on the novel rebuild; better to flag output
    if length(covered) == 1
      # choose the only fully covered allele
      allele_label = locus2alleles[covered[1].allele]
      allele_votes  = covered[1].votes
      warning = allele_votes < 0 ? " Call with negative votes: $allele_votes." : ""
      return LocusPresentCall("$allele_label", covered[1].depth,  "Called allele $allele_label.$warning")
    end
    if length(covered) > 1
      # choose the most voted, but indicate that there are more fully covered that might be possible:
      allele_label = locus2alleles[covered[1].allele]
      multiple_alleles = join([locus2alleles[x.allele] for x in covered],", ")
      multiple_votes = join([x.votes for x in covered],", ")
      multiple_depth = join([x.depth for x in covered],", ")
      report_txt = "Multiple possible alleles:$multiple_alleles with depth $multiple_depth and votes $multiple_votes. Most voted ($allele_label) is chosen on call file."
      alleles_to_check = [("$(locus2alleles[x.allele])", allele_seqs[x.allele], "Multiple alleles present for this locus. Depth:$(x.depth) Votes:$(x.votes)") for x in covered]
      return MultiplePresentCall("$allele_label", covered[1].depth, report_txt, alleles_to_check)
    end
    # got here, no allele is covered; if the uncovered part is small enough, try novel allele
    # Try to find novel; sort by # of uncovered kmers, smallest to largest
    sorted_allele_coverage = sort(allele_coverage, by=x->x.uncovered_kmers)
    uncovered_kmers = sorted_allele_coverage[1].uncovered_kmers
    covered_kmers = sorted_allele_coverage[1].covered_kmers
    coverage = covered_kmers / (covered_kmers + uncovered_kmers)
    # each mutation gives around k uncovered kmers, so if smallest_uncovered > k * max_mutations, then declare missing;
    if uncovered_kmers > k * max_mutations # consider not present; TODO: better criteria?
      allele = sorted_allele_coverage[1].allele
      allele_label = "$(locus2alleles[allele])"
      depth = sorted_allele_coverage[1].depth
      report_txt = "Not present; allele $allele_label is the best covered but below threshold with $uncovered_kmers/$(covered_kmers+uncovered_kmers) missing kmers."
      return UncoveredAlleleCall(report_txt, [(allele_label, allele_seqs[allele], "Best covered allele, but declared not present, with $uncovered_kmers uncovered kmers and coverage $coverage.")])
    end

    # If gap is small, try to find a novel allele using an existing allele as template;
    # and filling the 'gaps' of coverage searching for mutations.
    # sort candidates by smaller # of gaps, then larger vote:
    sorted_by_gap_allele_coverage = sort(allele_coverage, by=x->(length(x.gap_list),-x.votes))
    best = sorted_by_gap_allele_coverage[1]
    # get all templates with same # of gaps:
    allele_templates = [x for x in sorted_by_gap_allele_coverage if length(x.gap_list) == length(best.gap_list)]
    candidate_list = [] # list with all candidate novel;
    for al_cov in allele_templates # al_cov is a struct AlleleCoverage
      template_seq = allele_seqs[al_cov.allele]
      novel_allele = correct_template(k, allele_seqs[al_cov.allele], al_cov.gap_list, kmer_count, kmer_thr, max_mutations)
      novel_allele.template_allele = locus2alleles[al_cov.allele]
      push!(candidate_list, (novel_allele, al_cov))
    end
    # sort by no gaps, then minimum mutations; then more votes
    # sorted_candidates = sort(candidate_list,by=(novel,cov)->(length(novel.uncorrected_gaps),novel.n_mutations,-cov.votes))
    sorted_candidates = sort(candidate_list,by=x->(length(x[1].uncorrected_gaps),x[1].n_mutations,-x[2].votes))
    # get best
    novel_allele, al_cov = sorted_candidates[1]
    template_allele_label = "$(novel_allele.template_allele)"

    if length(novel_allele.uncorrected_gaps) == 0
      # report:
      mutations_txt = novel_allele.n_mutations > 1 ? "$(novel_allele.n_mutations) mutations" : "$(novel_allele.n_mutations) mutation"
      mutation_desc = join([describe_mutation(ev) for ev in novel_allele.mutations_list], ", ")
      # save, closest and novel:
      alleles_to_check =[(template_allele_label, novel_allele.sequence, "Template for novel allele."),
        ("Novel", novel_allele.sequence, "$mutations_txt from allele $template_allele_label: $mutation_desc")]
      report_txt = "Novel, $mutations_txt from allele $template_allele_label: $mutation_desc"
      return NovelAlleleCall(report_txt, novel_allele, alleles_to_check)

    else # TODO: if all gaps have not been covered, just call directly uncovered? (no novel at all)
      # did not a fully corrected novel; possibly an uncovered allele, maybe partially corrected; output allele/N
      report_txt = "Partially covered alelle or novel allele; Best allele $template_allele_label has $uncovered_kmers/$(covered_kmers+uncovered_kmers) missing kmers, and no novel was found. Gaps on positions: $(join(al_cov.gap_list, ','))"
      alleles_to_check =[("$(novel_allele.template_allele)", novel_allele.sequence, "Best allele $template_allele_label, but partially covered.")]
      return PartiallyCoveredAllele("$template_allele_label", round(coverage,4), novel_allele.depth, report_txt, alleles_to_check)
    end
  end
  # Function to call all alleles by vote only:
  function call_by_vote(loci, votes, loci_votes, loci2alleles)
    ties = DefaultDict{String,Vector{Int16}}(() -> Int16[])
    vote_log = []
    best_voted_alleles = []
    for (idx,locus) in enumerate(loci)
      # if locus has no votes (no kmer present), output a zero:
      if loci_votes[idx] == 0
        push!(best_voted_alleles, "0")
        continue
      end
      # if it has votes, sort and pick the best;
      sorted_vote = sort(collect(votes[idx]), by=x->-x[2])
      push!(best_voted_alleles, loci2alleles[idx][sorted_vote[1][1]])
      votes_txt = join(["$(loci2alleles[idx][al])($votes)" for (al,votes) in sorted_vote[1:min(20,end)]],",")
      push!(vote_log, (loci_votes[idx], votes_txt))
      # find if there was a tie:
      if length(sorted_vote) > 1 && sorted_vote[1][2] == sorted_vote[2][2] # there is a tie;
        for i in 1:length(sorted_vote)
          if sorted_vote[i][2] < sorted_vote[1][2]
            break
          end
          push!(ties[locus], sorted_vote[i][1])
        end
      end
    end
    return best_voted_alleles, vote_log, ties
  end

  # get all calls:
  allele_calls = [call_allele(votes[idx], loci_votes[idx], loci[idx], loci2alleles[idx], fasta_files[idx]) for (idx, locus) in enumerate(loci)]
# allele_calls = pmap(x->call_allele(x...), [(k, kmer_count, votes[idx], loci_votes[idx], loci[idx], loci2alleles[idx], fasta_files[idx], kmer_thr, max_mutations) for (idx, locus) in enumerate(loci)])

  # Also do voting output?
  if output_votes
    voting_result = call_by_vote(loci, votes, loci_votes, loci2alleles)
  else
    voting_result = nothing
  end

  return allele_calls, voting_result
end


function write_calls(loci2alleles, allele_calls, loci, voting_result, sample, filename, profile, output_special_cases)
  # write the main call:
  allele_label_calls = ["$(call.allele)$(call.flag)" for call in allele_calls]
  st, clonal_complex = _find_profile(allele_label_calls, profile)
  open(filename, "w") do f
    header = join(vcat(["Sample"], loci, ["ST", "clonal_complex"]), "\t")
    write(f,  "$header\n")
    calls = join(vcat([sample], allele_label_calls, ["$st", clonal_complex]), "\t")
    write(f, "$calls\n")
  end
  # write the call report:
  open("$filename.coverage.txt", "w") do f
    write(f, "Locus\tCoverage\tMinKmerDepth\tCall\n")
    for (locus, call) in zip(loci, allele_calls)
      write(f, "$locus\t$(call.coverage)\t$(call.depth)\t$(call.report_txt)\n")
    end
  end
  # write the alleles with novel, missing, or multiple calls:
  alleles_to_check = [(locus, check...) for (locus, call) in zip(loci, allele_calls) for check in call.alleles_to_check if length(call.alleles_to_check) > 0]
  if output_special_cases && length(alleles_to_check) > 0
    open("$filename.special_cases.fa", "w") do f
      for (locus, al, seq, desc) in alleles_to_check
        write(f,">$(locus)_$al $desc\n$seq\n")
      end
    end
  end
  # write novel alleles:
  novel_alleles = [(locus, call.novel_allele) for (locus, call) in zip(loci, allele_calls) if call.novel_allele.template_allele != -1]
  if length(novel_alleles) > 0
    open("$filename.novel.fa", "w") do fasta
      open("$filename.novel.txt", "w") do text
        write(text, "Loci\tMinKmerDepth\tNmut\tDesc\n")
        for (locus, novel_allele) in novel_alleles
          write(fasta, ">$locus\n$(novel_allele.sequence)\n")
          mutation_desc = join([describe_mutation(ev) for ev in novel_allele.mutations_list], ", ")
          write(text, "$locus\t$(novel_allele.depth)\t$(novel_allele.n_mutations)\tFrom allele $(novel_allele.template_allele), $mutation_desc.\n")
        end
      end
    end
  end
  # write also the votes if we got them:
  if voting_result != nothing
    best_voted_alleles, vote_log, ties = voting_result
    st, clonal_complex = _find_profile(best_voted_alleles, profile)
    open("$filename.byvote", "w") do f
      header = join(vcat(["Sample"], loci, ["ST", "clonal_complex"]), "\t")
      write(f,  "$header\n")
      calls = join(vcat([sample], best_voted_alleles, ["$st", clonal_complex]), "\t")
      write(f, "$calls\n")
    end
    # write the detailed votes:
    open("$filename.votes.txt", "w") do f
      write(f, "Locus\tTotalLocusVotes\tAllele(votes),...\n")
      for (locus, data) in zip(loci, vote_log)
        write(f, "$locus\t$(join(data,'\t'))\n")
      end
    end

    # write ties:
    if length(ties) > 0
      open("$filename.ties.txt", "w") do f
        for (locus, tied_alleles) in sort(collect(ties), by=x->x[1])
          ties_txt = join(["$t" for t in tied_alleles],", ")
          write(f, "$locus\t$ties_txt\n")
        end
      end
    end
  end

end

function count_kmers_in_db_only{k}(::Type{DNAKmer{k}}, files, kmer_db)
  # Count kmers:
  kmer_count = DefaultDict{DNAKmer{k},Int}(0)
  for f in files
    istream = fastq_open(f)
    while (fq = fastq_read(istream))!=false
      for (pos, kmer) in each(DNAKmer{k}, DNASequence(fq.sequence.seq), 1)
        kmer = canonical(kmer)
        if haskey(kmer_db, kmer)
          kmer_count[kmer] += 1
        end
      end
    end
  end
  return kmer_count
end

function count_kmers{k}(::Type{DNAKmer{k}}, files)
  # Count kmers:
  kmer_count = DefaultDict{DNAKmer{k},Int}(0)
  for f in files
    istream = fastq_open(f)
    while (fq = fastq_read(istream))!=false
      for (pos, kmer) in each(DNAKmer{k}, DNASequence(fq.sequence.seq), 1)
        kmer_count[canonical(kmer)] += 1
      end
    end
  end
  return kmer_count
end

function count_votes(kmer_count, kmer_db, loci2alleles)
# function count_kmers_and_vote{k}(::Type{DNAKmer{k}}, files, kmer_db, loci2alleles)
  # Now count votes:
  # 0 votes for all alleles everyone at the start:
  votes = Dict(locus_idx => Dict{Int16, Int}(i => 0 for i in 1:length(alleles)) for (locus_idx,alleles) in loci2alleles)
  # votes per locus, to decide presence/absence:
  loci_votes = DefaultDict{Int16, Int}(0)
  for (kmer, count) in kmer_count
    if haskey(kmer_db, kmer)
      for (locus, weight, alleles) in kmer_db[kmer]
        v = weight * count
        loci_votes[locus] += abs(v)
        for allele in alleles
          votes[locus][allele] += v
        end
      end
    end
  end
  return votes, loci_votes
end

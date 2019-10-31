using BioSequences
using DataStructures: DefaultDict, OrderedDict
using TextWrap: wrap
import JLD: load, save
import Blosc
import JSON
import GZip
### Main calling function:
function run_calling_pipeline(args)
  # get samples/fastq files from command line parameters:
  sample_files = build_sample_files(args["sample_input_file"], args["1"], args["2"])

  # check if DB and fast(Q/A) files exist:
  check_files([args["db"];[fastq for fq_files in values(sample_files) for fastq in fq_files]])
  # array for saving each sample result:
  sample_results = []
  @info("Opening kmer database ... ")
  kmer_db, loci, loci2alleles, db_coverage, k, profile, build_args = open_db(args["db"])

  # prepend base path of the kmer database on to the relative paths for allele fasta files stored in the db
  for (index, file) in enumerate(build_args["fasta_files"])
    build_args["fasta_files"][index] = joinpath(dirname(args["db"]), file)
  end
  # check if scheme fasta files exist
    check_files(build_args["fasta_files"])
  # process each sample:
  for (sample, fq_files) in sample_files
    @info("Sample: $sample. Opening fastq file(s) and counting kmers ... ")
    kmer_count = count_kmers(DNAKmer{k}, fq_files, args["fasta"])
    @info("Voting for alleles ... ")
    votes, loci_votes = count_votes(kmer_count, kmer_db, loci2alleles, db_coverage)
    @info("Calling alleles and novel alleles ...")
    # If fasta are given as input, set kt as 1:
    kt = args["fasta"] ? 1 : args["kt"]
    allele_calls, voting_result = call_alleles(DNAKmer{k}, kmer_count, votes, loci_votes, loci, loci2alleles, build_args["fasta_files"], kt, args["mutation_threshold"], args["output_votes"])
    push!(sample_results, (sample, allele_calls, voting_result))
  end
  @info("Writing output ...")
  write_calls(sample_results, loci, loci2alleles, args["o"], profile, args["output_special"], args["output_votes"])

  @info("Done.")
end


## Helper structs:
struct AlleleCoverage
  allele::Int16
  votes::Int64
  depth::Int16
  covered_kmers::Int16
  uncovered_kmers::Int16
  gap_list::Vector{Tuple{Int16,Int16}}
end

mutable struct NovelAllele
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

struct AlleleCall
  allele::String
  flag::String
  coverage::Float16
  depth::Int
  report_txt::String
  novel_allele::NovelAllele
  alleles_to_check::Vector{Tuple{String,String, String}} # tuples of allele label, sequence, report.
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
  d = load("$filename")
  @info("Finished the JLD load, building alleles list...")
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
  @info("Decompressing weight list...")
  build_args = JSON.parse(d["args"])
  k = build_args["k"]
  loci = d["loci"]
  loci_list = Blosc.decompress(Int32, d["loci_list"])
  # weight_list = Blosc.decompress(Int8, d["weight_list"])
  weight_list = Blosc.decompress(Int16, d["weight_list"])
  allele_ids_per_locus = Blosc.decompress(Int, d["allele_ids_per_locus"])
  kmer_str = d["kmer_list"]
  # build a dict to transform allele idx (1,2,...) to original allele ids:
  # TODO: loci2alleles can be a vector, since indexes are 1..length(loci)
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
  @info("Building kmer index ...")
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
  # db_coverage: [i][j] -> number of kmers that cover allele j in locus i
  db_coverage = [Int16[] for (idx,locus) in enumerate(loci)]
  current_idx = 1
  if haskey(d, "allele_coverages") # might not exist for older dbs;
    allele_coverages_list = Blosc.decompress(Int16,d["allele_coverages"])
    for (idx,locus) in enumerate(loci)
      n_alleles = length(loci2alleles[idx])
      db_coverage[idx] = allele_coverages_list[current_idx:current_idx+n_alleles-1]
      current_idx += n_alleles
    end
  end
  return kmer_classification, loci, loci2alleles, db_coverage, k, profile, build_args
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
  fh = open(fastafile, "r")
  reader = FASTA.Reader(fh)
  record = FASTA.Record()
  while !eof(reader)
    read!(reader, record)
    if in(idx, ids)
      alleles[idx] = FASTA.sequence(record)
    end
    idx += 1
  end
  close(fh)
  return alleles
end

function find_allele(fastafile, wanted_seq)
  idx = 1
  reader = FASTA.Reader(open(fastafile, "r"))
  record = FASTA.Record()
  while !eof(reader)
    read!(reader, record)
    if FASTA.sequence(fastafile) == wanted_seq
      return idx
    end
    idx += 1
  end
  return nothing
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

N_SEARCH = 40 # number of alleles to search, in order of voting;
function call_alleles(::Type{DNAKmer{k}}, kmer_count, votes, loci_votes, loci, loci2alleles, fasta_files, kmer_thr, max_mutations, output_votes) where {k}

  # Function to call for each allele
  function call_allele(votes_per_allele, locus_votes, locus, locus2alleles, fasta_file)
    # given a subset of the voted alleles, find the allele call:
    function call_from_selected_allele_votes(selected_allele_votes)
      allele_set = Set([al for (al, votes) in selected_allele_votes])
      allele_seqs = read_alleles(fasta_file, allele_set)
      allele_coverage = [AlleleCoverage(al, votes, sequence_coverage(DNAKmer{k}, allele_seqs[al], kmer_count, kmer_thr)...) for (al, votes) in selected_allele_votes]
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
        novel_allele = correct_template(DNAKmer{k}, allele_seqs[al_cov.allele], al_cov.gap_list, kmer_count, kmer_thr, max_mutations)
        novel_allele.template_allele = locus2alleles[al_cov.allele]
        push!(candidate_list, (novel_allele, al_cov))
      end
      # sort by no gaps, then minimum mutations; then more votes
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
        report_txt = "Partially covered allele or novel allele; Best allele $template_allele_label has $uncovered_kmers/$(covered_kmers+uncovered_kmers) missing kmers, and no novel was found. Gaps on positions: $(join(al_cov.gap_list, ','))"
        alleles_to_check =[("$(novel_allele.template_allele)", novel_allele.sequence, "Best allele $template_allele_label, but partially covered.")]
        return PartiallyCoveredAllele("$template_allele_label", round(coverage,digits=4), 0, report_txt, alleles_to_check) # if it's not 100% covered, then depth=0
      end
    end

    # if locus_votes == 0   # 1st - if no locus votes, no presence:
    #   return LocusNotPresentCall()
    # end
    # Votes: Dict{Int16,Dict{Int16,Int64}}; locus_idx -> { allele_idx -> votes};

    sorted_voted_alleles = sort(collect(votes_per_allele), by=x->-x[2]) # sort by max votes for this locus;
    tie_idx = 1
    while tie_idx < length(sorted_voted_alleles) && sorted_voted_alleles[tie_idx+1][2] == sorted_voted_alleles[1][2]
      tie_idx += 1
    end
    # tied_best = [x for x in sorted_voted_alleles if x[2] == sorted_voted_alleles[1][2]]
    tie_idx +=10 # add a bit more
    tied_best = sorted_voted_alleles[1:min(tie_idx,end)]
    allele_call = call_from_selected_allele_votes(tied_best)

    # if it is novel; check if it already exists in DB:
    # if allele_call.allele == "N"
    #   allele_idx = find_allele(fastafile, wanted_seq)
    #   if exists != nothing

    # If I did not get a full covered call, or novel, try looking at ALL alleles
    if (allele_call.depth == 0 || allele_call.allele == "N") && (tie_idx <= length(sorted_voted_alleles))
      second_call = call_from_selected_allele_votes(sorted_voted_alleles[tie_idx:end])
      if second_call.depth == 0
        return allele_call
      else
        return second_call
      end
    else
      return allele_call
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
        push!(vote_log, ("0","No kmers matches on this locus."))
        continue
      end
      # if it has votes, sort and pick the best;
      sorted_vote = sort(collect(votes[idx]), by=x->-x[2])
      push!(best_voted_alleles, loci2alleles[idx][sorted_vote[1][1]])
      votes_txt = join(["$(loci2alleles[idx][al])($votes)" for (al,votes) in sorted_vote[1:min(N_SEARCH,end)]],",")
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

  # Also do voting output?
  if output_votes
    voting_result = call_by_vote(loci, votes, loci_votes, loci2alleles)
  else
    voting_result = nothing
  end

  return allele_calls, voting_result
end


function write_calls(sample_results, loci, loci2alleles, filename, profile, output_special_cases, output_votes)
# function write_calls(loci2alleles, allele_calls, loci, voting_result, sample, filename, profile, output_special_cases)
  # write the main call:
  open(filename, "w") do f
    header = join(vcat(["Sample"], loci, ["ST", "clonal_complex"]), "\t")
    write(f,  "$header\n")
    for (sample, allele_calls, voting_result) in sample_results
      allele_label_calls = ["$(call.allele)$(call.flag)" for call in allele_calls]
      st, clonal_complex = _find_profile(allele_label_calls, profile)
      calls = join(vcat([sample], allele_label_calls, ["$st", clonal_complex]), "\t")
      write(f, "$calls\n")
    end
  end
  # write the call report:
  open("$filename.coverage.txt", "w") do f
    write(f, "Sample\tLocus\tCoverage\tMinKmerDepth\tCall\n")
    for (sample, allele_calls, voting_result) in sample_results
      for (locus, call) in zip(loci, allele_calls)
        write(f, "$sample\t$locus\t$(call.coverage)\t$(call.depth)\t$(call.report_txt)\n")
      end
    end
  end
  # write the alleles with novel, missing, or multiple calls:
  alleles_to_check = [(locus, allele_to_check...) for (sample, allele_calls, voting_result) in sample_results for (locus, call) in zip(loci, allele_calls) for allele_to_check in call.alleles_to_check if length(call.alleles_to_check) > 0]
  if output_special_cases && length(alleles_to_check) > 0
    open("$filename.special_cases.fa", "w") do f
      # TODO: deal with repeated alleles;
      for (locus, al, seq, desc) in alleles_to_check
        write(f,">$(locus)_$al $desc\n$seq\n")
      end
    end
  end
  # write novel alleles
  # seen_novel = DefaultDict{String, Set{String}}(() -> Set{String}()) # locus to set of seen novel sequences
  times_seen = DefaultDict{String, Int}(0) # times that a novel DNA sequence was seen; (DNA -> count)
  novel_id = Dict{String, Int}() # novel id of a DNA seq (DNA -> id)
  current_novel_id = DefaultDict{String, Int}(0) # current novel id for a given locus. (locus -> current_id)
  novel_to_fasta = []
  open("$filename.novel.txt", "w") do text
    write(text, "Sample\tLocus\tNovel_id\tMinKmerDepth\tNmut\tDesc\n")
    for (sample, allele_calls, voting_result) in sample_results
      # println(novel_to_fasta)
      # println("------------")
      novel_alleles = [(locus, call.novel_allele) for (locus, call) in zip(loci, allele_calls) if call.novel_allele.template_allele != -1]
      # print("novel_allele: ")
      # println(novel_alleles)
      for (locus, novel_allele) in novel_alleles
        mutation_desc = join([describe_mutation(ev) for ev in novel_allele.mutations_list], ", ")
        # check if this novel was seen before:
        if (times_seen[novel_allele.sequence] == 0)
          # println("new novel!")
          current_novel_id[locus] += 1
          # println(locus)
          # println(current_novel_id[locus])
          push!(novel_to_fasta, (locus, novel_allele.sequence))
          novel_id[novel_allele.sequence] = current_novel_id[locus]
          write(text, "$sample\t$locus\tN$(novel_id[novel_allele.sequence])\t$(novel_allele.depth)\t$(novel_allele.n_mutations)\tFrom allele $(novel_allele.template_allele), $mutation_desc.\n")
        end
        times_seen[novel_allele.sequence] += 1
        # println(">>> Times we seen this allele")
        # println(times_seen[novel_allele.sequence])
        # println("--------")
        # println(novel_allele.sequence)
        # println(novel_id[novel_allele.sequence])

      end
    end
  end
  # now write the fasta:
  # println("<><><><><><><><><><><><><><><><><><>")
  # println(novel_to_fasta)
  open("$filename.novel.fa", "w") do f
    for (locus, seq) in novel_to_fasta
      desc = "Seen in $(times_seen[seq]) sample(s)"
      write(f, ">$(locus)_N$(novel_id[seq]) $desc.\n$(wrap(seq, width=120))\n")
    end
  end

  # TODO: with multiple samples, we have to check how many
  # novel_alleles = [(locus, call.novel_allele) for (locus, call) in zip(loci, allele_calls) if call.novel_allele.template_allele != -1]
  # if length(novel_alleles) > 0
  #   open("$filename.novel.fa", "w") do fasta
  #       for (locus, novel_allele) in novel_alleles
  #         write(fasta, ">$locus\n$(novel_allele.sequence)\n")
  #         mutation_desc = join([describe_mutation(ev) for ev in novel_allele.mutations_list], ", ")
  #         write(text, "$locus\t$(novel_allele.depth)\t$(novel_allele.n_mutations)\tFrom allele $(novel_allele.template_allele), $mutation_desc.\n")
  #       end
  #     end
  #   end
  # end
  # write also the votes if we got them:
  if output_votes
    f_vote_call = open("$filename.byvote", "w")
    f_details = open("$filename.votes.txt", "w")
    f_ties = open("$filename.ties.txt", "w")
    # headers:
    write(f_vote_call,  "$(join(vcat(["Sample"], loci, ["ST", "clonal_complex"]), "\t"))\n")
    write(f_details, "Sample\tLocus\tTotal locus votes\tAllele(relative votes),...\n")
    write(f_ties, "Sample\tLocus\tTied Alleles\n")
    # loop per sample:
    for (sample, allele_calls, voting_result) in sample_results
      # write the voting call:
      best_voted_alleles, vote_log, ties = voting_result
      st, clonal_complex = _find_profile(best_voted_alleles, profile)
      calls = join(vcat([sample], best_voted_alleles, ["$st", clonal_complex]), "\t")
      write(f_vote_call, "$calls\n")
      # write the detailed votes:
      for (locus, data) in zip(loci, vote_log)
        write(f_details, "$sample\t$locus\t$(join(data,'\t'))\n")
      end
      # write ties:
      if length(ties) > 0
        for (locus, tied_alleles) in sort(collect(ties), by=x->x[1])
          ties_txt = join(["$t" for t in sort(tied_alleles)],", ")
          write(f_ties, "$sample\t$locus\t$ties_txt\n")
        end
      end
    end
    close(f_vote_call)
    close(f_details)
    close(f_ties)
  end
end

function count_kmers_in_fasta_file(::Type{DNAKmer{k}}, fasta_file, kmer_count) where {k}
  reader = FASTA.Reader(open(fasta_file, "r"))
	record = FASTA.Record()
  while !eof(reader)
      read!(reader, record)
      for (pos, kmer) in each(DNAKmer{k}, FASTA.sequence(record))
        kmer_count[canonical(kmer)] += 1
      end
    end
end

function count_kmers_in_fastq_file(::Type{DNAKmer{k}}, fastq_file, kmer_count) where {k}
  reader = FASTQ.Reader(endswith(fastq_file, ".gz") ? GZip.open(fastq_file, "r") : open(fastq_file, "r"))
  record = FASTQ.Record()
  while !eof(reader)
      read!(reader, record)
      for (pos, kmer) in each(DNAKmer{k}, FASTQ.sequence(record))
        kmer_count[canonical(kmer)] += 1
      end
    end
end

function count_kmers(::Type{DNAKmer{k}}, files, is_fasta=false) where {k}
  # Count kmers:
  kmer_count = DefaultDict{DNAKmer{k},Int}(0)
  for f in files
    # check type:
    if is_fasta
      count_kmers_in_fasta_file(DNAKmer{k}, f, kmer_count)
    else
      count_kmers_in_fastq_file(DNAKmer{k}, f, kmer_count)
    end
  end
  return kmer_count
end

function count_votes(kmer_count, kmer_db, loci2alleles, db_coverage)
# function count_kmers_and_vote{k}(::Type{DNAKmer{k}}, files, kmer_db, loci2alleles)
  # Now count votes:
  # 0 votes for all alleles everyone at the start:
  # TODO: also list, instead of DICT?
  votes = Dict(locus_idx => Dict{Int16, Int}(i => 0 for i in 1:length(alleles)) for (locus_idx,alleles) in loci2alleles)
  # votes per locus, to decide presence/absence:
  loci_votes = DefaultDict{Int16, Int}(0)
  # for (locus,) in loci2alleles
  #     loci_votes[locus] = 0
  # end
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
  # adjust for allele coverage:
  # for (locus, allele_votes) in votes
  #   for al in keys(allele_votes)
  #     allele_votes[al] = Int(round(allele_votes[al]/db_coverage[locus][al]))
  #   end
  # end
  return votes, loci_votes
end

# Calling:
function sequence_coverage(::Type{DNAKmer{k}}, sequence, kmer_count, kmer_thr=6) where {k}
  # TODO: not using the gap list anymore, just for debug. Eventually delete!
  # find the minimum coverage depth of all kmers of the given variant, and
  # how many kmers are covered/uncovered
  gap_list = Tuple{Int,Int}[]
  smallest_coverage::Int16 = 32766
  covered = 0
  uncovered = 0

  in_gap = false
  gap_start = 0

  for (pos, kmer) in each(DNAKmer{k}, DNASequence(sequence))
    can_kmer = canonical(kmer)
    cov = get(kmer_count, can_kmer, 0)
    if cov >= kmer_thr
      covered += 1
      # gap:
      if in_gap
        push!(gap_list, (gap_start, pos-1))
        in_gap = false
      end
    else
      uncovered += 1
      # gap
      if !in_gap
        gap_start = pos
        in_gap = true
      end
    end
    if cov < smallest_coverage
      smallest_coverage = cov
    end
  end
  if in_gap
    push!(gap_list, (gap_start, length(sequence)-k+1))
  end
  # merge too close gaps: if by chance we get a kmer match in the middle of a gap, this will break a gap in two, even
  # with just one mutation, giving errors. TODO: How to prevent that?
  # simple idea here, merge gaps (s1,e1) and (s2,e2) if s1 + k >= s2;
  merged_gap_list = []
  i = 1
  while i <= length(gap_list)
    # if too close to the next, merge:
    if i<length(gap_list) && gap_list[i][1] + k >= gap_list[i+1][1]
      push!(merged_gap_list, (gap_list[i][1], gap_list[i+1][2]))
      i += 2
    else
      push!(merged_gap_list, gap_list[i])
      i += 1
    end
  end
  return smallest_coverage, covered, uncovered, merged_gap_list
end

function cover_sequence_gap(::Type{DNAKmer{k}}, sequence, kmer_count, kmer_thr=8, max_mutations=4) where {k}
  gap_cover_list = []
  found = Set{String}()
  function scan_variant(n_mut, sequence, events, start_pos=1)
    if n_mut > max_mutations
      return
    end
    last_kmer_present = false
    for (pos, kmer) in each(DNAKmer{k}, DNASequence(sequence))
      if pos < start_pos
        continue
      end
      can_kmer = canonical(kmer)
      if get(kmer_count, can_kmer, 0) >= kmer_thr
        if pos > start_pos +1 && (!last_kmer_present)
          # found a kmer and the last was not there; check if I can find kmer variations for the previous position
          km = "$kmer"
          for base in ["A","C","T","G"]
            if get(kmer_count, canonical(DNAKmer{k}("$base$(km[1:end-1])")), 0) >= kmer_thr
              # found a kmer variant;
              # Test one more kmer before adding to the pool:
              # substitution at position pos-1
              substitution_seq = sequence[1:pos-2] * base * sequence[pos:end]
              push!(candidate_list, (n_mut+1, substitution_seq, [events; [("S",pos-1,"$(sequence[pos-1])->$base")]], start_pos))
              # insertion:
              insertion_seq = sequence[1:pos-1] * base * sequence[pos:end]
              push!(candidate_list, (n_mut+1, insertion_seq, [events; [("I",pos-1,base)]], start_pos))
              # deletion:
              for i in 0:2
                if pos-2-i < 1
                  break
                end
                if sequence[pos-2-i:pos-2-i] == base
                  deletion_seq = sequence[1:pos-2-i] * sequence[pos:end]
                  push!(candidate_list, (n_mut+1+i, deletion_seq, [events; [("D",pos-1,"$(i+1)")]], start_pos))
                  break
                end
              end
              # this variant has mismatches, so just return;
              return
            end
          end
        end
        last_kmer_present = true
      else
        if last_kmer_present
          # last_kmer_present was present, but this not; check kmer variations for nest position
          km = "$kmer"
          for base in ["A","C","T","G"]
            if get(kmer_count, canonical(DNAKmer{k}("$(km[1:end-1])$base")), 0) >= kmer_thr
              # println("Found a $base at pos $(pos+k-1)")
              # substitution at position pos+k-1;
              substitution_seq = sequence[1:pos+k-2] * base * sequence[pos+k:end]
              push!(candidate_list, (n_mut+1, substitution_seq, [events; [("S",pos+k-1,"$(sequence[pos+k-1])->$base")]], pos-1))
              # insertion:
              insertion_seq = sequence[1:pos+k-2] * base * sequence[pos+k-1:end]
              push!(candidate_list, (n_mut+1, insertion_seq, [events; [("I",pos+k-1,base)]], pos-1))
              # deletion:
              for i in 0:2
                # println("trying $(pos+k+i), found $(sequence[pos+k+i]) (want $base) $(sequence[pos+k+i:pos+k+i] == base).")
                if pos+k+i > length(sequence)
                  break
                end
                if sequence[pos+k+i:pos+k+i] == base
                  # println("Found a possible deletion, pos $(pos+k-1) to $(pos+k+i-1)")
                  deletion_seq = sequence[1:pos+k-2] * sequence[pos+k+i:end]
                  push!(candidate_list, (n_mut+1+i, deletion_seq, [events; [("D",pos+k-1,"$(i+1)")]], pos-1))
                  break
                end
              end
              # the original sequence has mismatches, so just return;
              return
            end
          end
        end
        last_kmer_present = false
      end
    end
    # passed without suggesting mutations, test this variant:
    seq_cov = sequence_coverage(DNAKmer{k}, sequence, kmer_count)[1]
    if seq_cov[1] >= kmer_thr
      if !in(sequence, found)
        push!(gap_cover_list, (n_mut, sequence, events, seq_cov))
        push!(found, sequence)
        max_mutations = n_mut
      end
    end
  end
  candidate_list = [(0, sequence, [], 0)]
  idx = 1
  while idx <= length(candidate_list)
    n_mut, candidate_sequence, mutations, start_pos = candidate_list[idx]
    scan_variant(n_mut, candidate_sequence, mutations, start_pos)
    idx += 1
  end

  if length(gap_cover_list) > 0
    # sort and return the one with minimum mutations: TODO: tie?
    return sort(gap_cover_list, by=x->x[1])[1]
  else
    return nothing
  end

end
# struct CorrectedTemplate
#
# end
function correct_template(::Type{DNAKmer{k}}, template_seq, gap_list, kmer_count, kmer_thr, max_mutations) where {k}
  function find_next_gap(sequence, skip=1)
    in_gap = false
    gap_start = 0
    for (pos, kmer) in each(DNAKmer{k}, DNASequence(sequence))
      if pos < skip
        continue
      end
      can_kmer = canonical(kmer)
      cov = get(kmer_count, can_kmer, 0)
      if cov >= kmer_thr
        if in_gap
          return (gap_start, pos-1)
        end
      else
        if !in_gap
          gap_start = pos
          in_gap = true
        end
      end
    end # loop
    if in_gap
      return (gap_start, length(sequence)-k+1)
    end
    return nothing
  end # end auxiliar function find_next_gap

  # start correction:
  current_skip = 1
  uncorrected_gaps = []
  corrected_seq = template_seq
  total_mut = 0
  max_depth = 0
  mutations_list = []
  while current_skip < length(corrected_seq)
    gap = find_next_gap(corrected_seq, current_skip)
    if gap == nothing
      break
    end
    st_pos, end_pos = gap
    adj_start = max(st_pos-1, 1)
    adj_end = min(end_pos+k, length(corrected_seq))
    gap_seq = corrected_seq[adj_start:adj_end]
    gap_cover = cover_sequence_gap(DNAKmer{k}, gap_seq, kmer_count, kmer_thr, max_mutations)
    if gap_cover == nothing
      push!(uncorrected_gaps, (st_pos, end_pos))
      current_skip = end_pos + 1
    else
      n_mut, gap_cover_seq, mut_list, depth = gap_cover
      total_mut += n_mut
      if depth > max_depth
        max_depth = depth
      end
      for mut in mut_list
        # TODO: correct idx
        push!(mutations_list, (mut[1], mut[2] + adj_start - 1, mut[3]))
      end
      corrected_seq = corrected_seq[1:adj_start-1] * gap_cover_seq * corrected_seq[adj_end+1:end]
      current_skip = adj_start + length(gap_cover_seq) - k
    end
  end
  # 1st field, template_allele, is filled later:
  return NovelAllele(-1, corrected_seq, total_mut, mutations_list, max_depth, uncorrected_gaps)
end

function describe_mutation(mut)
  mut, pos, desc = mut
  if mut == "D"
    return "Del of len $desc at pos $pos"
  elseif mut == "I"
    return "Ins of base $desc at pos $pos"
  elseif mut == "S"
    return "Subst $desc at pos $pos"
  else
    return ""
  end
end


### input file helper calling_functions

function build_sample_files(input_file, forward_files, reverse_files)
  # First, sanity checks
  # at least one input option:
  if all([x == nothing for x in [input_file, forward_files, reverse_files]])
    exit_error("Please give at least one input file using the parameters -m, --1 and/or --2.")
  end
  #
  if length(reverse_files)>0 && length(forward_files) != length(reverse_files)
    exit_error("Forward and reverse input file does not match, got $(length(forward_files)) forward and $(length(reverse_files)) reverse.")
  end

  # Create an ordered dict to get all samples. Combinations of -m and -1,-2 are allowed.
  sample_files = OrderedDict{String, Vector{String}}()
  # Check
  if input_file != nothing
      check_files([input_file])
      open(input_file) do f
        for ln in eachline(f)
          sample, fastq = split(strip(ln))
          if !haskey(sample_files, sample)
            sample_files[sample] = String[]
          end
          push!(sample_files[sample], fastq)
        end
      end
  end
  if length(reverse_files) == 0 # single files only
    for fw_file in forward_files
      sample_files[remove_fastq_ext(fw_file)] = [fw_file]
    end
  end
  # build sample files for fw and rev:
  for (fw, rev) in zip(forward_files, reverse_files)
    sample = lcp([fw,rev]) # largest common prefix, removing special symbols at the end, will be the name of the sample.
    if sample == "" # if empty, fastq names dont match, probably wrong pairing.
      exit_error("No match between forward and reverse files $fw and $bw, please check the input options --1 and --2.")
    end
    sample_files[sample] = [fw,rev]
  end
  return sample_files
end

# longest common prefix
function lcp(str::Vector{String})
  r = IOBuffer()
  i = 1
  while all(i <= length(s) for s in str) && all(s == str[1][i] for s in getindex.(str, i))
    print(r, str[1][i])
    i += 1
  end
  return strip(basename(String(take!(r))),['_','.','-'])
end

# # Remove fastq
function remove_fastq_ext(str)
  return replace(str,  r"\.(fastq|fq)(\.gz)?" => "")
end

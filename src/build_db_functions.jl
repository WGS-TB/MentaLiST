using Suppressor
@suppress_err begin
using Bio.Seq
using DataStructures
import GZip
import JLD
import FileIO: File, @format_str
import Blosc
using OpenGene
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

  # +1 and -1 corresponding to the weight of the kmer:
  # weight_list = Int8[]
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
  # weight::Int8 = 0
  weight::Int16 = 0

  for (locus,kmer_class, allele_ids, kmer_weights) in results
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
      push!(kmer_list, "$kmer")
      push!(weight_list, weight)
      push!(alleles_list, length(allele_list))
      append!(alleles_list, allele_list)
      n_kmers_in_class += 1
    end
    # @show n_kmers_in_class
    push!(loci_list, n_kmers_in_class)
    push!(loci_list, l2int[locus])
    push!(allele_ids_per_locus, n_alleles)
    append!(allele_ids_per_locus, allele_ids)
  end
  return (loci_list, weight_list, alleles_list, kmer_list, allele_ids_per_locus)
end

function kmer_class_for_each_locus(k::Int8, files::Vector{String}, compress::Bool)
  results = []
  loci = String[]
  for file in files
    locus::String = splitext(basename(file))[1]
    push!(loci, locus)
    kmer_class, allele_ids, kmer_weights = kmer_class_for_locus(DNAKmer{k}, file, compress)
    push!(results, (locus,kmer_class, allele_ids, kmer_weights))
  end
  return results, loci
end

function kmer_class_for_locus{k}(::Type{DNAKmer{k}}, fastafile::String, compress::Bool)
  allowed_kmers = Dict()
  if compress
    # Find db graph contigs, and get 1st kmer of each:
    allowed_kmers = db_graph_contig_kmers(DNAKmer{k}, [fastafile])
  end
  record = FASTASeqRecord{BioSequence{DNAAlphabet{2}}}()
  kmer_class = DefaultDict{DNAKmer{k}, Vector{Int16}}(() -> Int16[])
  allele_ids = Int16[]
  allele_idx::Int16 = 1
  open(FASTAReader{BioSequence{DNAAlphabet{2}}}, fastafile) do reader
      while !eof(reader)
          read!(reader, record)
          for (pos, kmer) in each(DNAKmer{k}, record.seq)
            can_kmer = canonical(kmer)
            # if !compress || (can_kmer in allowed_kmers)
            if !compress || (haskey(allowed_kmers,can_kmer))
              push!(kmer_class[can_kmer], allele_idx)
            end
          end
          # update idx; the counter idx is incremental (1,2, ...) because we need the array sorted.
          # But this is not always sin the allele ordering, so we have to save the original id to restore it later;
          # find the separator; will assume that if I see a "_", that's it, otherwise try "-";
          separator = in('_', record.name) ? "_" : "-"
          allele_id = parse(Int16,split(record.name, separator)[end])
          push!(allele_ids, allele_id)
          allele_idx += 1
      end
  end
  return kmer_class, allele_ids, allowed_kmers
end

function save_db(k, kmer_db, loci, filename, profile)
  loci_list, weight_list, alleles_list, kmer_list, allele_ids_per_locus = kmer_db
  d = Dict(
    "loci_list"=> Blosc.compress(loci_list),
    "weight_list" => Blosc.compress(weight_list),
    "allele_ids_per_locus" => Blosc.compress(allele_ids_per_locus),
    "kmer_list" => join(kmer_list,""),
    "k"=>k,
    "loci"=>loci
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
  k = d["k"]
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
  kmer_classification = DefaultDict{DNAKmer{k}, Vector{Tuple{Int16, Int16, Vector{Int16}}}}(() -> Vector{Tuple{Int16, Int16, Vector{Int16}}}())
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
      push!(kmer_classification[kmer], (locus_idx, weight, current_allele_list))
    end
  end
  return kmer_classification, loci, loci2alleles, k, profile
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

function write_calls(votes, loci_votes, loci, loci2alleles, sample, filename, profile)
  # get best voted:
  best_voted_alleles = [sort(collect(votes[idx]), by=x->-x[2])[1][1] for (idx,locus) in enumerate(loci)]
  # translate alleles to original id:
  best_voted_alleles = [loci2alleles[locus_idx][best] for (locus_idx, best) in enumerate(best_voted_alleles)]
  # check if there is a type on the database:
  genotype, clonal_complex = _find_profile(best_voted_alleles, profile)
  # write:
  open(filename, "w") do f
    header = join(vcat(["Sample"], loci, ["ST", "clonal_complex"]), "\t")
    write(f,  "$header\n")
    calls = join(vcat([sample], best_voted_alleles, ["$genotype", clonal_complex]), "\t")
    write(f, "$calls\n")
  end
  # votes, find ties:
  ties = Dict{String,Vector{Int16}}()
  open("$filename.votes.txt", "w") do f
    # write(f, "Locus\tAllele(votes),...\n")
    write(f, "Locus\tTotalVotes\tTopVote\tAllele(votes),...\n")
    for (idx,locus) in enumerate(loci)
      max_idx = min(length(votes[idx]),10)
      sorted_vote = sort(collect(votes[idx]), by=x->-x[2])
      top_vote = sorted_vote[1][2] + (sorted_vote[end][2] < 0 ? abs(sorted_vote[end][2]) : 0)
      votes_txt = join(["$(loci2alleles[idx][a])($b)" for (a,b) in sorted_vote],",")
      # votes_txt = join(["$a($b)" for (a,b) in sorted_vote[1:max_idx]],",")
      write(f, "$locus\t$(loci_votes[idx])\t$top_vote\t$votes_txt\n")
      # ties:
      if length(sorted_vote) > 1 && sorted_vote[1][2] == sorted_vote[2][2]
        # find all ties:
        tied_val = sorted_vote[1][2]
        current = 1
        ties[locus] = Int16[]
        while (current <= length(sorted_vote) && sorted_vote[current][2] == tied_val)
          push!(ties[locus], sorted_vote[current][1])
          current += 1
        end
      end
    end
  end
  # write ties:
  open("$filename.ties.txt", "w") do f
    for (locus, tied_alleles) in sort(collect(ties), by=x->x[1])
      ties_txt = join(["$t" for t in tied_alleles],", ")
      write(f, "$locus\t$ties_txt\n")
    end
  end
end

function count_kmers_and_vote{k}(::Type{DNAKmer{k}}, files, kmer_db, loci2alleles)
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

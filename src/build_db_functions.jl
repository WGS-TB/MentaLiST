using Bio.Seq
using DataStructures
import GZip
import JLD
import Blosc

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
  kmer_list = String[]
  loci_list = Int[]
  weight_list = Int8[]
  alleles_list = Int16[]
  n_alleles_list = Int16[]
  # locus to int:
  l2int = Dict{String,Int16}(locus => idx for (idx, locus) in enumerate(loci))
  weight::Int8 = 0

  for (locus, kmer_class, idx_to_allele_id) in results
    n_alleles = length(idx_to_allele_id)
    # check if the relabeling is not sorted, which means we need to relabel the alleles;
    relabel = !issorted(idx_to_allele_id)
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
        weight = -1
      else
        weight = 1
      end
      # relabel the allele_list if needed:
      if relabel
        allele_list = [idx_to_allele_id[i] for i in allele_list]
      end
      push!(kmer_list, "$kmer")
      push!(weight_list, weight)
      push!(alleles_list, length(allele_list))
      append!(alleles_list, allele_list)
      n_kmers_in_class += 1
    end
    push!(loci_list, n_kmers_in_class)
    push!(loci_list, l2int[locus])
    push!(n_alleles_list, n_alleles)
  end
  return (loci_list, weight_list, alleles_list, kmer_list, n_alleles_list)
end

function kmer_class_for_each_locus(k::Int8, files::Vector{String})
  results = []
  loci = String[]
  for file in files
    locus::String = splitext(basename(file))[1]
    push!(loci, locus)
    kmer_class, idx_to_allele_id = kmer_class_for_locus(DNAKmer{k}, file)
    push!(results, (locus,kmer_class, idx_to_allele_id))
  end
  return results, loci
end

function kmer_class_for_locus{k}(::Type{DNAKmer{k}}, fastafile::String)
  record = FASTASeqRecord{BioSequence{DNAAlphabet{2}}}()
  kmer_class = DefaultDict{DNAKmer{k}, Vector{Int16}}(() -> Int16[])
  idx_to_allele_id = Int16[]
  allele_idx::Int16 = 1
  open(FASTAReader{BioSequence{DNAAlphabet{2}}}, fastafile) do reader
      while !eof(reader)
          read!(reader, record)
          for (pos, kmer) in each(DNAKmer{k}, record.seq)
            push!(kmer_class[canonical(kmer)], allele_idx)
          end
          # update idx; the counter idx is incremental (1,2, ...) because we need the array sorted.
          # But this is not always sin the allele ordering, so we have to save the original id to restore it later;
          allele_id = parse(Int16,split(record.name, "_")[2])
          push!(idx_to_allele_id, allele_id)
          allele_idx += 1
      end
  end
  return kmer_class, idx_to_allele_id
end

function save_db(k, kmer_db, loci, filename)
  loci_list, weight_list, alleles_list, kmer_list, n_alleles_list = kmer_db
  d = Dict(
    "k"=>k,
    "loci"=>loci,
    "loci_list"=> Blosc.compress(loci_list),
    "weight_list" => Blosc.compress(weight_list),
    "alleles_list" => Blosc.compress(alleles_list),
    "n_alleles_list" => Blosc.compress(n_alleles_list),
    "kmer_list" => join(kmer_list,"")
  )
  JLD.save("$filename.jld", d)
end
function open_db(filename)
  d = JLD.load("$filename.jld")
  k = d["k"]
  loci = d["loci"]
  alleles_list = Blosc.decompress(Int16, d["alleles_list"])
  loci_list = Blosc.decompress(Int, d["loci_list"])
  weight_list = Blosc.decompress(Int8, d["weight_list"])
  n_alleles_list = Blosc.decompress(Int16, d["n_alleles_list"])
  kmer_str = d["kmer_list"]
  # build the kmer db in the usual format:
  kmer_classification = DefaultDict{DNAKmer{k}, Vector{Tuple{Int16, Int8, Vector{Int16}}}}(() -> Vector{Tuple{Int16, Int8, Vector{Int16}}}())
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
  return kmer_classification, loci, n_alleles_list, k
end

function write_calls(votes, loci, sample, filename)
  best_votes_alleles = [sort(collect(votes[idx]), by=x->-x[2])[1][1] for (idx,locus) in enumerate(loci)]
  open(filename, "w") do f
    header = join(vcat(["Sample"], loci, ["ClonalComplex"]), "\t")
    write(f,  "$header\n")
    # todo: look in the database for the type, to assing a 'ClonalComplex'
    calls = join(vcat([sample], best_votes_alleles, ["0"]), "\t")
    write(f, "$calls\n")
  end
  # debug votes:
  open("$filename.votes.txt", "w") do f
    for (idx,locus) in enumerate(loci)
      sorted_vote = sort(collect(votes[idx]), by=x->-x[2])[1:10]
      write(f, "$locus:$sorted_vote\n")
    end
  end
end

function get_votes_for_sequence{k}(::Type{DNAKmer{k}}, seq, kmer_db, threshold=10, prefilter=false)
  if length(seq) < k
    return false, false
  end
  if prefilter  # stringMLST like pre_filter
    try
      half = div(length(seq)-k,2)
      testmer = DNAKmer{k}(seq[half:half+k-1])
      if !haskey(kmer_db, testmer)
        return false, false
      end
    catch LoadError
       # if it has an ambiguous base (N), throws an exception; I will not filter in this case
     end
  end
  # count all votes:
  votes = Dict()
  # count locus hits, for filtering
  locus_hits = DefaultDict{Int16,Int16}(0)
  for (pos, kmer) in each(DNAKmer{k}, DNASequence(seq), 1)
    kmer = canonical(kmer)
    if haskey(kmer_db, kmer)
      for (locus, val, alleles) in kmer_db[kmer]
        locus_hits[locus] += 1
        for allele in alleles
          if !haskey(votes, locus)
            votes[locus] = DefaultDict{Int16, Int16}(0)
          end
          votes[locus][allele] += val
        end
      end
    end
  end
  if isempty(votes)
    return false, false
  end
  # check if locus hits are above the threshold for returning votes:
  best_locus, n_hits = sort(collect(locus_hits), by=x->-x[2])[1]
  if n_hits > threshold
    return true, (best_locus, votes[best_locus])
  else
    return false, false
  end
end

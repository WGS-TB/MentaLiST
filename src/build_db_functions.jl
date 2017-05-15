using Bio.Seq
using DataStructures
import GZip
import JLD
import Blosc

function complement_alleles(orig_set, m)
  # array = collect(orig_set)
  # c_set = OrderedSet{Int16}()
  array = orig_set
  c_set = Int16[]

  expected::Int16 = 1
  i::Int16 = 1
  l::Int16 = length(array)
  while i <= l
    diff = array[i] - expected
    if diff < 0
      i += 1
      continue
    end
    while diff > 0
      push!(c_set, expected)
      expected += 1
      diff -= 1
    end
    i += 1
    expected +=1
  end
  while expected <= m
    push!(c_set, expected)
    expected += 1
  end
  return c_set
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

  for (locus, kmer_class, n_alleles) in results
    half = n_alleles/2
    n_kmers_in_class = 0
    for (kmer, allele_set) in kmer_class
      l_as = length(allele_set)
      if l_as == n_alleles
        continue
      elseif l_as > half
        allele_set = complement_alleles(allele_set, n_alleles)
        if isempty(allele_set)
          continue
        end
        weight = -1
      else
        weight = 1
      end
      push!(kmer_list, "$kmer")
      push!(weight_list, weight)
      push!(alleles_list, length(allele_set))
      append!(alleles_list, allele_set)
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
    kmer_class, n_alleles = kmer_class_for_locus(DNAKmer{k}, file)
    push!(results, (locus,kmer_class, n_alleles))
  end
  return results, loci
end

function kmer_class_for_locus{k}(::Type{DNAKmer{k}}, fastafile::String)
  record = FASTASeqRecord{BioSequence{DNAAlphabet{2}}}()
  # kmer_class = DefaultDict{DNAKmer{k}, OrderedSet{Int16}}(() -> OrderedSet{Int16}())
  kmer_class = DefaultDict{DNAKmer{k}, Vector{Int16}}(() -> Int16[])
  length_alleles::Int16 = 0
  open(FASTAReader{BioSequence{DNAAlphabet{2}}}, fastafile) do reader
    allele_idx::Int16 = 1
      while !eof(reader)
          read!(reader, record)
          allele_idx = parse(Int16,split(record.name, "_")[2])
          length_alleles += 1
          for (pos, kmer) in each(DNAKmer{k}, record.seq)
            push!(kmer_class[canonical(kmer)], allele_idx)
          end
      end
  end
  return kmer_class, length_alleles
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

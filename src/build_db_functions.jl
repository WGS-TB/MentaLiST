using Distributed
import GZip
import JLD: load, save
import FileIO: File, @format_str
import Blosc
import JSON
using Pkg
include("db_graph.jl")


# Complement a set; useful for 'compressing' large allele sets when building the scheme DB.
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

# Combine all loci kmer classifications in a format that is easy to serialize to disk.
function combine_loci_classification(::Type{DNAKmer{k}}, results, loci) where {k}
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

  # allele coverages: how many kmers cover each allele; same format as alleles_list;
  allele_coverage_per_locus = Int16[]

  for (locus,(kmer_class, allele_ids, kmer_weights, allele_coverages)) in zip(loci,results)
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
      dnakmer = DNAKmer{k}(kmer)
      push!(kmer_list, "$dnakmer")
      push!(weight_list, weight)
      push!(alleles_list, length(allele_list))
      append!(alleles_list, allele_list)
      n_kmers_in_class += 1
    end
    push!(loci_list, n_kmers_in_class)
    push!(loci_list, l2int[locus])
    push!(allele_ids_per_locus, n_alleles)
    append!(allele_ids_per_locus, allele_ids)
    append!(allele_coverage_per_locus, allele_coverages)
    # println("MIN COV:", m_val, "MAX VAR:", maximum(allele_coverages), "VAR>0: ", [x for x in allele_coverages if x > 0 ])
    # println("MIN COV:", m_val, " MAX VAR:", maximum(allele_coverages))
  end
  return (loci_list, weight_list, alleles_list, kmer_list, allele_ids_per_locus, allele_coverage_per_locus)

end

function kmer_class_for_each_locus(::Type{DNAKmer{k}}, files::Vector{String}, coverage) where {k}
  loci = [splitext(basename(file))[1] for file in files]
  if length(workers()) == 1 
    results = [build_db_graph(DNAKmer{k}, file, coverage) for file in files]
  else
    results = pmap(file->build_db_graph(DNAKmer{k}, file, coverage), files)
  end

  return results, loci
end

function save_db(::Type{DNAKmer{k}}, kmer_db, loci, filename, profile, args, version) where {k}

  loci_list, weight_list, alleles_list, kmer_list, allele_ids_per_locus, al_coverages = kmer_db
  d = Dict(
    "mentalist_version" => version,
    "loci_list"=> Blosc.compress(loci_list),
    "weight_list" => Blosc.compress(weight_list),
    "allele_ids_per_locus" => Blosc.compress(allele_ids_per_locus),
    "allele_coverages" => Blosc.compress(al_coverages),
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
  save(File(format"JLD", "$filename"), d)
  # Profile:
  if profile != nothing
    cp(profile, "$filename.profile", force=true)
  end
end

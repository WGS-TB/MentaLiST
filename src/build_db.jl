using ArgParse
using Bio.Seq
using DataStructures
using Logging

import GZip
function complement(orig_set, m)
  array = sort(collect(orig_set))
  c_set = Set{Int16}()
  expected::Int16 = 1
  i::Int16 = 1
  l::Int16 = length(array)
  while i <= l
    diff::Int16 = array[i] - expected
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

function save_db{k}(::Type{DNAKmer{k}}, kmer_db, loci, filename)
  GZip.open("$filename.gz", "w") do f
  # open(filename, "w") do f
    write(f, "$k\n")
    loci = join(loci,",")
    write(f, "$loci\n")
    for (kmer, v) in kmer_db
      for (locus, val, alleles) in v
        alleles = join(alleles,",")
        write(f, "$kmer\t$locus\t$val\t$alleles\n")
      end
    end
  end
end

function kmer_class_for_locus{k}(::Type{DNAKmer{k}}, fastafile::String)
  record = FASTASeqRecord{BioSequence{DNAAlphabet{2}}}()
  result = DefaultDict{DNAKmer{k}, Set{Int16}}(() -> Set{Int16}())
  length_alleles::Int16 = 0
  open(FASTAReader{BioSequence{DNAAlphabet{2}}}, fastafile) do reader
    allele_idx::Int16 = 1
      while !eof(reader)
          read!(reader, record)
          allele_idx = parse(Int16,split(record.name, "_")[2])
          length_alleles += 1
          for (pos, kmer) in each(DNAKmer{k}, record.seq)
            push!(result[canonical(kmer)], allele_idx)
          end
      end
  end
  # TODO: filter uninformative alleles: (present in all alleles)
  return result, length_alleles
end

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "-o"
          help = "Output file (kmer database)"
          arg_type = String
          required = true
        "-k"
          help = "Kmer size"
          required = true
          arg_type = Int8
        "files"
            nargs = '*'
            help = "Fasta files with the MLST scheme"
            required = true

    end
    return parse_args(s)
end


function main()
  # Logging.configure(level=DEBUG)
  Logging.configure(level=INFO)

  args = parse_commandline()
  k::Int8 = args["k"]
  results = []
  loci = Array{String,1}()
  info("Opening FASTA files ... ")
  for file in args["files"]
    locus::String = splitext(basename(file))[1]
    debug(" Parsing $locus ...")
    push!(loci, locus)
    kmer_class, n_alleles = kmer_class_for_locus(DNAKmer{k}, file)
    push!(results, (locus,kmer_class, n_alleles))
  end
  # Combine results:
  info("Combining results for each locus ...")
  kmer_classification = DefaultDict{DNAKmer{k}, Vector{Tuple{String, Int8, Set{Int16}}}}(() -> Vector{Tuple{String, Int8, Set{Int16}}}())
  for (locus, kmer_class, n_alleles) in results
    half = n_alleles/2
    for (kmer, allele_set) in kmer_class
      l_as = length(allele_set)
      if l_as == n_alleles
        continue
      elseif l_as > half
        allele_set = complement(allele_set, n_alleles)
        if length(allele_set) == 0
          println("Alelle set empty!!!!")
          exit()
        end
        push!(kmer_classification[kmer], (locus, -1, allele_set))
      else
        push!(kmer_classification[kmer], (locus, 1, allele_set))
      end
    end
  end
  info("Saving DB ...")
  save_db(DNAKmer{k}, kmer_classification, loci, args["o"])


  info("Done!")
  # @show kmer_alleles
end

main()

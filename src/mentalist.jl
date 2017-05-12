using OpenGene
using Logging
using ArgParse
using Bio.Seq
using DataStructures
import GZip

function open_db{k}(::Type{DNAKmer{k}}, filename)

  f = GZip.open("$filename")
  k_from_file = parse(Int8, readline(f))
  if k_from_file != k
    Logging.critical("Wrong kmer size, DB has $k !")
  end
  kmer_db = DefaultDict{DNAKmer{k}, Vector{Tuple{String, Int8, Vector{Int16}}}}(() -> Vector{Tuple{String, Int8, Vector{Int16}}}())
  loci = split(strip(readline(f)), ",")
  for ln in eachline(f)
    kmer, locus, val, alleles = split(ln, "\t")
    push!(kmer_db[DNAKmer{k}(kmer)], (locus, parse(Int8, val), map(x->parse(Int16,x), split(alleles,","))))
  end
  close(f)
  return kmer_db, loci
end

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "-k"
          help = "Kmer size"
          required = true
          arg_type = Int8
        "-o"
          help = "Output file with MLST call"
          arg_type = String
          required = true
        "--db"
          help = "Kmer database"
          required = true
          arg_type = String
        "-f"
          help = "FastQ input file"
          required = true
          arg_type = String
    end
    return parse_args(s)
end

function get_votes_for_sequence{k}(::Type{DNAKmer{k}}, seq, kmer_db, threshold=20)
  # stringMLST like pre_filter?
  # half = div(length(seq),2)
  # try
  #   testmer = DNAKmer{k}(seq[half:half+k-1])
  #   if !haskey(kmer_db, testmer)
  #     return false, false
  #   end
  # catch
  #   return false, false
  # end
  # count all votes:
  votes = Dict()
  # count locus hits, for filtering
  locus_hits = DefaultDict{String,Int16}(0)
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

function main()
  Logging.configure(level=INFO)
  args = parse_commandline()
  k::Int8 = args["k"]
  info("Opening kmer database ... ")
  kmer_db, loci = open_db(DNAKmer{k}, args["db"])

  votes = Dict(locus => DefaultDict{Int16, Int16}(0) for locus in loci)
  info("Opening fastq file(s) ... ")
  istream = fastq_open(args["f"])
  while (fq = fastq_read(istream))!=false
    good, locus_to_allele_votes = get_votes_for_sequence(DNAKmer{k}, fq.sequence.seq, kmer_db)
    if good
      locus, allele_votes = locus_to_allele_votes
      for (allele, val) in allele_votes
        votes[locus][allele] += val
      end
    end
  end
  println(join(loci,"\t"))
  println(join([sort(collect(votes[locus]), by=x->-x[2])[1][1] for locus in loci],"\t"))
  # for locus in loci
    # sorted_vote = sort(collect(votes[locus]), by=x->-x[2])[1:2]
    # println("$locus:$sorted_vote")
  # end
  info("Done.")
end

main()

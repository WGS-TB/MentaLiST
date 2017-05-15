include("build_db_functions.jl")
using OpenGene
using Logging
using ArgParse
using Bio.Seq
using DataStructures


function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "-o"
          help = "Output file with MLST call"
          arg_type = String
          required = true
        "-s"
          help = "Sample name"
          arg_type = String
          required = true
        "--db"
          help = "Kmer database"
          required = true
          arg_type = String
        "files"
          nargs = '*'
          help = "FastQ input files"
          required = true
          arg_type = String
    end
    return parse_args(s)
end

function get_votes_for_sequence{k}(::Type{DNAKmer{k}}, seq, kmer_db, threshold=10, prefilter=false)
  if prefilter
    # stringMLST like pre_filter
    half = div(length(seq),2)
    try
      testmer = DNAKmer{k}(seq[half:half+k-1])
      if !haskey(kmer_db, testmer)
        return false, false
      end
    catch
      return false, false
    end
  end
  # count all votes:
  votes = Dict()
  # count locus hits, for filtering
  # locus_hits = DefaultDict{String,Int16}(0)
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
  # @show locus_hits
  # @show best_locus
  # @show sort(collect(values(votes[best_locus])))
  if n_hits > threshold
    return true, (best_locus, votes[best_locus])
  else
    return false, false
  end
end

function main()
  Logging.configure(level=INFO)
  args = parse_commandline()
  info("Opening kmer database ... ")
  kmer_db, loci, n_alleles_list, k = open_db(args["db"])
  # 0 votes for all alleles everyone at the start:
  votes = Dict(idx => Dict{Int16, Int16}(i => 0 for i=1:n_alleles) for (idx,n_alleles) in enumerate(n_alleles_list))
  info("Opening fastq file(s) ... ")
  for f in args["files"]
    istream = fastq_open(f)
    while (fq = fastq_read(istream))!=false
      good, locus_to_allele_votes = get_votes_for_sequence(DNAKmer{k}, fq.sequence.seq, kmer_db)
      if good
        locus, allele_votes = locus_to_allele_votes
        for (allele, val) in allele_votes
          votes[locus][allele] += val
        end
      end
    end
  end
  info("Writing output ...")
  best_votes_alleles = [sort(collect(votes[idx]), by=x->-x[2])[1][1] for (idx,locus) in enumerate(loci)]
  open(args["o"], "w") do f
    header = join(vcat(["Sample"], loci, ["ClonalComplex"]), "\t")
    write(f,  "$header\n")
    # todo: look in the database for the type, to assing a 'ClonalComplex'
    calls = join(vcat([args["s"]], best_votes_alleles, ["0"]), "\t")
    write(f, "$calls\n")
  end
  # debug votes:
  # println(join(loci,"\t"))
  # println(join([sort(collect(votes[idx]), by=x->-x[2])[1][1] for (idx,locus) in enumerate(loci)],"\t"))

  for (idx,locus) in enumerate(loci)
    sorted_vote = sort(collect(votes[idx]), by=x->-x[2])[1:10]
    println("$locus:$sorted_vote")
  end
  info("Done.")
end

main()

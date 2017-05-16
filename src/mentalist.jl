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
        "-t"
          help = "A read is considered if it has at least this number of hits to the same locus. "
          arg_type = Int
          default = 10
        "-q"
          help = "Quick filter; if the first, middle and last kmers of a read are not in the kmer DB, the read is discarded. Disabled by default."
          action = :store_true

        "files"
          nargs = '*'
          help = "FastQ input files"
          required = true
          arg_type = String
    end
    return parse_args(s)
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
      good, locus_to_allele_votes = get_votes_for_sequence(DNAKmer{k}, fq.sequence.seq, kmer_db, args["t"], args["q"])
      if good
        locus, allele_votes = locus_to_allele_votes
        for (allele, val) in allele_votes
          votes[locus][allele] += val
        end
      end
    end
  end
  info("Writing output ...")
  write_calls(votes, loci, args["s"], args["o"])
  info("Done.")
end

main()

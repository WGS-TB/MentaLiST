module MentaLiST

include("build_db_functions.jl")
using OpenGene
using Logging
using ArgParse
using Bio.Seq
using DataStructures


function kmerize_kmc(files, k, threads=1)
  filepath = ""
  outpath, fh = mktemp()
  if length(files) > 1
  # create a tmp file with all files:
    filepath, f = mktemp()
    for file in files
      write(f, "$file\n")
    end
    filepath = "@$filepath"
    close(f)
  else
    filepath = files[1]
  end
  # now run:
  run(`/projects/pathogist/bin/kmc -k$k -t$threads -ci0 $filepath $outpath /tmp`)
  run(`/projects/pathogist/bin/kmc_tools transform $outpath dump $outpath`)
  return outpath
end

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
          help = "A read of length L is discarded if it has at less than (L - k) * t hits to the same locus in the kmer database, where k is the kmer length. 0 <= t <= 1"
          arg_type = Float64
          default = 0.3
          range_tester = x -> 0 <= x <= 1
        "-q"
          help = "Quick filter; if the first, middle and last kmers of a read are not in the kmer DB, the read is discarded. Disabled by default."
          action = :store_true
        "-e"
          help = "Use external kmc kmer counter. Disabled by default."
          action = :store_true
        "-j"
          help = "Skip length between consecutive k-mers. Defaults to 1."
          arg_type = Int
          default = 1
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
  kmer_db, loci, loci2alleles, k = open_db(args["db"])
  # 0 votes for all alleles everyone at the start:
  votes = Dict(locus_idx => Dict{Int16, Int}(i => 0 for i in 1:length(alleles)) for (locus_idx,alleles) in loci2alleles)
  if args["e"] # external kmer counter:
    info("Running kmc ... ")
    kmer_count_file = kmerize_kmc(args["files"], k)
    info("Counting votes from kmc kmers ... ")
    open(kmer_count_file) do f
      for ln in eachline(f)
        kmer, count = split(chomp(ln))
        kmer = DNAKmer{k}(kmer)
        if haskey(kmer_db, kmer)
          for (locus, val, alleles) in kmer_db[kmer]
            v = val * parse(Int,count)
            for allele in alleles
              if !haskey(votes, locus)
                votes[locus] = DefaultDict{Int16, Int16}(0)
              end
              votes[locus][allele] += v
            end
          end
        end
      end
    end
  else
    info("Opening fastq file(s) ... ")
    for f in args["files"]
      istream = fastq_open(f)
      while (fq = fastq_read(istream))!=false
        good, locus_to_allele_votes = get_votes_for_sequence(DNAKmer{k}, fq.sequence.seq, kmer_db, args["t"], args["q"], args["j"])
        if good
          locus, allele_votes = locus_to_allele_votes
          # @printf "%d to %s " minimum(collect(values(allele_votes))) maximum(collect(values(allele_votes)))
          # if haskey(allele_votes, 3)
          #   @show allele_votes[3]
          # end
          for (allele, val) in allele_votes
            votes[locus][allele] += val
          end
        end
      end
    end
  end
  info("Writing output ...")
  write_calls(votes, loci, loci2alleles, args["s"], args["o"])
  info("Done.")
end

main()

end

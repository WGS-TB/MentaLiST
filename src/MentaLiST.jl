#!/usr/bin/env julia
module MentaLiST

using Suppressor
@suppress_err begin
using Lumberjack
using ArgParse
end

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
  try
    run(`kmc -k$k -t$threads -ci0 $filepath $outpath /tmp`)
  catch e
    println("caught error $e")
    exit(1)
  end
  try
    run(`kmc_tools transform $outpath dump $outpath`)
  catch e
    println("caught error $e")
    exit(1)
  end
  return outpath
end

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
      "call"
        help = "MLST caller, given a sample and a k-mer database."
        action = :command
      "build_db"
        help = "Build a MLST k-mer database, given a list of FASTA files."
        action = :command
      "list_pubmlst"
        help = "List all available MLST schema from www.pubmlst.org. "
        action = :command
      "download_pubmlst"
        help = "Dowload a MLST scheme from pubmlst and build a MLST k-mer database."
        action = :command
      "list_cgmlst"
        help = "List all available cgMLST schema from www.cgmlst.org."
        action = :command
      "download_cgmlst"
          help = "Dowload a MLST scheme from cgmlst.org and build a MLST k-mer database."
          action = :command

    end
    # Calling MLST options:
    @add_arg_table s["call"] begin
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
          default = 0.2
          range_tester = x -> 0 <= x <= 1
        "-q"
          help = "Quick filter; if middle kmer of a read are not in the kmer DB, the read is discarded. Disabled by default."
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
    # Build DB from FASTA, options:
    @add_arg_table s["build_db"] begin
        "--db"
          help = "Output file (kmer database)"
          arg_type = String
          required = true
        "-k"
          help = "Kmer size"
          required = true
          arg_type = Int8
        "-f", "--fasta_files"
            nargs = '+'
            arg_type = String
            help = "Fasta files with the MLST scheme"
            required = true
        "-p", "--profile"
            arg_type = String
            help = "Profile file for known genotypes."
    end
    @add_arg_table s["list_pubmlst"] begin
      "-p", "--prefix"
      help = "Only list schema that starts with this prefix."
      arg_type = String
    end

    @add_arg_table s["list_cgmlst"] begin
      "-p", "--prefix"
      help = "Only list schema that starts with this prefix."
      arg_type = String
    end

    @add_arg_table s["download_pubmlst"] begin
      "-o", "--output"
        help = "Output folder for the schema files."
        arg_type = String
        required = true
      "-s", "--scheme"
        help = "Species name or ID of the scheme."
        arg_type = String
        required = true
      "-k"
        help = "K-mer size"
        required = true
        arg_type = Int8
      "--db"
        help = "Output file for the kmer database."
        arg_type = String
        required = true
    end

    @add_arg_table s["download_cgmlst"] begin
      "-o", "--output"
        help = "Output folder for the schema files."
        arg_type = String
        required = true
      "-s", "--scheme"
        help = "Species name or ID of the scheme"
        arg_type = String
        required = true
      "-k"
        help = "K-mer size"
        required = true
        arg_type = Int8
      "--db"
        help = "Output file for the kmer database."
        arg_type = String
        required = true
    end

    return parse_args(s)
end


#### Main COMMAND functions:
function call_mlst(args)
  include("build_db_functions.jl")
  info("Opening kmer database ... ")
  kmer_db, loci, loci2alleles, k, profile = open_db(args["db"])
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
          for (allele, val) in allele_votes
            votes[locus][allele] += val
          end
        end
      end
    end
  end
  info("Writing output ...")
  write_calls(votes, loci, loci2alleles, args["s"], args["o"], profile)
  info("Done.")
end

function list_pubmlst(args)
  include("mlst_download_functions.jl")
  list_pubmlst_schema(args["prefix"])
end

function download_pubmlst(args)
  include("mlst_download_functions.jl")
  loci_files, profile_file = download_pubmlst_scheme(args["scheme"], args["output"])
  info("Building the k-mer database ...")
  args["fasta_files"] = loci_files
  args["profile"] = profile_file
  build_db(args)
end

function list_cgmlst(args)
  include("mlst_download_functions.jl")
  list_cgmlst_schema(args["prefix"])
end

function download_cgmlst(args)
  include("mlst_download_functions.jl")
  loci_files = download_cgmlst_scheme(args["scheme"], args["output"])
  info("Building the k-mer database ...")
  args["fasta_files"] = loci_files
  args["profile"] = nothing
  build_db(args)
end


# function download_pubmlst(args)
#   include("mlst_download_functions.jl")
#   loci_files, profile_file = download_pubmlst_scheme(args["scheme"], args["output"])
#   info("Building the k-mer database ...")
#   args["fasta_files"] = loci_files
#   args["profile"] = profile_file
#   build_db(args)
# end

function build_db(args)
  include("build_db_functions.jl")
  k::Int8 = args["k"]
  info("Opening FASTA files ... ")
  results, loci = kmer_class_for_each_locus(k, args["fasta_files"])
  # Combine results:
  info("Combining results for each locus ...")
  kmer_classification = combine_loci_classification(k, results, loci)

  info("Saving DB ...")
  save_db(k, kmer_classification, loci, args["db"], args["profile"])
  info("Done!")
end

##### Main function: just calls the appropriate commands, with arguments:
function main()
  args = parse_commandline()
  # determine command:
  if args["%COMMAND%"] == "call"
    call_mlst(args["call"])
  elseif args["%COMMAND%"] == "build_db"
    build_db(args["build_db"])
  elseif args["%COMMAND%"] == "list_pubmlst"
    list_pubmlst(args["list_pubmlst"])
  elseif args["%COMMAND%"] == "download_pubmlst"
    download_pubmlst(args["download_pubmlst"])
  elseif args["%COMMAND%"] == "list_cgmlst"
    list_cgmlst(args["list_cgmlst"])
  elseif args["%COMMAND%"] == "download_cgmlst"
    download_cgmlst(args["download_cgmlst"])
  end
end

main()

end

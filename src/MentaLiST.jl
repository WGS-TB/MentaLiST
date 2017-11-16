#!/usr/bin/env julia

using Suppressor
@suppress_err begin
using Lumberjack
using ArgParse
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
        help = "List all available MLST schemes from www.pubmlst.org. "
        action = :command
      "download_pubmlst"
        help = "Dowload a MLST scheme from pubmlst and build a MLST k-mer database."
        action = :command
      "list_cgmlst"
        help = "List all available cgMLST schemes from www.cgmlst.org."
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
        "-t", "--mutation_threshold"
          help = "Maximum edit distance (number of mutations) when looking for novel alleles."
          arg_type = Int
          default = 6
        "--kt"
          help = "Minimum # of times a kmer is seen, to be considered 'solid', meaning actually present in the sample."
          arg_type = Int
          default = 10
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
        "-c", "--disable_compression"
          help = "Disables the default compression of the database, that stores only the most informative kmers. Not recommended unless for debugging."
          action = :store_true
    end
    @add_arg_table s["list_pubmlst"] begin
      "-p", "--prefix"
      help = "Only list schemes that starts with this prefix."
      arg_type = String
    end

    @add_arg_table s["list_cgmlst"] begin
      "-p", "--prefix"
      help = "Only list schemes that start with this prefix."
      arg_type = String
    end

    @add_arg_table s["download_pubmlst"] begin
      "-o", "--output"
        help = "Output folder for the scheme files."
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
      "-c", "--disable_compression"
        help = "Disables the default compression of the database, that stores only the most informative kmers. Not recommended unless for debugging."
        action = :store_true
    end

    @add_arg_table s["download_cgmlst"] begin
      "-o", "--output"
        help = "Output folder for the scheme files."
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
      "-c", "--disable_compression"
        help = "Disables the default compression of the database, that stores only the most informative kmers. Not recommended unless for debugging."
        action = :store_true
    end

    return parse_args(s)
end

#### Main COMMAND functions:
function call_mlst(args)
  include("build_db_functions.jl")
  info("Opening kmer database ... ")
  kmer_db, loci, loci2alleles, k, profile, build_args = open_db(args["db"])

  info("Opening fastq file(s) and counting kmers ... ")
  kmer_count = count_kmers(DNAKmer{k}, args["files"], kmer_db)
  info("Voting for alleles ... ")
  votes, loci_votes = count_votes(kmer_count, kmer_db, loci2alleles)
  info("Calling alleles ...")
  allele_calls, novel_alleles, best_voted_alleles, report, vote_log, ties, alleles_to_check = call_alleles_from_votes(DNAKmer{k}, kmer_count, votes, loci_votes, loci, loci2alleles, build_args["fasta_files"], args["kt"], args["mutation_threshold"])
  info("Writing output ...")
  write_calls(DNAKmer{k}, allele_calls, novel_alleles, best_voted_alleles, report, vote_log, ties, alleles_to_check, loci, args["s"], args["o"], profile)

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

function build_db(args)
  include("build_db_functions.jl")
  k::Int8 = args["k"]
  db_file = args["db"]
  profile = args["profile"]
  info("Opening FASTA files ... ")
  results, loci = kmer_class_for_each_locus(k, args["fasta_files"], !args["disable_compression"])
  # Combine results:
  info("Combining results for each locus ...")
  kmer_classification = combine_loci_classification(k, results, loci)

  info("Saving DB ...")
  save_db(k, kmer_classification, loci, db_file, profile, args)
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

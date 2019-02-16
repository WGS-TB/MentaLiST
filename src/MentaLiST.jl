# Alternative shebang:
#!/bin/bash
# -*- mode: julia -*-
#=
exec julia --color=yes --startup-file=no "${BASH_SOURCE[0]}" "$@"
=#
using ArgParse
using Distributed
using Printf

VERSION = "1.0.0"

function parse_commandline()
    s = ArgParseSettings()
    s.epilog = "MentaLiST -- The MLST pipeline developed by the PathOGiST research group. https://github.com/WGS-TB/MentaLiST\n" *
      "To cite: Feijao P, Yao H, Fornika D, Gardy J, Hsiao W, Chauve C, Chindelevitch L. 10/01/2018. Microbial Genomics 4(2): doi:10.1099/mgen.0.000146\n"
    s.version = "MentaLiST $VERSION."
    s.allow_ambiguous_opts= true # to allow -1 and -2.
    s.preformatted_epilog = true
    @add_arg_table s begin
      "call"
        help = "MLST caller, given a sample and a k-mer database."
        action = :command
      "build_db"
        help = "Build a MLST k-mer database, given a list of FASTA files."
        action = :command
      "db_info"
        help = "Extract information from an existing MentaLiST k-mer database"
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
      "download_enterobase"
          help = "Dowload a MLST scheme from Enterobase (enterobase.warwick.ac.uk) and build a MLST k-mer database."
          action = :command
      "-v", "--version"
          action = :show_version
          help = "show version information and exit"
      end

    # Calling MLST options:
    s["call"].preformatted_epilog = true
    s["call"].epilog = "MentaLiST MLST calling function. Calls alleles on a given MLST database.\n" *
      "You can create a custom DB with 'create_db' or other MentaLiST functions that download schemes from pubmlst, cgmlst.org or Enterobase.\n\n" *
      "Examples:\n" *
      "mentalist call -o my_sample.mlst --db my_scheme.db -1 sample_1.fastq.gz -2 sample_2.fastq.gz # one paired-end sample.\n" *
      "mentalist call -o all_samples.mlst --db my_scheme.db -1 *.fastq.gz -2 *.fastq.gz # multiple paired-end samples."

    @add_arg_table s["call"] begin
        "-o"
          help = "Output file with MLST call"
          arg_type = String
          required = true
        "--db"
          help = "Kmer database"
          required = true
          arg_type = String
        "-t", "--mutation_threshold"
          help = "Maximum number of mutations when looking for novel alleles."
          arg_type = Int
          default = 6
        "--kt"
          help = "Minimum # of times a kmer is seen to be considered present in the sample (solid)."
          arg_type = Int
          default = 10
        "--output_votes"
          help = "Outputs the results for the original voting algorithm."
          action = :store_true
        "--output_special"
          help = "Outputs a FASTA file with the alleles from 'special cases' such as incomplete coverage, novel, and multiple alleles. "
          action = :store_true
        "-i", "--sample_input_file"
          help = "Input TXT file for multiple samples. First column has the sample name, second the FASTQ file. Repeat the sample name for samples with more than one file (paired reads, f.i.)"
        "-1"
          nargs = '*'
          help = "FastQ input files, one per sample, forward reads (or unpaired reads)."
          arg_type = String
        "-2"
          nargs = '*'
          help = "FastQ input files, one per sample, reverse reads."
          arg_type = String
        "--fasta"
          help = "Input files are in FASTA format, instead of the default FASTQs."
          action = :store_true
    end

    ## Common option for all db building commands:
    s_db = ArgParseSettings()
    @add_arg_table s_db begin
      "--db"
        help = "Output file (kmer database)"
        arg_type = String
        required = true
      "-k"
        help = "Kmer size"
        required = true
        arg_type = Int8
      "--threads" # disabled for julia 1.1, the current syntax does not seem to work.
        arg_type = Int
        default = 1
        help = "Number of threads used in parallel."
      "-c", "--allele_coverage"
       arg_type = Float64
       default = 1.0
       help = "Minimum percentage of allele coverage in number of kmers of each allele."
    end
    # Build DB from FASTA, options:
    import_settings(s["build_db"], s_db)
    @add_arg_table s["build_db"] begin
        "-f", "--fasta_files"
            nargs = '+'
            arg_type = String
            help = "Fasta files with the MLST scheme"
            required = true
        "-p", "--profile"
            arg_type = String
            help = "Profile file for known genotypes."
    end

    @add_arg_table s["db_info"] begin
      "--db"
        help = "MentaLiST kmer database"
        arg_type = String
        required = true
    end

    # Listing functions, common options:
    s_list = ArgParseSettings()
    @add_arg_table s_list begin
      "-p", "--prefix"
      help = "Only list schemes where the species name starts with this prefix."
      arg_type = String
    end

    # List pubmlst
    import_settings(s["list_pubmlst"], s_list)

    # List cgmlst
    import_settings(s["list_cgmlst"], s_list)

    # Common options for MLST download functions:
    s_db_download = ArgParseSettings()
    import_settings(s_db_download, s_db) # import build_db common options
    @add_arg_table s_db_download begin
      "-o", "--output"
        help = "Output folder for the scheme Fasta files."
        arg_type = String
        required = true
      "-s", "--scheme"
        help = "Species name or scheme ID."
        arg_type = String
        required = true
    end

    # Download pubmlst:
    import_settings(s["download_pubmlst"], s_db_download) # import common options

    # Download cgmlst:
    import_settings(s["download_cgmlst"], s_db_download)  # import common options

    # Download enterobase:
    import_settings(s["download_enterobase"], s_db_download) # import common options
    # add (and override) specific enterobase options:
    s["download_enterobase"].error_on_conflict = false  # do not error-out when trying to override an option
    @add_arg_table s["download_enterobase"] begin
      "-s", "--scheme"
        help = "Letter identifying which scheme: (S)almonella, (Y)ersinia, or (E)scherichia/Shigella."
        arg_type = String
        required = true
      "-t", "--type"
        help = "Choose the type: 'cg' or 'wg' for cgMLST or wgMLST scheme, respectively."
        arg_type = String
        required = true
    end

    return parse_args(s)
end

# Check if files exist:
function check_files(files)
  dont_exist = [file for file in files if !isfile(file)]
  if length(dont_exist) > 0
    exit_error("The following input file(s) could not be found: $(join(dont_exist,',')).")
  end
end

### Error FUNCTIONS
@everywhere function exit_error(msg)
  @error(msg)
  println("exiting ...")
  exit(-1)
end

#### Main COMMAND functions:
function call_mlst(args)
  run_calling_pipeline(args)
end

function list_pubmlst(args)
  list_pubmlst_schema(args["prefix"])
end

function download_pubmlst(args)
  loci_files, profile_file = download_pubmlst_scheme(args["scheme"], args["output"])
  @info "Building the k-mer database ..."
  args["fasta_files"] = loci_files
  args["profile"] = profile_file
  build_db(args)
end

function list_cgmlst(args)
  list_cgmlst_schema(args["prefix"])
end

function download_cgmlst(args)
  loci_files = download_cgmlst_scheme(args["scheme"], args["output"])
  @info "Building the k-mer database ..."
  args["fasta_files"] = loci_files
  args["profile"] = nothing
  build_db(args)
end

function download_enterobase(args)
  loci_files = download_enterobase_scheme(args["scheme"], args["type"], args["output"])
  @info "Building the k-mer database ..."
  args["fasta_files"] = loci_files
  args["profile"] = nothing
  build_db(args)
end

function build_db(args, version=VERSION)
  # check if we need Gurobi:
  if args["allele_coverage"] < 1
    ok = include_gurobi()
    if !ok
      exit_error("Could not load libraries. To use the allele coverage functionality (option -c), a proper installation of Gurobi (www.gurobi.com) and the Gurobi and JuMP julia packages are required.\nCheck https://github.com/JuliaOpt/Gurobi.jl for installation instructions.")
    end
  end

  # check if files exist:
  check_files(args["fasta_files"])
  # get arguments and call the kmer db builder for each locus:
  k::Int8 = args["k"]
  db_file = args["db"]
  profile = args["profile"]
  @info "Opening FASTA files ... "
  results, loci = kmer_class_for_each_locus(DNAKmer{k}, args["fasta_files"], args["allele_coverage"])
  # Combine results:
  @info "Combining results for each locus ..."
  kmer_classification = combine_loci_classification(DNAKmer{k}, results, loci)

  for (index, fasta_file) in enumerate(args["fasta_files"])
    scheme_fasta_directory = pop!(split(dirname(fasta_file), "/"))
    if scheme_fasta_directory == ""
      scheme_fasta_directory = "."
    end
    args["fasta_files"][index] = scheme_fasta_directory * "/" * basename(fasta_file)
  end
  @info "Saving DB ..."
  save_db(DNAKmer{k}, kmer_classification, loci, db_file, profile, args, version)
  @info "Done!"
end

function db_info(args)
  filename = args["db"]
  @info "Opening kmer database ... "
  # Compressed database, open and decompress/decode in memory:
  d = load("$filename")
  @info "Finished the JLD load."
  build_args = JSON.parse(d["args"])
  k = build_args["k"]
  loci = d["loci"]
  loci_list = Blosc.decompress(Int32, d["loci_list"])
  num_loci = length(loci_list)
  mentalist_version = try
    d["mentalist_version"]
  catch
    "unknown"
  end
  println("mentalist_version\t$mentalist_version")
  println("k\t$k")
  println("num_loci\t$num_loci")
end

## Gurobi ILP:

# Kmer coverage is based on Gurobi:
function has_gurobi()
  packages = keys(Pkg.installed())
  return (in("Gurobi", packages) & in("JuMP", packages))
end

function include_gurobi()
  if has_gurobi()
    try
      include("kmer_coverage.jl") 
    catch e
      if isa(e, LoadError)
        @error("Error trying to load Gurobi package.")
        return false
      end
    end
    return true
  else
    return false
  end
end

##### Main function: just calls the appropriate commands, with arguments:

args = parse_commandline()
# determine command:
cmd = args["%COMMAND%"]
if cmd == "call"
  include("calling_functions.jl")
  call_mlst(args[cmd])

elseif cmd == "build_db"
  if args[cmd]["threads"] > 1
    addprocs(args[cmd]["threads"])
  end
  include("build_db_functions.jl")
  build_db(args[cmd])

elseif cmd == "db_info"
  import JLD: load, JSON, Blosc
  db_info(args[cmd])

elseif cmd == "list_pubmlst"
  include("mlst_download_functions.jl")
  list_pubmlst(args[cmd])

elseif cmd == "download_pubmlst"
  include("mlst_download_functions.jl")
  addprocs(args[cmd]["threads"])
  include("build_db_functions.jl")
  download_pubmlst(args[cmd])

elseif cmd == "list_cgmlst"
  include("mlst_download_functions.jl")
  list_cgmlst(args[cmd])

elseif cmd == "download_cgmlst"
  if args[cmd]["threads"] > 1
    addprocs(args[cmd]["threads"])
  end
  include("mlst_download_functions.jl")
  include("build_db_functions.jl")
  download_cgmlst(args[cmd])

elseif cmd == "download_enterobase"
  if args[cmd]["threads"] > 1
    addprocs(args[cmd]["threads"])
  end
  include("mlst_download_functions.jl")
  include("build_db_functions.jl")
  download_enterobase(args[cmd])

end

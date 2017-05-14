include("build_db_functions.jl")
using ArgParse
using Logging

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
            arg_type = String
            help = "Fasta files with the MLST scheme"
            required = true

    end
    return parse_args(s)
end


function main()
  Logging.configure(level=INFO)
  args = parse_commandline()
  k::Int8 = args["k"]

  info("Opening FASTA files ... ")
  results, loci = kmer_class_for_each_locus(k, args["files"])
  # Combine results:
  info("Combining results for each locus ...")
  kmer_classification = combine_loci_classification(k, results, loci)

  info("Saving DB ...")
  save_db(k, kmer_classification, loci, args["o"])
  info("Done!")
end

main()

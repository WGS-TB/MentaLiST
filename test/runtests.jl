using Test
# using Base.Test
# using Base.@__FILE__

include("../src/mlst_download_functions.jl")
include("../src/build_db_functions.jl")
include("../src/calling_functions.jl")
include("../src/db_graph.jl")

TEST_DIR = (dirname(@__FILE__))
TMPDIR = mktempdir()
VERSION = "testing"
# encapsulate in a big test to get a summary of all in the end; if a testset fails, it will show the specific testset results,
# otherwise will show just the total # of test passes.
@testset "mentalist_tests" begin

  c_jejuni_pubmlst_dir = string(TMPDIR, "/", "c_jejuni_pubmlst")
  c_jejuni_pubmlst_db_file = string(TMPDIR, "/", "c_jejuni_pubmlst_31.db")
  c_jejuni_pubmlst_copy_dir = string(TMPDIR, "/", "c_jejuni_pubmlst_copy")
  c_jejuni_fastq_file = string("../data/SRR5824107_small.fastq.gz")
  l_pneumophila_cgmlst_dir = string(TMPDIR, "/", "l_pneumophila_cgmlst")
  l_pneumophila_cgmlst_db_file = string(TMPDIR,"/", "l_pneumophila_cgmlst_31.db")
  l_pneumophila_fastq_file = string("../data/SRR801793.fastq.gz")

  c_jejuni_pubmlst_loci_files = ""
  @testset "download_pubmlst_scheme" begin
    mkdir(c_jejuni_pubmlst_dir)
    c_jejuni_pubmlst_loci_files, c_jejuni_profile_file = download_pubmlst_scheme("Campylobacter jejuni", c_jejuni_pubmlst_dir)
    c_jejuni_allele_filenames = open(string(TEST_DIR, "/", "test_data", "/", "c_jejuni_pubmlst_allele_filenames.txt"))
    for c_jejuni_allele_filename in eachline(c_jejuni_allele_filenames)
       @test isfile(string(c_jejuni_pubmlst_dir, "/", chomp(c_jejuni_allele_filename)))
    end
  end

  l_pneumophila_cgmlst_loci_files = ""
  @testset "download_cgmlst_scheme" begin
    mkdir(l_pneumophila_cgmlst_dir)
    l_pneumophila_cgmlst_loci_files = download_cgmlst_scheme("Legionella pneumophila", l_pneumophila_cgmlst_dir)
    l_pneumophila_allele_filenames = open(string(TEST_DIR, "/", "test_data", "/", "l_pneumophila_cgmlst_allele_filenames.txt"))
    for l_pneumophila_allele_filename in eachline(l_pneumophila_allele_filenames)
      @test isfile(string(l_pneumophila_cgmlst_dir, "/", chomp(l_pneumophila_allele_filename)))
    end
  end

  @testset "build_l_pneumophila_db" begin
    K::Int8 = 31
    results, loci = kmer_class_for_each_locus(DNAKmer{K}, l_pneumophila_cgmlst_loci_files, 1)
    @test typeof(results) == Vector{Tuple{Dict{UInt64,Vector{Int16}},Vector{Int16},Dict{UInt64,Int64},Vector{Int64}}}
    @test typeof(loci) == Vector{String}

    kmer_classification = combine_loci_classification(DNAKmer{K}, results, loci)
    loci_list, weight_list, alleles_list, kmer_list, allele_ids_per_locus = kmer_classification
    @test typeof(loci_list) == Vector{Int32}
    @test typeof(weight_list) == Vector{Int16}
    @test typeof(alleles_list) == Vector{Int16}
    @test typeof(kmer_list) == Vector{String}
    @test typeof(allele_ids_per_locus) == Vector{Int64}
    @test length(weight_list) == length(kmer_list)
    @test all([x != 0 for x in weight_list])
    @test all([length(kmer) == K for kmer in kmer_list])

    profile = nothing
    args = Dict("k" => K, "fasta_files" => l_pneumophila_cgmlst_loci_files)
    save_db(DNAKmer{K}, kmer_classification, loci, l_pneumophila_cgmlst_db_file, profile, args, VERSION)
    @test isfile(l_pneumophila_cgmlst_db_file)
  end

  @testset "build_c_jejuni_db" begin
    K::Int8 = 31
    results, loci = kmer_class_for_each_locus(DNAKmer{K}, c_jejuni_pubmlst_loci_files, 1)
    @test typeof(results) == Vector{Tuple{Dict{UInt64,Vector{Int16}},Vector{Int16},Dict{UInt64,Int64},Vector{Int64}}}
    @test typeof(loci) == Vector{String}

    kmer_classification = combine_loci_classification(DNAKmer{K}, results, loci)
    loci_list, weight_list, alleles_list, kmer_list, allele_ids_per_locus = kmer_classification
    @test typeof(loci_list) == Vector{Int32}
    @test typeof(weight_list) == Vector{Int16}
    @test typeof(alleles_list) == Vector{Int16}
    @test typeof(kmer_list) == Vector{String}
    @test typeof(allele_ids_per_locus) == Vector{Int64}
    @test length(weight_list) == length(kmer_list)
    @test all([x != 0 for x in weight_list])
    @test all([length(kmer) == K for kmer in kmer_list])

    profile = nothing
    args = Dict("k" => K, "fasta_files" => c_jejuni_pubmlst_loci_files)
    save_db(DNAKmer{K}, kmer_classification, loci, c_jejuni_pubmlst_db_file, profile, args, VERSION)
    @test isfile(c_jejuni_pubmlst_db_file)
  end

  @testset "call_l_pneumophila" begin
    # open the db:
    kmer_db, loci, loci2alleles, coverage, k, profile, build_args = open_db(l_pneumophila_cgmlst_db_file)
    @test typeof(kmer_db) == Dict{Kmer{DNA,31},Vector{Tuple{Int16,Int16,Vector{Int16}}}}
    @test typeof(loci) == Vector{String}
    @test typeof(loci2alleles) == Dict{Int16,Vector{Int16}}

    # count kmers:
    kmer_count = count_kmers(DNAKmer{k}, [l_pneumophila_fastq_file])
    # test kmer counts:
    @test kmer_count[DNAKmer{k}("AACTAAAACACTGGCTTTACTGCCCCAAATG")] == 19
    @test kmer_count[DNAKmer{k}("TTATGCCTATGAAAATATGATTTTATCTGAA")] == 24
    @test kmer_count[DNAKmer{k}("AAAGGAATTACAATGAATATCTAAATAATTC")] == 14
    @test kmer_count[DNAKmer{k}("AGGAAGTTTTGGTTATATGCTATATGATGCG")] == 21
    @test kmer_count[DNAKmer{k}("AAAAGAGCCAATTTGATGAGCCACCACATCA")] == 20

    # call:
    votes, loci_votes = count_votes(kmer_count, kmer_db, loci2alleles, coverage)
    @test typeof(votes) == Dict{Int16,Dict{Int16,Int64}}
    @test typeof(loci_votes) == DefaultDict{Int16,Int64,Int64}

    # some parameters:
    kmer_thr, max_mutations, output_votes = 1, 5, true
    # call:
    allele_calls, voting_result = call_alleles(DNAKmer{k}, kmer_count, votes, loci_votes, loci, loci2alleles, l_pneumophila_cgmlst_loci_files, kmer_thr, max_mutations, output_votes)
    @test typeof(allele_calls) == Vector{AlleleCall}
    @test allele_calls[1].allele == "1"
    @test allele_calls[1].depth == 6
    @test allele_calls[2].allele == "1"
    @test allele_calls[2].depth == 9
    @test allele_calls[132].allele == "N" # novel_allele
    @test allele_calls[662].allele == "N" # novel_allele
    @test allele_calls[536].allele == "1"
    @test allele_calls[536].flag == "+" # multiple alleles possible
    @test allele_calls[652].allele == "1"
    @test allele_calls[652].flag == "+" # multiple alleles possible
  end

  @testset "call_c_jejuni" begin
    # open the db:
    kmer_db, loci, loci2alleles, coverage, k, profile, build_args = open_db(c_jejuni_pubmlst_db_file)
    @test typeof(kmer_db) == Dict{Kmer{DNA,31},Vector{Tuple{Int16,Int16,Vector{Int16}}}}
    @test typeof(loci) == Vector{String}
    @test typeof(loci2alleles) == Dict{Int16,Vector{Int16}}

    # count kmers:
    kmer_count = count_kmers(DNAKmer{k}, [c_jejuni_fastq_file])

    # test kmer counts:
    @test kmer_count[DNAKmer{k}("TCATTTAAGGACTTTTCAGTGATTAAAATCA")] == 2
    @test kmer_count[DNAKmer{k}("CACTCCAATTTTTTCAAATAAAGTAGCTAAG")] == 0
    @test kmer_count[DNAKmer{k}("ATTCTTTTACTCCTATTATCGGTTATACTAA")] == 0
    @test kmer_count[DNAKmer{k}("GAAAAAAGTAATCCAAGGTGCGCAAAAAGCA")] == 4
    @test kmer_count[DNAKmer{k}("AAATATAGTCAATAAATTATAAAAAAAACTT")] == 0

    # call:
    votes, loci_votes = count_votes(kmer_count, kmer_db, loci2alleles, coverage)
    @test typeof(votes) == Dict{Int16,Dict{Int16,Int64}}
    @test typeof(loci_votes) == DefaultDict{Int16,Int64,Int64}

    # some parameters:
    kmer_thr, max_mutations, output_votes = 2, 5, true
    # call:
    allele_calls, voting_result = call_alleles(DNAKmer{k}, kmer_count, votes, loci_votes, loci, loci2alleles, c_jejuni_pubmlst_loci_files, kmer_thr, max_mutations, output_votes)
    @test typeof(allele_calls) == Vector{AlleleCall}
    @test [ac.allele for ac in allele_calls] == ["2", "17", "2", "3", "2", "1", "5"]
  end

  @testset "call_c_jejuni_after_moving_db" begin
    mkdir(c_jejuni_pubmlst_copy_dir)
    cp(c_jejuni_pubmlst_db_file, string(c_jejuni_pubmlst_copy_dir, "/", basename(c_jejuni_pubmlst_db_file)))
    cp(c_jejuni_pubmlst_dir, string(c_jejuni_pubmlst_copy_dir, "/", basename(c_jejuni_pubmlst_dir)))
    # open the db:
    kmer_db, loci, loci2alleles, coverage, k, profile, build_args = open_db(string(dirname(c_jejuni_pubmlst_copy_dir), "/", basename(c_jejuni_pubmlst_db_file)))
    @test typeof(kmer_db) == Dict{Kmer{DNA,31},Vector{Tuple{Int16,Int16,Vector{Int16}}}}
    @test typeof(loci) == Vector{String}
    @test typeof(loci2alleles) == Dict{Int16,Vector{Int16}}

    # count kmers:
    kmer_count = count_kmers(DNAKmer{k}, [c_jejuni_fastq_file])

    # test kmer counts:
    @test kmer_count[DNAKmer{k}("TCATTTAAGGACTTTTCAGTGATTAAAATCA")] == 2
    @test kmer_count[DNAKmer{k}("CACTCCAATTTTTTCAAATAAAGTAGCTAAG")] == 0
    @test kmer_count[DNAKmer{k}("ATTCTTTTACTCCTATTATCGGTTATACTAA")] == 0
    @test kmer_count[DNAKmer{k}("GAAAAAAGTAATCCAAGGTGCGCAAAAAGCA")] == 4
    @test kmer_count[DNAKmer{k}("AAATATAGTCAATAAATTATAAAAAAAACTT")] == 0

    # call:
    votes, loci_votes = count_votes(kmer_count, kmer_db, loci2alleles, coverage)
    @test typeof(votes) == Dict{Int16,Dict{Int16,Int64}}
    @test typeof(loci_votes) == DefaultDict{Int16,Int64,Int64}

    # some parameters:
    kmer_thr, max_mutations, output_votes = 2, 5, true
    # call:
    allele_calls, voting_result = call_alleles(DNAKmer{k}, kmer_count, votes, loci_votes, loci, loci2alleles, c_jejuni_pubmlst_loci_files, kmer_thr, max_mutations, output_votes)
    @test typeof(allele_calls) == Vector{AlleleCall}
    @test [ac.allele for ac in allele_calls] == ["2", "17", "2", "3", "2", "1", "5"]
  end

  @testset "complement_alleles" begin
    @test complement_alleles([], 0) == Int16[]
    @test complement_alleles([], 1) == [Int16(1)]
    @test complement_alleles([-1], 0) == Int16[]
    @test complement_alleles([-1], 1) == [Int16(1)]
    @test complement_alleles([0], 0) == Int16[]
    @test complement_alleles([0], 1) == [Int16(1)]
    @test complement_alleles([1], 0) == Int16[]
    @test complement_alleles([1], 1) == Int16[]
    @test complement_alleles([0,1,3,5,7,9], 10) == Int16[2,4,6,8,10]
    @test complement_alleles([0,2,4,6,8,10], 10) == Int16[1,3,5,7,9]
  end

  @testset "twin" begin
    @test twin(DNAKmer("A")) == DNAKmer("T")
    @test twin(DNAKmer("C")) == DNAKmer("G")
    @test twin(DNAKmer("G")) == DNAKmer("C")
    @test twin(DNAKmer("T")) == DNAKmer("A")
    @test twin(DNAKmer("ACTG")) == DNAKmer("CAGT")
  end

  # rm(TMPDIR, recursive=true)
end

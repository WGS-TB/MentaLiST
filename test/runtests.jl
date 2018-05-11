using Base.Test
using Base.@__FILE__
using Lumberjack

include("../src/mlst_download_functions.jl")
include("../src/build_db_functions.jl")
include("../src/db_graph.jl")

TEST_DIR = (dirname(@__FILE__))
TMPDIR = mktempdir()
VERSION = "testing"

@testset "mentalist_tests" begin
  c_jejuni_pubmlst_dir = string(TMPDIR, "/", "c_jejuni_pubmlst")
  c_jejuni_pubmlst_db_file = string(TMPDIR, "/", "c_jejuni_pubmlst_31.db")
  c_jejuni_fastq_file = string("../data/SRR5824107_small.fastq.gz")
  l_pneumophila_cgmlst_dir = string(TMPDIR, "/", "l_pneumophila_cgmlst")
  s_enterica_enterobase_cgmlst_dir = string(TMPDIR, "/", "s_enterica_enterobase_cgmlst")
  s_enterica_enterobase_cgmlst_db_file = string(TMPDIR,"/", "s_enterica_enterobase_cgmlst_31.db")

  c_jejuni_pubmlst_loci_files = ""
  @testset "download_pubmlst_scheme" begin
    mkdir(c_jejuni_pubmlst_dir)
    c_jejuni_pubmlst_loci_files, c_jejuni_profile_file = download_pubmlst_scheme("Campylobacter jejuni", c_jejuni_pubmlst_dir)
    c_jejuni_allele_filenames = open(string(TEST_DIR, "/", "test_data", "/", "c_jejuni_pubmlst_allele_filenames.txt"))
    for c_jejuni_allele_filename in eachline(c_jejuni_allele_filenames)
      @test isfile(string(c_jejuni_pubmlst_dir, "/", chomp(c_jejuni_allele_filename)))
    end
  end

  @testset "build_c_jejuni_db" begin
    K::Int8 = 31
    results, loci = kmer_class_for_each_locus(K, c_jejuni_pubmlst_loci_files, true)
    @test typeof(results) == Vector{Any}
    @test typeof(loci) == Vector{String}

    kmer_classification = combine_loci_classification(K, results, loci)
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
    save_db(K, kmer_classification, loci, c_jejuni_pubmlst_db_file, profile, args)
    @test isfile(c_jejuni_pubmlst_db_file)
  end
  
  @testset "call_c_jejuni" begin
    # open the db:
    kmer_classification, loci, loci2alleles, k, profile = open_db(c_jejuni_pubmlst_db_file)
    # @test typeof(kmer_classification) == DataStructures.DefaultDict{Bio.Seq.Kmer{Bio.Seq.DNANucleotide,31},Vector{Tuple{Int16,Int16,Vector{Int16}}}}
    @test typeof(loci) == Vector{String}
    @test typeof(loci2alleles) == Dict{Int16,Vector{Int16}}
  
    # count kmers:
    votes, loci_votes = count_kmers_and_vote(DNAKmer{k}, [c_jejuni_fastq_file], kmer_classification, loci2alleles)
    
    @test typeof(loci_votes) == DataStructures.DefaultDict{Int16,Int64,Int64}
    @test loci_votes[29] == 0
    @test loci_votes[306] == 0
    @test loci_votes[1090] == 0
    @test loci_votes[1316] == 0
    @test loci_votes[1333] == 0
  
    
    @test typeof(votes) == Dict{Int16,Dict{Int16,Int64}}  
    @test votes[7][288] == 482
    @test votes[7][306] == 437
    @test votes[7][520] == 516
    @test votes[4][306] == 1533
    @test votes[4][29] == 484
    @test votes[4][413] == -775
    
  end


  # @testset "download_cgmlst_scheme" begin
  #   mkdir(l_pneumophila_cgmlst_dir)
  #   download_cgmlst_scheme("Legionella pneumophila", l_pneumophila_cgmlst_dir)
  #   l_pneumophila_allele_filenames = open(string(TEST_DIR, "/", "test_data", "/", "l_pneumophila_cgmlst_allele_filenames.txt"))
  #   for l_pneumophila_allele_filename in eachline(l_pneumophila_allele_filenames)
  #     @test isfile(string(l_pneumophila_cgmlst_dir, "/", chomp(l_pneumophila_allele_filename)))
  #   end
  # end

  # s_enterica_enterobase_cgmlst_locus_files = ""
  # @testset "download_s_enterica_enterobase_cgmlst_scheme" begin
  #   mkdir(s_enterica_enterobase_cgmlst_dir)
  #   s_enterica_enterobase_cgmlst_locus_files = download_enterobase_scheme("S", "cg", s_enterica_enterobase_cgmlst_dir)
  #   s_enterica_enterobase_cgmlst_allele_filenames = open(string(TEST_DIR, "/", "test_data", "/", "s_enterica_enterobase_cgmlst_allele_filenames.txt"))
  #   for s_enterica_enterobase_cgmlst_allele_filename in eachline(s_enterica_enterobase_cgmlst_allele_filenames)
  #     @test isfile(string(s_enterica_enterobase_cgmlst_dir, "/", chomp(s_enterica_enterobase_cgmlst_allele_filename)))
  #   end
  # end

  # @testset "build_s_enterica_enterobase_cgmlst_db" begin
  #   K::Int8 = 31
  #   results, loci = kmer_class_for_each_locus(K, s_enterica_enterobase_cgmlst_locus_files, true)
  #   @test typeof(results) == Array{Any,1}
  #   @test typeof(loci) == Vector{String}
  #   
  #   kmer_classification = combine_loci_classification(K, results, loci)
  #   loci_list, weight_list, alleles_list, kmer_list, allele_ids_per_locus = kmer_classification
  #   @test typeof(loci_list) == Vector{Int32}
  #   @test typeof(weight_list) == Vector{Int16}
  #   @test typeof(alleles_list) == Vector{Int16}
  #   @test typeof(kmer_list) == Vector{String}
  #   @test typeof(allele_ids_per_locus) == Vector{Int64}
  #   @test length(weight_list) == length(kmer_list)
  #   @test all([x != 0 for x in weight_list])
  #   @test all([length(kmer) == K for kmer in kmer_list])
  #   
  #   profile = nothing
  #   args = Dict("k" => K, "fasta_files" => s_enterica_enterobase_cgmlst_locus_files)
  #   save_db(K, kmer_classification, loci, s_enterica_enterobase_cgmlst_db_file, profile, VERSION)
  #   @test isfile(s_enterica_enterobase_cgmlst_db_file)
  # end

  @testset "complement_alleles" begin
    @test complement_alleles([], 0) == Int16[]
    @test complement_alleles([], 1) == [Int16(1)]
    @test complement_alleles([-1], 0) == Int16[]
    @test complement_alleles([-1], 1) == [Int16(1)]
    @test complement_alleles([0], 0) == Int16[]
    @test complement_alleles([0], 1) == [Int16(1)]
    @test complement_alleles([1], 0) == Int16[]
    @test complement_alleles([1], 1) == Int16[]
    @test complement_alleles([], 0) == Int16[]
  end

  @testset "twin" begin
    @test twin(DNAKmer{0}, DNAKmer("")) == DNAKmer("")
    @test twin(DNAKmer{1}, DNAKmer("A")) == DNAKmer("T")
    @test twin(DNAKmer{1}, DNAKmer("C")) == DNAKmer("G")
    @test twin(DNAKmer{1}, DNAKmer("G")) == DNAKmer("C")
    @test twin(DNAKmer{1}, DNAKmer("T")) == DNAKmer("A")
  end

end

# rm(TMPDIR, recursive=true)

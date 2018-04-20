using Base.Test
using Base.@__FILE__

include("../src/mlst_download_functions.jl")
include("../src/build_db_functions.jl")
include("../src/db_graph.jl")

TEST_DIR = (dirname(@__FILE__))
TMPDIR = mktempdir()
VERSION = "testing"

c_jejuni_pubmlst_dir = string(TMPDIR, "/", "c_jejuni_pubmlst")
l_pneumophila_cgmlst_dir = string(TMPDIR, "/", "l_pneumophila_cgmlst")
s_enterica_enterobase_cgmlst_dir = string(TMPDIR, "/", "s_enterica_enterobase_cgmlst")
s_enterica_enterobase_cgmlst_db_file = string(TMPDIR,"/", "s_enterica_enterobase_cgmlst_31.db")

@testset "download_pubmlst_scheme" begin
  mkdir(c_jejuni_pubmlst_dir)
  download_pubmlst_scheme("Campylobacter jejuni", c_jejuni_pubmlst_dir)
  c_jejuni_allele_filenames = open(string(TEST_DIR, "/", "test_data", "/", "c_jejuni_pubmlst_allele_filenames.txt"))
  for c_jejuni_allele_filename in eachline(c_jejuni_allele_filenames)
     @test isfile(string(c_jejuni_pubmlst_dir, "/", chomp(c_jejuni_allele_filename)))
  end
end

@testset "download_cgmlst_scheme" begin
  mkdir(l_pneumophila_cgmlst_dir)
  download_cgmlst_scheme("Legionella pneumophila", l_pneumophila_cgmlst_dir)
  l_pneumophila_allele_filenames = open(string(TEST_DIR, "/", "test_data", "/", "l_pneumophila_cgmlst_allele_filenames.txt"))
  for l_pneumophila_allele_filename in eachline(l_pneumophila_allele_filenames)
    @test isfile(string(l_pneumophila_cgmlst_dir, "/", chomp(l_pneumophila_allele_filename)))
  end
end

s_enterica_enterobase_cgmlst_locus_files = ""
@testset "downloads_enterica_enterobase_cgmlst_scheme" begin
  mkdir(s_enterica_enterobase_cgmlst_dir)
  s_enterica_enterobase_cgmlst_locus_files = download_enterobase_scheme("S", "cg", s_enterica_enterobase_cgmlst_dir)
  s_enterica_enterobase_cgmlst_allele_filenames = open(string(TEST_DIR, "/", "test_data", "/", "s_enterica_enterobase_cgmlst_allele_filenames.txt"))
  for s_enterica_enterobase_cgmlst_allele_filename in eachline(s_enterica_enterobase_cgmlst_allele_filenames)
    @test isfile(string(s_enterica_enterobase_cgmlst_dir, "/", chomp(s_enterica_enterobase_cgmlst_allele_filename)))
  end
end

@testset "build_s_enterica_enterobase_cgmlst_db" begin
    K::Int8 = 31
    results, loci = kmer_class_for_each_locus(K, s_enterica_enterobase_cgmlst_locus_files, true)
    @test typeof(results) == Vector{Tuple{Dict{UInt64,Vector{Int16}},Vector{Int16},Dict{UInt64,Int64}}}
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
    args = Dict("k" => K, "fasta_files" => s_enterica_enterobase_cgmlst_locus_files)
    save_db(K, kmer_classification, loci, s_enterica_enterobase_cgmlst_db_file, profile, args, VERSION)
    @test isfile(l_pneumophila_cgmlst_db_file)
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
  @test complement_alleles([], 0) == Int16[]
end

@testset "twin" begin
  @test twin(DNAKmer{0}, DNAKmer("")) == DNAKmer("")
  @test twin(DNAKmer{1}, DNAKmer("A")) == DNAKmer("T")
  @test twin(DNAKmer{1}, DNAKmer("C")) == DNAKmer("G")
  @test twin(DNAKmer{1}, DNAKmer("G")) == DNAKmer("C")
  @test twin(DNAKmer{1}, DNAKmer("T")) == DNAKmer("A")
end

# rm(TMPDIR, recursive=true)

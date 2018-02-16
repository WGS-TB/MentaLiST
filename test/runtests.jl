using Base.Test
using Base.@__FILE__

include("../src/mlst_download_functions.jl")
include("../src/build_db_functions.jl")
include("../src/db_graph.jl")

TEST_DIR = (dirname(@__FILE__))
TMPDIR = mktempdir()

c_jejuni_pubmlst_dir = string(TMPDIR, "/", "c_jejuni_pubmlst")
l_pneumophila_cgmlst_dir = string(TMPDIR, "/", "l_pneumophila_cgmlst")

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

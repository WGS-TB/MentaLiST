using Base.Test
using MentaLiST

include("../src/build_db_functions.jl")

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



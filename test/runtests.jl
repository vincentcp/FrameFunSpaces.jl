
using FrameFunSpaces, Test
using FrameFunSpaces: DictionarySpace1d, DictionarySpace2d, DictionarySpaceNd


@testset "Simple functionality" begin
    F = Fourier(10)
    @test F isa DictionarySpace1d
    @test F[1] isa eltype(typeof(F))
    @test F[2](0)≈1
    @test evaluate.(points(F), F)≈[F[i](x) for x in points(F), i in 1:dimension(F) ]
    FN = Fourier((3,4))
    @test FN isa Dictionary
    @test domaindimension(FN)==2
    @test FN[1] isa eltype(typeof(FN))
    @test FN[1,2] isa eltype(typeof(FN))
    @test !(FN isa FrameFunSpaces.DictionarySpaceNd{3})
    @test 12==dimension(FN)
    @test (3,4)==dimensions(FN)
    @test FN isa FrameFunSpaces.DictionarySpace2d
    F = Fourier((3,4))
    x = points(F)
    A = [F[i](x...) for x in points(F), i in 1:dimension(F) ]
    @test A≈reshape(evaluate.(x, F), size(A))
end

@testset "OuterProduct" begin
    a1 = rand(9,4)
    a2 = rand(9,3)
    a = OuterProduct((a1,a2),(false,true))
    for i in 1:size(a,1)
        @test a[1,:,:]≈a1[1,:]*a2[1,:]'
    end
    dest = similar(a)
    @timed(copy!(dest,  a))[3] <= 200
    @test @timed(copy!(dest,  a))[3] <= 200
end

F = Fourier(10)
for f in (:domaintype, :rangetype, :dimension, :domaindimension, :domain, :points)
    @eval $f(F)
end
FN= Fourier((3,4))
for f in (:domaintype, :rangetype, :dimension, :domaindimension, :domain, :points)
    @eval $f(FN)
end

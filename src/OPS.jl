
for T in (:Fourier, :Taylor, :CosSpace, :SinSpace, :Laurent)
    @eval $T(n::Int, trailing...) =
        SubSpace($T(trailing...),1:n)
end

@reexport using ApproxFunOrthogonalPolynomials
import ApproxFunOrthogonalPolynomials:
    Chebyshev, NormalizedChebyshev, Hermite, NormalizedHermite,
    Jacobi, NormalizedJacobi, Legendre, NormalizedLegendre,
    Laguerre, NormalizedLaguerre, Ultraspherical, NormalizedUltraspherical, ContinuousSpace, ChebyshevDirichlet

for T in (:Chebyshev, :NormalizedChebyshev, :Hermite, :NormalizedHermite,
    :Jacobi, :NormalizedJacobi, :Legendre, :NormalizedLegendre,
    :Laguerre, :NormalizedLaguerre, :Ultraspherical, :NormalizedUltraspherical, :ContinuousSpace, :ChebyshevDirichlet)
    @eval $T(n::Int, trailing...) =
        SubSpace($T(trailing...),1:n)
end


for T in (:Chebyshev, :NormalizedChebyshev, :Hermite, :NormalizedHermite,
    :Jacobi, :NormalizedJacobi, :Legendre, :NormalizedLegendre,
    :Laguerre, :NormalizedLaguerre, :Ultraspherical, :NormalizedUltraspherical, :ContinuousSpace, :ChebyshevDirichlet)
    @eval $T(10)
end

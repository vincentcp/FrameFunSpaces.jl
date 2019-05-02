
@reexport using ApproxFunFourier
import ApproxFunFourier: Fourier, Taylor, CosSpace, SinSpace, Laurent

Fourier(n::Int, trailing...) =
    dictionary(SubSpace(Laurent(trailing...),1:n))

Fourier(n::NTuple{N,Int}) where N =
    dictionary(SubSpace(Laurent()^N,CartesianIndices(n)))

function Fourier(n::Int, a::Number, b::Number)
	T = float(promote_type(typeof(a),typeof(b)))
	Laurent(n, Interval(a, b))
end

instantiate(::Type{Fourier}, n::Int, ::Type{T}) where {T} = Fourier(n, zero(T), one(T))

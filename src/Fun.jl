struct DictContainer{DICT<:Dictionary}
    dict::DICT
end

abstract type DictEvaluation <: Base.Broadcast.BroadcastStyle
end

struct Univariate <: DictEvaluation
end

struct Multivariate <: DictEvaluation
end

abstract type DictionaryStyle <: Base.Broadcast.BroadcastStyle end

struct DictStyle <: DictionaryStyle
end

struct UVDictStyle <: DictionaryStyle
end

struct MVDictStyle <: DictionaryStyle
end

# Create something that supports axes, but do not give a range to the first dimension
Base.Broadcast.broadcastable(dict::Dictionary) = DictContainer(dict)
Base.axes(a::DictContainer{<:Dictionary}) = (FullSpace(rangetype(a.dict)),map(Base.OneTo, dimensions(a.dict))...)
# evaluation in a grid discretizes
Base.Broadcast._bcs1(a, b::FullSpace) = a
Base.Broadcast._bcs1(a::FullSpace, b) = b


Base.BroadcastStyle(::Type{<:DictContainer{<:Dictionary}}) = DictStyle()
Base.BroadcastStyle(::Type{<:DictContainer{<:SubSpaceDictionary{<:TensorSpace}}}) = MVDictStyle()

Base.Broadcast.result_style(::Base.Broadcast.DefaultArrayStyle, ::DictionaryStyle) =
    Univariate()
Base.Broadcast.result_style(::Base.Broadcast.DefaultArrayStyle, ::MVDictStyle) =
    Multivariate()
Base.Broadcast.result_style(::DictionaryStyle, ::Base.Broadcast.DefaultArrayStyle) =
    FunStyle()

Base.@propagate_inbounds Broadcast._getindex(args::Tuple{Any,<:DictContainer},I) =
    (Broadcast._broadcast_getindex(args[1], I), Broadcast._getindex(Base.tail(args), CartesianIndex(Base.tail(I.I)...))...)
Base.@propagate_inbounds Broadcast._broadcast_getindex(A::DictContainer, I) = A.dict[I.I...]
Base._return_type(::Core.Typeof(evaluate), ::Type{<:Tuple{FrameFunSpaces.SingletonDictionary{DS,IT,DD,RR},T}}) where {DS,IT,DD,RR,T} =
    promote_type(T,RR)
Broadcast.Base.broadcast(f, coefficients::AbstractArray, dict::Dictionary, pts::AbstractArray) =
    broadcast(*, Broadcast.Base.broadcasted(f, pts, dict), Broadcast.Base.broadcasted(coefficients))
Broadcast.broadcast(::Core.Typeof(evaluate), dict::Dictionary, x::AbstractArray) = Fun(dict, x)

Base.similar(bc::Base.Broadcast.Broadcasted{<:DictEvaluation}, ::Type{ElType}) where ElType =
    similar(Array{ElType}, axes(bc))

_univariate_elements(dc::DictContainer) =
    map(DictContainer, _univariate_elements(dc.dict))
_univariate_elements(a::SVector) = Tuple(a)
function _univariate_elements(a::AbstractArray{<:SVector}, i)
    A = Array{eltype(a[1])}(undef, size(a))
    for k in eachindex(a)
        A[k] = a[k][i]
    end
    A
end
_univariate_elements(t::Tuple{<:Array{<:SVector{N}},DictContainer{<:SubSpaceDictionary{<:TensorSpace}}}) where N =
    ntuple(k->(_univariate_elements(t[1], k),_univariate_elements(t[2], k)), Val(N))
_univariate_elements(bc::Broadcast.Broadcasted{<:Multivariate}) =
    map(x->Broadcast.Broadcasted{Univariate}(bc.f, x), _univariate_elements(bc.args))
Broadcast.materialize(bc::Broadcast.Broadcasted{<:Multivariate}) =
    OuterProduct(map(Broadcast.materialize, FrameFunSpaces._univariate_elements(bc)), ntuple(i-> i==1 ? false : true, Val(length(bc.args))))


struct OuterProduct{ARRAYS<:Tuple,BOOLS<:NTuple{M,Bool} where M, T,N,I} <: AbstractArray{T,N}
    arrays::ARRAYS
    outer::BOOLS
    size::NTuple{N,Int}
    index_scratch::I
    function OuterProduct(arrays::NTuple{N,AbstractArray}, bools::NTuple{N,Bool}=ntuple(true,Val(N))) where {N}
        _size(arrays, bools, i) = bools[i] ? map(x->size(x,i), arrays) : (size(arrays[1], i),)
        ElType = promote_type(map(eltype, arrays)...)
        n = 0
        for i in 1:N
            if bools[i]
                n += N
            else
                n1 = size(arrays[1], i)
                for a in arrays
                    @assert size(a, i) == n1
                end
                n +=1
            end
        end
        index_scratch = map( x-> Vector{Int}(undef, x ? N : 1), bools)
        new{typeof(arrays),typeof(bools),ElType,n,typeof(index_scratch)}(arrays, bools, Tuple(vcat(map(SVector, _size.(Ref(arrays),Ref(bools), 1:N))...)), index_scratch)
    end
end
Base.size(a::OuterProduct) = a.size
Base.IndexStyle() = Base.IndexCartesian()


Base.@propagate_inbounds function Base.getindex(a::OuterProduct{ARRAYS,<:NTuple{M,Bool},T,N} , I::Vararg{Int,N}) where {ARRAYS,M,T,N}
    Base.checkbounds(a, I...)
    Base.unsafe_getindex(a, I...)
end

function Base.unsafe_getindex(a::OuterProduct{ARRAYS,<:NTuple{M,Bool},T,N}, I::Vararg{Int,N}) where {ARRAYS,M,T,N}
    r = one(eltype(a))
    i = _splitI!(a,I...)
    for (j,ai) in enumerate(a.arrays)
        r *= ai[_getsplitindex(a, j)...]
    end
    r
end

function _splitI!(a::OuterProduct, I::Vararg)
    d = ndims(a.arrays[1])
    i = 1
    for expand in a.outer
        if expand
            copyto!(a.index_scratch[i], 1, I, i, d)
            i += d
        else
            a.index_scratch[i][1] = I[i]
            i += 1
        end
    end
    nothing
end

_getsplitindex(a::OuterProduct, i) =
    map( x->length(x)==1 ? x[1] : x[i], a.index_scratch)

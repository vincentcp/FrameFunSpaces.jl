using ApproxFunBase: SubSpace


export dofsize, doflength

"""
    const TruncatedSpace1D{DD,DD,RR}  = SubSpace{DS,<:UnitRange{Int},DD,RR} where {DS,DD,RR}

A space that only consists of the first `N` elements of a full space
"""
const TruncatedSpace1D{DS,DD,RR} = SubSpace{DS,<:UnitRange{Int},DD,RR} where {DS,DD,RR<:Number}

"""
    const TruncatedSpaceND{DD,DD,RR}  = SubSpace{DS,<:UnitRange{Int},DD,RR} where {DS,DD,RR}

A space that only consists of the first `N` elements of a full space
"""
const TruncatedSpaceND{DS,IT,DD,RR} = SubSpace{DS,IT,DD,RR} where {DS,IT<:CartesianIndices{N,NTuple{N,Base.OneTo{Int}}} where N,DD,RR}




"""
    abstract type Dictionary{SPACE<:SubSpace}

A space of a dictionary of spaces
"""
abstract type DictionarySpace{DD,RR}<:Space{DD,RR}
end
const Dictionary = DictionarySpace

# Useful abstraction for special cases
const DictionarySpace1d{T<:Number,RR} = DictionarySpace{<:Domain{T},RR} where {T,RR}
# Warning: not all Nd dictionaries have SVector type, they could have (S1,S2) type
const DictionarySpace2d{T<:Number,RR} = DictionarySpace{<:Domain{StaticArrays.SArray{Tuple{2},T,1,2}},RR} where {T,RR}
const DictionarySpace3d{T<:Number,RR} = DictionarySpace{<:Domain{StaticArrays.SArray{Tuple{3},T,1,3}},RR} where {T,RR}
const DictionarySpace4d{T<:Number,RR} = DictionarySpace{<:Domain{StaticArrays.SArray{Tuple{4},T,1,4}},RR} where {T,RR}
const DictionarySpaceNd{N,T<:Number,RR} = DictionarySpace{<:Domain{StaticArrays.SArray{Tuple{N},T,1,N}},RR} where {T,RR}

function Base.iterate(d::Dictionary)
    iter = eachindex(d)
    first_item, first_state = iterate(iter)
    (d[first_item], (iter, (first_item, first_state)))
end

function Base.iterate(d::Dictionary, state)
    iter, iter_tuple = state
    iter_item, iter_state = iter_tuple
    next_tuple = iterate(iter, iter_state)
    if next_tuple != nothing
        next_item, next_state = next_tuple
        (d[next_item], (iter,next_tuple))
    end
end

Base.eltype(::Type{DICT}) where {DICT<:Dictionary} = SingletonDictionary{DICT,Int,domaintype(DICT),rangetype(DICT)}

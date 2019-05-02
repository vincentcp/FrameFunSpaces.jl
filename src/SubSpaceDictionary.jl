struct SubSpaceDictionary{S,IT,DD,RR} <:DictionarySpace{DD,RR}
    space   ::  SubSpace{S,IT,DD,RR}
end

@forward SubSpaceDictionary.space dimension, domaindimension, domaintype,
    rangetype, canonicalspace, canonicaldomain, domain, points
evaluate(coefficients::AbstractVector, dict::SubSpaceDictionary, x) =
    evaluate(coefficients, space(dict), x)
coefficients(coefs::AbstractVector,s1::SubSpaceDictionary{S,IT,DD,RR}, s2::SubSpace{S,IT,DD,RR}) where {S,IT,DD,RR} =
    coefs
coefficients(coefs::AbstractVector,s1::SubSpace{S,IT,DD,RR}, s2::SubSpaceDictionary{S,IT,DD,RR}) where {S,IT,DD,RR} =
    coefs
coefficients(coefs::AbstractArray,s1::SubSpaceDictionary{S,IT,DD,RR}, s2::SubSpace{S,IT,DD,RR}) where {S,IT,DD,RR} =
    coefs
coefficients(coefs::AbstractArray,s1::SubSpace{S,IT,DD,RR}, s2::SubSpaceDictionary{S,IT,DD,RR}) where {S,IT,DD,RR} =
    coefs

for T in (:SumSpace,:PiecewiseSpace,:TensorSpace,:ConstantSpace,:Space) # Resolve conflict
    @eval begin
        coefficients(coefs::AbstractVector,s1::SubSpaceDictionary, s2::ApproxFunBase.$T) =
            coefficients(coefs, canonicalspace(s1), s2)
        coefficients(coefs::AbstractVector,s1::ApproxFunBase.$T, s2::SubSpaceDictionary) =
            coefficients(coefs, s1, canonicalspace(s2))
        coefficients(coefs::AbstractArray,s1::SubSpaceDictionary, s2::ApproxFunBase.$T) =
            coefficients(coefs, canonicalspace(s1), s2)
        coefficients(coefs::AbstractArray,s1::ApproxFunBase.$T, s2::SubSpaceDictionary) =
            coefficients(coefs, s1, canonicalspace(s2))
    end
end

dimensions(dict::SubSpaceDictionary) = size(space(dict).indexes)
dimensions(dict::SubSpaceDictionary) = size(space(dict).indexes)

"""
    dictionary(space::SubSpace)

Convert a subspace to a dictionary space.
"""
dictionary(space::SubSpace) = SubSpaceDictionary(space)

space(dict::SubSpaceDictionary) = dict.space

Base.getindex(dict::Dictionary, i) = dictionary(SubSpace(space(dict), i))
Base.getindex(dict::Dictionary, i::AbstractVector{Int}) = dictionary(SubSpace(space(dict), i))
Base.getindex(dict::Dictionary, i::Int...) = dictionary(SubSpace(space(dict), CartesianIndex(i)))


ApproxFunBase.SubSpace(sp::SubSpace,kr) = SubSpace(sp.space,ApproxFunBase.reindex(sp,sp.indexes,to_indexes(kr)))
ApproxFunBase.reindex(sp::SubSpace, indexes, i) = indexes[i]
to_indexes(a::AbstractArray{Int}) = a
to_indexes(a::AbstractArray{CartesianIndex}) = a
to_indexes(a::Int) = a
to_indexes(a::CartesianIndex) = a


SubSpace(space::TruncatedSpace1D, i::Int) = SubSpace(space.space, space.indexes[i])


Base.eltype(::Type{<:SubSpaceDictionary{S,IT,DD,RR}}) where {S,IT,DD,RR} = SingletonDictionary

_univariate_elements(a, i) =
    _univariate_elements(a)[i]
_univariate_elements(dict::Dictionary) =
    map(dictionary, _univariate_elements(space(dict)))
_univariate_elements(space::SubSpace{<:TensorSpace,<:CartesianIndices}) =
    map(SubSpace, _univariate_elements(space.space), _univariate_elements(space.indexes))
_univariate_elements(cart::CartesianIndices) = cart.indices
_univariate_elements(cart::CartesianIndices, i) = cart.indices[i]
_univariate_elements(space::TensorSpace) = space.spaces
_univariate_elements(space::TensorSpace, i) = space.spaces[i]

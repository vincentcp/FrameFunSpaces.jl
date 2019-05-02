const SingletonSubSpace{DS,DD,RR} = SubSpace{DS,IT,DD,RR} where {DS,IT<:Union{Int,CartesianIndex},DD,RR}

const SingletonDictionary{DS,IT,DD,RR} = SubSpaceDictionary{DS,IT,DD,RR} where {DS,IT<:Union{Int,CartesianIndex},DD,RR}

Base.getindex(dict::SingletonDictionary, i...) = throw(error("Singleton can not be indexed"))

evaluate(x, f::SingletonDictionary) = Fun(f,SVector(convert(rangetype(f),1)))(x...)

(f::SingletonDictionary)(x...) = evaluate(x, f)

evaluate(x, f::SingletonDictionary{<:TensorSpace}) = _tensorsingleton_evaluate(x, f)
evaluate(x::SVector, f::SingletonDictionary{<:TensorSpace}) = _tensorsingleton_evaluate(x, f)

_tensorsingleton_evaluate(x, f) = prod(map(evaluate, x, elements(f)))

elements(f::SingletonDictionary{<:TensorSpace}) =
    map(dictionary, map(SubSpace, space(f).space.spaces, space(f).indexes.I))

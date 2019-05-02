__precompile__()

module FrameFunSpaces

using Reexport, Base, StaticArrays, DomainSets
using DomainSets: FullSpace
using MacroTools: @forward
@reexport using ApproxFunBase

import ApproxFunBase: domaintype, rangetype, prectype, dimension, domaindimension,
    domain, points, canonicalspace, canonicaldomain, evaluate, hasfasttransform, coefficients,
    FunStyle, dimensions



export Dictionary, domaintype, rangetype, domaindimension, domain, dimension, dimensions, points,
    OuterProduct

# include("tensorspacesupport.jl")
include("Dictionary.jl")
include("SubSpaceDictionary.jl")
include("SingletonDict.jl")

include("Fun.jl")

include("Fourier.jl")



# export DictionaryOperator
# abstract type DictionaryOperator{T} <: Operator{T} end


end

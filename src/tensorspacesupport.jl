
function evaluate(f::AbstractArray,S::SubSpace{<:TensorSpace},x)
    csp=canonicalspace(S)
    if spacescompatible(csp,S)
        error("Override evaluate for " * string(typeof(csp)))
    else
        evaluate(coefficients(f,S,csp),csp,x...)
    end
end
function ApproxFunBase.subspace_coefficients(a::AbstractArray,sp::SubSpace,dropsp::SubSpace)
    if sp == dropsp
        a
    else
        coefficients(a,sp,canonicalspace(sp),dropsp)
    end
end


function ApproxFunBase.subspace_coefficients(a::AbstractArray,sp::Space,dropsp::SubSpace)
    n=length(a)
    if sp == dropsp.space
        a[dropsp.indexes]
        # Array{eltype(a)}(undef, size(d))
        # for k in dropsp.indexes
        #     if k > n
        #         return ret
        #     end
        #     push!(ret,v[k])
        # end
        # ret
    else
        coefficients(v,sp,canonicalspace(dropsp),dropsp)
    end
end

function ApproxFunBase.subspace_coefficients(a::AbstractArray,dropsp::SubSpace,sp::Space)
    if sp==dropsp.space
        isempty(v) && return v
        ret = zeros(eltype(v),Tuple(dropsp.indexes[length(a)]))
        for k = eachindex(a)
            ret[dropsp.indexes[k]] = a[k]
        end
        ret
    else
        coefficients(a,dropsp,canonicalspace(dropsp),sp)
    end
end

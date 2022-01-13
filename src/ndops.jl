#
# Date created: 2022-01-13
# Author: aradclif
#
#
############################################################################################
#### 2022-01-13: a common location for items falling under the category of
# N-dimensional operations; origin is ndadd!, and also ndpush!.
# Eventually, this might moved to a separate package, along with other indexing items.
################
function ndadd!(fs::Vector{Function}, A::Array{T, N}, v::T, ks::Vararg{Any, N}) where {N} where {T<:Real}
    # idxs::NTuple{N, Int} = ntuple(i -> fs[i](ks[i]), Val(N))
    # A[idxs...] += v
    idxs = ntuple(i -> fs[i](ks[i]), N)
    ndadd!(A, v, idxs)
    return A
end
function ndadd!(fs::Vector{Function}, A::Array{T, N}, v::T, ks::Vector) where {N} where {T<:Real}
    idxs = ntuple(i -> fs[i](ks[i]), N)
    ndadd!(A, v, idxs)
    return A
end

function ndadd1!(fs::Vector{Function}, A::Array{T, N}, ks::Vararg{Any, N}) where {N} where {T<:Real}
    idxs = ntuple(i -> fs[i](ks[i]), N)
    ndadd1!(A, idxs)
    return A
end
function ndadd1!(fs::Vector{Function}, A::Array{T, N}, ks::Vector) where {N} where {T<:Real}
    idxs = ntuple(i -> fs[i](ks[i]), N)
    ndadd1!(A, idxs)
    return A
end

@generated function ndadd!(A::Array{T, N}, v::T, idxs::NTuple{N, Int}) where {N} where {T<:Real}
    quote
        Base.Cartesian.@nextract $N i d -> idxs[d]
        (Base.Cartesian.@nref $N A i) += v
        return A
    end
end
function ndadd!(A::Array{T, N}, v::T, idxs::Vararg{S, N}) where {N} where {T<:Real} where {S<:Integer}
    ndadd!(A, v, idxs)
end

@generated function ndadd1!(A::Array{T, N}, idxs::NTuple{N, Int}) where {N} where {T<:Real}
    quote
        Base.Cartesian.@nextract $N i d -> idxs[d]
        (Base.Cartesian.@nref $N A i) += one(T)
        return A
    end
end
function ndadd1!(A::Array{T, N}, idxs::Vararg{S, N}) where {N} where {T<:Real} where {S<:Integer}
    ndadd1!(A, idxs)
end


@generated function ndadd1_v2!(fs::Vector{Function}, A::Array{T, N}, ks::Vector) where {N} where {T<:Real}
    quote
        Base.Cartesian.@nexprs $N i -> i_i = fs[i](ks[i])
        (Base.Cartesian.@nref $N A i) += one(T)
        return A
    end
end


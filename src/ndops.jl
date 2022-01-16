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
@generated function ndadd!(fs::Vector{Function}, A::Array{T, N}, v::T,
                           ks::Vararg{Any, N}) where {N} where {T<:Real}
    quote
        Base.Cartesian.@ncall $N ndadd! A v i -> fs[i](ks[i])
        return A
    end
end
@generated function ndadd!(fs::Vector{Function}, A::Array{T, N}, v::T,
                           ks::Vector) where {N} where {T<:Real}
    quote
        Base.Cartesian.@ncall $N ndadd! A v i -> fs[i](ks[i])
        return A
    end
end

@generated function ndadd1!(fs::Vector{Function}, A::Array{T, N},
                            ks::Vararg{Any, N}) where {N} where {T<:Real}
    quote
        Base.Cartesian.@ncall $N ndadd1! A i -> fs[i](ks[i])
        return A
    end
end
@generated function ndadd1!(fs::Vector{Function}, A::Array{T, N}, ks::Vector) where {T<:Real, N}
    # # Version 1
    # quote
    #     Base.Cartesian.@nexprs $N i -> i_i = fs[i](ks[i])
    #     (Base.Cartesian.@nref $N A i) += one(T)
    #     return A
    # end
    # # Version 2
    # quote
    #     ndadd1!(A, (Base.Cartesian.@ntuple $N i -> fs[i](ks[i])))
    # end
    # # Version 3: success
    quote
        Base.Cartesian.@ncall $N ndadd1! A i -> fs[i](ks[i])
        return A
    end
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

################
#### 2022-01-13: p. 639: push! strategy

@generated function ndpush!(fs::Vector{Function}, A::Array{Vector{Vector{T}}, N}, v::Vector{T},
                            ks::Vararg{S, N}) where {N} where {T} where {S}
    # idxs = ntuple(i -> fs[i](ks[i]), N)
    # ndpush!(A, v, idxs)
    # return A
    quote
        Base.Cartesian.@ncall $N ndpush! A v i -> fs[i](ks[i])
        return A
    end
end

@generated function ndpush!(fs::Vector{Function}, A::Array{Vector{Vector{T}}, N}, v::Vector{T},
                            ks::Vararg{Any, N}) where {N} where {T}
    # idxs = ntuple(i -> fs[i](ks[i]), N)
    # ndpush!(A, v, idxs)
    # return A
    quote
        Base.Cartesian.@ncall $N ndpush! A v i -> fs[i](ks[i])
        return A
    end
end

@generated function ndpush!(fs::Vector{Function}, A::Array{Vector{Vector{T}}, N}, v::Vector{T},
                            ks::Vector) where {N} where {T}
    # idxs = ntuple(i -> fs[i](ks[i]), N)
    # ndpush!(A, v, idxs)
    # return A
    quote
        Base.Cartesian.@ncall $N ndpush! A v i -> fs[i](ks[i])
        return A
    end
end

@generated function ndpush!(A::Array{Vector{Vector{T}}, N}, v::Vector{T},
                            idxs::NTuple{N, Int}) where {N} where {T}
    quote
        Base.Cartesian.@nextract $N i d -> idxs[d]
        push!((Base.Cartesian.@nref $N A i), v)
        return A
    end
end

function ndpush!(A::Array{Vector{Vector{T}}, N}, v::Vector{T},
                 idxs::Vararg{S, N}) where {N} where {T} where {S<:Integer}
    ndpush!(A, v, idxs)
end

# vv = [1,2,3]
# @code_warntype ndpush!(𝓃𝓂idxs, vv, (1, 2))
# fs = [x -> parse(Int, x), x -> 2]
# ks = ["1", 3]
# @code_warntype ndpush!(fs, 𝓃𝓂idxs, vv, ks...)
# @benchmark ndpush!(fs, 𝓃𝓂idxs, vv, ks)

# A = zeros(Int, 2, 2);
# @benchmark ndadd1!(fs, A, ks)

# @benchmark ndadd1!(A, (1, 2))

# @benchmark ndadd1_v2!(fs, A, ks)

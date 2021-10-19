#
# Date created: 2021-10-19
# Author: aradclif
#
#
############################################################################################
#### AbstractPathKey, AbstractPathKeys: p. 447-449, 456-458, 2021-09-15
# Type hierarchy of functors
abstract type AbstractPathKey{T} end

abstract type AbstractIndexedPathKey{T<:Function} <: AbstractPathKey{T} end

struct IndexedPathKey{T, U} <: AbstractIndexedPathKey{T}
    f::T
    i::U
    # IndexedPathKey{T, U}(f, i) where {T, U} = new(f, i)
end
# IndexedPathKey(f::T, i::U) where {T, U} = IndexedPathKey{T, U}(f, i)

#### Outer constructors for IndexedPathKey
IndexedPathKey(i::Int) = IndexedPathKey(identity, i)
IndexedPathKey(i::Vector{Int}) = IndexedPathKey(identity, i)

#### functor: default behavior for all IndexedPathKey
function (p::AbstractIndexedPathKey)(A)
    p.f(getindex(A, p.i))
end
################################################################
# Revised abstract type for AbstractPathKeys3
abstract type AbstractPathKeys{U<:Vector{T} where {T<:AbstractPathKey}} end
struct PathKeys{U} <: AbstractPathKeys{U}
    ftrs::U
    N::Int
    # PathKeys{T}(ftrs, N) where {T} = new(ftrs, N)
    # function PathKeys{U}(ftrs, N) where {U}
    #     new(copyto!(similar(ftrs), ftrs), N)
    # end
    # function PathKeys{U}(ftrs, N) where {U}
    #     let ftrs = copyto!(similar(ftrs), ftrs), N = N
    #         new(ftrs, N)
    #     end
    # end
end
# PathKeys(ftrs::T, N) where {T} = PathKeys{T}(ftrs, N)

#### Interface for AbstractPathKeys
Base.length(p::AbstractPathKeys) = p.N
Base.size(p::AbstractPathKeys) = (p.N,)

Base.iterate(p::AbstractPathKeys, state=1) =
    state > p.N ? nothing : (p.ftrs[state], state + 1)
Base.eltype(p::AbstractPathKeys) = eltype(p.ftrs)

Base.IndexStyle(::Type{<:AbstractPathKeys}) = IndexLinear()
Base.getindex(p::AbstractPathKeys, i::Int) = getindex(p.ftrs, i)
Base.getindex(p::AbstractPathKeys, I) = [p[i] for i in I]
Base.getindex(p::AbstractPathKeys, ::Colon) = p.ftrs
Base.firstindex(p::AbstractPathKeys) = 1
Base.lastindex(p::AbstractPathKeys) = p.N
Base.setindex!(p::AbstractPathKeys, v, i) = setindex!(p.ftrs, v, i)

function Base.isequal(p1::AbstractPathKeys, p2::AbstractPathKeys)
    p1 === p2 && return true
    length(p1) == length(p2) || return false
    for n = 1:length(p1)
        p1[n] == p2[n] || return false
    end
    return true
end
Base.:(==)(p1::AbstractPathKeys, p2::AbstractPathKeys) = isequal(p1, p2)

#### Outer constructors: PathKeys
PathKeys(ftrs) = PathKeys(ftrs, length(ftrs))

#### functor: PathKeys
function (p::AbstractPathKeys)(x, A)
    n = 1
    N = p.N
    fs = p.ftrs
    while n ≤ N
        @inbounds x[n] = fs[n](A)
        n += 1
    end
    return x
end

function (p::AbstractPathKeys)(A)
    x = Vector{Any}(undef, p.N)
    p(x, A)
end

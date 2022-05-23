#
# Date created: 2021-10-19
# Author: aradclif
#
#
############################################################################################
#### AbstractEdge, AbstractEdges: p. 447-449, 456-458, 2021-09-15
# Type hierarchy of functors
abstract type AbstractEdge{T} end

abstract type AbstractIndexedEdge{T<:Function} <: AbstractEdge{T} end

struct IndexedEdge{T, U} <: AbstractIndexedEdge{T}
    f::T
    i::U
    # IndexedEdge{T, U}(f, i) where {T, U} = new(f, i)
end
# IndexedEdge(f::T, i::U) where {T, U} = IndexedEdge{T, U}(f, i)

#### Interface for AbstractEdge
function Base.isequal(x::AbstractEdge, y::AbstractEdge)
    x === y && return true
    x.f == y.f || return false
    x.i == y.i || return false
    true
end
Base.:(==)(x::AbstractEdge, y::AbstractEdge) = isequal(x, y)

#### Outer constructors for IndexedEdge
IndexedEdge(i::Int) = IndexedEdge(identity, i)
IndexedEdge(i::Vector{Int}) = IndexedEdge(identity, i)
IndexedEdge(i::UnitRange{Int}) = IndexedEdge(identity, i)
IndexedEdge(i::CartesianIndex{N}) where {N} = IndexedEdge(identity, i)

#### functor: default behavior for all IndexedEdge
function (x::AbstractIndexedEdge)(A)
    @inbounds x.f(getindex(A, x.i))
end


############################################################################################
#### 2022-03-28: a revision to make N known at compile time
# As it happens, this makes effectively no difference in terms of performance based
# on tests of the relevant functor methods. Perhaps it would make a difference
# in other generated functions that resort to ::Val(p.N), but,
# testing reveals that performance is essentially the same.

abstract type AbstractEdges{U, N} end
struct Edges{U<:Tuple{Vararg{AbstractEdge, N}} where {N}, N} <: AbstractEdges{U, N}
    ftrs::U
end
Edges(ftrs::U) where {U<:Tuple{Vararg{AbstractEdge, N}}} where {N} = Edges{U,N}(ftrs)

Base.length(x::AbstractEdges{U, N}) where {U, N} = N

Base.iterate(x::AbstractEdges{U, N}, state=1) where {U, N} = state > N ? nothing : (x.ftrs[state], state + 1)
Base.eltype(x::AbstractEdges) = eltype(x.ftrs)

Base.IndexStyle(::Type{<:AbstractEdges}) = IndexLinear()
Base.getindex(x::AbstractEdges, i::Int) = getindex(x.ftrs, i)
Base.getindex(x::AbstractEdges, I) = Tuple(x[i] for i ∈ I)
Base.getindex(x::AbstractEdges, ::Colon) = x.ftrs
Base.firstindex(x::AbstractEdges) = 1
Base.lastindex(x::AbstractEdges{U, N}) where {U, N} = N

Base.keys(x::AbstractEdges) = keys(x.ftrs)

function Base.isequal(x::AbstractEdges, y::AbstractEdges)
    x === y && return true
    length(x) == length(y) || return false
    for i ∈ eachindex(x, y)
        x[i] == y[i] || return false
    end
    true
end
Base.:(==)(x::AbstractEdges, y::AbstractEdges) = isequal(x, y)

#### Outer constructors: Edges
Edges(gen::T) where {T<:Base.Generator} = Edges(Tuple(gen))
Edges(x::AbstractEdge) = Edges((x,))
Edges(xs::Vararg{AbstractEdge, N}) where {N} = Edges(xs)

#### functor: Edges
@generated function (x::AbstractEdges{U, N})(B, A) where {U, N}
    quote
        Base.Cartesian.@nexprs $N i -> B[i] = x[i](A)
        return B
    end
end

function (x::AbstractEdges{U, N})(A) where {U, N}
    B = Vector{Any}(undef, N)
    x(B, A)
end

#
# Date created: 2022-05-13
# Author: aradclif
#
#
############################################################################################
#### In essence, separate the concept of Edge and Label. While these can be implemented
# using the same functors, it is helpful for dispatch reasons (and conceptual clarity)
# to define separate types for each.

abstract type AbstractLabel{T} end

#### Interface for AbstractLabel
function Base.isequal(x::AbstractLabel, y::AbstractLabel)
    x === y && return true
    fieldcount(x) == fieldcount(y) || return false
    ns₁, ns₂ = fieldnames(x), fieldnames(y)
    for i ∈ eachindex(ns₁, ns₂)
        n₁, n₂ = ns₁[i], ns₂[i]
        n₁ === n₂ || return false
        getproperty(x, n₁) == getproperty(y, n₂) || return false
    end
    true
end
Base.:(==)(x::AbstractLabel, y::AbstractLabel) = isequal(x, y)

abstract type AbstractIndexedLabel{T<:Function} <: AbstractLabel{T} end

struct IndexedLabel{T, U} <: AbstractIndexedLabel{T}
    f::T
    i::U
    # IndexedLabel{T, U}(f, i) where {T, U} = new(f, i)
end
# IndexedLabel(f::T, i::U) where {T, U} = IndexedLabel{T, U}(f, i)

#### Outer constructors for IndexedLabel
IndexedLabel(i::Int) = IndexedLabel(identity, i)
IndexedLabel(i::Vector{Int}) = IndexedLabel(identity, i)
IndexedLabel(i::UnitRange{Int}) = IndexedLabel(identity, i)
IndexedLabel(i::CartesianIndex{N}) where {N} = IndexedLabel(identity, i)

#### functor: default behavior for all IndexedLabel
function (x::AbstractIndexedLabel)(A)
    @inbounds x.f(getindex(A, x.i))
end

#### For constructing singles, which can then be composed into a tuple

struct Label{T} <: AbstractLabel{T}
    f::T
end

Label() = Label(identity)

Base.isequal(x::Label, y::Label) = x.f === y.f

(x::AbstractLabel)(A) = x.f(A)

################################################################
#### For constructing tuples (or a single in the special case)

abstract type AbstractLabels{U, N} end
struct Labels{U<:Tuple{Vararg{T, N} where {T<:AbstractLabel}} where {N}, N} <: AbstractLabels{U, N}
    ftrs::U
end
Labels(ftrs::U) where {U<:Tuple{Vararg{T, N} where {T<:AbstractLabel}}} where {N} = Labels{U, N}(ftrs)

Base.length(x::AbstractLabels{U, N}) where {U, N} = N

Base.iterate(x::AbstractLabels{U, N}, state=1) where {U, N} = state > N ? nothing : (x.ftrs[state], state + 1)
Base.eltype(x::AbstractLabels) = eltype(x.ftrs)

Base.IndexStyle(::Type{<:AbstractLabels}) = IndexLinear()
Base.getindex(x::AbstractLabels, i::Int) = getindex(x.ftrs, i)
Base.getindex(x::AbstractLabels, I) = Tuple(x[i] for i ∈ I)
Base.getindex(x::AbstractLabels, ::Colon) = x.ftrs
Base.firstindex(x::AbstractLabels) = 1
Base.lastindex(x::AbstractLabels{U, N}) where {U, N} = N

Base.keys(x::AbstractLabels) = keys(x.ftrs)

function Base.isequal(x::AbstractLabels, y::AbstractLabels)
    x === y && return true
    length(x) == length(y) || return false
    for i ∈ eachindex(x, y)
        x[i] == y[i] || return false
    end
    true
end
Base.:(==)(x::AbstractLabels, y::AbstractLabels) = isequal(x, y)

#### Outer constructors: Labels
Labels(gen::T) where {T<:Base.Generator} = Labels(Tuple(gen))
Labels(x::AbstractLabel) = Labels((x,))
Labels(xs::Vararg{AbstractLabel, N}) where {N} = Labels(xs)

Labels(fs::Vararg{Function, N}) where {N} = Labels(ntuple(i -> Label(fs[i]), Val(N)))
Labels(fs::NTuple{N, Function}) where {N} = Labels(ntuple(i -> Label(fs[i]), Val(N)))

@generated function (x::AbstractLabels{U, N})(A) where {U, N}
    quote
        Base.Cartesian.@ntuple $N i -> x.ftrs[i](A)
    end
end

function (x::AbstractLabels{U, 1})(A) where {U}
    @inbounds x.ftrs[1](A)
end

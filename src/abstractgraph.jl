#
# Date created: 2021-10-26
# Author: aradclif
#
#
############################################################################################
abstract type AbstractGraph end

abstract type AbstractSimpleDiGraph <: AbstractGraph end

abstract type AbstractSimpleGraph <: AbstractGraph end

struct SimpleDiGraph <: AbstractSimpleDiGraph
    fadj::Dict{Any, SimpleDiGraph}
    data::Vector{Any}
end

struct SimpleGraph <: AbstractSimpleGraph
    fadj::Dict{Any, SimpleGraph}
    badj::Dict{Any, SimpleGraph}
    data::Vector{Any}
end

#### Interface: AbstractGraph
Base.iterate(g::A) where {A<:AbstractGraph} = iterate(g.fadj)
Base.iterate(g::A, state) where {A<:AbstractGraph} = iterate(g.fadj, state)

Base.length(g::A) where {A<:AbstractGraph} = length(g.fadj)

Base.eltype(::A) where {A<:AbstractGraph} = Pair{Any, A}

#### Opt-ins: AbstractGraph
Base.keys(g::A) where {A<:AbstractGraph} = keys(g.fadj)

Base.values(g::A) where {A<:AbstractGraph} = values(g.fadj)

Base.pairs(g::A) where {A<:AbstractGraph} = pairs(g.fadj)

# Other functions which can be feasible must return nothing in order for
# the same behavior to be achieved.
function Base.get(f::Function, g::A, k1)::Union{Nothing, A} where {A<:AbstractGraph}
    get(f, g.fadj, k1)
end

function Base.get(f::Function, g::A, k1, k2)::Union{Nothing, A} where {A<:AbstractGraph}
    tmp = get(f, g, k1);
    tmp === nothing ? nothing : get(f, tmp, k2)
end
function Base.get(f::Function, g::A, k1, k2, ks::Vararg{Any, N})::Union{Nothing, A} where {N} where {A<:AbstractGraph}
    tmp = get(f, g, k1)
    tmp === nothing ? (return nothing) : (tmp = get(f, tmp, k2))
    tmp === nothing && return nothing
    get(f, tmp, ks...)
    # Alternative 1
    # tmp = get(f, g, p)
    # tmp === nothing && return nothing
    # tmp = get(f, tmp, q)
    # tmp === nothing && return nothing
    # get(f, tmp, ps...)
end

# methods that supply () -> nothing; covers dispatches (1), (2) and (3)
_returnnothing() = nothing
Base.get(g::A, k1) where {A<:AbstractGraph} = get(_returnnothing, g, k1)
Base.get(g::A, k1, k2) where {A<:AbstractGraph} = get(_returnnothing, g, k1, k2)
Base.get(g::A, k1, k2, ks...) where {A<:AbstractGraph} = get(_returnnothing, g, k1, k2, ks...)

# See notes in abstractnode.jl on options and incomplete benchmarking
haspath(g::AbstractGraph, k1) = isa(get(_returnnothing, g, k1), AbstractGraph)
haspath(g::AbstractGraph, k1, ks...) = isa(get(_returnnothing, g, k1, ks...), AbstractGraph)

# Enforcing type-stability on the return value should not be necessary
function Base.get!(f::Function, x::A, k1)::A where {A<:AbstractNode}
    get!(f, x.link, k1)
end
function Base.get!(f::Function, x::A, k1, k2)::A where {A<:AbstractNode}
    get!(f, get!(f, x, k1), k2)
end
function Base.get!(f::Function, x::A, k1, k2, ks::Vararg{S, N})::A where {S,N} where {A<:AbstractNode}
    tmp = get!(f, x, k1, k2)
    for k in ks
        tmp = get!(f, tmp, k)
    end
    tmp
end
function Base.get!(f::Function, x::A, k1, k2, ks::Vararg{Any, N})::A where {N} where {A<:AbstractNode}
    get!(f, get!(f, get!(f, x, k1), k2), ks...)
end

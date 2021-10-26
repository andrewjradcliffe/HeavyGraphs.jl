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
# Despite the temptation to re-name keys≡edges and values≡vertices,
# it seems logical to avoid re-naming known functions.
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
function Base.get!(f::Function, g::A, k1) where {A<:AbstractGraph}
    get!(f, g.fadj, k1)
end
function Base.get!(f::Function, g::A, k1, k2) where {A<:AbstractGraph}
    get!(f, get!(f, g, k1), k2)
end
function Base.get!(f::Function, g::A, k1, k2, ks::Vararg{S, N}) where {S,N} where {A<:AbstractGraph}
    tmp = get!(f, g, k1, k2)
    for k in ks
        tmp = get!(f, tmp, k)
    end
    tmp
end
function Base.get!(f::Function, g::A, k1, k2, ks::Vararg{Any, N}) where {N} where {A<:AbstractGraph}
    get!(f, get!(f, get!(f, g, k1), k2), ks...)
end

# Likewise, type-stability enforcement by annotated return type should be unnecessary.
function Base.getindex(g::A, k1)::A where {A<:AbstractGraph}
    getindex(g.fadj, k1)
end
function Base.getindex(g::A, k1, k2)::A where {A<:AbstractGraph}
    getindex(getindex(g, k1), k2)
end
function Base.getindex(g::A, k1, k2, ks::Vararg{Any, N})::A where {N} where {A<:AbstractGraph}
    getindex(getindex(getindex(g, k1), k2), ks...)
end
function Base.getindex(g::A, k1, k2, ks::Vararg{S, N})::A where {S, N} where {A<:AbstractGraph}
    tmp = getindex(g, k1, k2)
    for k in ks
        tmp = getindex(tmp, k)
    end
    tmp
end

Base.setindex!(g::AbstractGraph, value, k1) = (setindex!(g.fadj, value, k1); g)
function Base.setindex!(g::AbstractGraph, value, k1, ks...)
    tmp = g[k1, ks[1:end-1]...]
    tmp[ks[end]] = value
    g
end

Base.push!(g::AbstractGraph, p::Pair) = setindex!(g, p.second, p.first)
Base.push!(g::AbstractGraph, p::Pair, q::Pair) = push!(push!(g, p), q)
Base.push!(g::AbstractGraph, p::Pair, q::Pair, r::Pair...) = push!(push!(push!(g, p), q), r...)

Base.sizehint!(g::AbstractGraph, newsz) = sizehint!(g.fadj, newsz)

function Base.empty!(g::AbstractGraph)
    empty!(g.fadj)
    empty!(g.data)
    return g
end
function Base.empty!(g::AbstractSimpleGraph)
    empty!(g.fadj)
    empty!(g.badj)
    empty!(g.data)
    return g
end

Base.pop!(g::AbstractGraph, key) = pop!(g.fadj, key)
Base.pop!(g::AbstractGraph, key, default) = pop!(g.fadj, key, default)
Base.pop!(g::AbstractGraph) = pop!(g.fadj)

Base.delete!(g::AbstractGraph, key) = (delete!(g.fadj, key); g)

Base.isempty(g::AbstractGraph) = isempty(g.fadj)
Base.isempty(g::AbstractSimpleGraph) = isempty(g.fadj) && isempty(g.badj)

Base.filter!(pred, g::AbstractGraph) = (filter!(pred, g.fadj); g)

function Base.merge!(a::AbstractGraph, bs::AbstractGraph...)
    merge!(a.fadj, getproperty.(bs, :fadj)...)
    for b in bs
        isempty(b.data) || append!(a.data, b.data)
    end
    a
end

function Base.merge!(a::AbstractSimpleGraph, bs::AbstractSimpleGraph...)
    merge!(a.fadj, getproperty.(bs, :fadj)...)
    merge!(a.badj, getproperty.(bs, :badj)...)
    for b in bs
        isempty(b.data) || append!(a.data, b.data)
    end
    a
end

function Base.merge(a::T, bs::T...) where {T<:AbstractGraph}
    tmp = T()
    merge!(tmp, a, bs...)
    tmp
end

################
function _size!(g::AbstractGraph, sz::Vector{Int}, C::Int)
    length(sz) < C && push!(sz, 0)
    @inbounds sz[C] += 1
    isempty(g) && return sz
    C̃ = C + 1
    for p in g
        _size!(p.second, sz, C̃)
    end
    return sz
end

Base.size(g::AbstractGraph) = tuple(_size!(g, Int[0], 1)...)

# type-assertion on return value should be unnecessary.
function _sizeat(g::AbstractGraph, N::Int, C::Int)#::Int
    C̃ = C + 1
    C̃ == N && return length(g)
    s = 0
    for p in g
        s += _sizeat(p.second, N, C̃)
    end
    return s
end

Base.size(g::AbstractGraph, d::Int) = _sizeat(g, d, 1)

################
function depth(g::AbstractGraph, C::Int=0)
    C̃ = C + 1
    isempty(g) && return C̃
    C̃ₘ = C̃
    for p in g
        C̃ₘ = max(C̃ₘ, depth(p.second, C̃))
    end
    return C̃ₘ
end

function maxbreadth(g::AbstractGraph)
    b::Int = length(g)
    isempty(g) && return b
    for p in g
        b = max(b, maxbreadth(p.second))
    end
    return b
end

function rlength(g::AbstractGraph, C::Int=0)
    C̃ = C + 1
    isempty(g) && return C̃
    for p in g
        C̃ += rlength(p.second, 0)
    end
    return C̃
end

################################################################
#### comparison operators
function Base.isequal(a::AbstractGraph, b::AbstractGraph)
    a === b && return true
    a.fadj == b.fadj && a.data == b.data
end
function Base.isequal(a::AbstractSimpleGraph, b::AbstractSimpleGraph)
    a === b && return true
    a.fadj == b.fadj && a.badj == b.badj && a.data == b.data
end
Base.:(==)(a::AbstractGraph, b::AbstractGraph) = Base.isequal(a, b)

################################################################
#### Outer constructors: SimpleDiGraph
SimpleDiGraph(fadj::Dict{Any,SimpleDiGraph}) = SimpleDiGraph(fadj, [])
SimpleDiGraph(data::Vector{Any}) = SimpleDiGraph(Dict{Any,SimpleDiGraph}(), data)
SimpleDiGraph() = SimpleDiGraph(Dict{Any,SimpleDiGraph}())

#### Outer constructors: SimpleGraph
SimpleGraph(fadj::Dict{Any,SimpleGraph}, badj::Dict{Any,SimpleGraph}) = SimpleGraph(fadj, badj, [])
SimpleGraph(fadj::Dict{Any,SimpleGraph}) = SimpleGraph(fadj, Dict{Any,SimpleGraph}())
SimpleGraph(data::Vector{Any}) = SimpleGraph(Dict{Any,SimpleGraph}(), Dict{Any,SimpleGraph}(), data)
SimpleGraph(fadj::Dict{Any,SimpleGraph}, data::Vector{Any}) =
    SimpleGraph(fadj, Dict{Any,SimpleGraph}(), data)
SimpleGraph() = SimpleGraph(Dict{Any,SimpleGraph}())

#### Conversion
# Tentatively defined, but far less complex than in node.jl and abstractnode.jl
Base.convert(::Type{T}, g::AbstractSimpleGraph) where {T<:AbstractSimpleGraph} = T(g)
Base.convert(::Type{T}, g::AbstractSimpleDiGraph) where {T<:AbstractSimpleDiGraph} = T(g)
Base.convert(::Type{T}, g::T) where {T<:AbstractSimpleGraph} = g
Base.convert(::Type{T}, g::T) where {T<:AbstractSimpleDiGraph} = g

# SimpleGraph to SimpleDiGraph is the most sensible conversion. From DiGraph to Graph is fuzzy.
Base.convert(::Type{T}, g::AbstractSimpleGraph) where {T<:AbstractSimpleDiGraph} = T(g)

#### Promotion
# Not currently well-defined, as conversion of a SimpleGraph to SimpleDiGraph is
# actually a simplification. In other words, one can convert SimpleGraph to SimpleDiGraph,
# but the converse is not true. This means that one should not define the conversion,
# as convert(SimpleGraph, convert(SimpleDiGraph, SimpleGraph)) will not result in the
# same SimpleGraph one started with.
# -- promotion becomes more meaningful once one has multiple subtypes under
# AbstractSimpleDiGraph and AbstractSimpleGraph

#### Outer constructors, within Graph types
SimpleDiGraph(g::AbstractSimpleDiGraph) = SimpleDiGraph(g.fadj, g.data)
SimpleGraph(g::AbstractSimpleGraph) = SimpleGraph(g.fadj, g.badj, g.data)

# between Graph types
SimpleDiGraph(g::AbstractSimpleGraph) = SimpleDiGraph(g.fadj, g.data)

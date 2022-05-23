#
# Date created: 2022-05-23
# Author: aradclif
#
#
############################################################################################
abstract type AbstractGraph{T} end

abstract type AbstractSimpleDiGraph{T} <: AbstractGraph{T} end

abstract type AbstractSimpleGraph{T} <: AbstractGraph{T} end

struct SimpleDiGraph{T} <: AbstractSimpleDiGraph{T}
    fadj::Dict{T, SimpleDiGraph{T}}
    data::Vector{Any}
end

struct SimpleGraph{T} <: AbstractSimpleGraph{T}
    fadj::Dict{T, SimpleGraph{T}}
    badj::Dict{T, SimpleGraph{T}}
    data::Vector{Any}
end

#### Interface: AbstractGraph
Base.iterate(g::A) where {A<:AbstractGraph} = iterate(g.fadj)
Base.iterate(g::A, state) where {A<:AbstractGraph} = iterate(g.fadj, state)

Base.length(g::A) where {A<:AbstractGraph} = length(g.fadj)

Base.eltype(::A) where {A<:AbstractGraph{T}} where {T} = Pair{T, A}

Base.hash(g::AbstractGraph, h::UInt) = hash(g.fadj, hash(g.data, h))

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
    tmp = get(f, g, k1)
    tmp === nothing ? nothing : get(f, tmp, k2)
end
function Base.get(f::Function, g::A, k1, k2, ks::Vararg{Any, N})::Union{Nothing, A} where {N} where {A<:AbstractGraph}
    tmp = get(f, g, k1)
    tmp === nothing ? (return nothing) : (tmp = get(f, tmp, k2))
    tmp === nothing && return nothing
    get(f, tmp, ks...)
    # Alternative 1
    # tmp = get(f, g, k1)
    # tmp === nothing && return nothing
    # tmp = get(f, tmp, k2)
    # tmp === nothing && return nothing
    # get(f, tmp, ks...)
end

# methods that supply () -> nothing; covers dispatches (1), (2) and (3)
_returnnothing() = nothing
Base.get(g::A, k1) where {A<:AbstractGraph} = get(_returnnothing, g, k1)
Base.get(g::A, k1, k2) where {A<:AbstractGraph} = get(_returnnothing, g, k1, k2)
Base.get(g::A, k1, k2, ks...) where {A<:AbstractGraph} = get(_returnnothing, g, k1, k2, ks...)

# See notes in abstractnode.jl on options and incomplete benchmarking.
# Technically, this should be named haswalk
haspath(g::A, k1) where {A<:AbstractGraph} = isa(get(_returnnothing, g, k1), AbstractGraph)
haspath(g::A, k1, k2) where {A<:AbstractGraph} = isa(get(_returnnothing, g, k1, k2), AbstractGraph)
haspath(g::A, k1, k2, ks...) where {A<:AbstractGraph} = isa(get(_returnnothing, g, k1, k2, ks...), AbstractGraph)


# Simplified methods, eliminating unnecessary 2-arg with Vararg version
@inline function Base.get!(f::Function, x::A, k₁) where {A<:AbstractGraph}
    get!(f, x.fadj, k₁)
end
@generated function Base.get!(f::Function, x::A, k₁,
                              ks::Vararg{T, N}) where {T, N} where {A<:AbstractGraph}
    quote
        tmp = get!(f, x, k₁)
        Base.Cartesian.@nexprs $N i -> tmp = get!(f, tmp, ks[i])
        return tmp
    end
end
@generated function Base.get!(f::Function, x::A, k₁,
                              ks::Vararg{Any, N}) where {N} where {A<:AbstractGraph}
    quote
        tmp = get!(f, x, k₁)
        Base.Cartesian.@nexprs $N i -> tmp = get!(f, tmp, ks[i])
        return tmp
    end
end

@inline function Base.getindex(x::A, k₁) where {A<:AbstractGraph}
    getindex(x.fadj, k₁)
end
@generated function Base.getindex(x::A, k₁, ks::Vararg{T, N}) where {T, N} where {A<:AbstractGraph}
    quote
        tmp = getindex(x, k₁)
        Base.Cartesian.@nexprs $N i -> tmp = getindex(tmp, ks[i])
        return tmp
    end
end
@generated function Base.getindex(x::A, k₁, ks::Vararg{Any, N}) where {N} where {A<:AbstractGraph}
    quote
        tmp = getindex(x, k₁)
        Base.Cartesian.@nexprs $N i -> tmp = getindex(tmp, ks[i])
        return tmp
    end
end

# setindex! is always lowered to (setindex!(A, v, k); v).
# In the event that these functions are called directly, the signatures below
# ensure that the correct return values are given in both cases:
# g[k] = v                     : returns v
# g[k, ks...] = v              : returns v
# setindex!(g, v, k)           : returns g
# setindex!(g, v, k, ks...)    : returns g
Base.setindex!(g::A, value, k1) where {A<:AbstractGraph} = (setindex!(g.fadj, value, k1); g)
function Base.setindex!(g::A, value, k1, ks...) where {A<:AbstractGraph}
    tmp = g[k1, ks[1:end-1]...]
    tmp[ks[end]] = value
    # value
    g
end

Base.push!(g::A, p::Pair) where {A<:AbstractGraph} = setindex!(g, p.second, p.first) # push!(g.fadj, p)
Base.push!(g::A, p::Pair, q::Pair) where {A<:AbstractGraph} = push!(push!(g, p), q)
Base.push!(g::A, p::Pair, q::Pair, r::Pair...) where {A<:AbstractGraph} = push!(push!(push!(g, p), q), r...)

Base.sizehint!(g::A, newsz) where {A<:AbstractGraph} = sizehint!(g.fadj, newsz)

function Base.empty!(g::A) where {A<:AbstractGraph}
    empty!(g.fadj)
    empty!(g.data)
    return g
end
function Base.empty!(g::A) where {A<:AbstractSimpleGraph}
    empty!(g.fadj)
    empty!(g.badj)
    empty!(g.data)
    return g
end

Base.pop!(g::A, key) where {A<:AbstractGraph} = pop!(g.fadj, key)
Base.pop!(g::A, key, default) where {A<:AbstractGraph} = pop!(g.fadj, key, default)
function Base.pop!(g::A)::Pair{Any, A} where {A<:AbstractGraph}
    pop!(g.fadj)
end


Base.delete!(g::A, key) where {A<:AbstractGraph} = (delete!(g.fadj, key); g)

Base.isempty(g::A) where {A<:AbstractGraph} = isempty(g.fadj)
# Base.isempty(g::A) where {A<:AbstractSimpleGraph} = isempty(g.fadj) && isempty(g.badj)

function Base.filter(pred, g::A) where {A<:AbstractSimpleDiGraph}
    A(filter(pred, g.fadj), g.data)
end

Base.filter!(pred, g::A) where {A<:AbstractGraph} = (filter!(pred, g.fadj); g)

# function Base.intersect(a::T, b::T) where {T<:AbstractGraph}
#     T(Dict(a.fadj ∩ b.fadj), a.data ∩ b.data)
# end

function Base.merge!(a::A, bs::A...) where {A<:AbstractGraph}
    merge!(a.fadj, getproperty.(bs, :fadj)...)
    for b ∈ bs
        isempty(b.data) || append!(a.data, b.data)
    end
    a
end

function Base.merge!(a::AbstractSimpleGraph, bs::AbstractSimpleGraph...)
    merge!(a.fadj, getproperty.(bs, :fadj)...)
    merge!(a.badj, getproperty.(bs, :badj)...)
    for b ∈ bs
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
    for p ∈ g
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
    for p ∈ g
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
    for p ∈ g
        C̃ₘ = max(C̃ₘ, depth(p.second, C̃))
    end
    return C̃ₘ
end

function maxbreadth(g::AbstractGraph)
    b = length(g)
    isempty(g) && return b
    for p ∈ g
        b = max(b, maxbreadth(p.second))
    end
    return b
end

function findmaxbreadth(g::AbstractGraph, C::Int=0)
    C̃ = C + 1
    b = length(g)
    isempty(g) && return (b, C̃)
    C̃ₘ = C̃
    for p ∈ g
        (b̂, Ĉ) = findmaxbreadth(p.second, C̃)
        (b, C̃ₘ) = b̂ > b ? (b̂, Ĉ) : (b, C̃ₘ)
    end
    return (b, C̃ₘ)
end

function rlength(g::AbstractGraph, C::Int=0)
    C̃ = C + 1
    isempty(g) && return C̃
    for p ∈ g
        C̃ += rlength(p.second, 0)
    end
    return C̃
end

# Sometimes faster
function rlength2(g::AbstractGraph)
    C̃ = 1
    isempty(g) && return C̃
    for p ∈ g
        C̃ += rlength2(p.second)
    end
    return C̃
end

"""
    maxsize(i::Int, D::Int, N::Vector{Int})

Compute the maximum number of vertices for a graph with `D` recursively expanded levels,
the number of unique edges per which are given by `N`, starting from level `i`.

# Arguments
- `i::Int`: index of current level. On initial call, the level at which recursion begins.
- `D::Int`: total levels. On initial call, the total recursive depth.

# Examples
```julia-repl
julia> N = [3800, 440, 154];
julia> maxsize(1, 3, N)
259163800
julia> maxsize(2, 3, N)
68200
julia> maxsize(1, 2, N)
1675800
```
"""
function maxsize(i::Int, D::Int, N::Vector{Int})
    # 1 ≤ i ≤ D - 1 ? begin return N[i] * nodeenum(i + 1, D, N) + N[i] end : begin
    #     return N[i] end
    1 ≤ i ≤ D - 1 ? N[i] * maxsize(i + 1, D, N) + N[i] : N[i]
end

"""
    maxsize(N::Vector{Int})

Compute the maximun number of vertices for a graph with recursively expanded levels,
the number of unique edges per which are given by `N`.
"""
maxsize(N::Vector{Int}) = maxsize(1, length(N), N)
################################################################
#### comparison operators
function Base.isequal(a::AbstractGraph, b::AbstractGraph)
    a === b && return true
    # a.data == b.data || return false
    a.fadj == b.fadj && a.data == b.data
    # isequal(a.fadj, b.fadj) && isequal(a.data, b.data)
end
# Actually, given the circular definition which results from using both
# forward and backward adjacency lists, SimpleGraph can be shown to be equal
# if only the forward adjacency lists are checked.
# function Base.isequal(a::AbstractSimpleGraph, b::AbstractSimpleGraph)
#     a === b && return true
#     a.fadj == b.fadj && a.badj == b.badj && a.data == b.data
# end
Base.:(==)(a::AbstractGraph, b::AbstractGraph) = Base.isequal(a, b)

#### Outer constructors: SimpleDiGraph
SimpleDiGraph(fadj::Dict{T,SimpleDiGraph{T}}) where {T} = SimpleDiGraph(fadj, [])
SimpleDiGraph{T}(data::Vector{Any}) where {T} = SimpleDiGraph(Dict{T,SimpleDiGraph{T}}(), data)
SimpleDiGraph(data::Vector{Any}) = SimpleDiGraph{Any}(data)

SimpleDiGraph{T}() where {T} = SimpleDiGraph(Dict{T, SimpleDiGraph{T}}())
SimpleDiGraph() = SimpleDiGraph{Any}()
SimpleDiGraph(::Type{T}) where {T} = SimpleDiGraph{T}()

# @benchmark SimpleDiGraph{UInt}()

# fadj = Dict{UInt, SimpleDiGraph{UInt}}()
# @timev SimpleDiGraph(fadj, [])
# @timev SimpleDiGraph(Dict{UInt, SimpleDiGraph{UInt}}(), [])
# x = SimpleDiGraph(fadj, [])
# @code_warntype keys(x)
# data = []
# @code_warntype SimpleDiGraph(fadj, data)

# usdg() = SimpleDiGraph{UInt}()
# @timev get!(usdg, usdg(), 0x1, 0x2, 0x3)
# @code_warntype get!(usdg, x, 0x1, 0x2, 0x3)
# @code_warntype get!(sdg, sdg(), 0x1, 0x2, 0x3)
# @timev get!(sdg, sdg(), 0x1, 0x2, 0x3)

#### Outer constructors: SimpleGraph
SimpleGraph(fadj::Dict{T,SimpleGraph{T}}, badj::Dict{T,SimpleGraph{T}}) where {T} = SimpleGraph(fadj, badj, [])
SimpleGraph(fadj::Dict{T,SimpleGraph{T}}) where {T} = SimpleGraph(fadj, Dict{T,SimpleGraph{T}}())
SimpleGraph{T}(data::Vector{Any}) where {T} = SimpleGraph(Dict{T,SimpleGraph{T}}(), Dict{T,SimpleGraph{T}}(), data)
SimpleGraph(data::Vector{Any}) = SimpleGraph{Any}(data)

SimpleGraph(fadj::Dict{T,SimpleGraph{T}}, data::Vector{Any}) where {T} =
    SimpleGraph(fadj, Dict{T,SimpleGraph{T}}(), data)

SimpleGraph{T}() where {T} = SimpleGraph(Dict{T, SimpleGraph{T}}())
SimpleGraph() = SimpleGraph{Any}()
SimpleGraph(::Type{T}) where {T} = SimpleGraph{T}()

#### Conversion
# Tentatively defined, but far less complex than in node.jl and abstractnode.jl
# Base.convert(::Type{T}, g::AbstractSimpleGraph) where {T<:AbstractSimpleGraph} = T(g)
# Base.convert(::Type{T}, g::AbstractSimpleDiGraph) where {T<:AbstractSimpleDiGraph} = T(g)
# Base.convert(::Type{T}, g::T) where {T<:AbstractSimpleGraph} = g
# Base.convert(::Type{T}, g::T) where {T<:AbstractSimpleDiGraph} = g

# SimpleGraph to SimpleDiGraph is the most sensible conversion. From DiGraph to Graph is fuzzy.
# Actually, they are both well-defined.
# Base.convert(::Type{T}, g::AbstractSimpleGraph) where {T<:AbstractSimpleDiGraph} = T(g)
# Base.convert(::Type{T}, g::AbstractSimpleDiGraph) where {T<:AbstractSimpleGraph} = T(g)

#### Conversions in light of types
# SimpleDiGraph
Base.convert(::Type{T}, x::T) where {T<:AbstractGraph} = x
Base.convert(::Type{T}, x::AbstractSimpleDiGraph{S}) where {U, T<:AbstractSimpleDiGraph{U}, S} = T(x)

SimpleDiGraph{T}(x::SimpleDiGraph{T}) where {T} = x
SimpleDiGraph{T}(x::SimpleDiGraph{U}) where {T, U} =
    SimpleDiGraph(convert(Dict{T, SimpleDiGraph{T}}, x.fadj), x.data)

Base.promote_rule(::Type{T₁}, ::Type{T₂}) where {U₁, U₂, T₁<:AbstractSimpleDiGraph{U₁}, T₂<:AbstractSimpleDiGraph{U₂}} = SimpleDiGraph{promote_type(U₁, U₂)}

# SimpleGraph
Base.convert(::Type{T}, x::AbstractSimpleGraph{S}) where {U, T<:AbstractSimpleGraph{U}, S} = T(x)
SimpleGraph{T}(x::SimpleGraph{T}) where {T} = x
SimpleGraph{T}(x::SimpleGraph{U}) where {T, U} =
    SimpleGraph(convert(Dict{T, SimpleGraph{T}}, x.fadj), convert(Dict{T, SimpleGraph{T}}, x.badj), x.data)

Base.convert(::Type{T}, x::AbstractSimpleDiGraph{S}) where {U, T<:AbstractSimpleGraph{U}, S} = T(x)

SimpleGraph{T}(x::SimpleDiGraph{T}) where {T} = SimpleGraph(x.fadj, x.data)
SimpleGraph{T}(x::SimpleDiGraph{U}) where {T, U} =
    SimpleGraph(convert(Dict{T, SimpleGraph{T}}, x.fadj), x.data)

Base.promote_rule(::Type{T₁}, ::Type{T₂}) where {U₁, U₂, T₁<:AbstractSimpleGraph{U₁}, T₂<:AbstractSimpleGraph{U₂}} = SimpleGraph{promote_type(U₁, U₂)}
#### Promotion
# Not currently well-defined, as conversion of a SimpleGraph to SimpleDiGraph is
# actually a simplification. In other words, one can convert SimpleGraph to SimpleDiGraph,
# but the converse is not true. This means that one should not define the conversion,
# as convert(SimpleGraph, convert(SimpleDiGraph, SimpleGraph)) will not result in the
# same SimpleGraph one started with.
# -- promotion becomes more meaningful once one has multiple subtypes under
# # AbstractSimpleDiGraph and AbstractSimpleGraph
# function Base.promote_rule(::Type{T}, ::Type{U}) where {T<:AbstractSimpleGraph} where {U<:AbstractSimpleDiGraph}
#     T
# end


#### Outer constructors, within Graph types
SimpleDiGraph(g::AbstractSimpleDiGraph) = SimpleDiGraph(g.fadj, g.data)
SimpleGraph(g::AbstractSimpleGraph) = SimpleGraph(g.fadj, g.badj, g.data)

# between Graph types
SimpleDiGraph(g::AbstractSimpleGraph) = SimpleDiGraph(g.fadj, g.data)
# identical to:
# SimpleDiGraph(g::AbstractSimpleGraph) =
#     SimpleDiGraph(convert(Dict{Any, SimpleDiGraph}, g.fadj), g.data)

# see p. 553-557
function SimpleGraph(g::SimpleDiGraph{T}) where {T}
    sg = SimpleGraph{T}(g.data)
    fadj = sg.fadj
    badj = sg.badj
    for p ∈ g
        fadj[p.first] = SimpleGraph(p.second)
        fadj[p.first].badj[p.first] = sg
    end
    sg
end

#### Convenience shorthand constructors
sdg() = SimpleDiGraph()
sg() = SimpleGraph()

################################################################
#### Opt-ins: AbstractSimpleGraph
# bi-directional: p. 460-463, 2021-10-06
function isbidirectional(V₁::AbstractSimpleGraph, V₂::AbstractSimpleGraph, k)
    # Original
    # (k => V₂) ∈ V₁.fadj && (k => V₁) ∈ V₂.badj && return true
    # (k => V₂) ∈ V₁.badj && (k => V₁) ∈ V₂.fadj && return true
    # false
    # More efficient
    p₁ = (k => V₂)
    p₂ = (k => V₁)
    p₁ ∈ V₁.fadj && p₂ ∈ V₂.badj && return true
    p₁ ∈ V₁.badj && p₂ ∈ V₂.fadj && return true
    false
end

# this should actually be forwardget (abbrev. fget), and a backwardget should exist
# In fact, providing the dispatch to get! on AbstractSimpleGraph should be enough
function Base.get!(f::Function, V₁::AbstractSimpleGraph, k)
    V₂ = get!(f, V₁.fadj, k)
    setindex!(V₂.badj, V₁, k)
    V₂
end
function bget!(f::Function, V₁::AbstractSimpleGraph, k)
    # V₂ = get!(f, V₁, k) # get!(f, V₁.fadj, k)
    V₂ = get!(f, V₁.fadj, k)
    setindex!(V₂.badj, V₁, k)
    V₂
end
bget!(f::Function, V₁::AbstractSimpleGraph, k1, k2) = bget!(f, bget!(f, V₁, k1), k2)
bget!(f::Function, V₁::AbstractSimpleGraph, k1, k2, ks::Vararg{Any, N}) where {N} =
    bget!(bget!(f, bget!(f, V₁, k1), k2), ks...)
function bget!(f::Function, V₁::AbstractSimpleGraph, k1, k2, ks::Vararg{S, N}) where {S, N}
    tmp = bget!(f, V₁, k1, k2)
    for k ∈ ks
        tmp = bget!(f, tmp, k)
    end
    tmp
end

# get needs to be specialized, but it should be possible to handle it entirely with
# the dispatch on get for AbstractSimpleGraph. Interestingly,
# the b- methods are ≈ 1.2x faster than dispatching onto the generic methods.
# This may be due to the overhead of dispatch itself, but this is a surprising difference.
# function Base.get(f::Function, V₁::A, k)::Union{Nothing, A} where {A<:AbstractSimpleGraph}
#     V₂ = get(f, V₁.fadj, k)
#     # V₂ === nothing && return nothing
#     # isbidirectional(V₁, V₂, k) || return nothing
#     (V₂ === nothing || !isbidirectional(V₁, V₂, k)) && return nothing
#     V₂
# end
# function Base.get(f::Function, V₁::A, k1, k2)::Union{Nothing, A} where {A<:AbstractSimpleGraph}
#     V₂ = get(f, V₁, k1)
#     V₂ === nothing && return nothing
#     V₃ = get(f, V₂, k2)
# end

# function Base.get(f::Function, V₁::A, k1, k2, ks::Vararg{Any, N})::Union{Nothing, A} where {N} where {A<:AbstractSimpleGraph}
#     V₂ = get(f, V₁, k1)
#     V₂ === nothing && return nothing
#     V₃ = get(f, V₂, k2)
#     V₃ === nothing && return nothing
#     get(f, V₃, ks...)
# end

function bget(f::Function, V₁::A, k)::Union{Nothing, A} where {A<:AbstractSimpleGraph}
    V₂ = get(f, V₁.fadj, k)
    # V₂ = get(f, V₁, k)
    # V₂ === nothing && return nothing
    # isbidirectional(V₁, V₂, k) || return nothing
    (V₂ === nothing || !isbidirectional(V₁, V₂, k)) && return nothing
    V₂
end
function bget(f::Function, V₁::A, k1, k2)::Union{Nothing, A} where {A<:AbstractSimpleGraph}
    V₂ = bget(f, V₁, k1)
    V₂ === nothing && return nothing
    V₃ = bget(f, V₂, k2)
end

function bget(f::Function, V₁::A, k1, k2, ks::Vararg{Any, N})::Union{Nothing, A} where {N} where {A<:AbstractSimpleGraph}
    V₂ = bget(f, V₁, k1)
    V₂ === nothing && return nothing
    V₃ = bget(f, V₂, k2)
    V₃ === nothing && return nothing
    bget(f, V₃, ks...)
end

function hasbipath(V₁::AbstractSimpleGraph, k)
    isa(bget(_returnnothing, V₁, k), AbstractSimpleGraph)
end
function hasbipath(V₁::AbstractSimpleGraph, k1, k2)
    isa(bget(_returnnothing, V₁, k1, k2), AbstractSimpleGraph)
end
function hasbipath(V₁::AbstractSimpleGraph, k1, k2, ks...)
    isa(bget(_returnnothing, V₁, k1, k2, ks...), AbstractSimpleGraph)
end

# 2021-10-27, p. 554.
# See the comment on setindex!(g::AbstractGraph, v, k) -- this maintains ensures
# that the return value is consistent
function Base.setindex!(g::AbstractSimpleGraph, v, k)
    setindex!(g.fadj, v, k)
    v.badj[k] = g
    # v
end

function bgetindex(g::AbstractSimpleGraph, k)
    g.badj[k]
end
function bsetindex!(g::AbstractSimpleGraph, v, k)
    g.badj[k] = v
    g
end

############################################################################################
#### Methods of get_:(:op):(:field)!
# Code generation to cover the various function definitions.
# If used, one should then provide function prototypes to make docstrings visible.
# for op = (:push!, :append!), field = (:data)
#     fname = Symbol(:get_, field, op)
#     eval(:(function $fname(f::Function, g::AbstractGraph, v, k1)
#                tmp = get!(f, g, k1)
#                $op(tmp.$field, v)
#                tmp
#            end))
#     eval(:(function $fname(f::Function, g::AbstractGraph, v, k1, k2)
#                tmp = get!(f, g, k1, k2)
#                $op(tmp.$field, v)
#                tmp
#            end))
#     eval(:(function $fname(f::Function, g::AbstractGraph, v, k1, k2, ks::Vararg{Any, N}) where {N}
#                tmp = get!(f, g, k1, k2, ks...)
#                $op(tmp.$field, v)
#                tmp
#            end))
# end

# Direct approach: example of some code which would be generated
function get_datapush!(f::Function, g::AbstractGraph, v, k1)
    tmp = get!(f, g, k1)
    push!(tmp.data, v)
    tmp
end

function get_datapush!(f::Function, g::AbstractGraph, v, k1, ks...)
    tmp = get!(f, g, k1, ks...)
    push!(tmp.data, v)
    tmp
end

############################################################################################
#### Methods of grow_(:field)! : see p. 443-446, 451-455, 2021-09-14/15
function grow!(f::Function, g::AbstractGraph, p::AbstractEdges)
    for x ∈ p
        get!(f, g, x...)
    end
    return g
end

# Replaces splatting method
@generated function grow!(f::Function, g::AbstractGraph, p::AbstractEdges{U, N}, itr) where {U, N}
    quote
        for item ∈ itr
            tmp = g
            Base.Cartesian.@nexprs $N i -> tmp = get!(f, tmp, p[i](item))
        end
        return g
    end
end

# Convenience wrappers, but useful nonetheless
grow(f::Function, p::AbstractEdges) = grow!(f, f(), p)
grow(f::Function, p::AbstractEdges, itr) = grow!(f, f(), p, itr)

function datagrow!(f::Function, g::AbstractGraph, v::T, p::AbstractEdges) where {T<:Union{AbstractLabel, AbstractLabels}}
    for x ∈ p
        get_datapush!(f, g, v(x), x...)
    end
    return g
end

@generated function datagrow!(f::Function, g::AbstractGraph, v::T,
                              p::AbstractEdges{U, N}, itr) where {U, N} where {T<:Union{AbstractLabel, AbstractLabels}}
    quote
        for item ∈ itr
            tmp = g
            Base.Cartesian.@nexprs $N i -> tmp = get!(f, tmp, p[i](item))
            push!(tmp.data, v(item))
        end
        return g
    end
end


# Convenience wrappers, but useful nonetheless
datagrow(f::Function, v, p) = datagrow!(f, f(), v, p)
datagrow(f::Function, v, p, itr) = datagrow!(f, f(), v, p, itr)

# 2022-01-10: multi-step growth
function datagrow!(f::Function, x::AbstractGraph, vs::Vector{<:AbstractEdge},
                   ps::Vector{<:AbstractEdges}, itrs::Vector)
    for i ∈ eachindex(vs, ps, itrs)
        datagrow!(f, x, vs[i], ps[i], itrs[i])
    end
    return x
end
datagrow(f::Function, vs::Vector{<:AbstractEdge}, ps::Vector{<:AbstractEdges}, itrs::Vector) =
    datagrow!(f, f(), vs, ps, itrs)

# 2022-01-11: alternative meta
# @generated function datagrow!(f::Function, x::AbstractGraph,
#                               vs::Tuple{Vararg{T, N} where {T<:Union{AbstractLabel, AbstractLabels}}},
#                               ps::Tuple{Vararg{S, N} where {S<:AbstractEdges}}, itrs) where {N}
#     quote
#         Base.Cartesian.@nexprs $N i -> datagrow!(f, x, vs[i], ps[i], itrs[i])
#     end
# end
@generated function datagrow!(f::Function, x::AbstractGraph,
                              vs::Tuple{Vararg{Union{AbstractLabel, AbstractLabels}, N}},
                              ps::Tuple{Vararg{AbstractEdges, N}}, itrs) where {N}
    quote
        Base.Cartesian.@nexprs $N i -> datagrow!(f, x, vs[i], ps[i], itrs[i])
    end
end
# Convenience wrapper
# function datagrow(f::Function,
#                   vs::Tuple{Vararg{T, N} where {T<:Union{AbstractLabel, AbstractLabels}}},
#                   ps::Tuple{Vararg{S, N} where {S<:AbstractEdges}}, itrs) where {N}
#     datagrow!(f, f(), vs, ps, itrs)
# end
function datagrow(f::Function,
                  vs::Tuple{Vararg{Union{AbstractLabel, AbstractLabels}, N}},
                  ps::Tuple{Vararg{AbstractEdges, N}}, itrs) where {N}
    datagrow!(f, f(), vs, ps, itrs)
end


################################
# Growth from non-flat sources - p. 475, 2021-09-22
function datagrow!(f::Function, g::AbstractGraph, v::T, p::AbstractEdges, vitr, pitr) where {T<:Union{AbstractLabel, AbstractLabels}}
    x = Vector{Any}(undef, p.N)
    for (a, b) ∈ zip(vitr, pitr)
        p(x, b)
        get_datapush!(f, g, v(a), x...)
    end
    return g
end

# Possible parallel growth method
function tdatagrow!(f::Function, g::AbstractGraph, v::T, p::AbstractEdges,
                    itrsource::AbstractDict) where {T<:Union{AbstractLabel, AbstractLabels}}
    @sync for pₜ ∈ t
        Threads.@spawn datagrow!(f, pₜ.second, v, p, eachcol(itrsource[pₜ.first]))
        # Alternative:
        # let itr = eachcol(itrsource[p.first])
        #     Threads.@spawn datagrow!(f, p.second, v, p, itr)
        # end
    end
    return g
end

################
# 2021-11-10: map constructors for vectors of iterables
function mapgrow!(f::Function, gs::Vector{A}, p::AbstractEdges, itrs::Vector{T}) where {T} where {A<:AbstractGraph}
    @inbounds for i ∈ eachindex(gs, itrs)
        grow!(f, gs[i], p, itrs[i])
    end
    return gs
end
function mapgrow(f::Function, p::AbstractEdges, itrs::Vector{T}) where {T}
    gs = [f() for _ = 1:length(itrs)]
    mapgrow!(f, gs, p, itrs)
end

function tmapgrow!(f::Function, gs::Vector{A}, p::AbstractEdges, itrs::Vector{T}) where {T} where {A<:AbstractGraph}
    @sync @inbounds for i ∈ eachindex(gs, itrs)
        Threads.@spawn grow!(f, gs[i], p, itrs[i])
    end
    return gs
end

function tmapgrow(f::Function, p::AbstractEdges, itrs::Vector{T}) where {T}
    gs = [f() for _ = 1:length(itrs)]
    tmapgrow!(f, gs, p, itrs)
end

function tʳmapgrow!(f::Function, gs::Vector{A}, p::AbstractEdges, itrs::Vector{T}, M::Int=Threads.nthreads()) where {T} where {A<:AbstractGraph}
    ranges = equalranges(length(itrs), M)
    @sync for r ∈ ranges
        Threads.@spawn mapgrow!(f, gs[r], p, itrs[r])
    end
    return gs
end

function tʳmapgrow(f::Function, p::AbstractEdges, itrs::Vector{T}, M::Int=Threads.nthreads()) where {T}
    gs = [f() for _ = 1:length(itrs)]
    tʳmapgrow!(f, gs, p, itrs, M)
end

####
function mapdatagrow!(f::Function, gs::Vector{A}, v::U, p::AbstractEdges, itrs::Vector{T}) where {T} where {A<:AbstractGraph} where {U<:Union{AbstractLabel, AbstractLabels}}
    @inbounds for i ∈ eachindex(gs, itrs)
        datagrow!(f, gs[i], v, p, itrs[i])
    end
    return gs
end
function mapdatagrow(f::Function, v::U, p::AbstractEdges, itrs::Vector{T}) where {T} where {U<:Union{AbstractLabel, AbstractLabels}}
    gs = [f() for _ = 1:length(itrs)]
    mapdatagrow!(f, gs, v, p, itrs)
end

function tmapdatagrow!(f::Function, gs::Vector{A}, v::U, p::AbstractEdges, itrs::Vector{T}) where {T} where {A<:AbstractGraph} where {U<:Union{AbstractLabel, AbstractLabels}}
    @sync @inbounds for i ∈ eachindex(gs, itrs)
        Threads.@spawn datagrow!(f, gs[i], v, p, itrs[i])
    end
    return gs
end

function tmapdatagrow(f::Function, v::U, p::AbstractEdges, itrs::Vector{T}) where {T} where {U<:Union{AbstractLabel, AbstractLabels}}
    gs = [f() for _ = 1:length(itrs)]
    tmapdatagrow!(f, gs, v, p, itrs)
end

function tʳmapdatagrow!(f::Function, gs::Vector{A}, v::U, p::AbstractEdges, itrs::Vector{T}, M::Int=Threads.nthreads()) where {T} where {A<:AbstractGraph} where {U<:Union{AbstractLabel, AbstractLabels}}
    ranges = equalranges(length(itrs), M)
    @sync for r ∈ ranges
        Threads.@spawn mapdatagrow!(f, gs[r], v, p, itrs[r])
    end
    return gs
end

function tʳmapdatagrow(f::Function, v::U, p::AbstractEdges, itrs::Vector{T}, M::Int=Threads.nthreads()) where {T} where {U<:Union{AbstractLabel, AbstractLabels}}
    gs = [f() for _ = 1:length(itrs)]
    tʳmapdatagrow!(f, gs, v, p, itrs, M)
end


# function mapdatagrow!(f::Function, gs::Vector{A}, v::AbstractEdge, p::AbstractEdges, itrs::Vector{T}) where {T} where {A<:AbstractGraph}
#     for i ∈ eachindex(gs, itrs)
#         datagrow!(f, gs[i], v, p, itrs[i])
#     end
#     return gs
# end
# function mapdatagrow(f::Function, v::AbstractEdge, p::AbstractEdges, itrs::Vector{T}) where {T}
#     gs = Vector{typeof(f())}(undef, length(itrs))
#     # @inbounds for n ∈ eachindex(itrs)
#     #     gs[n] = datagrow(f, v, p, eachcol(itrs[n]))
#     # end
#     mapdatagrow!(f, gs, v, p, itrs)
#     return gs
# end

# function tmapdatagrow!(f::Function, gs::Vector{A}, v::AbstractEdge, p::AbstractEdges, itrs::Vector{T}, M::Int=Threads.nthreads()) where {T} where {A<:AbstractGraph}
#     N = length(itrs)
#     ranges = equalranges(N, M)
#     @inbounds @sync for m ∈ eachindex(ranges)
#         Threads.@spawn mapdatagrow!(f, gs[ranges[m]], v, p, itrs[ranges[m]])
#     end
#     return gs
# end

# function tmapdatagrow(f::Function, v::AbstractEdge, p::AbstractEdges, itrs::Vector{T},
#                       M::Int=Threads.nthreads()) where {T}
#     N = length(itrs)
#     ranges = equalranges(N, M)
#     gs = Vector{Vector{typeof(f())}}(undef, M)
#     # Threads.@threads for m ∈ 1:M #eachindex(ranges)
#     @inbounds @sync for m ∈ eachindex(ranges)
#         Threads.@spawn gs[m] = mapdatagrow(f, v, p, itrs[ranges[m]])
#     end
#     return vcat(gs...)
# end

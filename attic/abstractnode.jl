#
# Date created: 2021-10-25
# Author: aradclif
#
#
############################################################################################
#### Experiment in condensing the abstract type hierarchy
################################################################
abstract type AbstractNode end

# Simple Node
abstract type AbstractSimpleNode{T<:AbstractDict, U<:AbstractVector} <: AbstractNode end

struct SimpleNode{T, U} <: AbstractSimpleNode{T, U}
    link::T
    val::U
    SimpleNode{T, U}(link, val) where {T, U} = new(link, val)
end
SimpleNode(link::T, val::U) where {T, U} = SimpleNode{T, U}(link, val)

# Complex Node
abstract type AbstractComplexNode{T<:AbstractDict, U<:AbstractVector, V<:AbstractDict, W<:AbstractVector} <: AbstractNode end

struct ComplexNode{T, U, V, W} <: AbstractComplexNode{T, U, V, W}
    link::T
    val::U
    spec::V
    sval::W
    ComplexNode{T, U, V, W}(link, val, spec, sval) where {T, U, V, W} = new(link, val, spec, sval)
end
ComplexNode(link::T, val::U, spec::V, sval::W) where {T, U, V, W} =
    ComplexNode{T, U, V, W}(link, val, spec, sval)

#### Interface: AbstractNode
Base.iterate(x::A) where {A<:AbstractNode} = iterate(x.link)
Base.iterate(x::A, state) where {A<:AbstractNode} = iterate(x.link, state)
# Base.iterate(x::A) where {A<:AbstractNode{T, U}} where {T, U} = iterate(x.link)
# Base.iterate(x::A, state) where {A<:AbstractNode{T, U}} where {T, U} = iterate(x.link, state)
# function Base.iterate(t::A)::Union{Nothing, Tuple{Pair{Any, A}, Int}} where {T,U,A<:AbstractNode{T,U}}
#     iterate(t.link)
# end
# function Base.iterate(t::A, state)::Union{Nothing, Tuple{Pair{Any, A}, Int}} where {T,U,A<:AbstractNode{T,U}}
#     iterate(t.link, state)
# end

Base.length(x::A) where {A<:AbstractNode} = length(x.link)

Base.eltype(::A) where {A<:AbstractNode} = Pair{Any, A}
#### Opt-ins: AbstractNode
Base.keys(x::A) where {A<:AbstractNode} = keys(x.link)

Base.values(x::A) where {A<:AbstractNode} = values(x.link)

Base.pairs(x::A) where {A<:AbstractNode} = pairs(x.link)

# Base.get(x::A where {A<:AbstractNode{T, U}}, key, default) where {T, U} = get(x.link, key, default)

# f must be () -> nothing for this to work properly.
# Other functions which can be feasible must return nothing in order for
# the same behavior to be achieved.
function Base.get(f::Function, t::A, k1)::Union{Nothing, A} where {A<:AbstractNode}
    get(f, x.link, k1)
end

function Base.get(f::Function, t::A, k1, k2)::Union{Nothing, A} where {A<:AbstractNode}
    tmp = get(f, t, k1);
    tmp === nothing ? nothing : get(f, tmp, k2)
end
function Base.get(f::Function, t::A, k1, k2, ks::Vararg{Any, N})::Union{Nothing, A} where {N} where {A<:AbstractNode}
    tmp = get(f, t, k1)
    tmp === nothing ? (return nothing) : (tmp = get(f, tmp, k2))
    tmp === nothing && return nothing
    get(f, tmp, ks...)
    # Alternative 1
    # tmp = get(f, t, p)
    # tmp === nothing && return nothing
    # tmp = get(f, tmp, q)
    # tmp === nothing && return nothing
    # get(f, tmp, ps...)
end
# methods that supply () -> nothing; covers dispatches (1), (2) and (3)
# In order for this to work, one must opt out of get(...) on line 50.
# Most likely, that is a worthwhile thing, as the definition on line 50 is easily
# replicated using a direct call of the definition on line 55.
# The alternative is to simply supply () -> nothing to each call to get(f, t, ...)
_returnnothing() = nothing
Base.get(t::A, p) where {A<:AbstractNode} = get(_returnnothing, t, p)
Base.get(t::A, p, q) where {A<:AbstractNode} = get(_returnnothing, t, p, q)
Base.get(t::A, p, q, ps...) where {A<:AbstractNode} =
    get(_returnnothing, t, p, q, ps...)
# actually, the following should achieve the same result
# Base.get(t::A where {A<:AbstractNode{T, U}}, ps...) = get(() -> nothing, t, ps...)

# Convenient definition of haskey based on get
# haspath(t::AbstractNode, p) = get(_returnnothing, t, p) === nothing ? false : true
haspath(t::AbstractNode, p) = isa(get(_returnnothing, t, p), AbstractNode)
# better, likely more performant definition; benchmark to determine which to use,
# then change (2), (3). Even if the above is faster, consider get(...) !== nothing,
# when one can tolerate receiving something other than a AbstractNode.
# haspath(t::AbstractNode, p) = isa(get(_returnnothing, t, p), AbstractNode)
# haspath(t::AbstractNode, p, q) = isa(get(f, t, p, q), AbstractNode)
# haspath(t::AbstractNode, p, q, ps...) = isa(get(f, t, p, q, ps...), AbstractNode)
# Follow pattern: Varargs dispatch to underlying function where possible, rather
# than dispatch on (1), (2), (3) for each method. Single arg form is likely still useful.
haspath(t::AbstractNode, p, ps...) = isa(get(_returnnothing, t, p, ps...), AbstractNode)

# default should be a sub-type of AbstractNode
Base.get!(x::A, key, default) where {A<:AbstractNode} =
    get!(x.link, key, default)

# # f must be some function which returns a sub-type of AbstractNode
# Base.get!(f::Function, x::A where {A<:AbstractNode{T, U}}, key) where {T, U} = get!(f, x.link, key)
# Base.get!(f::Function, x::A where {A<:AbstractNode{T, U}}, k1, k2) where {T, U} =
#     get!(f, get!(f, x, k1), k2)
# # An interesting note: in the method signature, replacing ks... with Vararg{Any, N} where {N}
# # results in an approximate 1.177 speedup (15%) and 4% reduction in memory due to
# # specialization.
# Base.get!(f::Function, x::A where {A<:AbstractNode{T, U}}, k1, k2, ks::Vararg{Any, N}) where {N} where {T, U} =
#     get!(f, get!(f, get!(f, x, k1), k2), ks...)
# # An interesting note: get!.(() -> SimpleNode(), Ref(x), [1:100;])
# # does exactly what one would expect.

## 2021-10-20: Experiments with type stability
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

function Base.getindex(x::A, key)::A where {A<:AbstractNode}
    getindex(x.link, key)
end
function Base.getindex(x::A, k1, k2)::A where {A<:AbstractNode}
    getindex(getindex(x, k1), k2)
end
function Base.getindex(x::A, k1, k2, ks::Vararg{Any, N})::A where {N} where {A<:AbstractNode}
    getindex(getindex(getindex(x, k1), k2), ks...)
end
function Base.getindex(x::A, k1, k2, ks::Vararg{S, N})::A where {S, N} where {A<:AbstractNode}
    tmp = getindex(x, k1, k2)
    for k in ks
        tmp = getindex(tmp, k)
    end
    tmp
end
# function Base.getindex(x::A, ::Colon)::Vector{A} where {T, U, A<:AbstractNode{T, U}}
#     # convert(Vector{A}, collect(values(x)))
#     collect(values(x))
# end
# function getindex3(x::A, key) where {T, U, A<:AbstractNode{T, U}}
#     convert(A, getindex(x.link, key))
# end
# Base.getindex(x::AbstractNode, k1, k2) = getindex(getindex(x, k1), k2)
# Base.getindex(x::AbstractNode, k1, k2, ks::Vararg{Any, N}) where {N} = getindex(getindex(getindex(x, k1), k2), ks...)

Base.setindex!(x::AbstractNode, value, key) = (setindex!(x.link, value, key); x)
function Base.setindex!(x::AbstractNode, value, k1, ks...)
    tmp = x[k1, ks[1:end-1]...]
    tmp[ks[end]] = value
    x
end

Base.push!(x::AbstractNode, p::Pair) = setindex!(x, p.second, p.first)
Base.push!(x::AbstractNode, p::Pair, q::Pair) = push!(push!(x, p), q)
Base.push!(x::AbstractNode, p::Pair, q::Pair, r::Pair...) = push!(push!(push!(x, p), q), r...)

Base.sizehint!(x::AbstractNode, newsz) = sizehint!(x.link, newsz)

function Base.empty!(x::AbstractNode)
    empty!(x.link)
    empty!(x.val)
    return x
end

Base.pop!(x::AbstractNode, key) = pop!(x.link, key)
Base.pop!(x::AbstractNode, key, default) = pop!(x.link, key, default)
Base.pop!(x::AbstractNode) = pop!(x.link)
# function _pop!(x::A)::Pair{Any,A} where {T,U,A<:AbstractNode{T,U}}
#     pop!(x.link)
# end

Base.delete!(x::AbstractNode, key) = (delete!(x.link, key); x)

Base.isempty(x::AbstractNode) = isempty(x.link)

Base.filter!(pred, x::AbstractNode) = (filter!(pred, x.link); x)

function Base.merge(a::AbstractSimpleNode, bs::AbstractSimpleNode...)
    ltmp = merge(a.link, getproperty.(bs, :link)...)
    vtmp = Vector{promote_type(eltype(a.val), eltype.(getproperty.(bs, :val))...)}()
    isempty(a.val) || append!(vtmp, deepcopy(a.val))
    for b in bs
        isempty(b.val) || append!(vtmp, deepcopy(b.val))
    end
    SimpleNode(ltmp, vtmp)
end

function Base.merge!(a::AbstractSimpleNode, bs::AbstractSimpleNode...)
    merge!(a.link, getproperty.(bs, :link)...)
    for b in bs
        isempty(b.val) || append!(a.val, b.val)
    end
    a
end

# p. 528-529, 2021-10-15
function maxdepth(t::AbstractNode, C::Int=0)
    C̃ = C + 1
    isempty(t) && return C̃
    C̃ₘ = C̃
    for p in t
        C̃ₘ = max(C̃ₘ, maxdepth(p.second, C̃))
    end
    return C̃ₘ
end

function maxbreadth(t::AbstractNode)
    b::Int = length(t)
    isempty(t) && return b
    for p in t
        b = max(b, maxbreadth(p.second))
    end
    return b
end

# Most likely rename to `rlength`
function rlength(t::AbstractNode, C::Int=0)
    C̃ = C + 1
    isempty(t) && return C̃
    for p in t
        C̃ += rlength(p.second, 0)
    end
    return C̃
end

################
function _size!(t::AbstractNode, sz::Vector{Int}, C::Int)
    length(sz) < C && push!(sz, 0)
    @inbounds sz[C] += 1
    isempty(t) && return sz
    C̃ = C + 1
    for p in t
        _size!(p.second, sz, C̃)
    end
    return sz
end

Base.size(t::AbstractNode) = tuple(_size!(t, Int[0], 1)...)

function _sizeat(t::AbstractNode, N::Int, C::Int)::Int
    C̃ = C + 1
    C̃ == N && return length(t)
    s = 0
    for p in t
        s += _sizeat(p.second, N, C̃)
    end
    return s
end

Base.size(t::AbstractNode, d::Int) = _sizeat(t, d, 1)

################################################################


#### comparison operators
# function Base.isequal(a::A where {A<:AbstractNode{T₁, U₁}},
#                       b::B where {B<:AbstractNode{T₂, U₂}}) where {T₁, U₁} where {T₂, U₂}
#     a.link == b.link && a.val == b.val
# end
# Base.:(==)(a::A where {A<:AbstractNode{T₁, U₁}},
#            b::B where {B<:AbstractNode{T₂, U₂}}) where {T₁, U₁} where {T₂, U₂} =
#                Base.isequal(a, b)
## Identical definitions
# function Base.isequal(a::AbstractNode{T₁, U₁},
#                       b::AbstractNode{T₂, U₂}) where {T₁, U₁} where {T₂, U₂}
#     a.link == b.link && a.val == b.val
# end
# Base.:(==)(a::AbstractNode{T₁, U₁},
#            b::AbstractNode{T₂, U₂}) where {T₁, U₁} where {T₂, U₂} =
#                Base.isequal(a, b)
## Another identical definition
function Base.isequal(a::AbstractNode, b::AbstractNode)
    a === b && return true
    a.link == b.link && a.val == b.val
end
Base.:(==)(a::AbstractNode, b::AbstractNode) = Base.isequal(a, b)

SimpleNode(link::T) where {T<:AbstractDict} = SimpleNode(link, [])
SimpleNode(val::U) where {U<:AbstractVector} = SimpleNode(Dict(), val)
SimpleNode() = SimpleNode(Dict())

#### Outer constructors: ComplexNode
ComplexNode(link::T, val::U, spec::V) where {T<:AbstractDict, U<:AbstractVector, V<:AbstractDict} =
    ComplexNode(link, val, spec, [])
ComplexNode(link::T, val::U, sval::W) where {T<:AbstractDict, U<:AbstractVector, W<:AbstractVector} =
    ComplexNode(link, val, Dict(), sval)
ComplexNode(link::T, val::U) where {T<:AbstractDict, U<:AbstractVector} =
    ComplexNode(link, val, Dict())
ComplexNode(link::T) where {T<:AbstractDict} = ComplexNode(link, [])
ComplexNode(val::U) where {U<:AbstractVector} = ComplexNode(Dict(), val)
ComplexNode() = ComplexNode(Dict())

#### Conversion
Base.convert(::Type{T}, x::AbstractSimpleNode) where {T<:AbstractComplexNode} = T(x)
Base.convert(::Type{T}, x::AbstractComplexNode) where {T<:AbstractSimpleNode} = T(x)
Base.convert(::Type{T}, x::T) where {T<:AbstractNode} = x
# or, if needed:
Base.convert(::Type{T}, x::T) where {T<:AbstractSimpleNode} = x
Base.convert(::Type{T}, x::T) where {T<:AbstractComplexNode} = x
# necessary:
Base.convert(::Type{T}, x::AbstractSimpleNode) where {T<:AbstractSimpleNode} = T(x)
Base.convert(::Type{T}, x::AbstractComplexNode) where {T<:AbstractComplexNode} = T(x)

#### Promotion
Base.promote_rule(::Type{S₁} where {S₁<:AbstractSimpleNode{T₁, U₁}},
                  ::Type{S₂} where {S₂<:AbstractSimpleNode{T₂, U₂}}) where {T₁, U₁, T₂, U₂} =
                      SimpleNode{promote_type(T₁, T₂), promote_type(U₁, U₂)}
Base.promote_rule(::Type{C₁} where {C₁<:AbstractComplexNode{T₁, U₁, V₁, W₁}},
                  ::Type{C₂} where {C₂<:AbstractComplexNode{T₂, U₂, V₂, W₂}}) where {T₁, U₁, V₁, W₁, T₂, U₂, V₂, W₂} =
                      ComplexNode{promote_type(T₁, T₂), promote_type(U₁, U₂), promote_type(V₁, V₂), promote_type(W₁, W₂)}
Base.promote_rule(::Type{S} where {S<:AbstractSimpleNode{T₁, U₁}},
                  ::Type{C} where {C<:AbstractComplexNode{T₂, U₂, V, W}}) where {T₁, U₁, T₂, U₂, V, W} =
                      ComplexNode{promote_type(T₁, T₂), promote_type(U₁, U₂), V, W}
# An option to facilitate promotion is to call eltype on the types, using
# the result to specify the resultant types as:
# Dict{promote_type(eltype(T₁, T₂))}, Vector{promote_type(eltype(U₁, U₂))}
# This is a useful option for forcing the type to be Any rather than an ambiguous type
# such as Dict{K, Any} where K, or, Vector{T} where T.

#### Outer constructors, between Node types
SimpleNode(x::AbstractSimpleNode) = SimpleNode(deepcopy(x.link), copyto!(similar(x.val), x.val))
SimpleNode(x::AbstractComplexNode) = SimpleNode(deepcopy(x.link), copyto!(similar(x.val), x.val))
# or, technically always possible to convert to a SimpleNode, but leave inactive for now.
# SimpleNode(x::AbstractNode) = SimpleNode(deepcopy(x.link), copyto!(similar(x.val), x.val))

ComplexNode(x::AbstractComplexNode) = ComplexNode(deepcopy(x.link), copyto!(similar(x.val), x.val),
                                                  deepcopy(x.spec), copyto!(similar(x.sval), x.sval))
ComplexNode(x::AbstractSimpleNode) = ComplexNode(deepcopy(x.link), copyto!(similar(x.val), x.val))

#### More elaborate outer constructors in light of promotion rules
# For both Node types, the ternary operator can be replaced by deepcopy(convert(T₁, x.link))
# as a full copy is necessary in both cases.
# Originally, the first term in the arguments to the constructors was:
# T₂ === T₁ ? deepcopy(x.link) : deepcopy(convert(T₁, x.link))
# However, the ternary operator can be replaced by deepcopy(convert(T₁, x.link))
# as a full copy is necessary in both cases of the ternary
SimpleNode{T₁, U₁}(x::AbstractComplexNode{T₂, U₂, V, W}) where {T₁, U₁} where {T₂, U₂, V, W} =
    # SimpleNode(T₂ === T₁ ? deepcopy(x.link) : deepcopy(convert(T₁, x.link)),
    #            copyto!(similar(T₁, x.val), x.val))
    SimpleNode(deepcopy(convert(T₁, x.link)), copyto!(similar(x.val, U₁), x.val))

ComplexNode{T₁, U₁, V, W}(x::AbstractSimpleNode{T₂, U₂})  where {T₁, U₁, V, W} where {T₂, U₂} =
    # ComplexNode(T₂ === T₁ ? deepcopy(x.link) : deepcopy(convert(T₁, x.link)),
    #             copyto!(similar(T₁, x.val), x.val))
    ComplexNode(deepcopy(convert(T₁, x.link)), copyto!(similar(x.val, U₁), x.val))

#### More outer constructors which may be needed. Second definition is perhaps not as safe?
# Interestingly, the first definition seems to be lower performance by a factor of ≈ 2.18.
# For now, I will prefer the second definition as it fits the intention exactly, in addition
# to being the original definition.
SimpleNode{T, U}(x::AbstractSimpleNode) where {T, U} =
    # SimpleNode(deepcopy(convert(T, x.link)), copyto!(similar(x.val, eltype(U)), x.val))
    SimpleNode(deepcopy(convert(T, x.link)), copyto!(similar(U, axes(x.val)), x.val))

ComplexNode{T, U, V, W}(x::AbstractComplexNode) where {T, U, V, W} =
    # ComplexNode(deepcopy(convert(T, x.link)), copyto!(similar(x.val, eltype(U)), x.val),
    #             deepcopy(convert(V, x.spec)), copyto!(similar(x.sval, eltype(W)), x.sval))
    ComplexNode(deepcopy(convert(T, x.link)), copyto!(similar(U, axes(x.val)), x.val),
                deepcopy(convert(V, x.spec)), copyto!(similar(W, axes(x.sval)), x.sval))

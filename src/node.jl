#
# Date created: 2021-09-09
# Author: aradclif
#
#
############################################################################################
#### Draft of AbstractNode, SimpleNode, ComplexNode; see p. 406-414, 2021-09-09
################################################################
abstract type AbstractNode{T<:AbstractDict, U<:AbstractArray} end

# Simple Node
abstract type AbstractSimpleNode{T, U<:AbstractVector} <: AbstractNode{T, U} end

struct SimpleNode{T, U} <: AbstractSimpleNode{T, U}
    link::T
    val::U
    SimpleNode{T, U}(link, val) where {T, U} = new(link, val)
end
SimpleNode(link::T, val::U) where {T, U} = SimpleNode{T, U}(link, val)

# Complex Node
abstract type AbstractComplexNode{T, U<:AbstractVector, V<:AbstractDict, W<:AbstractVector} <: AbstractNode{T, U} end

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
Base.iterate(x::A where {A<:AbstractNode{T, U}}) where {T, U} = iterate(x.link)
Base.iterate(x::A where {A<:AbstractNode{T, U}}, state) where {T, U} = iterate(x.link, state)
# Base.iterate(x::A) where {A<:AbstractNode{T, U}} where {T, U} = iterate(x.link)
# Base.iterate(x::A, state) where {A<:AbstractNode{T, U}} where {T, U} = iterate(x.link, state)

Base.length(x::A where {A<:AbstractNode{T, U}}) where {T, U} = length(x.link)

Base.eltype(::A) where {A<:AbstractNode{T, U}} where {T, U} = Pair{Any, A}
#### Opt-ins: AbstractNode
Base.keys(x::A where {A<:AbstractNode{T, U}}) where {T, U} = keys(x.link)

Base.values(x::A where {A<:AbstractNode{T, U}}) where {T, U} = values(x.link)

Base.pairs(x::A where {A<:AbstractNode{T, U}}) where {T, U} = pairs(x.link)

# Base.get(x::A where {A<:AbstractNode{T, U}}, key, default) where {T, U} = get(x.link, key, default)

# f must be () -> nothing for this to work properly.
# Other functions which can be feasible must return nothing in order for
# the same behavior to be achieved.
Base.get(f::Function, x::A where {A<:AbstractNode{T, U}}, key) where {T, U} = get(f, x.link, key)
Base.get(f::Function, t::A where {A<:AbstractNode{T, U}}, p, q) where {T, U} =
    (tmp = get(f, t, p); tmp === nothing ? nothing : get(f, tmp, q))
function Base.get(f::Function, t::A where {A<:AbstractNode{T, U}}, p, q, ps::Vararg{Any, N}) where {N} where {T, U}
    tmp = get(f, t, p)
    tmp === nothing ? (return nothing) : (tmp = get(f, tmp, q))
    tmp === nothing && return nothing
    get(f, tmp, ps...)
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
Base.get(t::A where {A<:AbstractNode{T, U}}, p) where {T, U} = get(_returnnothing, t, p)
Base.get(t::A where {A<:AbstractNode{T, U}}, p, q) where {T, U} = get(_returnnothing, t, p, q)
Base.get(t::A where {A<:AbstractNode{T, U}}, p, q, ps...) where {T, U} =
    get(_returnnothing, t, p, q, ps...)
# actually, the following should achieve the same result
# Base.get(t::A where {A<:AbstractNode{T, U}}, ps...) = get(() -> nothing, t, ps...)

# Convenient definition of haskey based on get
# haspath(t::AbstractNode, p) = get(_returnnothing, t, p) === nothing ? false : true
haspath(t::AbstractNode, p) = isa(get(_returnnothing, t, p), AbstractNode)
# better, likely more performant definition; benchmark to determine which to use,
# then change (2), (3). Even if the above is faster, consider get(...) !== nothing,
# when one can tolerate receiving something other than a AbstractNode.
haspath(t::AbstractNode, p) = isa(get(_returnnothing, t, p), AbstractNode)
# haspath(t::AbstractNode, p, q) = isa(get(f, t, p, q), AbstractNode)
# haspath(t::AbstractNode, p, q, ps...) = isa(get(f, t, p, q, ps...), AbstractNode)
# Follow pattern: Varargs dispatch to underlying function where possible, rather
# than dispatch on (1), (2), (3) for each method. Single arg form is likely still useful.
haspath(t::AbstractNode, p, ps...) = isa(get(_returnnothing, t, p, ps...), AbstractNode)

# default should be a sub-type of AbstractNode
Base.get!(x::A where {A<:AbstractNode{T, U}}, key, default) where {T, U} =
    get!(x.link, key, default)

# f must be some function which returns a sub-type of AbstractNode
Base.get!(f::Function, x::A where {A<:AbstractNode{T, U}}, key) where {T, U} = get!(f, x.link, key)
Base.get!(f::Function, x::A where {A<:AbstractNode{T, U}}, k1, k2) where {T, U} =
    get!(f, get!(f, x, k1), k2)
# An interesting note: in the method signature, replacing ks... with Vararg{Any, N} where {N}
# results in an approximate 1.177 speedup (15%) and 4% reduction in memory due to
# specialization.
Base.get!(f::Function, x::A where {A<:AbstractNode{T, U}}, k1, k2, ks::Vararg{Any, N}) where {N} where {T, U} =
    get!(f, get!(f, get!(f, x, k1), k2), ks...)
# An interesting note: get!.(() -> SimpleNode(), Ref(x), [1:100;])
# does exactly what one would expect.

## 2021-10-20: Experiments with type stability
function get5!(f::Function, x::A, k1)::A where {T, U, A<:AbstractNode{T, U}}
    get!(f, x.link, k1)
end
function get5!(f::Function, x::A, k1, k2)::A where {T,U,A<:AbstractNode{T,U}}
    get5!(f, get5!(f, x, k1), k2)
end
function get5!(f::Function, x::A, k1, k2, ks::Vararg{S, N})::A where {S,N} where {T,U,A<:AbstractNode{T,U}}
    tmp = get5!(f, x, k1, k2)
    for k in ks
        tmp = get5!(f, tmp, k)
    end
    tmp
end
function get5!(f::Function, x::A, k1, k2, ks::Vararg{Any, N})::A where {N} where {T,U,A<:AbstractNode{T,U}}
    get5!(f, get5!(f, get5!(f, x, k1), k2), ks...)
end

# Result: Test on 10x20000
# 217.104ms for original, 147.63ms for Vararg handled by loop, 148.221ms for
# Vararg handled by recursion. It is not entirely clear which is better.
@benchmark grow!(gf, SimpleNode(), pni6, eachcol(qmat)) seconds=30
@benchmark grow3!(gf, SimpleNode(), pni6, eachcol(qmat)) seconds=30
@benchmark grow4!(gf, SimpleNode(), pni6, eachcol(qmat)) seconds=30
# Result: Test on 200x1000
# 866.170ms for original, 113.30ms for Vararg handled by loop, 725.592ms for
# Vararg handled by recursion. It is apparent that for deeper trees, loop is superior.
@benchmark grow!(gf, SimpleNode(), pni6, eachcol(qmat2)) seconds=30
@benchmark grow3!(gf, SimpleNode(), pni6, eachcol(qmat2)) seconds=30
@benchmark grow4!(gf, SimpleNode(), pni6, eachcol(qmat2)) seconds=30
# Result: Test on 20x10000
# 222.434ms for original, 127.251ms for Vararg handled by loop, 135.743ms for
# Vararg handled by recursion. At 2x the depth of the original, loop begins to pull away
@benchmark grow!(gf, SimpleNode(), pni2, eachcol(pmat)) seconds=30
@benchmark grow3!(gf, SimpleNode(), pni2, eachcol(pmat)) seconds=30
@benchmark grow4!(gf, SimpleNode(), pni2, eachcol(pmat)) seconds=30
# Result: Test on 4x50000
# 254.646ms for original, 198.503ms for Vararg handled by loop, 202.396ms for
# Vararg handled by recursion. Even at small depths, there is no effective difference.
@benchmark grow!(gf, SimpleNode(), pni, eachcol(mat)) seconds=30
@benchmark grow3!(gf, SimpleNode(), pni, eachcol(mat)) seconds=30
@benchmark grow4!(gf, SimpleNode(), pni, eachcol(mat)) seconds=30
# Result: Test on 50x4000
# 310.738 for original, 118.202 for Vararg handled by loop, 192.994 for
# Vararg handled by recursion. The loop over homogeneous keys shines again.
@benchmark grow!(gf, SimpleNode(), pni3, eachcol(qmat3)) seconds=30
@benchmark grow3!(gf, SimpleNode(), pni3, eachcol(qmat3)) seconds=30
@benchmark grow4!(gf, SimpleNode(), pni3, eachcol(qmat3)) seconds=30
@benchmark grow5!(gf, SimpleNode(), pni3, eachcol(qmat3)) seconds=30
# However, for heterogenous ks, the recursion seems to be superior.
# At depth 10, times are: 225ms, 207ms, 168ms
# At depth 200, the loop is still better, but it is 4x slower than with homogeneous ks.
# At depth 20, the slurp/splat recursion is fastest: abs(198.7 - 142.8)/198.7
# Or, about 30% faster.
# At depth 4, the recursion is fastest: abs(234 - 222)/234
##
# Essentially for most trees of reasonable depth, the slurp/splat recursion
# will be superior when heterogeneous keys are present. Only at trees of
# depth 200 is the loop over heterogeneous keys superior.
# Alas, it is undeniable: with homogeneous keys, the loop really shines.
## Multiple dispatch to the rescue!
# Use Vararg{S, N} where {S,N} to loop on homogeneous keys
# Use Vararg{Any,N} where {S,N} to slurp/splat recursion on heterogeneous keys
# Significantly, as the slurp/splat recursion chips away at ks,
# this leaves open the possibility of initiating a loop, since
# get!(f, get!(f, get!(f, x, k1), k2), ks...) has the chance to dispatch again.
# ---- true beauty at work.



Base.getindex(x::AbstractNode, key) = getindex(x.link, key)
## Solves type-instability at cost of 50% more time and 32 bytes allocation
# function getindex2(x::A, key)::A where {T, U, A<:AbstractNode{T, U}}
#     getindex(x.link, key)
# end
# function getindex3(x::A, key) where {T, U, A<:AbstractNode{T, U}}
#     convert(A, getindex(x.link, key))
# end
Base.getindex(x::AbstractNode, k1, k2) = getindex(getindex(x, k1), k2)
Base.getindex(x::AbstractNode, k1, k2, ks...) = getindex(getindex(getindex(x, k1), k2), ks...)

Base.setindex!(x::AbstractNode, value, key) = (setindex!(x.link, value, key); x)

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

Base.delete!(x::AbstractNode, key) = (delete!(x.link, key); x)

Base.isempty(x::AbstractNode) = isempty(x.link)

Base.filter!(pred, x::AbstractNode) = (filter!(pred, x.link); x)

# p. 485-487, 2021-10-04
# No reason to handle 3-arg case uniquely, as merge/merge! both direct to the same function.
# Hence, it is more efficient to just direct to the varargs.
# function Base.merge(a::AbstractSimpleNode, b::AbstractSimpleNode)
#     ltmp = merge(a.link, b.link)
#     vtmp = Vector{promote_type(eltype(a.val), eltype(b.val))}()
#     isempty(a.val) || append!(vtmp, deepcopy(a.val))
#     isempty(b.val) || append!(vtmp, deepcopy(b.val))
#     AbstractSimpleNode(ltmp, vtmp)
# end
# function Base.merge(a::AbstractSimpleNode, b::AbstractSimpleNode, c::AbstractSimpleNode)
#     ltmp = merge(a.link, b.link, c.link)
#     vtmp = Vector{promote_type(eltype(a.val), eltype(b.val), eltype(c.val))}()
#     isempty(a.val) || append!(vtmp, deepcopy(a.val))
#     isempty(b.val) || append!(vtmp, deepcopy(b.val))
#     isempty(c.val) || append!(vtmp, deepcopy(c.val))
#     AbstractSimpleNode(ltmp, vtmp)
# end
# function Base.merge(a::AbstractSimpleNode, b::AbstractSimpleNode,
#                     c::AbstractSimpleNode, ds::AbstractSimpleNode...)
#     ltmp = merge(a.link, b.link, c.link, getproperty.(ds, :link)...)
#     vtmp = Vector{promote_type(eltype(a.val), eltype(b.val), eltype(c.val),
#                                eltype.(getproperty.(ds, :val))...)}()
#     isempty(a.val) || append!(vtmp, deepcopy(a.val))
#     isempty(b.val) || append!(vtmp, deepcopy(b.val))
#     isempty(c.val) || append!(vtmp, deepcopy(c.val))
#     for d in ds
#         isempty(d.val) || append!(vtmp, deepcopy(d.val))
#     end
#     AbstractSimpleNode(ltmp, vtmp)
# end

# function Base.merge!(a::AbstractSimpleNode, b::AbstractSimpleNode)
#     merge!(a.link, b.link)
#     isempty(b.val) || append!(a.val, b.val)
#     a
# end
# function Base.merge!(a::AbstractSimpleNode, b::AbstractSimpleNode, c::AbstractSimpleNode)
#     merge!(a.link, b.link, c.link)
#     isempty(b.val) || append!(a.val, b.val)
#     isempty(c.val) || append!(a.val, c.val)
#     a
# end
# Base.merge!(a::AbstractSimpleNode, b::AbstractSimpleNode,
#             c::AbstractSimpleNode, ds::AbstractSimpleNode...) =
#     merge!(merge!(merge!(a, b), c), ds...)

function Base.merge(a::AbstractSimpleNode, bs::AbstractSimpleNode...)
    ltmp = merge(a.link, getproperty.(bs, :link)...)
    vtmp = Vector{promote_type(eltype(a.val), eltype.(getproperty.(bs, :val))...)}()
    isempty(a.val) || append!(vtmp, deepcopy(a.val))
    for b in bs
        isempty(b.val) || append!(vtmp, deepcopy(b.val))
    end
    AbstractSimpleNode(ltmp, vtmp)
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


#### Opt-ins: AbstractComplexNode
function Base.isequal(a::AbstractComplexNode, b::AbstractComplexNode)
    a.link == b.link && a.val == b.val && a.spec == b.spec && a.sval == b.sval
end
Base.:(==)(a::AbstractComplexNode, b::AbstractComplexNode) = Base.isequal(a, b)

# bi-directional: p. 460-463, 2021-10-06
function isbidirectional(V₁::AbstractComplexNode, V₂::AbstractComplexNode, k)
    (K => V₂) ∈ V₁.flink && (K => V₁) ∈ V₂.blink && return true
    (K => V₂) ∈ V₁.blink && (K => V₁) ∈ V₂.flink && return true
    false
end

# this should actually be forwardget (abbrev. fget), and a backwardget should exist
function bget!(f::Function, V₁::AbstractComplexNode, k)
    V₂ = get!(f, V₁, k)
    setindex!(V₂.blink, V₁, k)
    V₂
end
bget!(f::Function, V₁::AbstractComplexNode, k1, k2) = bget!(f, bget!(f, V₁, k1), k1)
bget!(f::Function, V₁::AbstractComplexNode, k1, k2, ks...) =
    bget!(bget!(f, bget!(f, V₁, k1), k1), ks...)

function bget(f::Function, V₁::AbstractComplexNode, k)
    V₂ = get(f, V₁, k)
    # V₂ === nothing && return nothing
    # isbidirectional(V₁, V₂, k) || return nothing
    V₂ === nothing || !isbidirectional(V₁, V₂, k) && return nothing
    V₂
end
function bget(f::Function, V₁::AbstractComplexNode, k1, k2)
    V₂ = bget(f, V₁, k)
    V₂ === nothing && return nothing
    V₃ = bget(f, V₂, k)
end

function bget(f::Function, V₁::AbstractComplexNode, k1, k2, ks...)
    V₂ = bget(f, V₁, k)
    V₂ === nothing && return nothing
    V₃ = bget(f, V₂, k)
    V₃ === nothing && return nothing
    bget(f, V₃, ks...)
end

function hasbipath(V₁::AbstractComplexNode, k)
    isa(bget(_returnnothing, V₁, k), AbstractComplexNode)
end
function hasbipath(V₁::AbstractComplexNode, k1, k2)
    isa(bget(_returnnothing, V₁, k1, k2), AbstractComplexNode)
end
function hasbipath(V₁::AbstractComplexNode, k1, k2, ks...)
    isa(bget(_returnnothing, V₁, k1, k2, ks...), AbstractComplexNode)
end

#### END Opt-ins: AbstractComplexNode

#### Constructors, Conversion and Promotion -- initial experiment
# ComplexNode(x::AbstractSimpleNode) = ComplexNode(copy(x.link), copyto!(similar(x.val), x.val))

# Base.convert(::Type{A}, x::A) where {A<:AbstractNode} = x
# Base.convert(::Type{A} where {A<:AbstractSimpleNode}, x::AbstractComplexNode) =
#     A(copy(x.link), copyto!(similar(x.val), x.val))

#### Outer constructors: SimpleNode
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

################ Tests: will be moved to package tests.jl when ready
using BenchmarkTools
using Test
# Constructors
# Conversion -- SimpleNode-SimpleNode
x = SimpleNode(Dict(), [])
@test convert(SimpleNode{Dict{Any, Any}, Vector{Any}}, x) === x
x = SimpleNode(Dict{String, Any}(), Int[])
y = convert(SimpleNode{Dict{Any, Any}, Vector{Any}}, x)
@test x == y
@test x !== y
@test typeof(x) !== typeof(y)
x = SimpleNode(Dict{String, Int}(), Int[])
y = convert(SimpleNode{Dict{String, Any}, Vector{Float64}}, x)
@test x == y
@test x !== y
@test typeof(x) !== typeof(y)
x = SimpleNode(Dict{Any, Int}(), Int[])
y = convert(SimpleNode{Dict{String, Any}, Vector{Float64}}, x)
@test x == y
@test x !== y
@test typeof(x) !== typeof(y)
# Conversion -- ComplexNode-ComplexNode
x = ComplexNode(Dict(), [], Dict(), [])
@test convert(ComplexNode{Dict{Any, Any}, Vector{Any}, Dict{Any, Any}, Vector{Any}}, x) === x
x = ComplexNode(Dict{String, Any}(), Int[], Dict(), [])
y = convert(ComplexNode{Dict{Any, Any}, Vector{Any}, Dict{Any, Any}, Vector{Any}}, x)
@test x == y
@test x !== y
@test typeof(x) !== typeof(y)
# Conversion -- SimpleNode-ComplexNode
x = SimpleNode(Dict{String, Any}(), Int[])
y = convert(ComplexNode{Dict{Any, Any}, Vector{Any}, Dict{Any, Any}, Vector{Any}}, x)
@test x == y
@test x !== y
@test typeof(x) !== typeof(y)
# Conversion -- ComplexNode-SimpleNode
x = ComplexNode(Dict{String, Any}(), Int[], Dict(), [])
y = convert(SimpleNode{Dict{Any, Any}, Vector{Any}}, x)
@test x == y
@test x !== y
@test typeof(x) !== typeof(y)
# Promotion -- SimpleNode-SimpleNode
x = SimpleNode(Dict{String, Int}(), Int[])
y = SimpleNode(Dict{Int, Float64}(), Float64[])
promote(x, y)
x = SimpleNode{Dict{String, Any}, Vector{Any}}
y = SimpleNode{Dict{Int, Any}, Vector{Any}}
z = SimpleNode{Dict{Any, Any}, Vector{Any}}
@test promote_type(x, z) !== z # previously === returned false
@test promote_type(x, y) !== z # previously === returned false
a = SimpleNode{Dict{Any, Any}, Vector{Int}}
b = SimpleNode{Dict{Any, Any}, Vector{String}}
c = z
@test promote_type(a, b) !== c # previously === returned false
# Promotion -- ComplexNode-ComplexNode
xc = ComplexNode(Dict{String, Int}(), Int[], Dict{Float64, Any}(), Int[])
yc = ComplexNode(Dict{String, Int}(), Float64[], Dict{Float64, Any}(), Real[])
promote(xc, yc)
xc = ComplexNode{Dict{String, Any}, Vector{Any}, Dict{Any, Any}, Vector{Any}}
yc = ComplexNode{Dict{Int, Any}, Vector{Any}, Dict{Any, Any}, Vector{Any}}
zc = ComplexNode{Dict{Any, Any}, Vector{Any}, Dict{Any, Any}, Vector{Any}}
@test promote_type(xc, yc) !== zc # previously === returned false -- and used x,y,z
ac = ComplexNode{Dict{Any, Any}, Vector{Int}, Dict{Any, Any}, Vector{Any}}
bc = ComplexNode{Dict{Any, Any}, Vector{String}, Dict{Any, Any}, Vector{Any}}
cc = zc
@test promote_type(ac, bc) !== cc # previously === returned false
# Promotion -- ComplexNode-SimpleNode
@test promote_type(xc, y) <: AbstractComplexNode
@test promote_type(xc, y) !== zc # previously === returned false
@test promote_type(xc, x) === xc
@test promote_type(ac, a) === ac
@test promote_type(ac, b) <: AbstractComplexNode
@test promote_type(ac, b) !== zc # previously === returned false
@test promote_type(bc, b) === bc
############################################################################################
#### Methods of get_:(:op):(:field)!
# One has the signatures
# (1) -- (f::Function, t::AbstractNode, v, p)
# (2) -- (f::Function, t::AbstractNode, v, p, q)
# (3) -- (f::Function, t::AbstractNode, v, p, q, ps...)
# One must cover fields
# :val -- AbstractNode
# :spec -- AbstractComplexNode
# :sval -- AbstractComplexNode
# 3 function definitions per (:op, :field), 1 :field for AbstractNode, 2 for AbstractComplexNode
# yields 9 function definitions for just one :op. Hence, the metaprogramming.
# One should then provide function prototypes to make docstrings visible.
####
# # original version
# function get_pushval!(f::Function, t::AbstractNode, v, p)
#     tmp = get!(f, t, p)
#     push!(tmp.val, v)
#     tmp
# end
# function get_pushval!(f::Function, t::AbstractNode, v, p, q)
#     tmp = get!(f, t, p, q)
#     push!(tmp.val, v)
#     tmp
# end
# function get_pushval!(f::Function, t::AbstractNode, v, p, q, ps...)
#     tmp = get!(f, t, p, q, ps...)
#     push!(tmp.val, v)
#     tmp
# end
# An alternate way: code generation -- needs another look, but seems to work.
# for op = (:push!, :append!, :empty!), field = (:val)
#     fname = Symbol(:get_, field, op)
#     eval(:(function $fname(f::Function, t::AbstractNode, v, p)
#                tmp = get!(f, t, p)
#                $op(tmp.$field, v)
#                tmp
#            end))
#     eval(:(function $fname(f::Function, t::AbstractNode, v, p, q)
#                tmp = get!(f, t, p, q)
#                $op(tmp.$field, v)
#                tmp
#            end))
#     eval(:(function $fname(f::Function, t::AbstractNode, v, p, q, ps...)
#                tmp = get!(f, t, p, q, ps...)
#                $op(tmp.$field, v)
#                tmp
#            end))
# end
# macro getfieldop1arg(field, op)
#     :(function $(Symbol(:get_, field, op))(f::Function, t::AbstractNode, v, p)
#           tmp = get!(f, t, p)
#           $op(tmp.$field, v)
#           tmp
#       end)
# end
# macro getfieldop2arg(field, op)
#     :(function $(Symbol(:get_, field, op))(f::Function, t::AbstractNode, v, p, q)
#           tmp = get!(f, t, p, q)
#           $op(tmp.$field, v)
#           tmp
#       end)
# end
# macro getfieldopvararg(field, op)
#     :(function $(Symbol(:get_, field, op))(f::Function, t::AbstractNode, v, p, q, ps...)
#           tmp = get!(f, t, p, q, ps...)
#           $op(tmp.$field, v)
#           tmp
#       end)
# end
# # Does not work as field op are consumed directly as symbols
# for field = (:val,), op = (:empty!,)
#     @getfieldop1arg field op
#     @getfieldop2arg field op
#     @getfieldopvararg field op
# end


# Best approach: just use dispatch on get! to limit number of possible dispatches
function get_valpush!(f::Function, t::AbstractNode, v, p)
    tmp = get!(f, t, p)
    push!(tmp.val, v)
    tmp
end
function get_valpush!(f::Function, t::AbstractNode, v, p, ps::Vararg{Any, N}) where {N}
    tmp = get!(f, t, p, ps...)
    push!(tmp.val, v)
    tmp
end
function get_specpush!(f::Function, t::AbstractNode, v::Pair, p)
    tmp = get!(f, t, p)
    push!(tmp.spec, v)
    tmp
end
function get_specpush!(f::Function, t::AbstractNode, v::Pair, p, ps::Vararg{Any, N}) where {N}
    tmp = get!(f, t, p, ps...)
    push!(tmp.spec, v)
    tmp
end
function get_svalpush!(f::Function, t::AbstractNode, v, p)
    tmp = get!(f, t, p)
    push!(tmp.sval, v)
    tmp
end
function get_svalpush!(f::Function, t::AbstractNode, v, p, ps::Vararg{Any, N}) where {N}
    tmp = get!(f, t, p, ps...)
    push!(tmp.sval, v)
    tmp
end

#### Methods of grow_(:field)! : see p. 443-446, 451-455, 2021-09-14/15
# grow! -> get! only has 2 signatures:
# (f::Function, t::AbstractNode, p::AbstractPathKeys)
# (f::Function, t::AbstractNode, p::AbstractPathKeys, itr)
# grow_val! -> get_pushval! only has 2 signatures:
# (f::Function, t::AbstractNode, v::AbstractPathKeys, p::AbstractPathKeys)
# (f::Function, t::AbstractNode, v::AbstractPathKeys, p::AbstractPathKeys, itr)
# grow_spec! -> get_pushspec! has 2 signatures, same as grow_val!
# grow_sval! -> get_pushsval! has 2 signatures, same as grow_val!

function grow!(f::Function, t::AbstractNode, p::AbstractPathKeys)
    for x in p
        get!(f, t, x...)
    end
    return t
end
function grow!(f::Function, t::AbstractNode, p::AbstractPathKeys, itr)
    x = Vector{Any}(undef, p.N)
    for item in itr
        p(x, item)
        get!(f, t, x...)
    end
    return t
end

# Alternately, just remove the type on v

# function _valgrow!(f::Function, t::AbstractNode, v, p::AbstractPathKeys)
#     for x in p
#         get_valpush!(f, t, v(x), x...)
#     end
#     return t
# end
# function _valgrow!(f::Function, t::AbstractNode, v, p::AbstractPathKeys, itr)
#     x = Vector{Any}(undef, p.N)
#     for item in itr
#         p(x, item)
#         get_valpush!(f, t, v(item), x...)
#     end
#     return t
# end
# valgrow!(f::Function, t::AbstractNode, v::AbstractPathKeys, p::AbstractPathKeys) =
#     _valgrow!(f, t, v, p)
# valgrow!(f::Function, t::AbstractNode, v::AbstractPathKey, p::AbstractPathKeys) =
#     _valgrow!(f, t, v, p)
# valgrow!(f::Function, t::AbstractNode, v::AbstractPathKeys, p::AbstractPathKeys, itr) =
#     _valgrow!(f, t, v, p, itr)
# valgrow!(f::Function, t::AbstractNode, v::AbstractPathKey, p::AbstractPathKeys, itr) =
#     _valgrow!(f, t, v, p, itr)

function valgrow!(f::Function, t::AbstractNode, v::AbstractPathKeys, p::AbstractPathKeys)
    for x in p
        get_valpush!(f, t, v(x), x...)
    end
    return t
end
function valgrow!(f::Function, t::AbstractNode, v::AbstractPathKeys, p::AbstractPathKeys, itr)
    x = Vector{Any}(undef, p.N)
    for item in itr
        p(x, item)
        get_valpush!(f, t, v(item), x...)
    end
    return t
end
function valgrow!(f::Function, t::AbstractNode, v::AbstractPathKey, p::AbstractPathKeys)
    for x in p
        get_valpush!(f, t, v(x), x...)
    end
    return t
end
function valgrow!(f::Function, t::AbstractNode, v::AbstractPathKey, p::AbstractPathKeys, itr)
    x = Vector{Any}(undef, p.N)
    for item in itr
        p(x, item)
        get_valpush!(f, t, v(item), x...)
    end
    return t
end
# Growth from non-flat sources - p. 475, 2021-09-22
function valgrow!(f::Function, t::AbstractNode, v::AbstractPathKeys, p::AbstractPathKeys,
                  vitr, pitr)
    x = Vector{Any}(undef, p.N)
    for (a, b) in zip(vitr, pitr)
        p(x, b)
        get_valpush!(f, t, v(a), x...)
    end
    return t
end
function valgrow!(f::Function, t::AbstractNode, v::AbstractPathKey, p::AbstractPathKeys,
                  vitr, pitr)
    x = Vector{Any}(undef, p.N)
    for (a, b) in zip(vitr, pitr)
        p(x, b)
        get_valpush!(f, t, v(a), x...)
    end
    return t
end
# Possible parallel growth method
function tvalgrow!(f::Function, t::AbstractNode, v::AbstractPathKeys, p::AbstractPathKeys,
                   itrsource::AbstractDict)
    @sync for p in t
        Threads.@spawn valgrow!(f, p.second, v, p, eachcol(itrsource[p.first]))
        # Alternative:
        # let itr = eachcol(itrsource[p.first])
        #     Threads.@spawn valgrow!(f, p.second, v, p, itr)
        # end
    end
    t
end

function specgrow!(f::Function, t::AbstractNode, v::AbstractPathKeys, p::AbstractPathKeys)
    for x in p
        get_specpush!(f, t, v(x), x...)
    end
    return t
end
function specgrow!(f::Function, t::AbstractNode, v::AbstractPathKeys, p::AbstractPathKeys, itr)
    x = Vector{Any}(undef, p.N)
    for item in itr
        p(x, item)
        get_specpush!(f, t, v(item), x...)
    end
    return t
end
function specgrow!(f::Function, t::AbstractNode, v::AbstractPathKey, p::AbstractPathKeys)
    for x in p
        get_specpush!(f, t, v(x), x...)
    end
    return t
end
function specgrow!(f::Function, t::AbstractNode, v::AbstractPathKey, p::AbstractPathKeys, itr)
    x = Vector{Any}(undef, p.N)
    for item in itr
        p(x, item)
        get_specpush!(f, t, v(item), x...)
    end
    return t
end

function svalgrow!(f::Function, t::AbstractNode, v::AbstractPathKeys, p::AbstractPathKeys)
    for x in p
        get_svalpush!(f, t, v(x), x...)
    end
    return t
end
function svalgrow!(f::Function, t::AbstractNode, v::AbstractPathKeys, p::AbstractPathKeys, itr)
    x = Vector{Any}(undef, p.N)
    for item in itr
        p(x, item)
        get_svalpush!(f, t, v(item), x...)
    end
    return t
end
function svalgrow!(f::Function, t::AbstractNode, v::AbstractPathKey, p::AbstractPathKeys)
    for x in p
        get_svalpush!(f, t, v(x), x...)
    end
    return t
end
function svalgrow!(f::Function, t::AbstractNode, v::AbstractPathKey, p::AbstractPathKeys, itr)
    x = Vector{Any}(undef, p.N)
    for item in itr
        p(x, item)
        get_svalpush!(f, t, v(item), x...)
    end
    return t
end

############################################################################################
#### AbstractPathKey, AbstractPathKeys: p. 447-449, 456-458, 2021-09-15
# # Type hierarchy of functors
# abstract type AbstractPathKey{T} end
# # abstract type AbstractPathKey1{Function} end

# abstract type AbstractIndexedPathKey{T<:Function} <: AbstractPathKey{T} end

# # struct LinearPathKey{T} <: AbstractIndexedPathKey{T}
# #     f::T
# #     i::Int
# #     # LinearPathKey{T, U}(f, i) where {T, U} = new(f, i)
# # end
# # # LinearPathKey(f::T, i::U) where {T, U} = LinearPathKey{T, U}(f, i)

# # struct MultipleLinearPathKey{T} <: AbstractIndexedPathKey{T}
# #     f::T
# #     i::Vector{Int}
# # end

# struct IndexedPathKey{T, U} <: AbstractIndexedPathKey{T}
#     f::T
#     i::U
#     # IndexedPathKey{T, U}(f, i) where {T, U} = new(f, i)
# end
# # IndexedPathKey(f::T, i::U) where {T, U} = IndexedPathKey{T, U}(f, i)

# #### Outer constructors for IndexedPathKey
# IndexedPathKey(i::Int) = IndexedPathKey(identity, i)
# IndexedPathKey(i::Vector{Int}) = IndexedPathKey(identity, i)

# #### functor: default behavior for all IndexedPathKey
# function (p::AbstractIndexedPathKey)(A)
#     p.f(getindex(A, p.i))
# end

# # ipk1 = IndexedPathKey(identity, 2)
# # ipk2 = IndexedPathKey(x -> "2" .* x, [1, 2])
# # @code_warntype ipk1(a)
# # vi = [ipk1, ipk2]
# # pi = PathKeys(vi, 2)
# # @code_warntype pi(y, a)
# # #### Attempt 1
# # abstract type AbstractIndexedPathKey1{T} <: AbstractPathKey1{T} end
# # struct LinearPathKey1{T} <: AbstractIndexedPathKey1{T}
# #     f::T
# #     i::Int
# # end
# # LinearPathKey1(i::Int) = LinearPathKey1(identity, i)
# # (p::LinearPathKey1)(A) = p.f(getindex(A, p.i))
# # ii1 = LinearPathKey1(identity, 1)
# # ii3 = LinearPathKey1(identity, 3)
# # ii4 = LinearPathKey1(csv, 4)
# # ii2 = LinearPathKey1(x -> 25, 2)
# # ii1 = LinearPathKey1(x -> 1, 2)
# # ii3 = LinearPathKey1(x -> 2, 2)
# # ii4 = LinearPathKey1(x -> 2, 2)
# # @code_warntype ii2(a)

# # abstract type AbstractPathKeys4{U<:Tuple{Vararg{T, N} where {T<:AbstractPathKey1}} where {N}} end
# # struct PathKeys4{U} <: AbstractPathKeys4{U}
# #     ftrs::U
# #     N::Int
# # end

# # p4 = PathKeys4((ii1, ii3, ii4, ii2), 4)
# # @code_warntype p4(a)
# # @code_warntype p4(y, a)
# # @benchmark p4(y, a)
# # @benchmark ii1(a)

# # #### Attempt 2
# # abstract type AbstractPathKeys7{U<:Vector{T} where {T<:LinearPathKey1}} end
# # struct PathKeys7{U} <: AbstractPathKeys7{U}
# #     ftrs::U
# #     N::Int
# # end

# # p7 = PathKeys7([ii1, ii3, ii4, ii2], 4)
# # @code_warntype p7(a)
# # @code_warntype p7(y, a)
# # @benchmark p7(y, a)
# # @benchmark ii1(a)
# # @code_warntype ii2(a)

# # abstract type AbstractPathKeys{T<:NTuple{<:M, <:AbstractPathKey}} end
# # abstract type AbstractPathKeys{U<:Tuple{Vararg{T, N} where {T<:AbstractPathKey}} where {N}} end
# # abstract type AbstractPathKeys1{U<:Tuple{Vararg{T, N} where N where {T<:AbstractPathKey}}} end
# # abstract type AbstractPathKeys2{U<:Tuple{Vararg{T, N} where N} where {T<:AbstractPathKey}} end
# # abstract type AbstractPathKeys3{U<:Tuple{Vararg{T, N}} where N where {T<:AbstractPathKey}} end
# # Tuple{Vararg{T, N} where {T<:AbstractArray}} where N
# # Tuple{Vararg{T, N} where N where {T<:AbstractArray}}
# # isconcretetype(ans) # false
# # isabstracttype(ans) # false
# # ans{Vector, Vector}
# # Tuple{Vararg{T, N} where N} where {T<:AbstractArray}
# # ans{Vector{Int}} # Tuple{Vararg{Vector{Int64}}}
# # isconcretetype(ans) # false
# # isabstracttype(ans) # false
# # Tuple{Vararg{T, N}} where N where {T<:AbstractArray}
# # ans{Vector{Int}} # Tuple{Vararg{Vector{Int64}, N}} where N :same as NTuple{N, Vector{Int}} where N
# # ans{Vector{Int}, 2} # Tuple{Vector{Int64}, Vector{Int64}}
# # isconcretetype(ans) # true
# # const U2 = Tuple{Vararg{T, N}} where {N} where {T<:AbstractPathKey}
# # const U3 = Tuple{Vararg{T, N} where {T<:AbstractPathKey, N}}
# # const U4 = Tuple{Vararg{T, N} where {T<:AbstractPathKey}} where {N}
# # const U5 = Tuple{Vararg{T} where {T<:AbstractPathKey}}
# # abstract type AbstractPathKeys3{U<:Tuple{Vararg{T, N} where {T<:AbstractPathKey, N}}} end
# # struct PathKeys3{U} <: AbstractPathKeys3{U}
# #     ftrs::U
# #     N::Int
# #     # PathKeys{T}(ftrs, N) where {T} = new(ftrs, N)
# #     # function PathKeys{U}(ftrs, N) where {U}
# #     #     new(copyto!(similar(ftrs), ftrs), N)
# #     # end
# # end
# # abstract type AbstractStaticPathKeys{U<:SVector{S, T} where {S, T<:AbstractPathKey}} end
# # struct StaticPathKeys{U} <: AbstractStaticPathKeys{U}
# #     ftrs::U
# #     N::Int
# # end

# # Revised abstract type for AbstractPathKeys
# abstract type AbstractPathKeys{U<:Vector{T} where {T<:AbstractPathKey}} end
# struct PathKeys{U} <: AbstractPathKeys{U}
#     ftrs::U
#     N::Int
#     # PathKeys{T}(ftrs, N) where {T} = new(ftrs, N)
#     # function PathKeys{U}(ftrs, N) where {U}
#     #     new(copyto!(similar(ftrs), ftrs), N)
#     # end
#     # function PathKeys{U}(ftrs, N) where {U}
#     #     let ftrs = copyto!(similar(ftrs), ftrs), N = N
#     #         new(ftrs, N)
#     #     end
#     # end
# end
# # PathKeys(ftrs::T, N) where {T} = PathKeys{T}(ftrs, N)


# #### Interface for AbstractPathKeys
# Base.length(p::AbstractPathKeys) = p.N
# Base.size(p::AbstractPathKeys) = (p.N,)

# Base.iterate(p::AbstractPathKeys, state=1) =
#     state > p.N ? nothing : (p.ftrs[state], state + 1)
# Base.eltype(p::AbstractPathKeys) = eltype(p.ftrs)

# Base.IndexStyle(::Type{<:AbstractPathKeys}) = IndexLinear()
# Base.getindex(p::AbstractPathKeys, i::Int) = getindex(p.ftrs, i)
# Base.getindex(p::AbstractPathKeys, I) = [p[i] for i in I]
# Base.getindex(p::AbstractPathKeys, ::Colon) = p.ftrs
# Base.firstindex(p::AbstractPathKeys) = 1
# Base.lastindex(p::AbstractPathKeys) = p.N
# Base.setindex!(p::AbstractPathKeys, v, i) = setindex!(p.ftrs, v, i)

# function Base.isequal(p1::AbstractPathKeys, p2::AbstractPathKeys)
#     p1 === p2 && return true
#     length(p1) == length(p2) || return false
#     for n = 1:length(p1)
#         p1[n] == p2[n] || return false
#     end
#     return true
# end
# Base.:(==)(p1::AbstractPathKeys, p2::AbstractPathKeys) = isequal(p1, p2)

# #### Outer constructors: PathKeys
# PathKeys(ftrs) = PathKeys(ftrs, length(ftrs))

# #### functor: PathKeys
# # both internal functions yield no performance gain
# # function _pcall(p::AbstractPathKey, A)
# #     p(A)
# # end
# # function _pfunctor(fs::Vector{T} where {T<:AbstractPathKey}, x, A, N)
# #     n = 1
# #     while n ≤ N
# #         x[n] = fs[n](A)
# #         n += 1
# #     end
# #     return x
# # end
# function (p::AbstractPathKeys)(x, A)
#     n = 1
#     N = p.N
#     # N = length(p)
#     fs = p.ftrs
#     while n ≤ N#p.N
#         # x[n] = p.ftrs[n](A)
#         # x[n] = _pcall(p.ftrs[n], A) # still type-unstable
#         # Realistic options: use local fs, or, access via getindex
#         x[n] = fs[n](A)
#         # x[n] = p[n](A)
#         n += 1
#     end
#     return x
#     # _pfunctor(p.ftrs, x, A, p.N)
#     # # Same, despite performance annotations
#     # fs = p.ftrs
#     # @simd for n = 1:p.N
#     #     @inbounds x[n] = fs[n](A)
#     # end
#     # return x
#     # Meta 1 -- much slower
#     # fs = p.ftrs
#     # N = p.N
#     # eval(quote
#     #          Base.Cartesian.@nexprs $N i -> x[i] = fs[i](A)
#     #      end)
#     # return x
#     # Meta 2 -- @generated does not work
#     # N = p.N
#     # fs = p.ftrs
#     # ex = quote
#     #     Base.Cartesian.@nexprs $N i -> x[i] = fs[i](A)
#     # end
#     # return :($ex; x)
# end
# function (p::AbstractPathKeys)(A)
#     x = Vector{Any}(undef, p.N)
#     p(x, A)
# end

# struct PathKeysLet{U} <: AbstractPathKeys{U}
#     ftrs::U
#     N::Int
#     function PathKeysLet{U}(ftrs, N) where {U}
#         let fts = ftrs, M = N
#             new(fts, M)
#         end
#     end
# end
# PathKeysLet(ftrs::U, N) where {U} = PathKeys{U}(ftrs, N)

# abstract type AbstractStaticPathKeys{U<:SVector{S, T} where {S, T<:AbstractPathKey}} end
# struct StaticPathKeys{U} <: AbstractStaticPathKeys{U}
#     ftrs::U
#     N::Int
# end

# struct PathKeys1{U} <: AbstractPathKeys1{U}
#     ftrs::U
#     N::Int
# end
# struct PathKeys2{U} <: AbstractPathKeys2{U}
#     ftrs::U
#     N::Int
# end
# struct PathKeys3{U} <: AbstractPathKeys3{U}
#     ftrs::U
#     N::Int
# end

# mutable struct MPathKeys{U} <: AbstractPathKeys{U}
#     ftrs::U
#     N::Int
# end

# p = PathKeys((i1, i3), 2)
# p1 = PathKeys1((i1, i3), 2)
# p2 = PathKeys2((i1, i3), 2)
# p3 = PathKeys3((i1, i3), 2)
# p = PathKeys((i1, i3, i4), 3)
# p1 = PathKeys1((i1, i3, i4), 3)
# p2 = PathKeys2((i1, i3, i4), 3)
# p3 = PathKeys3((i1, i3, i4), 3)
# PathKeys((1, "this"), 2)
# PathKeys1((1, "this"), 2)

# after revision to type of ftrs: Vector
# csv(x) = x * ".csv"
# a = ["a", "b", "c", "d"]
# y = Vector{Any}(undef, 4)
# i1 = LinearPathKey(1)
# i3 = LinearPathKey(3)
# i2 = LinearPathKey(x -> 25, 2)
# i4 = LinearPathKey(csv, 4)
# lpks = [i1, i2, i3, i4]
# p = PathKeys(lpks, 4)

# imult = MultipleLinearPathKey(identity, [1, 2])
# imult2 = MultipleLinearPathKey(identity, [1, 2])
# v = [i1, i4, imult, imult2]
# v2 = [i1, i2, imult, imult2]
# pv = PathKeys(v, 4)
# pv2 = PathKeys(v2)

# @code_warntype p(y, a)
# @benchmark p(y, a)
# @code_warntype pv(y, a)
# @benchmark pv(y, a)
# @code_warntype p(a)
# @code_warntype pv(a)

# # Testing to see if let is worthwhile -- no difference
# p2 = PathKeysLet(lpks, 4)
# @code_warntype p2(y, a)
# @benchmark p2(y, a)

# # Testing whether a direct copy is worthwhile -- no difference
# p3 = PathKeys([LinearPathKey(1), LinearPathKey(x -> 25, 2),
#                LinearPathKey(3), LinearPathKey(csv, 4)], 4)
# @code_warntype p3(y, a)
# @benchmark p3(y, a)

# # Testing whether StaticArrays is worthwhile -- slower
# p4 = StaticPathKeys(SVector{4}(v), 4)
# @code_warntype p4(y, a)
# @benchmark p4(y, a)

# # Testing whether mutability matters -- it has no effect, which, in reality, might be convenient
# mp = MPathKeys(lpks, 4)

# mpv = MPathKeys(v, 4)

# @code_warntype mp(y, a)
# @benchmark mp(y, a)
# @code_warntype mpv(y, a)
# @benchmark mpv(y, a)

# # Comparison to direct benchmark
# ff(x) = 25
# function test!(x, a)
#     x[1] = identity(getindex(a, 1))
#     x[2] = ff(getindex(a, 2))
#     x[3] = identity(getindex(a, 3))
#     x[4] = csv(getindex(a, 4))
#     return x
# end
# @code_warntype test!(y, a)
# @benchmark test!(y, a)
# Comparison to direct benchmark using Vector{Function}, Vector{Int}
# Conclusion: same time and allocation as using the functor
# function vtest!(x, a, fs::Vector{<:Function}, idxs, N::Int)
#     n = 1
#     while n ≤ N
#         x[n] = fs[n](getindex(a, idxs[n]))
#         n += 1
#     end
#     return x
# end
# vfs = [identity, ff, identity, csv]
# indexes = [1, 2, 3, 4]
# @code_warntype vtest!(y, a, vfs, indexes, 4)
# @benchmark vtest!(y, a, vfs, indexes, 4)

# #### Comparison vs. simplification
# ipk1 = IndexedPathKey(identity, 1)
# ipk2 = IndexedPathKey(ff, 2)
# ipk3 = IndexedPathKey(identity, 3)
# ipk4 = IndexedPathKey(csv, 4)
# vi = [ipk1, ipk2, ipk3, ipk4]
# @code_warntype ipk1(a)
# pi = PathKeys(vi, length(vi))
# pi2 = PathKeys([IndexedPathKey(identity, 1), IndexedPathKey(ff, 2), IndexedPathKey(identity, 3), IndexedPathKey(csv, 4)], 4)
# # using let statement in inner constructor -- same speed
# pi2_2 = PathKeys2([IndexedPathKey(identity, 1), IndexedPathKey(ff, 2), IndexedPathKey(identity, 3), IndexedPathKey(csv, 4)], 4)
# # using StaticArrays -- slower than Vector
# pi_s = StaticPathKeys(SVector{4}([IndexedPathKey(identity, 1), IndexedPathKey(ff, 2), IndexedPathKey(identity, 3), IndexedPathKey(csv, 4)]), 4)
# @code_warntype pi(y, a)
# @code_warntype pi2(y, a)
# @code_warntype pi2_2(y, a)
# @code_warntype pi_s(y, a)
# @benchmark pi(y, a)
# @benchmark pi2(y, a)
# @benchmark pi2_2(y, a)
# @benchmark pi_s(y, a)
# # #
# imult_pk1 = IndexedPathKey(identity, [1, 2])
# imult_pk2 = IndexedPathKey(identity, [1, 2])
# iv = [ipk1, ipk4, imult_pk1, imult_pk2]
# @code_warntype imult_pk1(a)
# pvi = PathKeys(iv, 4)
# @code_warntype pvi(y, a)
# @benchmark pvi(y, a)
# @benchmark pvi(a)

# # PathKeysTuple test
# pkt = PathKeysTuple(iv, 4)
# @code_warntype pkt(a)
# @benchmark pkt(a)

#### Comparison of Tuple vs. Vector arg -- Tuple is still worse.
# tvi = (ipk1, ipk2, ipk3, ipk4)
# tpi = PathKeys3(tvi, length(tvi))
# @benchmark tpi(y, a)
# @code_warntype tpi(y, a)
# #
# tiv = (ipk1, ipk4, imult_pk1, imult_pk2)
# tpvi = PathKeys3(tiv, 4)
# @code_warntype tpvi(y, a)
# @benchmark tpvi(y, a)

# #### Outer constructors: LinearPathKey
# LinearPathKey(i::Int) = LinearPathKey(identity, i)

#### functor: LinearPathKey
# function (p::AbstractIndexedPathKey)(A)
#     p.f(getindex(A, p.i))
# end

#### Outer constructors: PathKeys
# PathKeys(ftrs) = PathKeys(ftrs, length(ftrs))

#### functor: PathKeys
# both internal functions yield no performance gain
# function _pcall(p::AbstractPathKey, A)
#     p(A)
# end
# function _pfunctor(fs::Vector{T} where {T<:AbstractPathKey}, x, A, N)
#     n = 1
#     while n ≤ N
#         x[n] = fs[n](A)
#         n += 1
#     end
#     return x
# end
# function (p::AbstractPathKeys)(x, A)
#     n = 1
#     N = p.N
#     # N = length(p)
#     fs = p.ftrs
#     while n ≤ N#p.N
#         # x[n] = p.ftrs[n](A)
#         # x[n] = _pcall(p.ftrs[n], A) # still type-unstable
#         # Realistic options: use local fs, or, access via getindex
#         x[n] = fs[n](A)
#         # x[n] = p[n](A)
#         n += 1
#     end
#     return x
#     # _pfunctor(p.ftrs, x, A, p.N)
# end
# function (p::AbstractPathKeys)(A)
#     x = Vector{Any}(undef, p.N)
#     p(x, A)
# end

#### special functor: PathKeys3
# function _pk3functorcall(x, A, fs::Vararg{Any, N}) where {N}
#     n = 1
#     while n ≤ N
#         x[n] = fs[n](A)
#         n += 1
#     end
#     return x
# end
# function _pk3functorcall2(x, A, fs::NTuple{N, Any}) where {N}
#     n = 1
#     while n ≤ N
#         x[n] = fs[n](A)
#         n += 1
#     end
#     return x
# end
# function (p::AbstractPathKeys3)(x, A)
#     _pk3functorcall2(x, A, p.ftrs)
# end


# struct PathKeysTuple{T} <: AbstractPathKeys{T}
#     ftrs::T
#     N::Int
# end

#### Outer constructors: PathKeysTuple
# PathKeysTuple(ftrs) = PathKeysTuple(ftrs, length(ftrs))

#### functor: PathKeysTuple
# function (p::PathKeysTuple)(A)
#     Tuple(p.ftrs[n](A) for n = 1:p.N)
# end

# struct PathKeysScalar{T} <: AbstractPathKeys{T}
#     ftrs::T
#     N::Int
# end

# #### Outer constructors: PathKeysScalar
# PathKeysScalar(ftrs) = PathKeysScalar(ftrs, length(ftrs))

# #### functor: PathKeysScalar
# function (p::PathKeysScalar)(A)
#     p.ftrs[p.N](A)
# end
# psc = PathKeysScalar([ipk1], 1)
# @benchmark ipk1(a)
# @benchmark psc(a)
################ Some random testing
# Interesting topic: fg = SimpleNode benchmarks at about 1.33 times slower
# It is because fg is simply a reference to the type, which, when called as fg() dispatches;
# gf is a declared function with a specific method.
gf() = SimpleNode()
mat = reshape([1:200000;], (4, 50000));
# pv2 = PathKeys(v2)
# pn = PathKeys([i1, i3, i2, i3])
# pni = PathKeys([ipk1, ipk3, ipk2, ipk3], 4)
pni = PathKeys([IndexedPathKey(i) for i = 1:4]);
# @benchmark grow!(gf, SimpleNode(), pn, eachcol(mat))
# @benchmark grow!(gf, SimpleNode(), pv2, eachcol(mat))
@benchmark grow!(gf, SimpleNode(), pni, eachcol(mat)) seconds=30
# mat[1:2, :] .= 1;
# @benchmark grow!(gf, SimpleNode(), pn, eachcol(mat))
# @benchmark grow!(gf, SimpleNode(), pv2, eachcol(mat))
# very long trees -- stack overflow
pmat = reshape([1:200000;], (20, 10000));
# pn2 = PathKeys([LinearPathKey(i) for i = 1:20]);
pni2 = PathKeys([IndexedPathKey(i) for i = 1:20]);
# @benchmark grow!(gf, SimpleNode(), pn2, eachcol(pmat))
@benchmark grow!(gf, SimpleNode(), pni2, eachcol(pmat)) seconds=30
#### Another test
pm = PathKeys([[IndexedPathKey(i) for i = 1:3]; IndexedPathKey(string, 4)]);
@benchmark grow!(gf, SimpleNode(), pm, eachcol(mat)) seconds=30
pm2 = PathKeys([[IndexedPathKey(i) for i = 1:10]; IndexedPathKey(string, 11); [IndexedPathKey(i) for i = 12:20]]);
@benchmark grow!(gf, SimpleNode(), pm2, eachcol(pmat)) seconds=30
#### Much longer trees
qmat = reshape([1:200000;], (50, 4000));
pni3 = PathKeys([IndexedPathKey(i) for i = 1:50]);
pm3 = PathKeys([[IndexedPathKey(i) for i = 1:10]; IndexedPathKey(string, 11); [IndexedPathKey(i) for i = 12:50]]);
@benchmark grow!(gf, SimpleNode(), pni3, eachcol(qmat)) seconds=30
@benchmark grow!(gf, SimpleNode(), pm3, eachcol(qmat)) seconds=30
#
qmat = reshape([1:200000;], (100, 2000));
pni4 = PathKeys([IndexedPathKey(i) for i = 1:100]);
pm4 = PathKeys([[IndexedPathKey(i) for i = 1:10]; IndexedPathKey(string, 11); [IndexedPathKey(i) for i = 12:100]]);
@benchmark grow!(gf, SimpleNode(), pni4, eachcol(qmat)) seconds=30
@benchmark grow!(gf, SimpleNode(), pm4, eachcol(qmat)) seconds=30
#
qmat = reshape([1:200000;], (200, 1000));
pni5 = PathKeys([IndexedPathKey(i) for i = 1:200]);
pm5 = PathKeys([[IndexedPathKey(i) for i = 1:10]; IndexedPathKey(string, 11); [IndexedPathKey(i) for i = 12:200]]);
@benchmark grow!(gf, SimpleNode(), pni5, eachcol(qmat)) seconds=30
@benchmark grow!(gf, SimpleNode(), pm5, eachcol(qmat)) seconds=30
#
qmat = reshape([1:200000;], (10, 20000));
qmatc = [qmat[:, j] for j = 1:size(qmat, 2)];
pni6 = PathKeys([IndexedPathKey(i) for i = 1:10]);
pm6 = PathKeys([[IndexedPathKey(i) for i = 1:9]; IndexedPathKey(string, 10)]);
@benchmark grow!(gf, SimpleNode(), pni6, eachcol(qmat)) seconds=30
@benchmark grow!(gf, SimpleNode(), pni6, qmatc) seconds=30
@benchmark grow!(gf, SimpleNode(), pm6, eachcol(qmat)) seconds=30
@benchmark grow!(gf, SimpleNode(), pm6, qmatc) seconds=30
#### Interesting test comparison of ks... vs. Vararg{Any, N} where {N}
# 4x50000: (281.766, 300.994) .- (280.666, 300.635)
# 20x10000: (324.17, 335.412) .- (251.17, 248.56)
# 50x4000: (477.709, 472.685) .- (387.73, 389.24)
# 100x2000: (670.685, 647.508) .- (577.006, 546.230)
# 200x1000: (1.201, 1.209) .- (1.029, 1.015)
# 10x20000: (287.992, 319.813) .- (223.464, 252.375)
# The gain seems to be approximately 70-90ms
A = [1 20; 1 50; 1 100; 1 200; 1 10]; b = [73, 89.89, 93.67, 172, 64.52]; # b ≡ ks... - Vararg
A \ b # [71.935, 0.2397] if first 3, [57.213, 0.5397] if first 4, [57.992, 0.5345] if 5
# For 4 arguments, the gain is negligible. The overall performance gain increases
# with the number of arguments slurped. Linear estimate: intercept 57.992ms, slope 0.53ms/arg.
# The relative gains are highest at low number of args, but, it is still appreciable even
# at 200 args. Expected use cases involve up to 10:20 args, at which point it is
# a relative performance gain of 22%.
# - Conclusion: worthwhile based on generic testing, but tests on a real use case will
#   provide useful insight.
# - Second conclusion: the trends below demonstrate that specialization of ks... as
#   Vararg{Any, N} where {N} is worthwhile.
# - Third conclusion: specialization of get_valpush!s ps... as a Vararg{Any, N} where {N}
#   also leads to performance enhancements.
#   Direct comparison to the previous constructors yields a reduction in memory of
#   the data structure by 50%, and memory allocated during construction by about 25%;
#   performance is about 200ms faster. However, the new constructors are far
#   easier to wield, and, they allow for optimized computation of elements in
#   the stored data itself: pre-computation of parse_split("RU_CRU") and parse_split("bank_chain")
#   Surprisingly, the parse_split_first / parse_split_last calls cost ≈ 12% of
#   the overall time for "RU_CRU", and ≈ 15% of the overall time for "bank_chain".
#   Pre-computing these would save ≈ 30% overall.
#   Furthermore, they create a vast number of temporaries which do not need to be generated
#   otherwise. This constitutes 896 bytes per column, which, obviously, quickly adds up.
#   Consider the approximation: 1215165cols/218MiB = 5.574151376146789e6 col/GiB
#   At 896 bytes per column, for the ≈ 24GiB of data for 4299 wafers, this
#   corresponds to 119.866GiB of temporaries just to parse the strings to Int.
# - Fourth conclusion: Pre-computing the "RU_CRU" and "bank_chain" parsing,
#   then storing in Arrow file reduces construction time by (8.051 - 5.724) = 2.327 sec
#   = 2327ms. This is ≈ 29% reduction in construction time.
#   For the test case of 1673512 nodes, this amounts to 937.16672 MiB which are needed
#   hence, the total 2.16GiB allocated during construction means that approximately
#   1222.8339999999998 is spent during construction.
#   Consider that previously, even under the new scheme, 3113MiB was spent during construction,
#   of which 1088MiB (1215165 * 896) was spent purely on parsing.
#   Comparatively, the original scheme, 3952MiB was spent during construction.
#   This is a reduction of (3113 - 2116) MiB, or a 32% reduction w.r.t. the non-optimized
#   Arrow format under the new scheme. With respect to the original scheme,
#   These changes represent a (3952 - 2116) MiB reduction, or 46% reduction overall.
################
prefix = "/nfs/site/home/aradclif/my_cit_scratch1/diagnosis/1276/x76se/gt_g1m/2021-09-15";
pni7 = PathKeys([IndexedPathKey(i) for i = 1:10]);
# 10x9478: (63.419,) .- (40.976,)
suffix = "D127YKQ0_696.arrow"; # 1.8MiB
file = joinpath(prefix, suffix);
dmat = permutedims(hcat([convert(Vector{Any}, copy(col)) for col in Arrow.Table(file)]...));
@benchmark grow!(gf, SimpleNode(), pni7, eachcol(dmat)) seconds=30
# 10x20210: (142.174,) .- (93.131,)
suffix2 = "D105EK60_242.arrow" # 3.6MiB
file2 = joinpath(prefix, suffix2);
dmat2 = permutedims(hcat([convert(Vector{Any}, copy(col)) for col in Arrow.Table(file2)]...));
@benchmark grow!(gf, SimpleNode(), pni7, eachcol(dmat2)) seconds=30
# 10x40383: (316.092,) .- (205.188,)
suffix3 = "D105EK60_241.arrow" # 7.2MiB
file3 = joinpath(prefix, suffix3);
dmat3 = permutedims(hcat([convert(Vector{Any}, copy(col)) for col in Arrow.Table(file3)]...));
@benchmark grow!(gf, SimpleNode(), pni7, eachcol(dmat3)) seconds=30
# 10x82459: (732.377,) .- (548.165,)
suffix4 = "D108E4F0_939.arrow" # 15MiB
file4 = joinpath(prefix, suffix4);
dmat4 = permutedims(hcat([convert(Vector{Any}, copy(col)) for col in Arrow.Table(file4)]...));
@benchmark grow!(gf, SimpleNode(), pni7, eachcol(dmat4)) seconds=30
# 10x165617: (1551,) .- (1128,)
suffix5 = "D112EF4A_689.arrow" # 30MiB
file5 = joinpath(prefix, suffix5);
dmat5 = permutedims(hcat([convert(Vector{Any}, copy(col)) for col in Arrow.Table(file5)]...));
@benchmark grow!(gf, SimpleNode(), pni7, eachcol(dmat5)) seconds=30
# 10x336706: (3189,) .- (2402,)
suffix6 = "D116E1D0_189.arrow" # 61MiB
file6 = joinpath(prefix, suffix6);
dmat6 = permutedims(hcat([convert(Vector{Any}, copy(col)) for col in Arrow.Table(file6)]...));
@benchmark grow!(gf, SimpleNode(), pni7, eachcol(dmat6)) seconds=60
# 10x623040: (6115,) .- (4646,)
suffix7 = "D126YK60_751.arrow" # 123MiB
file7 = joinpath(prefix, suffix7);
dmat7 = permutedims(hcat([convert(Vector{Any}, copy(col)) for col in Arrow.Table(file7)]...));
@benchmark grow!(gf, SimpleNode(), pni7, eachcol(dmat7)) seconds=120
# 10x1215165: (11242,) .- (9226,)
suffix8 = "D117E4C0_36.arrow" # 218MiB
file8 = joinpath(prefix, suffix8);
dmat8 = permutedims(hcat([convert(Vector{Any}, copy(col)) for col in Arrow.Table(file8)]...));
@benchmark grow!(gf, SimpleNode(), pni7, eachcol(dmat8)) seconds=240
#### Comparison to extant methods
pni8 = PathKeys([IndexedPathKey(x -> Tuple(x), [5,6]), IndexedPathKey(parse_split_first, 1),
                 IndexedPathKey(parse_split_last, 1),
                 IndexedPathKey(parse_split_first, 2),
                 IndexedPathKey(parse_split_last, 2),
                 IndexedPathKey(7), IndexedPathKey(8)])
vi8 = PathKeys([IndexedPathKey(10), IndexedPathKey(pathprefix, 9)])
@benchmark valgrow!(gf, SimpleNode(), vi8, pni8, eachcol(dmat8)) seconds=240
pni9 = PathKeys([IndexedPathKey(i) for i = 1:7]);
vi9 = PathKeys([IndexedPathKey(i) for i = 9:10]);
cellinst(x) = @inbounds (x[1], pathprefix(x[2]))
vi9 = IndexedPathKey(cellinst, [10, 9])
@benchmark valgrow!(gf, SimpleNode(), vi9, pni8, eachcol(dmat8)) seconds=240
#### Comparison given optimized Arrow files
prefix10 = "/nfs/site/home/aradclif/my_cit_scratch1/diagnosis/1276/x76se/gt_g1m/2021-09-24";
suffix10 = "D117E4C0_36.arrow";
file10 = joinpath(prefix10, suffix10);
dmat10 = permutedims(hcat([convert(Vector{Any}, copy(col)) for col in Arrow.Table(file10)]...));
pni10 = PathKeys([IndexedPathKey(x -> Tuple(x), [1, 2]); [IndexedPathKey(i) for i = 3:8]])
vi10 = IndexedPathKey(cellinst, [10, 9])
@benchmark valgrow!(gf, SimpleNode(), vi10, pni10, eachcol(dmat10)) seconds=240
############################################################################################
#### Iterators: see p. 421-430, 2021-09-13
# - foreach variants
# - findall variants
# - forall variants
# - count variants
# - countall variants
#### Design concept
# "at" means apply f to the nodes at a given level
# "from" means apply f to the nodes at a given level, then recursively applying f
# to all nodes linked to each node below.
# "through" means apply f to all nodes up to on levels up to and including the
# given level.
################
# # modified 2021-10-03; see p. 480-481
# function foreachat(f::Function, t::AbstractNode, N::Int, C::Int)
#     # if C < N - 1
#     #     if !isempty(t)
#     #         for p in t
#     #             foreachat!(f, p.second, N, C + 1)
#     #         end
#     #     end
#     # elseif C == N - 1
#     #     if !isempty(t)
#     #         for p in t
#     #             f(p)
#     #         end
#     #     end
#     # end
#     # nothing
#     # A terse definition
#     isempty(t) && return nothing
#     C̃ = C + 1
#     if C̃ < N
#         for p in t
#             foreachat(f, p.second, N, C̃)
#         end
#     elseif C̃ == N
#         # for p in t
#         #     f(p)
#         # end
#         # Alternative that may lead to more inlining:
#         foreach(f, t)
#     end
#     nothing
# end
# function foreachat(f::Function, t::AbstractNode, N::Int)
#     foreachat(f, t, N, 1)
# end

# function foreachats!(f::Function, x::AbstractNode, Ns::NTuple{M, Int}) where {M}
#     states = C + 1 .== Ns
#     if any(states) && !isempty(x)
#         for i in eachindex(states)
#             if states[i]
#                 for v in values(x)
#                     f(v)
#                 end
#             end
#         end
#     else
#         if !isempty(x)
#             for v in values(x)
#                 foreachats!(f, v, Ns, C + 1)
#             end
#         end
#     end
#     nothing
#     # A terse definition
#     # isempty(x) && return nothing
#     # C̃ = C + 1
#     # for N in Ns
#     #     if C̃ == N
#     #         for p in x
#     #             f(p.second)
#     #         end
#     #         # or, foreach(p -> f(p.second), x)
#     #     end
#     # end
#     # for p in x
#     #     foreachats!(f, p.second, Ns, C̃)
#     # end
#     # return nothing
# end

# # Not exactly the most value, since foreachat covers this, but, perhaps, there
# # is minor value in providing a convenience iterator that applies f directly
# # to values rather than to pairs of (K => V)
# function foreachat_values!(f::Function, x::AbstractNode, N::Int, C::Int)
#     if C < N - 1
#         if !isempty(x)
#             for v in values(x)
#                 foreachat!(f, v, N, C + 1)
#             end
#         end
#     elseif C == N - 1
#         if !isempty(x)
#             for v in values(x)
#                 f(v)
#             end
#         end
#     end
#     nothing
#     # A terse definition
#     # isempty(x) && return nothing
#     # C̃ = C + 1
#     # if C̃ < N
#     #     for p in x
#     #         foreachat!(f, p.second, N, C̃)
#     #     end
#     # elseif C̃ == N
#     #     for p in x
#     #         f(p.second)
#     #     end
#     # end
#     # return nothing
# end
# function foreachat_values!(f::Function, x::AbstractNode, N::Int)
#     C = 1
#     foreachat_values!(f, x, N, C)
# end

# # Cover both data and functional filters; return an iterable of keys.
# function subsetgen(ks, x::AbstractNode)
#     return keys(x) ∩ ks
# end
# function subsetgen(f::Function, x::AbstractNode)
#     return f(x)
# end

# function foreachat_filter!(f::Function, x::AbstractNode, N::Int, C::Int, subsets)
#     if C < N - 1
#         if !isempty(x)
#             dKs = subsetgen(subsets[C], x)
#             for k in dKs
#                 foreachat_filter!(f, x[k], N, C + 1, subsets)
#             end
#         end
#     elseif C == N - 1
#         if !isempty(x)
#             dKs = subsetgen(subsets[C], x)
#             for k in dKs
#                 f(x[k])
#             end
#         end
#     end
#     nothing
#     # Terse version
#     # isempty(x) && return nothing
#     # C̃ = C + 1
#     # dKs = subsetgen(subsets[C], x)
#     # if C̃ < N
#     #     for k in dKs
#     #         foreachat_filter!(f, x[k], N, C̃, subsets)
#     #     end
#     # elseif C̃ == N
#     #     for k in dKs
#     #         f(x[k])
#     #     end
#     # end
#     # return nothing
# end

# function foreachat_filter!(f::Function, x::AbstractNode, subsets)
#     N = length(subsets)
#     C = 1
#     foreachat_filter!(f, x, N, C, subsets)
#     # foreachat_filter!(f, x, length(subsets), 1, subsets)
# end

# ####
# # modified 2021-10-03; see p. 480-481
# function findallpathsat!(f::Function, A::Vector{<:Vector}, ks::Vector, t::AbstractNode,
#                          N::Int, C::Int)
#     # if C < N - 1
#     #     if !isempty(t)
#     #         for p in t
#     #             findallat!(f, A, setindex!(ks, p.first, C), p.second, N, C + 1)
#     #         end
#     #     end
#     # elseif C == N - 1
#     #     if !isempty(t)
#     #         for p in t
#     #             setindex!(ks, p.first, C)
#     #             f(p) && push!(A, copyto!(similar(ks), ks))
#     #         end
#     #     end
#     # end
#     # return A
#     # Terse form
#     isempty(t) && return A
#     C̃ = C + 1 #tmp = N - 1
#     if C̃ < N # equivalent to C < N - 1
#         for p in t
#             @inbounds setindex!(ks, p.first, C)
#             findallat!(f, A, ks, p.second, N, C̃)
#         end
#     elseif C̃ == N
#         for p in t
#             # @inbounds setindex!(ks, p.first, C)
#             # f(p) && push!(A, copyto!(similar(ks), ks))
#             f(p) && (setindex!(ks, p.first, C); push!(A, copyto!(similar(ks), ks)))
#         end
#     end
#     return A
# end

# function findallpathsat(f::Function, t::AbstractNode, N::Int)
#     ks = Vector{Any}(undef, N - 1)
#     A = Vector{Vector{Any}}()
#     C = 1
#     findallpathsat!(f, A, ks, t, N, C)
# end

# ####
# function countat!(f::Function, A::Vector{Int}, t::AbstractNode, N::Int, C::Int)
#     # if C < N - 1
#     #     if !isempty(t)
#     #         for v in values(t)
#     #             countat!(f, A, v, N, C + 1)
#     #         end
#     #     end
#     # elseif C == N - 1
#     #     if !isempty(t)
#     #         for p in t
#     #             A[1] += f(p)
#     #         end
#     #     end
#     # end
#     # return A
#     # Terse form
#     isempty(t) && return A
#     C̃ = C + 1
#     if C̃ < N
#         for p in t
#             countat!(f, A, p.second, N, C̃)
#         end
#     elseif C̃ == N
#         for p in t
#             A[1] += f(p)
#         end
#     end
#     return A
# end

# function countat(f::Function, t::AbstractNode, N::Int)
#     A = Int[0]
#     C = 1
#     countat!(f, A, t, N, C)
#     A[1]
# end

# ####
# function forall_depthfirst!(f::Function, t::AbstractNode)
#     if !isempty(t)
#         for p in t
#             forall_depthfirst!(f, p.second)
#         end
#     end
#     f(t)
#     return nothing
#     # Terse form
#     # isempty(t) && (f(t); return nothing)
#     # for p in t
#     #     forall_depthfirst!(f, p.second)
#     # end
#     # f(t)
#     # return nothing
# end

# function forall_breadthfirst!(f::Function, t::AbstractNode)
#     f(t)
#     if !isempty(t)
#         for p in t
#             forall_breadthfirst!(f, p.second)
#         end
#     end
#     return t
#     # Terse form
#     # f(t)
#     # isempty(t) && return nothing
#     # for p in t
#     #     forall_breadthfirst!(f, p.second)
#     # end
#     # return nothing
# end

# function forallfrom!(f::Function, t::AbstractNode, N::Int, C::Int)
#     # if C < N - 1
#     #     if !isempty(t)
#     #         for p in t
#     #             forallfrom!(f, p.second, N, C + 1)
#     #         end
#     #     end
#     # elseif C == N - 1
#     #     if !isempty(t)
#     #         for p in t
#     #             forall_breadthfirst!(f, p.second)
#     #         end
#     #     end
#     # end
#     # nothing
#     # Terse form
#     isempty(t) && return nothing
#     C̃ = C + 1
#     if C̃ < N
#         for p in t
#             forallfrom!(f, p.second, N, C̃)
#         end
#     elseif C̃ == N
#         for p in t
#             forall_breadthfirst!(f, p.second)
#         end
#     end
#     nothing
# end

# function forallfrom!(f::Function, t::AbstractNode, N::Int)
#     forallfrom!(f, t, N, 1)
# end

# function forallthrough!(f::Function, t::AbstractNode, N::Int, C::Int)
#     f(t)
#     isempty(t) && return nothing
#     C̃ = C + 1
#     if C̃ < N
#         for p in t
#             forallthrough!(f, p.second, N, C̃)
#         end
#     elseif C̃ == N
#         for p in t
#             f(p.second)
#         end
#     end
#     nothing
#     # Works, but enters function for each node at last level. Might be less efficient.
#     # f(t)
#     # isempty(t) && return nothing
#     # C̃ = C + 1
#     # if C̃ ≤ N
#     #     for p in t
#     #         forallthrough!(f, p.second, N, C̃)
#     #     end
#     # end
#     # return nothing
# end

# function forallthrough(f::Function, t::AbstractNode, N::Int)
#     forallthrough!(f, t, N, 1)
# end

# ####
# function countall!(f::Function, A::Vector{Int}, t::AbstractNode)
#     # if !isempty(t)
#     #     for p in t
#     #         countall!(f, A, p.second)
#     #     end
#     # end
#     # A[1] += f(t)
#     # return A
#     # Terse form -- 2021-09-24: approtimately 18% faster
#     A[1] += f(t)
#     isempty(t) && return A
#     for p in t
#         countall!(f, A, p.second)
#     end
#     return A
# end

# function countall(f::Function, t::AbstractNode)
#     A = Int[0]
#     countall!(f, A, t)
#     A[1]
# end

# function countallfrom!(f::Function, A::Vector{Int}, t::AbstractNode, N::Int, C::Int)
#     # if C < N - 1
#     #     if !isempty(t)
#     #         for p in t
#     #             countallfrom!(f, A, p.second, N, C + 1)
#     #         end
#     #     end
#     # elseif C == N - 1
#     #     if !isempty(t)
#     #         for p in t
#     #             countall!(f, A, p.second)
#     #         end
#     #     end
#     # end
#     # return A
#     # Terse form
#     # A[1] += f(t) # Not necessary as countall! will perform this, even if empty
#     isempty(t) && return A
#     C̃ = C + 1
#     if C̃ < N
#         for p in t
#             countallfrom!(f, A, p.second, N, C̃)
#         end
#     elseif C̃ == N
#         for p in t
#             countall!(f, A, p.second)
#         end
#     end
#     return A
# end

# function countallfrom(f::Function, t::AbstractNode, N::Int)
#     A = Int[0]
#     C = 1
#     countallfrom!(f, A, t, N, C)
#     A[1]
# end

# function countallthrough!(f::Function, A::Vector{Int}, t::AbstractNode, N::Int, C::Int)
#     A[1] += f(t)
#     isempty(t) && return A
#     C̃ = C + 1
#     if C̃ < N
#         for p in t
#             countallthrough!(f, A, p.second, N, C̃)
#         end
#     elseif C̃ == N
#         for p in t
#             A[1] += f(p.second)
#         end
#     end
#     return A
#     # Works, but enters function for each node at last level. Might be less efficient.
#     # A[1] += f(t)
#     # isempty(t) && return A
#     # C̃ = C + 1
#     # if C̃ ≤ N
#     #     for p in t
#     #         countallthrough!(f, A, p.second, N, C̃)
#     #     end
#     # end
#     # return A
# end

# function countallthrough(f::Function, t::AbstractNode, N::Int)
#     A = Int[0]
#     C = 1
#     countallthrough!(f, A, t, N, C)
#     A[1]
# end
# ################################################################
# #### 2021-10-03: see p. 479-484
# function Base.foreach(f::Function, t::AbstractNode)
#     for p in t
#         f(p)
#     end
#     nothing
# end

# function Base.findall(f::Function, t::AbstractNode)
#     idxs = Any[] # or, idxs = Vector{eltype(keys(t))}()
#     for p in t
#         f(p) && (push!(idxs, p.first))
#     end
#     idxs
#     # or:
#     collect(first(p) for p in t if f(p))
# end

# function Base.count(f::Function, t::AbstractNode)
#     cnt = 0
#     for p in t
#         f(p) && (cnt += 1)
#     end
#     cnt
# end

# function findallpaths!(f::Function, A::Vector{<:Vector}, ks::Vector, t::AbstractNode)
#     if !isempty(t)
#         push!(ks, nothing)
#         C = lastindex(ks)
#         for p in t
#             f(p) && (setindex!(ks, p.first, C); push!(A, copyto!(similar(ks), ks)))
#         end
#         # then proceed recursively
#         for p in t
#             setindex!(ks, p.first, C)
#             findallpaths!(f, A, ks, p.second)
#         end
#         pop!(ks) # or, deleteat!(ks, C)
#     end
#     return A
# end

# function findallpaths(f::Function, t::AbstractNode)
#     A = Vector{Vector{Any}}()
#     ks = Any[]
#     findallpaths!(f, A, ks, t)
# end

# #### A draft of analyzeat!, p. 483, 2021-10-03
# function analyzeat!(f::Function, dest::Tuple{Vararg{Array{T}} where T}, ks::Vector,
#                     t::AbstractNode, N::Int, C::Int)
#     isempty(t) && return dest
#     C̃ = C + 1
#     if C̃ < N
#         for p in t
#             setindex!(ks, p.first, C)
#             analyzeat!(f, dest, ks, p.second, N, C̃)
#         end
#     elseif C̃ == N
#         for p in t
#             setindex!(ks, p.first, C)
#             f(dest, p, ks)
#         end
#     end
#     return dest
# end
# function analyzeat!(f::Function, dest::Tuple{Vararg{Array{T}} where T}, t::AbstractNode, N::Int)
#     ks = Vector{Any}(undef, N - 1)
#     analyzeat!(f, dest, ks, t, N, 1)
# end
# #### 2021-10-05: p. 1-4
# function analyzeat!(f::Function, dest::Tuple{Vararg{Array{T}} where T},
#                     ts::Vector{<:AbstractNode}, N::Int)
#     ks = Vector{Any}(undef, N - 1)
#     C = 1
#     for t in ts
#         analyzeat!(f, dest, ks, t, C, N)
#     end
#     return dest
# end
# function analyzeat(f::Function, dims::Tuple{Vararg{NTuple{S, Int} where S}},
#                    ts::Vector{<:AbstractNode}, N::Int)
#     dest = [zeros(Int, d) for d in dims]
#     analyzeat!(f, dest, ts, N)
# end
# #
# function tanalyzeat(f::Function, dims::Tuple{Vararg{NTuple{S, Int} where S}},
#                     gs::Vector{<:AbstractNode}, L::Int)
#     N = length(gs)
#     M = Threads.threads()
#     ranges = equalranges(N, M)
#     A = Vector{Tuple{(Array{Int, length(d)} for d in dims)...}}(undef, M)
#     Threads.@threads for m = 1:M
#         A[m] = analyzeat(f, dims, gs[ranges[m]], L)
#     end
#     return A
# end

############################################################################################
# #### p. 1-5, 2021-10-05
# """
#     mthpair(q::Int, r::Int, M::Int, m::Int)

# Return the `m`ᵗʰ pair of (lower,upper) bounds for division of a length, N = `q``M` + `r`,
# into `M` parts, where the remainder of Euclidean division, `r` is distributed equally
# across the first `r` parts.

# # Notes
# Formally, the expression is:
# (l, u) = ((m - 1)q + 1 + (m - 1), mq + m) if m ≤ r
#          ((m - 1)q + 1 + r, mq + r)       if r < m ≤ M
# This has a variety of algebraic simplifications, salient examples
# of which are:
# (l, u) = ((m - 1)q + m, mq + m)           if m ≤ r
#          (mq + r - q + 1, mq + r)         if r < m ≤ M
# (l, u) = (mq - q + m, mq + m)             if m ≤ r
#          (mq + r - q + 1, mq + r)         if r < m ≤ M
# (l, u) = (m(q + 1) - q, m(q + 1))         if m ≤ r
#          (mq + r - q + 1, mq + r)         if r < m ≤ M

# # Examples
# ```jldoctest
# julia> M = 7; N = 22; q, r = divrem(N, M)
# (3, 1)

# julia> mthpair(q, r, M, 1)
# (1, 4)

# julia> mthpair(q, r, M, 2)
# (5, 7)
# ```
# """
# function mthpair(q::Int, r::Int, M::Int, m::Int)
#     if m ≤ r
#         u = m * q + m
#         l = u - q
#         return (l, u)
#     elseif r < m ≤ M
#         u = m * q + r
#         l = u - q + 1
#         return (l, u)
#     end
# end

# """
#     mthpair(N::Int, M::Int, m::Int)

# Return the `m`ᵗʰ pair of (lower, upper) bounds for division of a length `N` into `M` parts,
# where the remainder (if it is nonzero) is equally distributed across the parts.
# """
# function mthpair(N::Int, M::Int, m::Int)
#     q, r = divrem(N, M)
#     mthpair(q, r, M, m)
# end

# """
#     mthrange(q::Int, r::Int, M::Int, m::Int)

# Return the `m`ᵗʰ range of (lower,upper) bounds for division of a length, N = `q``M` + `r`,
# into `M` parts, where the remainder of Euclidean division, `r` is distributed equally
# across the first `r` parts.

# # Examples
# ```jldoctest
# julia> M = 7; N = 22; q, r = divrem(N, M)
# (3, 1)

# julia> mthrange(q, r, M, 1)
# 1:4

# julia> mthrange(q, r, M, 2)
# 5:7
# ```
# """
# function mthrange(q::Int, r::Int, M::Int, m::Int)
#     if m ≤ r
#         u = m * q + m
#         l = u - q
#         return l:u
#     elseif r < m ≤ M
#         u = m * q + r
#         l = u - q + 1
#         return l:u
#     end
# end

# """
#     mthrange(N::Int, M::Int, m::Int)

# Return the `m`ᵗʰ range of (lower, upper) bounds for division of a length `N` into `M` parts,
# where the remainder (if it is nonzero) is equally distributed across the parts.
# """
# function mthrange(N::Int, M::Int, m::Int)
#     q, r = divrem(N, M)
#     mthrange(q, r, M, m)
# end

# """
#     equalpairs(q::Int, r::Int, M::Int)

# Return the `M` pairs of (lower, upper) bounds which span a length of N = `q``M` + `r`.
# More efficient than `[mthpair(q, r, M, m) for m = 1:M]`.

# # Examples
# ```jldoctest
# julia> M = 7; N = 22; q, r = divrem(N, M)
# (3, 1)

# julia> equalpairs(q, r, M)
# 7-element Vector{Tuple{Int64, Int64}}:
#  (1, 4)
#  (5, 7)
#  (8, 10)
#  (11, 13)
#  (14, 16)
#  (17, 19)
#  (20, 22)
# ```
# """
# function equalpairs(q::Int, r::Int, M::Int)
#     ps = Vector{Tuple{Int,Int}}(undef, M)
#     qp1 = q + 1
#     u = qp1
#     l = 1
#     ps[1] = (l, u)
#     m = 2
#     while m ≤ M
#         if m ≤ r
#             u += qp1
#             l = u - q
#             ps[m] = (l, u)
#         elseif r < m ≤ M
#             u += q
#             l = u - q + 1
#             ps[m] = (l, u)
#         end
#         m +=1
#     end
#     return ps
# end

# """
#     equalpairs(N::Int, M::Int)

# Return the `M` pairs of (lower, upper) bounds which divide a length `N` into `M` parts,
# where the remainder (if it is nonzero) is equally distributed across the parts.
# """
# function equalpairs(N::Int, M::Int)
#     q, r = divrem(N, M)
#     equalpairs(q, r, M)
# end

# """
#     equalranges(q::Int, r::Int, M::Int)

# Return the `M` ranges of (lower, upper) bounds which span a length of N = `q``M` + `r`.
# More efficient than `[mthrange(q, r, M, m) for m = 1:M]`.

# # Examples
# ```jldoctest
# julia> M = 7; N = 22; q, r = divrem(N, M)
# (3, 1)

# julia> equalranges(q, r, M)
# 7-element Vector{UnitRange{Int64}}:
#  1:4
#  5:7
#  8:10
#  11:13
#  14:16
#  17:19
#  20:22
# ```
# """
# function equalranges(q::Int, r::Int, M::Int)
#     ps = Vector{UnitRange{Int}}(undef, M)
#     qp1 = q + 1
#     u = qp1
#     l = 1
#     ps[1] = l:u
#     m = 2
#     while m ≤ M
#         if m ≤ r
#             u += qp1
#             l = u - q
#             ps[m] = l:u
#         elseif r < m ≤ M
#             u += q
#             l = u - q + 1
#             ps[m] = l:u
#         end
#         m +=1
#     end
#     return ps
# end

# """
#     equalranges(N::Int, M::Int)

# Return the `M` ranges of (lower, upper) bounds which divide a length `N` into `M` parts,
# where the remainder (if it is nonzero) is equally distributed across the parts.
# """
# function equalranges(N::Int, M::Int)
#     q, r = divrem(N, M)
#     equalranges(q, r, M)
# end
# # Clearly less efficient.
# # Invokes 1 multiply, 1(2) add, 1 subtract per mthpair call, which
# # can be reduced to 1(2) add, 1 subtract per pair.
# function mpairs(q::Int, r::Int, M::Int)
#     ps = Vector{Tuple{Int,Int}}(undef, M)
#     m = 1
#     while m ≤ M
#         ps[m] = mthpair(q, r, M, m)
#         m += 1
#     end
#     return ps
# end
# function mranges(q::Int, r::Int, M::Int)
#     ps = Vector{UnitRange{Int}}(undef, M)
#     m = 1
#     while m ≤ M
#         ps[m] = mthrange(q, r, M, m)
#         m += 1
#     end
#     return ps
# end
# M = 7
# N = 22
# q, r = divrem(N, M)
# easy = [mthpair(q, r, M, m) for m = 1:M]
# fast = equalpairs(q, r, M)
# decent = mpairs(q, r, M)
# easy == fast
# @benchmark equalpairs(q, r, M) # 45.46ns
# @benchmark mpairs(q, r, M) # 49.64ns
# easyr = [mthrange(q, r, M, m) for m = 1:M]
# fastr = equalranges(q, r, M)
# decentr = mranges(q, r, M)
# easyr == fastr
# @benchmark equalranges(q, r, M) # 47.35
# @benchmark mranges(q, r, M) # 53.00

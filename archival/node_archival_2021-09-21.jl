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
function Base.get(f::Function, t::A where {A<:AbstractNode{T, U}}, p, q, ps...) where {T, U}
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
Base.get(t::A where {A<:AbstractNode{T, U}}, p) where {T, U} = get(() -> nothing, t, p)
Base.get(t::A where {A<:AbstractNode{T, U}}, p, q) where {T, U} = get(() -> nothing, t, p, q)
Base.get(t::A where {A<:AbstractNode{T, U}}, p, q, ps...) where {T, U} =
    get(() -> nothing, t, p, q, ps...)
# actually, the following should achieve the same result
# Base.get(t::A where {A<:AbstractNode{T, U}}, ps...) = get(() -> nothing, t, ps...)

# Convenient definition of haskey based on get
haspath(t::AbstractNode, p) = get(f, t, p) === nothing ? false : true
# better, likely more performant definition; benchmark to determine which to use,
# then change (2), (3). Even if the above is faster, consider get(...) !== nothing,
# when one can tolerate receiving something other than a AbstractNode.
haspath(t::AbstractNode, p) = isa(get(f, t, p), AbstractNode)
# haspath(t::AbstractNode, p, q) = isa(get(f, t, p, q), AbstractNode)
# haspath(t::AbstractNode, p, q, ps...) = isa(get(f, t, p, q, ps...), AbstractNode)
# Follow pattern: Varargs dispatch to underlying function where possible, rather
# than dispatch on (1), (2), (3) for each method. Single arg form is likely still useful.
haspath(t::AbstractNode, p, ps...) = isa(get(f, t, p, ps...), AbstractNode)

# default should be a sub-type of AbstractNode
Base.get!(x::A where {A<:AbstractNode{T, U}}, key, default) where {T, U} = get!(x.link, key, default)

# f must be some function which returns a sub-type of AbstractNode
Base.get!(f::Function, x::A where {A<:AbstractNode{T, U}}, key) where {T, U} = get!(f, x.link, key)
Base.get!(f::Function, x::A where {A<:AbstractNode{T, U}}, k1, k2) where {T, U} =
    get!(f, get!(f, x, k1), k2)
Base.get!(f::Function, x::A where {A<:AbstractNode{T, U}}, k1, k2, ks...) where {T, U} =
    get!(f, get!(f, get!(f, x, k1), k2), ks...)

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
@test promote_type(x, y) !== z # previously === returned false
ac = ComplexNode{Dict{Any, Any}, Vector{Int}, Dict{Any, Any}, Vector{Any}}
bc = ComplexNode{Dict{Any, Any}, Vector{String}, Dict{Any, Any}, Vector{Any}}
cc = zc
@test promote_type(a, b) !== c # previously === returned false
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
# Best approach: just use dispatch on get! to limit number of possible dispatches
function get_valpush!(f::Function, t::AbstractNode, v, p)
    tmp = get!(f, t, p)
    push!(tmp.val, v)
    tmp
end
function get_valpush!(f::Function, t::AbstractNode, v, p, ps...)
    tmp = get!(f, t, p, ps...)
    push!(tmp.val, v)
    tmp
end
function get_specpush!(f::Function, t::AbstractNode, v::Pair, p)
    tmp = get!(f, t, p)
    push!(tmp.spec, v)
    tmp
end
function get_specpush!(f::Function, t::AbstractNode, v::Pair, p, ps...)
    tmp = get!(f, t, p, ps...)
    push!(tmp.spec, v)
    tmp
end
function get_svalpush!(f::Function, t::AbstractNode, v, p)
    tmp = get!(f, t, p)
    push!(tmp.sval, v)
    tmp
end
function get_svalpush!(f::Function, t::AbstractNode, v, p, ps...)
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
# Type hierarchy of functors
abstract type AbstractPathKey{T} end
# abstract type AbstractPathKey1{Function} end

abstract type AbstractIndexedPathKey{T<:Function} <: AbstractPathKey{T} end

struct LinearPathKey{T} <: AbstractIndexedPathKey{T}
    f::T
    i::Int
    # LinearPathKey{T, U}(f, i) where {T, U} = new(f, i)
end
# LinearPathKey(f::T, i::U) where {T, U} = LinearPathKey{T, U}(f, i)

struct MultipleLinearPathKey{T} <: AbstractIndexedPathKey{T}
    f::T
    i::Vector{Int}
end

struct IndexedPathKey{T, U} <: AbstractIndexedPathKey{T}
    f::T
    i::U
end
IndexedPathKey(i::Int) = IndexedPathKey(identity, i)
ipk1 = IndexedPathKey(identity, 2)
ipk2 = IndexedPathKey(x -> "2" .* x, [1, 2])
@code_warntype ipk1(a)
vi = [ipk1, ipk2]
pi = PathKeys(vi, 2)
@code_warntype pi(y, a)
# #### Attempt 1
# abstract type AbstractIndexedPathKey1{T} <: AbstractPathKey1{T} end
# struct LinearPathKey1{T} <: AbstractIndexedPathKey1{T}
#     f::T
#     i::Int
# end
# LinearPathKey1(i::Int) = LinearPathKey1(identity, i)
# (p::LinearPathKey1)(A) = p.f(getindex(A, p.i))
# ii1 = LinearPathKey1(identity, 1)
# ii3 = LinearPathKey1(identity, 3)
# ii4 = LinearPathKey1(csv, 4)
# ii2 = LinearPathKey1(x -> 25, 2)
# ii1 = LinearPathKey1(x -> 1, 2)
# ii3 = LinearPathKey1(x -> 2, 2)
# ii4 = LinearPathKey1(x -> 2, 2)
# @code_warntype ii2(a)

# abstract type AbstractPathKeys4{U<:Tuple{Vararg{T, N} where {T<:AbstractPathKey1}} where {N}} end
# struct PathKeys4{U} <: AbstractPathKeys4{U}
#     ftrs::U
#     N::Int
# end

# p4 = PathKeys4((ii1, ii3, ii4, ii2), 4)
# @code_warntype p4(a)
# @code_warntype p4(y, a)
# @benchmark p4(y, a)
# @benchmark ii1(a)

# #### Attempt 2
# abstract type AbstractPathKeys7{U<:Vector{T} where {T<:LinearPathKey1}} end
# struct PathKeys7{U} <: AbstractPathKeys7{U}
#     ftrs::U
#     N::Int
# end

# p7 = PathKeys7([ii1, ii3, ii4, ii2], 4)
# @code_warntype p7(a)
# @code_warntype p7(y, a)
# @benchmark p7(y, a)
# @benchmark ii1(a)
# @code_warntype ii2(a)

# abstract type AbstractPathKeys{T<:NTuple{<:M, <:AbstractPathKey}} end
# abstract type AbstractPathKeys{U<:Tuple{Vararg{T, N} where {T<:AbstractPathKey}} where {N}} end
# abstract type AbstractPathKeys1{U<:Tuple{Vararg{T, N} where N where {T<:AbstractPathKey}}} end
# abstract type AbstractPathKeys2{U<:Tuple{Vararg{T, N} where N} where {T<:AbstractPathKey}} end
# abstract type AbstractPathKeys3{U<:Tuple{Vararg{T, N}} where N where {T<:AbstractPathKey}} end
# Tuple{Vararg{T, N} where {T<:AbstractArray}} where N
# Tuple{Vararg{T, N} where N where {T<:AbstractArray}}
# isconcretetype(ans) # false
# isabstracttype(ans) # false
# ans{Vector, Vector}
# Tuple{Vararg{T, N} where N} where {T<:AbstractArray}
# ans{Vector{Int}} # Tuple{Vararg{Vector{Int64}}}
# isconcretetype(ans) # false
# isabstracttype(ans) # false
# Tuple{Vararg{T, N}} where N where {T<:AbstractArray}
# ans{Vector{Int}} # Tuple{Vararg{Vector{Int64}, N}} where N :same as NTuple{N, Vector{Int}} where N
# ans{Vector{Int}, 2} # Tuple{Vector{Int64}, Vector{Int64}}
# isconcretetype(ans) # true

# Revised abstract type for AbstractPathKeys
abstract type AbstractPathKeys{U<:Vector{T} where {T<:AbstractPathKey}} end
struct PathKeys{U} <: AbstractPathKeys{U}
    ftrs::U
    N::Int
    # PathKeys{T}(ftrs, N) where {T} = new(ftrs, N)
    # function PathKeys{U}(ftrs, N) where {U}
    #     new(copyto!(similar(ftrs), ftrs), N)
    # end
end
# PathKeys(ftrs::T, N) where {T} = PathKeys{T}(ftrs, N)

Base.length(p::AbstractPathKeys) = p.N
Base.size(p::AbstractPathKeys) = (p.N,)

Base.iterate(p::AbstractPathKeys, state=1) =
    state > p.N ? nothing : (p.ftrs[state], state + 1)
Base.eltype(p::AbstractPathKeys) = eltype(p.ftrs)

Base.IndexStyle(::Type{<:AbstractPathKeys}) = IndexLinear()
Base.getindex(p::AbstractPathKeys, i::Int) = getindex(p.ftrs, i)
Base.getindex(p::AbstractPathKeys, I) = [p[i] for i in I]
Base.firstindex(p::AbstractPathKeys) = 1
Base.lastindex(p::AbstractPathKeys) = p.N

function Base.isequal(p1::AbstractPathKeys, p2::AbstractPathKeys)
    p1 === p2 && return true
    length(p1) == length(p2) || return false
    for n = 1:length(p1)
        p1[n] == p2[n] || return false
    end
    return true
end
Base.:(==)(p1::AbstractPathKeys, p2::AbstractPathKeys) = isequal(p1, p2)

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

mutable struct MPathKeys{U} <: AbstractPathKeys{U}
    ftrs::U
    N::Int
end

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
i1 = LinearPathKey(1)
i3 = LinearPathKey(3)
csv(x) = x * ".csv"
i2 = LinearPathKey(x -> 25, 2)
i4 = LinearPathKey(csv, 4)
a = ["a", "b", "c", "d"]
y = Vector{Any}(undef, 4)
lpks = [i1, i2, i3, i4]
p = PathKeys(lpks, 4)

imult = MultipleLinearPathKey(identity, [1, 2])
imult2 = MultipleLinearPathKey(identity, [1, 2])
v = [i1, i4, imult, imult2]
v2 = [i1, i2, imult, imult2]
pv = PathKeys(v, 4)
pv2 = PathKeys(v2)

@code_warntype p(y, a)
@benchmark p(y, a)
@code_warntype pv(y, a)
@benchmark pv(y, a)
@code_warntype p(a)
@code_warntype pv(a)

# Testing to see if let is worthwhile -- no difference
p2 = PathKeysLet(lpks, 4)
@code_warntype p2(y, a)
@benchmark p2(y, a)

# Testing whether a direct copy is worthwhile -- no difference
p3 = PathKeys([LinearPathKey(1), LinearPathKey(x -> 25, 2),
               LinearPathKey(3), LinearPathKey(csv, 4)], 4)
@code_warntype p3(y, a)
@benchmark p3(y, a)

# Testing whether StaticArrays is worthwhile -- slower
p4 = StaticPathKeys(SVector{4}(v), 4)
@code_warntype p4(y, a)
@benchmark p4(y, a)

# Testing whether mutability matters -- it has no effect, which, in reality, might be convenient
mp = MPathKeys(lpks, 4)

mpv = MPathKeys(v, 4)

@code_warntype mp(y, a)
@benchmark mp(y, a)
@code_warntype mpv(y, a)
@benchmark mpv(y, a)

# Comparison to direct benchmark
ff(x) = 25
function test!(x, a)
    x[1] = identity(getindex(a, 1))
    x[2] = ff(getindex(a, 2))
    x[3] = identity(getindex(a, 3))
    x[4] = csv(getindex(a, 4))
    return x
end
@code_warntype test!(y, a)
@benchmark test!(y, a)

#### Comparison vs. simplification
ipk1 = IndexedPathKey(identity, 1)
ipk2 = IndexedPathKey(x -> 25, 2)
ipk3 = IndexedPathKey(identity, 3)
ipk4 = IndexedPathKey(csv, 4)
vi = [ipk1, ipk2, ipk3, ipk4]
@code_warntype ipk1(a)
pi = PathKeys(vi, length(vi))
@code_warntype pi(y, a)
@benchmark pi(y, a)
#
imult_pk1 = IndexedPathKey(identity, [1, 2])
imult_pk2 = IndexedPathKey(identity, [1, 2])
iv = [ipk1, ipk4, imult_pk1, imult_pk2]
@code_warntype imult_pk1(a)
pvi = PathKeys(iv, 4)
@code_warntype pvi(y, a)
@benchmark pvi(y, a)
@benchmark pvi(a)

# PathKeysTuple test
pkt = PathKeysTuple(iv, 4)
@code_warntype pkt(a)
@benchmark pkt(a)

#### Outer constructors: LinearPathKey
LinearPathKey(i::Int) = LinearPathKey(identity, i)

#### functor: LinearPathKey
function (p::AbstractIndexedPathKey)(A)
    p.f(getindex(A, p.i))
end

#### Outer constructors: PathKeys
PathKeys(ftrs) = PathKeys(ftrs, length(ftrs))

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
function (p::AbstractPathKeys)(x, A)
    n = 1
    N = p.N
    # N = length(p)
    fs = p.ftrs
    while n ≤ N#p.N
        # x[n] = p.ftrs[n](A)
        # x[n] = _pcall(p.ftrs[n], A) # still type-unstable
        # Realistic options: use local fs, or, access via getindex
        x[n] = fs[n](A)
        # x[n] = p[n](A)
        n += 1
    end
    return x
    # _pfunctor(p.ftrs, x, A, p.N)
end
function (p::AbstractPathKeys)(A)
    x = Vector{Any}(undef, p.N)
    p(x, A)
end

struct PathKeysTuple{T} <: AbstractPathKeys{T}
    ftrs::T
    N::Int
end

#### Outer constructors: PathKeysTuple
PathKeysTuple(ftrs) = PathKeysTuple(ftrs, length(ftrs))

#### functor: PathKeysTuple
function (p::PathKeysTuple)(A)
    Tuple(p.ftrs[n](A) for n = 1:p.N)
end

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
mat = reshape([1:200000;], (4, 50000));
pv2 = PathKeys(v2)
pn = PathKeys([i1, i3, i2, i3])
pni = PathKeys([ipk1, ipk3, ipk2, ipk3], 4)
@benchmark grow!(gf, SimpleNode(), pn, eachcol(mat))
@benchmark grow!(gf, SimpleNode(), pv2, eachcol(mat))
@benchmark grow!(gf, SimpleNode(), pni, eachcol(mat))
mat[1:2, :] .= 1;
@benchmark grow!(gf, SimpleNode(), pn, eachcol(mat))
@benchmark grow!(gf, SimpleNode(), pv2, eachcol(mat))
# very long trees -- stack overflow
pmat = reshape([1:200000;], (20, 10000));
pn2 = PathKeys([LinearPathKey(i) for i = 1:20]);
pni2 = PathKeys([IndexedPathKey(i) for i = 1:20]);
@benchmark grow!(gf, SimpleNode(), pn2, eachcol(pmat))
@benchmark grow!(gf, SimpleNode(), pni2, eachcol(pmat))

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
function foreachat!(f::Function, x::AbstractNode, N::Int, C::Int)
    if C < N - 1
        if !isempty(x)
            for p in x
                foreachat!(f, p.second, N, C + 1)
            end
        end
    elseif C == N - 1
        if !isempty(x)
            for p in x
                f(p)
            end
        end
    end
    nothing
    # A terse definition
    # isempty(x) && return nothing
    # C̃ = C + 1
    # if C̃ < N
    #     for p in x
    #         foreachat!(f, p.second, N, C̃)
    #     end
    # elseif C̃ == N
    #     for p in x
    #         f(p)
    #     end
    # end
    # nothing

end
function foreachat!(f::Function, x::AbstractNode, N::Int)
    C = 1
    foreachat!(f, x, N, C)
end

function foreachats!(f::Function, x::AbstractNode, Ns::NTuple{M, Int}) where {M}
    states = C + 1 .== Ns
    if any(states) && !isempty(x)
        for i in eachindex(states)
            if states[i]
                for v in values(x)
                    f(v)
                end
            end
        end
    else
        if !isempty(x)
            for v in values(x)
                foreachats!(f, v, Ns, C + 1)
            end
        end
    end
    nothing
    # A terse definition
    # isempty(x) && return nothing
    # C̃ = C + 1
    # for N in Ns
    #     if C̃ == N
    #         for p in x
    #             f(p.second)
    #         end
    #         # or, foreach(p -> f(p.second), x)
    #     end
    # end
    # for p in x
    #     foreachats!(f, p.second, Ns, C̃)
    # end
    # return nothing
end

function foreachat_values!(f::Function, x::AbstractNode, N::Int, C::Int)
    if C < N - 1
        if !isempty(x)
            for v in values(x)
                foreachat!(f, v, N, C + 1)
            end
        end
    elseif C == N - 1
        if !isempty(x)
            for v in values(x)
                f(v)
            end
        end
    end
    nothing
    # A terse definition
    # isempty(x) && return nothing
    # C̃ = C + 1
    # if C̃ < N
    #     for p in x
    #         foreachat!(f, p.second, N, C̃)
    #     end
    # elseif C̃ == N
    #     for p in x
    #         f(p.second)
    #     end
    # end
    # return nothing
end
function foreachat_values!(f::Function, x::AbstractNode, N::Int)
    C = 1
    foreachat_values!(f, x, N, C)
end

# Cover both data and functional filters; return an iterable of keys.
function subsetgen(ks, x::AbstractNode)
    return keys(x) ∩ ks
end
function subsetgen(f::Function, x::AbstractNode)
    return f(x)
end

function foreachat_filter!(f::Function, x::AbstractNode, N::Int, C::Int, subsets)
    if C < N - 1
        if !isempty(x)
            dKs = subsetgen(subsets[C], x)
            for k in dKs
                foreachat_filter!(f, x[k], N, C + 1, subsets)
            end
        end
    elseif C == N - 1
        if !isempty(x)
            dKs = subsetgen(subsets[C], x)
            for k in dKs
                f(x[k])
            end
        end
    end
    nothing
    # Terse version
    # isempty(x) && return nothing
    # C̃ = C + 1
    # dKs = subsetgen(subsets[C], x)
    # if C̃ < N
    #     for k in dKs
    #         foreachat_filter!(f, x[k], N, C̃, subsets)
    #     end
    # elseif C̃ == N
    #     for k in dKs
    #         f(x[k])
    #     end
    # end
    # return nothing
end

function foreachat_filter!(f::Function, x::AbstractNode, subsets)
    N = length(subsets)
    C = 1
    foreachat_filter!(f, x, N, C, subsets)
    # foreachat_filter!(f, x, length(subsets), 1, subsets)
end

####
function findallat!(f::Function, A::Vector{<:Vector}, ks::Vector, x::AbstractNode, N::Int, C::Int)
    if C < N - 1
        if !isempty(x)
            for p in x
                findallat!(f, A, setindex!(ks, p.first, C), p.second, N, C + 1)
            end
        end
    elseif C == N - 1
        if !isempty(x)
            for p in x
                setindex!(ks, p.first, C)
                f(p) && push!(A, copyto!(similar(ks), ks))
            end
        end
    end
    return A
    # Terse form
    # isempty(x) && return A
    # C̃ = C + 1 #tmp = N - 1
    # if C̃ < N # equivalent to C < N - 1
    #     for p in x
    #         @inbounds setindex!(ks, p.first, C)
    #         findallat!(f, A, ks, p.second, N, C̃)
    #     end
    # elseif C̃ == N
    #     for p in x
    #         @inbounds setindex!(ks, p.first, C)
    #         f(p) && push!(A, copyto!(similar(ks), ks))
    #     end
    # end
    # return A
end

function findallat!(f::Function, x::AbstractNode, N::Int)
    ks = Vector{Any}(undef, N - 1)
    A = Vector{Vector{Any}}()
    C = 1
    findallat!(f, A, ks, x, N, C)
end

####
function countat!(f::Function, A::Vector{Int}, x::AbstractNode, N::Int, C::Int)
    if C < N - 1
        if !isempty(x)
            for v in values(x)
                countat!(f, A, v, N, C + 1)
            end
        end
    elseif C == N - 1
        if !isempty(x)
            for p in x
                A[1] += f(p)
            end
        end
    end
    return A
    # Terse form
    # isempty(x) && return A
    # C̃ = C + 1
    # if C̃ < N
    #     for p in x
    #         countat!(f, A, p.second, N, C̃)
    #     end
    # elseif C̃ == N
    #     for p in x
    #         A[1] += f(p)
    #     end
    # end
    # return A
end

function countat(f::Function, x::AbstractNode, N::Int)
    A = Int[0]
    C = 1
    countat!(f, A, x, N, C)
end

####
function forall_depthfirst!(f::Function, x::AbstractNode)
    if !isempty(x)
        for p in x
            forall_depthfirst!(f, p.second)
        end
    end
    f(x)
    return nothing
    # Terse form
    # isempty(x) && (f(x); return nothing)
    # for p in x
    #     forall_depthfirst!(f, p.second)
    # end
    # f(x)
    # return nothing
end

function forall_breadthfirst!(f::Function, x::AbstractNode)
    f(x)
    if !isempty(x)
        for p in x
            forall_breadthfirst!(f, p.second)
        end
    end
    return x
    # Terse form
    # f(x)
    # isempty(x) && return nothing
    # for p in x
    #     forall_breadthfirst!(f, p.second)
    # end
    # return nothing
end

function forallfrom!(f::Function, x::AbstractNode, N::Int, C::Int)
    if C < N - 1
        if !isempty(x)
            for p in x
                forallfrom!(f, p.second, N, C + 1)
            end
        end
    elseif C == N - 1
        if !isempty(x)
            for p in x
                forall_breadthfirst!(f, p.second)
            end
        end
    end
    nothing
    # Terse form
    # isempty(x) && return nothing
    # C̃ = C + 1
    # if C̃ < N
    #     for p in x
    #         forallfrom!(f, p.second, N, C̃)
    #     end
    # elseif C̃ == N
    #     for p in x
    #         forall_breadthfirst!(f, p.second)
    #     end
    # end
    # return nothing
end

function forallfrom!(f::Function, x::AbstractNode, N::Int)
    C = 1
    forallfrom!(f, x, N, C)
end

function forallthrough!(f::Function, x::AbstractNode, N::Int, C::Int)
    f(x)
    isempty(x) && return nothing
    C̃ = C + 1
    if C̃ < N
        for p in x
            forallthrough!(f, p.second, N, C̃)
        end
    elseif C̃ == N
        for p in x
            f(p.second)
        end
    end
    return nothing
    # Works, but enters function for each node at last level. Might be less efficient.
    # f(x)
    # isempty(x) && return nothing
    # C̃ = C + 1
    # if C̃ ≤ N
    #     for p in x
    #         forallthrough!(f, p.second, N, C̃)
    #     end
    # end
    # return nothing
end

function forallthrough(f::Function, x::AbstractNode, N::Int)
    C = 1
    forallthrough!(f, x, N, C)
end


####
function countall!(f::Function, A::Vector{Int}, x::AbstractNode)
    if !isempty(x)
        for p in x
            countall!(f, A, p.second)
        end
    end
    A[1] += f(x)
    return A
    # Terse form
    # A[1] += f(x)
    # isempty(x) && return A
    # for p in x
    #     countall!(f, A, p.second)
    # end
    # return A
end

function countall(f::Function, x::AbstractNode)
    A = Int[0]
    countall!(f, A, x)
end

function countallfrom!(f::Function, A::Vector{Int}, x::AbstractNode, N::Int, C::Int)
    if C < N - 1
        if !isempty(x)
            for p in x
                countallfrom!(f, A, p.second, N, C + 1)
            end
        end
    elseif C == N - 1
        if !isempty(x)
            for p in x
                countall!(f, A, p.second)
            end
        end
    end
    return A
    # Terse form
    # # A[1] += f(x)
    # isempty(x) && return A
    # C̃ = C + 1
    # if C̃ < N
    #     for p in x
    #         countallfrom!(f, A, p.second, N, C̃)
    #     end
    # elseif C̃ == N
    #     for p in x
    #         countall!(f, A, p.second)
    #     end
    # end
    # return A
end

function countallfrom(f::Function, x::AbstractNode, N::Int)
    A = Int[0]
    C = 1
    countallfrom!(f, A, x, N, C)
end

function countallthrough!(f::Function, A::Vector{Int}, x::AbstractNode, N::Int, C::Int)
    A[1] += f(x)
    isempty(x) && return A
    C̃ = C + 1
    if C̃ < N
        for p in x
            countallthrough!(f, A, p.second, N, C̃)
        end
    elseif C̃ == N
        for p in x
            A[1] += f(p.second)
        end
    end
    return A
    # Works, but enters function for each node at last level. Might be less efficient.
    # A[1] += f(x)
    # isempty(x) && return A
    # C̃ = C + 1
    # if C̃ ≤ N
    #     for p in x
    #         countallthrough!(f, A, p.second, N, C̃)
    #     end
    # end
    # return A
end

function countallthrough(f::Function, x::AbstractNode, N::Int)
    A = Int[0]
    C = 1
    countallthrough!(f, A, x, N, C)
end

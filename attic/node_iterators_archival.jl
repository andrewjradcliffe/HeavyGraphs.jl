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

Base.get(x::A where {A<:AbstractNode{T, U}}, key, default) where {T, U} = get(x.link, key, default)
Base.get(f::Function, x::A where {A<:AbstractNode{T, U}}, key) where {T, U} = get(f, x.link, key)
Base.get!(x::A where {A<:AbstractNode{T, U}}, key, default) where {T, U} = get!(x.link, key, default)
Base.get!(f::Function, x::A where {A<:AbstractNode{T, U}}, key) where {T, U} = get!(f, x.link, key)

Base.getindex(x::AbstractNode, key) = getindex(x.link, key)
## Solves type-instability at cost of 50% more time and 32 bytes allocation
# function getindex2(x::A, key)::A where {T, U, A<:AbstractNode{T, U}}
#     getindex(x.link, key)
# end
# function getindex3(x::A, key) where {T, U, A<:AbstractNode{T, U}}
#     convert(A, getindex(x.link, key))
# end

Base.setindex!(x::AbstractNode, value, key) = setindex!(x.link, value, key)

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
ComplexNode(link::T, val::U, spec::V) where {T<:AbstractDict, U<:AbstractVector, V::AbstractDict} =
    ComplexNode(link, val, spec, [])
ComplexNode(link::T, val::U, sval::W) where {T<:AbstractDict, U<:AbstractVector, W::AbstractVector} =
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
# Base.convert(::Type{T}, x::T) where {T<:AbstractSimpleNode} = x
# Base.convert(::Type{T}, x::T) where {T<:AbstractComplexNode} = x

#### Promotion
promote_rule(::Type{S₁} where {S₁<:AbstractSimpleNode{T₁, U₁}},
             ::Type{S₂} where {S₂<:AbstractSimpleNode{T₂, U₂}}) where {T₁, U₁, T₂, U₂} =
                 SimpleNode{promote_type(T₁, T₂), promote_type(U₁, U₂)}
promote_rule(::Type{C₁} where {C₁<:AbstractComplexNode{T₁, U₁, V₁, W₁}},
             ::Type{C₂} where {C₂<:AbstractComplexNode{T₂, U₂, V₂, W₂}}) where {T₁, U₁, V₁, W₁, T₂, U₂, V₂, W₂} =
                 ComplexNode{promote_type(T₁, T₂), promote_type(U₁, U₂), promote_type(V₁, V₂), promote_type(W₁, W₂)}
promote_rule(::Type{S} where {S<:AbstractSimpleNode{T₁, U₁}},
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
    SimpleNode(deepcopy(convert(T₁, x.link)), copyto!(similar(T₁, x.val), x.val))

ComplexNode{T₁, U₁, V, W}(x::AbstractSimpleNode{T₂, U₂})  where {T₁, U₁, V, W} where {T₂, U₂} =
    # ComplexNode(T₂ === T₁ ? deepcopy(x.link) : deepcopy(convert(T₁, x.link)),
    #             copyto!(similar(T₁, x.val), x.val))
    ComplexNode(deepcopy(convert(T₁, x.link)), copyto!(similar(T₁, x.val), x.val))

#### More outer constructors which may be needed
SimpleNode{T, U}(x::AbstractSimpleNode) where {T, U} =
    SimpleNode(deepcopy(convert(T, x.link)), copyto!(similar(U, x.val), x.val))

ComplexNode{T, U, V, W}(x::AbstractComplexNode) where {T, U, V, W} =
    ComplexNode(deepcopy(convert(T, x.link)), copyto!(similar(U, x.val), x.val),
                deepcopy(convert(V, x.link)), copyto!(similar(W, x.val), x.val))

################ Tests: will be moved to package tests.jl when ready

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

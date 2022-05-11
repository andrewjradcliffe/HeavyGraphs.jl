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
function Base.isequal(e₁::AbstractEdge, e₂::AbstractEdge)
    e₁ === e₂ && return true
    e₁.f == e₂.f || return false
    e₁.i == e₂.i || return false
    true
end
Base.:(==)(e₁::AbstractEdge, e₂::AbstractEdge) = isequal(e₁, e₂)

#### Outer constructors for IndexedEdge
IndexedEdge(i::Int) = IndexedEdge(identity, i)
IndexedEdge(i::Vector{Int}) = IndexedEdge(identity, i)
IndexedEdge(i::UnitRange{Int}) = IndexedEdge(identity, i)
IndexedEdge(i::CartesianIndex{N}) where {N} = IndexedEdge(identity, i)

#### functor: default behavior for all IndexedEdge
function (p::AbstractIndexedEdge)(A)
    @inbounds p.f(getindex(A, p.i))
end

# ################################################################
# # Revised abstract type for AbstractEdges3
# # abstract type AbstractEdges{U<:Vector{T} where {T<:AbstractEdge}} end
# abstract type AbstractEdges{U<:Tuple{Vararg{T, N} where {T<:AbstractEdge}} where {N}} end
# struct Edges{U} <: AbstractEdges{U}
#     ftrs::U
#     N::Int
#     # Edges{T}(ftrs, N) where {T} = new(ftrs, N)
#     # function Edges{U}(ftrs, N) where {U}
#     #     new(copyto!(similar(ftrs), ftrs), N)
#     # end
#     # function Edges{U}(ftrs, N) where {U}
#     #     let ftrs = copyto!(similar(ftrs), ftrs), N = N
#     #         new(ftrs, N)
#     #     end
#     # end
# end
# # Edges(ftrs::T, N) where {T} = Edges{T}(ftrs, N)


# #### Interface for AbstractEdges
# Base.length(p::AbstractEdges) = p.N
# Base.size(p::AbstractEdges) = (p.N,)

# Base.iterate(p::AbstractEdges, state=1) = state > p.N ? nothing : (p.ftrs[state], state + 1)
# Base.eltype(p::AbstractEdges) = eltype(p.ftrs)

# Base.IndexStyle(::Type{<:AbstractEdges}) = IndexLinear()
# Base.getindex(p::AbstractEdges, i::Int) = getindex(p.ftrs, i)
# Base.getindex(p::AbstractEdges, I) = Tuple(p[i] for i ∈ I)
# Base.getindex(p::AbstractEdges, ::Colon) = p.ftrs
# Base.firstindex(p::AbstractEdges) = 1
# Base.lastindex(p::AbstractEdges) = p.N

# function Base.isequal(p₁::AbstractEdges, p₂::AbstractEdges)
#     p₁ === p₂ && return true
#     length(p₁) == length(p₂) || return false
#     for n = 1:length(p₁)
#         p₁[n] == p₂[n] || return false
#     end
#     return true
# end
# Base.:(==)(p₁::AbstractEdges, p₂::AbstractEdges) = isequal(p₁, p₂)

# #### Outer constructors: Edges
# Edges(ftrs) = Edges(ftrs, length(ftrs))
# Edges(gen::T) where {T<:Base.Generator} = Edges(Tuple(gen))
# Edges(e::AbstractEdge) = Edges((e,))
# Edges(es::Vararg{AbstractEdge, N}) where {N} = Edges(es)

# #### functor: Edges
# function (p::AbstractEdges)(x, A)
#     # n = 1
#     # N = p.N
#     # fs = p.ftrs
#     # while n ≤ N
#     #     @inbounds x[n] = fs[n](A)
#     #     n += 1
#     # end
#     # 2021-10-19 -- Interesting: the code below is slower, but it revealed an interesting
#     # aspect: if A is a vector rather than a view, the code is non-allocating,
#     # and faster by a factor of 2. This suggests that all iterables passed
#     # to the growth algorithms should actually just be Vector{Vector}.
#     # Consider! for 1.2M columns, the 480bytes per column adds up to 0.583 GB
#     # which could be eliminated. In fact, the Arrow format lends itself to this
#     # as the native form is not a matrix but a vector of vectors. Overall,
#     # this represents another 27% reduction in memory for the worst-case.
#     # 2021-10-20 -- Addendum: testing on Ints, with 1 out of 10 being converted to string,
#     # the speed gain is approximately abs(209.414 - 219.528) / 219.528 ≡ 0.0461
#     # and memory reduction is: abs(139.10 - 148.25) / 148.25 ≡ 0.0617
#     # Thus, expected to be about 5% speed and 6% memory per path key.
#     # Still, it appears to be worthwhile, especially given the fact that
#     # the natural way to store the flat files as Arrow is to store/load as
#     # a Vector{Vector}, thus, one can save on memory allocations there.
#     # setindex!.(Ref(x), Ref(A) .|> p.ftrs, 1:p.N)
#     _apk!(x, A, p.ftrs)
#     return x
# end

# function (p::AbstractEdges)(A)
#     x = Vector{Any}(undef, p.N)
#     p(x, A)
# end

# @generated function _apk!(x, A, ftrs::Tuple{Vararg{T, N} where {T<:AbstractEdge}}) where {N}
#     quote
#         Base.Cartesian.@nexprs $N i -> x[i] = ftrs[i](A)
#     end
# end

############################################################################################
############################################################################################
############################################################################################

# #### 2022-01-11: an experiment which embodied the original idea of the AbstractEdges
# # functor, finally realized via correct metaprogramming. Efficiency gains
# # for an N=4 AbstractEdges are: ≈15% speed, ≈20-25% memory.
# # The only change that needs to occur is that setindex! must be forbidden
# abstract type AbstractEdges4{U<:Tuple{Vararg{T, N} where {T<:AbstractEdge}} where {N}} end
# # abstract type AbstractEdges10{U<:Tuple{Vararg{T, N} where {T<:AbstractEdge}} where {N}, N} end
# # struct Edges10{U, N} <: AbstractEdges10{U, N}
# #     ftrs::U
# # end
# # Edges9(ftrs::Tuple{Vararg{U, N}}) where {U<:AbstractEdge, N} = Edges9{U, N}(ftrs, N)
# # Edges10(ftrs::Tuple{Vararg{U, N}}) where {N} where {U<:AbstractEdge} =
# #     Edges10{typeof(ftrs), N}(ftrs)
# # Edges10(ftrs::Tuple{Vararg{U, N} where {U<:AbstractEdge}}) where {N} =
# #     Edges10{typeof(ftrs), N}(ftrs)
# struct Edges4{U} <: AbstractEdges4{U}
#     ftrs::U
#     N::Int
#     # Edges{T}(ftrs, N) where {T} = new(ftrs, N)
#     # function Edges{U}(ftrs, N) where {U}
#     #     new(copyto!(similar(ftrs), ftrs), N)
#     # end
#     # function Edges{U}(ftrs, N) where {U}
#     #     let ftrs = copyto!(similar(ftrs), ftrs), N = N
#     #         new(ftrs, N)
#     #     end
#     # end
# end
# # Edges(ftrs::T, N) where {T} = Edges{T}(ftrs, N)

# #### Interface for AbstractEdges4
# Base.length(p::AbstractEdges4) = p.N
# Base.size(p::AbstractEdges4) = (p.N,)

# Base.iterate(p::AbstractEdges4, state=1) = state > p.N ? nothing : (p.ftrs[state], state + 1)
# Base.eltype(p::AbstractEdges4) = eltype(p.ftrs)

# Base.IndexStyle(::Type{<:AbstractEdges4}) = IndexLinear()
# Base.getindex(p::AbstractEdges4, i::Int) = getindex(p.ftrs, i)
# Base.getindex(p::AbstractEdges4, I) = [p[i] for i ∈ I]
# Base.getindex(p::AbstractEdges4, ::Colon) = p.ftrs
# Base.firstindex(p::AbstractEdges4) = 1
# Base.lastindex(p::AbstractEdges4) = p.N
# # Base.setindex!(p::AbstractEdges4, v, i) = setindex!(p.ftrs, v, i)

# function Base.isequal(p₁::AbstractEdges4, p₂::AbstractEdges4)
#     p₁ === p₂ && return true
#     length(p₁) == length(p₂) || return false
#     for n = 1:length(p₁)
#         p₁[n] == p₂[n] || return false
#     end
#     return true
# end
# Base.:(==)(p₁::AbstractEdges4, p₂::AbstractEdges4) = isequal(p₁, p₂)

# #### Outer constructors: Edges
# Edges4(ftrs) = Edges4(ftrs, length(ftrs))

# #### functor: Edges
# function (p::AbstractEdges4)(x, A)
#     # n = 1
#     # N = p.N
#     # fs = p.ftrs
#     # while n ≤ N
#     #     @inbounds x[n] = fs[n](A)
#     #     n += 1
#     # end
#     # 2021-10-19 -- Interesting: the code below is slower, but it revealed an interesting
#     # aspect: if A is a vector rather than a view, the code is non-allocating,
#     # and faster by a factor of 2. This suggests that all iterables passed
#     # to the growth algorithms should actually just be Vector{Vector}.
#     # Consider! for 1.2M columns, the 480bytes per column adds up to 0.583 GB
#     # which could be eliminated. In fact, the Arrow format lends itself to this
#     # as the native form is not a matrix but a vector of vectors. Overall,
#     # this represents another 27% reduction in memory for the worst-case.
#     # 2021-10-20 -- Addendum: testing on Ints, with 1 out of 10 being converted to string,
#     # the speed gain is approximately abs(209.414 - 219.528) / 219.528 ≡ 0.0461
#     # and memory reduction is: abs(139.10 - 148.25) / 148.25 ≡ 0.0617
#     # Thus, expected to be about 5% speed and 6% memory per path key.
#     # Still, it appears to be worthwhile, especially given the fact that
#     # the natural way to store the flat files as Arrow is to store/load as
#     # a Vector{Vector}, thus, one can save on memory allocations there.
#     # setindex!.(Ref(x), Ref(A) .|> p.ftrs, 1:p.N)
#     _apk!(x, A, p.ftrs)
#     return x
# end

# function (p::AbstractEdges4)(A)
#     x = Vector{Any}(undef, p.N)
#     p(x, A)
# end

# @generated function _apk!(x, A, ftrs::Tuple{Vararg{T, N} where {T<:AbstractEdge}}) where {N}
#     quote
#         Base.Cartesian.@nexprs $N i -> x[i] = ftrs[i](A)
#     end
# end

############################################################################################
#### 2022-03-28: a revision to make N known at compile time
# As it happens, this makes effectively no difference in terms of performance based
# on tests of the relevant functor methods. Perhaps it would make a difference
# in other generated functions that resort to ::Val(p.N), but,
# testing reveals that performance is essentially the same.

abstract type AbstractEdges{U, N} end
struct Edges{U<:Tuple{Vararg{T, N} where {T<:AbstractEdge}} where {N}, N} <: AbstractEdges{U, N}
    ftrs::U
end
Edges(ftrs::U) where {U<:Tuple{Vararg{T, N} where {T<:AbstractEdge}}} where {N} =
    Edges{U,N}(ftrs)

Base.length(p::AbstractEdges{U, N}) where {U, N} = N

Base.iterate(p::AbstractEdges{U, N}, state=1) where {U, N} = state > N ? nothing : (p.ftrs[state], state + 1)
Base.eltype(p::AbstractEdges) = eltype(p.ftrs)

Base.IndexStyle(::Type{<:AbstractEdges}) = IndexLinear()
Base.getindex(p::AbstractEdges, i::Int) = getindex(p.ftrs, i)
Base.getindex(p::AbstractEdges, I) = Tuple(p[i] for i ∈ I)
Base.getindex(p::AbstractEdges, ::Colon) = p.ftrs
Base.firstindex(p::AbstractEdges) = 1
Base.lastindex(p::AbstractEdges{U, N}) where {U, N} = N

function Base.isequal(p₁::AbstractEdges, p₂::AbstractEdges)
    p₁ === p₂ && return true
    length(p₁) == length(p₂) || return false
    for n = 1:length(p₁)
        p₁[n] == p₂[n] || return false
    end
    return true
end
Base.:(==)(p₁::AbstractEdges, p₂::AbstractEdges) = isequal(p₁, p₂)

#### Outer constructors: Edges
Edges(gen::T) where {T<:Base.Generator} = Edges(Tuple(gen))
Edges(e::AbstractEdge) = Edges((e,))
Edges(es::Vararg{AbstractEdge, N}) where {N} = Edges(es)

#### functor: Edges
@generated function (p::AbstractEdges{U, N})(x, A) where {U, N}
    quote
        Base.Cartesian.@nexprs $N i -> x[i] = p[i](A)
        return x
    end
end

function (p::AbstractEdges{U, N})(A) where {U, N}
    x = Vector{Any}(undef, N)
    p(x, A)
end

# #### respective datagrow variants
# @generated function datagrow2!(f::F, g::AbstractGraph, v::AbstractEdge, p::AbstractEdges2{U, N}, itr) where {F, U, N}
#     quote
#         for item ∈ itr
#             tmp = g
#             Base.Cartesian.@nexprs $N i -> tmp = get!(f, tmp, p[i](item))
#             push!(tmp.data, v(item))
#         end
#         return g
#     end
# end
# datagrow2(f::Function, v::AbstractEdge, p::AbstractEdges2) = datagrow2!(f, f(), v, p)
# datagrow2(f::Function, v::AbstractEdge, p::AbstractEdges2, itr) = datagrow2!(f, f(), v, p, itr)

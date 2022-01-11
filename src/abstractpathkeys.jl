#
# Date created: 2021-10-19
# Author: aradclif
#
#
############################################################################################
#### AbstractPathKey, AbstractPathKeys: p. 447-449, 456-458, 2021-09-15
# Type hierarchy of functors
abstract type AbstractPathKey{T} end

abstract type AbstractIndexedPathKey{T<:Function} <: AbstractPathKey{T} end

struct IndexedPathKey{T, U} <: AbstractIndexedPathKey{T}
    f::T
    i::U
    # IndexedPathKey{T, U}(f, i) where {T, U} = new(f, i)
end
# IndexedPathKey(f::T, i::U) where {T, U} = IndexedPathKey{T, U}(f, i)

#### Outer constructors for IndexedPathKey
IndexedPathKey(i::Int) = IndexedPathKey(identity, i)
IndexedPathKey(i::Vector{Int}) = IndexedPathKey(identity, i)
IndexedPathKey(i::UnitRange{Int}) = IndexedPathKey(identity, i)

#### functor: default behavior for all IndexedPathKey
function (p::AbstractIndexedPathKey)(A)
    @inbounds p.f(getindex(A, p.i))
end
################################################################
# Revised abstract type for AbstractPathKeys3
# abstract type AbstractPathKeys{U<:Vector{T} where {T<:AbstractPathKey}} end
abstract type AbstractPathKeys{U<:Tuple{Vararg{T, N} where {T<:AbstractPathKey}} where {N}} end
struct PathKeys{U} <: AbstractPathKeys{U}
    ftrs::U
    N::Int
    # PathKeys{T}(ftrs, N) where {T} = new(ftrs, N)
    # function PathKeys{U}(ftrs, N) where {U}
    #     new(copyto!(similar(ftrs), ftrs), N)
    # end
    # function PathKeys{U}(ftrs, N) where {U}
    #     let ftrs = copyto!(similar(ftrs), ftrs), N = N
    #         new(ftrs, N)
    #     end
    # end
end
# PathKeys(ftrs::T, N) where {T} = PathKeys{T}(ftrs, N)

#### Interface for AbstractPathKeys
Base.length(p::AbstractPathKeys) = p.N
Base.size(p::AbstractPathKeys) = (p.N,)

Base.iterate(p::AbstractPathKeys, state=1) = state > p.N ? nothing : (p.ftrs[state], state + 1)
Base.eltype(p::AbstractPathKeys) = eltype(p.ftrs)

Base.IndexStyle(::Type{<:AbstractPathKeys}) = IndexLinear()
Base.getindex(p::AbstractPathKeys, i::Int) = getindex(p.ftrs, i)
Base.getindex(p::AbstractPathKeys, I) = Tuple(p[i] for i ∈ I)#[p[i] for i ∈ I]
Base.getindex(p::AbstractPathKeys, ::Colon) = p.ftrs
Base.firstindex(p::AbstractPathKeys) = 1
Base.lastindex(p::AbstractPathKeys) = p.N
# Base.setindex!(p::AbstractPathKeys, v, i) = setindex!(p.ftrs, v, i)

function Base.isequal(p₁::AbstractPathKeys, p₂::AbstractPathKeys)
    p₁ === p₂ && return true
    length(p₁) == length(p₂) || return false
    for n = 1:length(p₁)
        p₁[n] == p₂[n] || return false
    end
    return true
end
Base.:(==)(p₁::AbstractPathKeys, p₂::AbstractPathKeys) = isequal(p₁, p₂)

#### Outer constructors: PathKeys
PathKeys(ftrs) = PathKeys(ftrs, length(ftrs))
PathKeys(gen::T) where {T<:Base.Generator} = PathKeys(Tuple(gen))

#### functor: PathKeys
function (p::AbstractPathKeys)(x, A)
    # n = 1
    # N = p.N
    # fs = p.ftrs
    # while n ≤ N
    #     @inbounds x[n] = fs[n](A)
    #     n += 1
    # end
    # 2021-10-19 -- Interesting: the code below is slower, but it revealed an interesting
    # aspect: if A is a vector rather than a view, the code is non-allocating,
    # and faster by a factor of 2. This suggests that all iterables passed
    # to the growth algorithms should actually just be Vector{Vector}.
    # Consider! for 1.2M columns, the 480bytes per column adds up to 0.583 GB
    # which could be eliminated. In fact, the Arrow format lends itself to this
    # as the native form is not a matrix but a vector of vectors. Overall,
    # this represents another 27% reduction in memory for the worst-case.
    # 2021-10-20 -- Addendum: testing on Ints, with 1 out of 10 being converted to string,
    # the speed gain is approximately abs(209.414 - 219.528) / 219.528 ≡ 0.0461
    # and memory reduction is: abs(139.10 - 148.25) / 148.25 ≡ 0.0617
    # Thus, expected to be about 5% speed and 6% memory per path key.
    # Still, it appears to be worthwhile, especially given the fact that
    # the natural way to store the flat files as Arrow is to store/load as
    # a Vector{Vector}, thus, one can save on memory allocations there.
    # setindex!.(Ref(x), Ref(A) .|> p.ftrs, 1:p.N)
    _apk!(x, A, p.ftrs)
    return x
end

function (p::AbstractPathKeys)(A)
    x = Vector{Any}(undef, p.N)
    p(x, A)
end

@generated function _apk!(x, A, ftrs::Tuple{Vararg{T, N} where {T<:AbstractPathKey}}) where {N}
    quote
        Base.Cartesian.@nexprs $N i -> x[i] = ftrs[i](A)
    end
end

############################################################################################
############################################################################################
############################################################################################

# #### 2022-01-11: an experiment which embodied the original idea of the AbstractPathKeys
# # functor, finally realized via correct metaprogramming. Efficiency gains
# # for an N=4 AbstractPathKeys are: ≈15% speed, ≈20-25% memory.
# # The only change that needs to occur is that setindex! must be forbidden
# abstract type AbstractPathKeys4{U<:Tuple{Vararg{T, N} where {T<:AbstractPathKey}} where {N}} end
# # abstract type AbstractPathKeys10{U<:Tuple{Vararg{T, N} where {T<:AbstractPathKey}} where {N}, N} end
# # struct PathKeys10{U, N} <: AbstractPathKeys10{U, N}
# #     ftrs::U
# # end
# # PathKeys9(ftrs::Tuple{Vararg{U, N}}) where {U<:AbstractPathKey, N} = PathKeys9{U, N}(ftrs, N)
# # PathKeys10(ftrs::Tuple{Vararg{U, N}}) where {N} where {U<:AbstractPathKey} =
# #     PathKeys10{typeof(ftrs), N}(ftrs)
# # PathKeys10(ftrs::Tuple{Vararg{U, N} where {U<:AbstractPathKey}}) where {N} =
# #     PathKeys10{typeof(ftrs), N}(ftrs)
# struct PathKeys4{U} <: AbstractPathKeys4{U}
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

# #### Interface for AbstractPathKeys4
# Base.length(p::AbstractPathKeys4) = p.N
# Base.size(p::AbstractPathKeys4) = (p.N,)

# Base.iterate(p::AbstractPathKeys4, state=1) = state > p.N ? nothing : (p.ftrs[state], state + 1)
# Base.eltype(p::AbstractPathKeys4) = eltype(p.ftrs)

# Base.IndexStyle(::Type{<:AbstractPathKeys4}) = IndexLinear()
# Base.getindex(p::AbstractPathKeys4, i::Int) = getindex(p.ftrs, i)
# Base.getindex(p::AbstractPathKeys4, I) = [p[i] for i ∈ I]
# Base.getindex(p::AbstractPathKeys4, ::Colon) = p.ftrs
# Base.firstindex(p::AbstractPathKeys4) = 1
# Base.lastindex(p::AbstractPathKeys4) = p.N
# # Base.setindex!(p::AbstractPathKeys4, v, i) = setindex!(p.ftrs, v, i)

# function Base.isequal(p₁::AbstractPathKeys4, p₂::AbstractPathKeys4)
#     p₁ === p₂ && return true
#     length(p₁) == length(p₂) || return false
#     for n = 1:length(p₁)
#         p₁[n] == p₂[n] || return false
#     end
#     return true
# end
# Base.:(==)(p₁::AbstractPathKeys4, p₂::AbstractPathKeys4) = isequal(p₁, p₂)

# #### Outer constructors: PathKeys
# PathKeys4(ftrs) = PathKeys4(ftrs, length(ftrs))

# #### functor: PathKeys
# function (p::AbstractPathKeys4)(x, A)
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

# function (p::AbstractPathKeys4)(A)
#     x = Vector{Any}(undef, p.N)
#     p(x, A)
# end

# @generated function _apk!(x, A, ftrs::Tuple{Vararg{T, N} where {T<:AbstractPathKey}}) where {N}
#     quote
#         Base.Cartesian.@nexprs $N i -> x[i] = ftrs[i](A)
#     end
# end

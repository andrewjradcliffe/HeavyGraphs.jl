#
# Date created: 2021-10-14
# Author: aradclif
#
#
############################################################################################
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
# modified 2021-10-03; see p. 480-481
"""
    foreachat(f::Function, t::AbstractNode, N::Int, C::Int)

Apply `f` to all nodes at a given level, `N`. Caller is responsible for `C + 1` ≤ `N`.
"""
function foreachat(f::Function, t::AbstractNode, N::Int, C::Int)
    isempty(t) && return nothing
    C̃ = C + 1
    if C̃ < N
        for p in t
            foreachat(f, p.second, N, C̃)
        end
    else#if C̃ == N
        for p in t
            f(p.second)
        end
    end
    return nothing
end

"""
    foreachat(f::Function, t::AbstractNode, N::Int)

Apply `f` to all nodes at a given level, `N`.
"""
foreachat(f::Function, t::AbstractNode, N::Int) = foreachat(f, t, N, 1)

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

"""
    foreachfilterat(f::Function, fs::Vector{Function}, t::AbstractNode, N::Int, C::Int)

Apply `f` to all nodes at a given level, `N`, using a filtered traversal, with
filter specified at each level by the respective element of `fs`. Caller is
responsible for `C + 1` ≤ `N`.

Call signature of `f` is: `f(t::AbstractNode)`.
Call signature of `fs[C]` is: `fs[C](p::Pair)`
"""
function foreachfilterat(f::Function, fs::Vector{Function}, t::AbstractNode, N::Int, C::Int)
    isempty(t) && return nothing
    C̃ = C + 1
    g = fs[C]
    if C̃ < N
        for p in t
            g(p) && foreachfilterat(f, fs, p.second, N, C̃)
        end
    else#if C̃ == N
        for p in t
            g(p) && f(p.second)
        end
    end
    nothing
end

foreachfilterat(f::Function, fs::Vector{Function}, t::AbstractNode, N::Int) =
    foreachfilterat(f, fs, t, N, 1)

foreachfilterat(f::Function, fs::Vector{Function}, t::AbstractNode) =
    foreachfilterat(f, fs, t, length(fs), 1)

####
# modified 2021-10-03; see p. 480-481
function findallpathsat!(f::Function, A::Vector{<:Vector}, ks::Vector, t::AbstractNode,
                         N::Int, C::Int)
    isempty(t) && return A
    C̃ = C + 1
    if C̃ < N
        for p in t
            @inbounds setindex!(ks, p.first, C)
            findallpathsat!(f, A, ks, p.second, N, C̃)
        end
    elseif C̃ == N
        for p in t
            f(p) && (setindex!(ks, p.first, C); push!(A, copyto!(similar(ks), ks)))
        end
    end
    return A
end

function findallpathsat(f::Function, t::AbstractNode, N::Int)
    ks = Vector{Any}(undef, N - 1)
    A = Vector{Vector{Any}}()
    findallpathsat!(f, A, ks, t, N, 1)
end

####


####
function forall_depthfirst(f::Function, t::AbstractNode)
    if !isempty(t)
        for p in t
            forall_depthfirst(f, p.second)
        end
    end
    f(t)
    return nothing
    # Terse form
    # isempty(t) && (f(t); return nothing)
    # for p in t
    #     forall_depthfirst!(f, p.second)
    # end
    # f(t)
    # return nothing
end

function forall_breadthfirst(f::Function, t::AbstractNode)
    f(t)
    isempty(t) && return nothing
    for p in t
        forall_breadthfirst!(f, p.second)
    end
    return nothing
end

"""
    forallfrom(f::Function, t::AbstractNode, N::Int, C::Int)

Apply `f` to all nodes starting at given level, `N`, then proceeding
recursively through any graph/trees which are present.
Caller is responsible for `C + 1` ≤ `N`.
"""
function forallfrom(f::Function, t::AbstractNode, N::Int, C::Int)
    isempty(t) && return nothing
    C̃ = C + 1
    if C̃ < N
        for p in t
            forallfrom(f, p.second, N, C̃)
        end
    else#if C̃ == N
        for p in t
            forall_breadthfirst(f, p.second)
        end
    end
    return nothing
end

"""
    forallfrom(f::Function, t::AbstractNode, N::Int)

Apply `f` to all nodes starting at given level, `N`, then proceeding
recursively through any graph/trees which are present.
"""
forallfrom(f::Function, t::AbstractNode, N::Int) = forallfrom!(f, t, N, 1)

"""
    forallthrough(f::Function, t::AbstractNode, N::Int, C::Int)

Apply `f` to all nodes up to and including those at given level, `N`, which
is the stopping point of the recursion.
Caller is responsible for `C + 1` ≤ `N`.
"""
function forallthrough(f::Function, t::AbstractNode, N::Int, C::Int)
    # f(t)
    # isempty(t) && return nothing
    # C̃ = C + 1
    # if C̃ < N
    #     for p in t
    #         forallthrough!(f, p.second, N, C̃)
    #     end
    # elseif C̃ == N
    #     for p in t
    #         f(p.second)
    #     end
    # end
    # return nothing
    # Works, but enters function for each node at last level. Might be less efficient.
    f(t)
    isempty(t) && return nothing
    C̃ = C + 1
    if C̃ ≤ N
        for p in t
            forallthrough!(f, p.second, N, C̃)
        end
    end
    return nothing
end

"""
    forallthrough(f::Function, t::AbstractNode, N::Int)

Apply `f` to all nodes up to and including those at given level, `N`, which
is the stopping point of the recursion.
"""
forallthrough(f::Function, t::AbstractNode, N::Int) = forallthrough!(f, t, N, 1)


################################################################
"""
    countat(f::Function, t::AbstractNode, N::Int, C::Int)

Count the number of elements at the given level `N`, starting at the level `C`,
for which the function `f` returns true. Caller is responsible for `C + 1` ≤ `N`.

Call signature of `f` is: `f(t::AbstractNode)`.
"""
function countat(f::Function, t::AbstractNode, N::Int, C::Int)
    s = 0
    isempty(t) && return s
    C̃ = C + 1
    if C̃ < N
        for p in t
            s += countat(f, p.second, N, C̃)
        end
    else#if C̃ == N
        for p in t
            f(p.second) && (s += 1)
        end
    end
    return s
end

"""
    countat(f::Function, t::AbstractNode, N::Int)

Count the number of elements at the given level `N`, starting at the level `C=1`,
for which the function `f` returns true.

Call signature of `f` is: `f(t::AbstractNode)`.
"""
countat(f::Function, t::AbstractNode, N::Int) = countat(f, t, N, 1)

# Much more efficient -- 16 bytes vs. 112
function countall(f::Function, t::AbstractNode)
    s = 0
    f(t) && (s +=1)
    isempty(t) && return s
    # isempty(t) && return f(t) ? 1 : 0
    for p in t
        s += countall(f, p.second)
    end
    return s
    # # Alternate phrasing to enforce interface: f(p::Pair{Any,AbstractNode})
    # s = 0
    # isempty(t) && return s
    # for p in t
    #     f(p) && (s += 1)
    #     s += countall(f, p.second)
    # end
    # return s
end
# @benchmark countall(_rettrue, t)

function countfrom(f::Function, t::AbstractNode, N::Int, C::Int)
    s = 0
    isempty(t) && return s
    C̃ = C + 1
    if C̃ < N
        for p in t
            s += countfrom(f, p.second, N, C̃)
        end
    else#if C̃ == N
        for p in t
            s += countall(f, p.second)
        end
    end
    return s
end

countfrom(f::Function, t::AbstractNode, N::Int) = countfrom(f, t, N, 1)
# @benchmark countfrom(_rettrue, t, 10)

function countthrough(f::Function, t::AbstractNode, N::Int, C::Int)
    s = 0
    f(t) && (s += 1)
    isempty(t) && return s
    C̃ = C + 1
    if C̃ ≤ N
        for p in t
            s += countthrough(f, p.second, N, C̃)
        end
    end
    return s
end

countthrough(f::Function, t::AbstractNode, N::Int) = countthrough(f, t, N, 1)
# @benchmark countthrough(_rettrue, t, 2)

function countupto(f::Function, t::AbstractNode, N::Int, C::Int)
    s = 0
    # C < N || return s # necessary to ensure safety, but technically optional
    f(t) && (s += 1)
    # s = f(t) ? 1 : 0
    isempty(t) && return s
    C̃ = C + 1
    if C̃ < N
        for p in t
            s += countupto(f, p.second, N, C̃)
        end
    end
    # # Variant 2
    # isempty(t) && return s
    # C̃ = C + 1
    # C̃ < N || return s
    # for p in t
    #     s += countupto(f, p.second, N, C̃)
    # end
    # # Variant 3
    # C̃ = C + 1
    # (C̃ ≥ N || isempty(t)) && return s
    # for p in t
    #     s += countupto(f, p.second, N, C̃)
    # end
    return s
end

countupto(f::Function, t::AbstractNode, N::Int) = countupto(f, t, N, 1)
# @benchmark countupto(_rettrue, t, 2)

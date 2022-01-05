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
# - foreach variants
# - count variants
# - countall variants
################ Design concept
#### Suffixes
# "at"      : proceed to N-1ᵗʰ level, act on Nᵗʰ
# "from"    : proceed to N-1ᵗʰ level, begin recursion which acts on all levels,
# the first of which is Nᵗʰ
# "through" : begin recursion which acts on all levels up to and including N
# "upto"    : begin recursion which acts on all levels up to but not including N.
# same as   : begin recursion which acts on all levels up to and including N-1.
# "all"     : begin recursion which acts on all levels
#### Prefixes
# "count"     : return number of nodes for which f(node)::Bool returns true
# "sum"       : return total computed by evaluating f(node) on all nodes visited by traversal
# "foreach"   : return nothing. Evaluate f(node) on all nodes visited by traversal
# "findpaths" : return Cartesian indices of nodes for which f(node)::Bool return true
################
# modified 2021-10-03; see p. 480-481
"""
    foreachat(f::Function, t::AbstractGraph, N::Int, C::Int)

Apply `f` to all nodes at a given level, `N`. Caller is responsible for `C + 1` ≤ `N`.
"""
function foreachat(f::Function, t::AbstractGraph, N::Int, C::Int)
    isempty(t) && return nothing
    C̃ = C + 1
    if C̃ < N
        for p ∈ t
            foreachat(f, p.second, N, C̃)
        end
    else#if C̃ == N
        for p ∈ t
            f(p.second)
        end
    end
    return nothing
end

"""
    foreachat(f::Function, t::AbstractGraph, N::Int)

Apply `f` to all nodes at a given level, `N`.
"""
foreachat(f::Function, t::AbstractGraph, N::Int) = foreachat(f, t, N, 1)

# function foreachats!(f::Function, x::AbstractGraph, Ns::NTuple{M, Int}) where {M}
#     states = C + 1 .== Ns
#     if any(states) && !isempty(x)
#         for i ∈ eachindex(states)
#             if states[i]
#                 for v ∈ values(x)
#                     f(v)
#                 end
#             end
#         end
#     else
#         if !isempty(x)
#             for v ∈ values(x)
#                 foreachats!(f, v, Ns, C + 1)
#             end
#         end
#     end
#     nothing
#     # A terse definition
#     # isempty(x) && return nothing
#     # C̃ = C + 1
#     # for N ∈ Ns
#     #     if C̃ == N
#     #         for p ∈ x
#     #             f(p.second)
#     #         end
#     #         # or, foreach(p -> f(p.second), x)
#     #     end
#     # end
#     # for p ∈ x
#     #     foreachats!(f, p.second, Ns, C̃)
#     # end
#     # return nothing
# end

"""
    foreachfilterat(f::Function, fs::Vector{Function}, t::AbstractGraph, N::Int, C::Int)

Apply `f` to all nodes at a given level, `N`, using a filtered traversal, with
filter specified at each level by the respective element of `fs`. Caller is
responsible for `C + 1` ≤ `N`.

Call signature of `f` is: `f(t::AbstractGraph)`.
Call signature of `fs[C]` is: `fs[C](p::Pair)`
"""
function foreachfilterat(f::Function, fs::Vector{Function}, t::AbstractGraph, N::Int, C::Int)
    isempty(t) && return nothing
    C̃ = C + 1
    g = fs[C]
    if C̃ < N
        for p ∈ t
            g(p) && foreachfilterat(f, fs, p.second, N, C̃)
        end
    else#if C̃ == N
        for p ∈ t
            g(p) && f(p.second)
        end
    end
    nothing
end

foreachfilterat(f::Function, fs::Vector{Function}, t::AbstractGraph, N::Int) =
    foreachfilterat(f, fs, t, N, 1)

foreachfilterat(f::Function, fs::Vector{Function}, t::AbstractGraph) =
    foreachfilterat(f, fs, t, length(fs), 1)

####
function foreach_depthfirst(f::Function, t::AbstractGraph)
    if !isempty(t)
        for p ∈ t
            foreach_depthfirst(f, p.second)
        end
    end
    f(t)
    return nothing
    # Terse form
    # isempty(t) && (f(t); return nothing)
    # for p ∈ t
    #     foreach_depthfirst!(f, p.second)
    # end
    # f(t)
    # return nothing
end

function foreach_breadthfirst(f::Function, t::AbstractGraph)
    f(t)
    isempty(t) && return nothing
    for p ∈ t
        foreach_breadthfirst(f, p.second)
    end
    return nothing
end

"""
    foreachfrom(f::Function, t::AbstractGraph, N::Int, C::Int)

Apply `f` to all nodes starting at given level, `N`, then proceeding
recursively through any graph/trees which are present.
Caller is responsible for `C + 1` ≤ `N`.
"""
function foreachfrom(f::Function, t::AbstractGraph, N::Int, C::Int)
    isempty(t) && return nothing
    C̃ = C + 1
    if C̃ < N
        for p ∈ t
            foreachfrom(f, p.second, N, C̃)
        end
    else#if C̃ == N
        for p ∈ t
            foreach_breadthfirst(f, p.second)
        end
    end
    return nothing
end

"""
    foreachfrom(f::Function, t::AbstractGraph, N::Int)

Apply `f` to all nodes starting at given level, `N`, then proceeding
recursively through any graph/trees which are present.
"""
foreachfrom(f::Function, t::AbstractGraph, N::Int) = foreachfrom!(f, t, N, 1)

"""
    foreachthrough(f::Function, t::AbstractGraph, N::Int, C::Int)

Apply `f` to all nodes up to and including those at given level, `N`, which
is the stopping point of the recursion.
Caller is responsible for `C + 1` ≤ `N`.
"""
function foreachthrough(f::Function, t::AbstractGraph, N::Int, C::Int)
    # f(t)
    # isempty(t) && return nothing
    # C̃ = C + 1
    # if C̃ < N
    #     for p ∈ t
    #         foreachthrough!(f, p.second, N, C̃)
    #     end
    # elseif C̃ == N
    #     for p ∈ t
    #         f(p.second)
    #     end
    # end
    # return nothing
    # Works, but enters function for each node at last level. Might be less efficient.
    f(t)
    isempty(t) && return nothing
    C̃ = C + 1
    if C̃ ≤ N
        for p ∈ t
            foreachthrough(f, p.second, N, C̃)
        end
    end
    return nothing
end

"""
    foreachthrough(f::Function, t::AbstractGraph, N::Int)

Apply `f` to all nodes up to and including those at given level, `N`, which
is the stopping point of the recursion.
"""
foreachthrough(f::Function, t::AbstractGraph, N::Int) = foreachthrough(f, t, N, 1)

"""
    foreachupto(f::Function, t::AbstractGraph, N::Int, C::Int)

Apply `f` to all nodes up to and including those at given level, `N`, which
is the stopping point of the recursion. Caller is responsible for `C + 1` ≤ `N`.
"""
foreachupto(f::Function, t::AbstractGraph, N::Int, C::Int) =
    foreachthrough(f, t, N - 1, C)

"""
    foreachupto(f::Function, t::AbstractGraph, N::Int)

Apply `f` to all nodes up to and including those at given level, `N`, which
is the stopping point of the recursion.
"""
foreachupto(f::Function, t::AbstractGraph, N::Int) =
    foreachthrough(f, t, N - 1)

################################################################
"""
    countat(f::Function, t::AbstractGraph, N::Int, C::Int)

Count the number of elements at the given level `N`, starting at the level `C`,
for which the function `f` returns true. Caller is responsible for `C + 1` ≤ `N`.

Call signature of `f` is: `f(t::AbstractGraph)`.
"""
function countat(f::Function, t::AbstractGraph, N::Int, C::Int)
    s = 0
    isempty(t) && return s
    C̃ = C + 1
    if C̃ < N
        for p ∈ t
            s += countat(f, p.second, N, C̃)
        end
    else#if C̃ == N
        for p ∈ t
            f(p.second) && (s += 1)
        end
    end
    return s
end

"""
    countat(f::Function, t::AbstractGraph, N::Int)

Count the number of elements at the given level `N`, starting at the level `C=1`,
for which the function `f` returns true.

Call signature of `f` is: `f(t::AbstractGraph)`.
"""
countat(f::Function, t::AbstractGraph, N::Int) = countat(f, t, N, 1)

# Much more efficient -- 16 bytes vs. 112
"""
    countall(f::Function, t::AbstractGraph)

Count the number of elements for which `f` returns true. Recursion runs over
all nodes on all levels.

Call signature of `f` is: `f(t::AbstractGraph)`.
"""
function countall(f::Function, t::AbstractGraph)
    s = 0
    f(t) && (s += 1)
    isempty(t) && return s
    # isempty(t) && return f(t) ? 1 : 0
    for p ∈ t
        s += countall(f, p.second)
    end
    return s
    # # Alternate phrasing to enforce interface: f(p::Pair{Any,AbstractGraph})
    # s = 0
    # isempty(t) && return s
    # for p ∈ t
    #     f(p) && (s += 1)
    #     s += countall(f, p.second)
    # end
    # return s
end

"""
    countfrom(f::Function, t::AbstractGraph, N::Int, C::Int)

Count the number of elements for which `f` returns true. Proceed to `N - 1`th level,
begin recursion which acts on all levels, the first of which is `N`.

Call signature of `f` is: `f(t::AbstractGraph)`.
"""
function countfrom(f::Function, t::AbstractGraph, N::Int, C::Int)
    s = 0
    isempty(t) && return s
    C̃ = C + 1
    if C̃ < N
        for p ∈ t
            s += countfrom(f, p.second, N, C̃)
        end
    else#if C̃ == N
        for p ∈ t
            s += countall(f, p.second)
        end
    end
    return s
end

countfrom(f::Function, t::AbstractGraph, N::Int) = countfrom(f, t, N, 1)

"""
    countthrough(f::Function, t::AbstractGraph, N::Int, C::Int)

Count the number of elements for which `f` returns true. Recursion runs on
all nodes on all levels up to and including `N`.

Call signature of `f` is: `f(t::AbstractGraph)`.
"""
function countthrough(f::Function, t::AbstractGraph, N::Int, C::Int)
    s = 0
    f(t) && (s += 1)
    isempty(t) && return s
    C̃ = C + 1
    if C̃ ≤ N
        for p ∈ t
            s += countthrough(f, p.second, N, C̃)
        end
    end
    return s
end

countthrough(f::Function, t::AbstractGraph, N::Int) = countthrough(f, t, N, 1)

"""
    countupto(f::Function, t::AbstractGraph, N::Int, C::Int)

Count the number of elements for which `f` returns true. Recursion runs on
all nodes on all levels up to _but not_ including `N` (i.e. up to and including `N - 1`).

Call signature of `f` is: `f(t::AbstractGraph)`.
"""
countupto(f::Function, t::AbstractGraph, N::Int, C::Int) = countthrough(f, t, N - 1, C)
# why N - 1 ? its the only difference in code.

# function countupto(f::Function, t::AbstractGraph, N::Int, C::Int)
#     s = 0
#     # C < N || return s # necessary to ensure safety, but technically optional
#     f(t) && (s += 1)
#     isempty(t) && return s
#     C̃ = C + 1
#     if C̃ < N
#         for p ∈ t
#             s += countupto(f, p.second, N, C̃)
#         end
#     end
#     # # Variant 2
#     # isempty(t) && return s
#     # C̃ = C + 1
#     # C̃ < N || return s
#     # for p ∈ t
#     #     s += countupto(f, p.second, N, C̃)
#     # end
#     # # Variant 3
#     # C̃ = C + 1
#     # (C̃ ≥ N || isempty(t)) && return s
#     # for p ∈ t
#     #     s += countupto(f, p.second, N, C̃)
#     # end
#     return s
# end

countupto(f::Function, t::AbstractGraph, N::Int) = countthrough(f, t, N - 1) #countupto(f, t, N, 1)

################################################################
####
# modified 2021-10-03; see p. 480-481
# Re-named from findallpathsat -> findpathsat. Obeys naming convention

"""
    findpathsat!(f::Function, A::Vector{<:Vector}, ks::Vector, t::AbstractGraph, N::Int, C::Int)

`push!` into `A` the Cartesian indices of elements at the given level `N`, starting at the
level `C`, for which `f` returns true. Call signature of `f` is: `f(t::AbstractGraph)`.

See also: [`findpathsat`](@ref), [`findpathsall!`](@ref), [`findpathsfrom!`](@ref),
[`findpathsthrough!`](@ref), [`findpathsupto!`](@ref)
"""
function findpathsat!(f::Function, A::Vector{<:Vector}, ks::Vector, t::AbstractGraph, N::Int, C::Int)
    isempty(t) && return A
    C̃ = C + 1
    if C̃ < N
        for p ∈ t
            @inbounds setindex!(ks, p.first, C)
            findpathsat!(f, A, ks, p.second, N, C̃)
        end
    elseif C̃ == N
        for p ∈ t
            f(p.second) && (setindex!(ks, p.first, C); push!(A, copyto!(similar(ks), ks)))
        end
    end
    return A
end

"""
    findpathsat(f::Function, t::AbstractGraph, N::Int)

Return the Cartesian indices of elements at the given level `N`,  starting at the
level `C=1`, for which `f` returns true. Call signature of `f` is: `f(t::AbstractGraph)`.

See also: [`findpathsat!`](@ref)
"""
function findpathsat(f::Function, t::AbstractGraph, N::Int)
    findpathsat!(f, Vector{Vector{Any}}(), Vector{Any}(undef, N - 1), t, N, 1)
end

####
# Re-named from findallpaths -> findpathsall. Obeys naming convention
"""
    findpathsall!(f::Function, A::Vector{<:Vector}, ks::Vector, t::AbstractGraph)

`push!` into `A` the Cartesian indices of elements for which `f` returns true.
Recursion runs over all nodes on all levels. Call signature of `f` is: `f(t::AbstractGraph)`.

See also: [`findpathsall`](@ref), [`findpathsat!`](@ref), [`findpathsfrom!`](@ref),
[`findpathsthrough!`](@ref), [`findpathsupto!`](@ref)
"""
function findpathsall!(f::Function, A::Vector{<:Vector}, ks::Vector, t::AbstractGraph)
    if !isempty(t)
        push!(ks, nothing)
        C = lastindex(ks)
        for p ∈ t
            f(p.second) && (setindex!(ks, p.first, C); push!(A, copyto!(similar(ks), ks)))
        end
        # then proceed recursively
        for p ∈ t
            setindex!(ks, p.first, C)
            findpathsall!(f, A, ks, p.second)
        end
        pop!(ks) # or, deleteat!(ks, C)
    end
    return A
end

"""
    findpathsall(f::Function, t::AbstractGraph)

Return a vector of the Cartesian indices of elements for which `f` returns true.
Recursion runs over all nodes on all levels. Call signature of `f` is: `f(t::AbstractGraph)`.

See also: [`findpathsall!`](@ref)
"""
function findpathsall(f::Function, t::AbstractGraph)
    findpathsall!(f, Vector{Vector{Any}}(), Any[], t)
end

"""
    findpathsfrom!(f::Function, A::Vector{<:Vector}, t::AbstractGraph, N::Int, C::Int)

`push!` into `A` the Cartesian indices of elements for which `f` returns true.
Proceed to `N - 1`th level, begin recursion which acts on all levels, the first of which is `N`.
Call signature of `f` is: `f(t::AbstractGraph)`.
"""
function findpathsfrom!(f::Function, A::Vector{<:Vector}, t::AbstractGraph, N::Int, C::Int)
    isempty(t) && return A
    C̃ = C + 1
    if C̃ < N
        for p ∈ t
            findpathsfrom!(f, A, p.second, N, C̃)
        end
    else#if C̃ == N
        for p ∈ t
            tmp = findpathsall(f, p.second)
            isempty(tmp[1]) || append!(A, tmp)
        end
    end
    return A
end

"""
    findpathsfrom(f::Function, t::AbstractGraph, N::Int)

Return a vector of the Cartesian indices of elements for which `f` returns true.
Proceed to `N - 1`th level, begin recursion which acts on all levels, the first of which is `N`.
Call signature of `f` is: `f(t::AbstractGraph)`.
"""
function findpathsfrom(f::Function, t::AbstractGraph, N::Int)
    findpathsfrom!(f, Vector{Vector{Any}}(), t, N, 1)
end

"""
    findpathsthrough!(f::Function, A::Vector{<:Vector}, ks::Vector, t::AbstractGraph, N::Int)

`push!` into `A` the Cartesian indices of elements for which `f` returns true.
Recursion runs over all nodes on all levels up to and including `N`.
Call signature of `f` is: `f(t::AbstractGraph)`.
"""
function findpathsthrough!(f::Function, A::Vector{<:Vector}, ks::Vector, t::AbstractGraph, N::Int)
    if !isempty(t) && length(ks) < N
        push!(ks, nothing)
        C = lastindex(ks)
        for p ∈ t
            f(p.second) && (setindex!(ks, p.first, C); push!(A, copyto!(similar(ks), ks)))
        end
        # then proceed recursively
        for p ∈ t
            setindex!(ks, p.first, C)
            findpathsthrough!(f, A, ks, p.second)
        end
        pop!(ks) # or, deleteat!(ks, C)
    end
    return A
end

"""
    findpathsthrough(f::Function, t::AbstractGraph, N::Int)

Return a vector of the Cartesian indices of elements for which `f` returns true.
Recursion runs over all nodes on all levels up to and including `N`.
Call signature of `f` is: `f(t::AbstractGraph)`.
"""
function findpathsthrough(f::Function, t::AbstractGraph, N::Int)
    findpathsthrough!(f, Vector{Vector{Any}}(), Any[], t, N)
end

"""
    findpathsupto!(f::Function, A::Vector{<:Vector}, ks::Vector, t::AbstractGraph, N::Int)

`push!` into `A` the Cartesian indices of elements for which `f` returns true.
Recursion runs over all nodes on all levels up to _but not_
including `N` (i.e. up to and including `N - 1`). Call signature of `f` is: `f(t::AbstractGraph)`.
"""
findpathsupto!(f::Function, A::Vector{<:Vector}, ks::Vector, t::AbstractGraph, N::Int) =
    findpathsthrough!(f, A, ks, t, N - 1)
# why N - 1 ? its the only difference in code.
# Moreover, why re-compute N - 1 every time one enters the function -- just compute it once.

"""
    findpathsupto(f::Function, t::AbstractGraph, N::Int)

Return a vector of the Cartesian indices of elements for which `f` returns true.
Recursion runs over all nodes on all levels up to _but not_ including
`N` (i.e. up to and including `N - 1`). Call signature of `f` is: `f(t::AbstractGraph)`.
"""
findpathsupto(f::Function, t::AbstractGraph, N::Int) = findpathsthrough(f, t, N - 1)

################################################################
#### 2022-01-05: natural extension of extant iterators.
"""
    allat(f::Function, t::AbstractGraph, N::Int, C::Int)

Determine whether predicate `f` returns `true` for all elements at the given level `N`,
starting at the level `C`, returning `false` as soon as the first item in `t`
for which `f` returns `false` is encountered. Caller is responsible for `C + 1` ≤ `N`.

Call signature of `f` is: `f(t::AbstractGraph)`.
"""
function allat(f::Function, t::AbstractGraph, N::Int, C::Int)
    isempty(t) && return false
    # s = true
    C̃ = C + 1
    if C̃ < N
        # for p ∈ t
        #     s &= allat(f, p.second, N, C̃)
        #     s || break
        # end
        for p ∈ t
            allat(f, p.second, N, C̃) || return false
        end
    else#if C̃ == N
        # for p ∈ t
        #     s &= f(p.second)
        #     s || break
        # end
        for p ∈ t
            f(p.second) || return false
        end
    end
    # return s
    return true
end

"""
    allat(f::Function, t::AbstractGraph, N::Int)

Determine whether predicate `f` returns `true` for all elements at the given level `N`,
starting at the level `C`, returning `false` as soon as the first item in `t`
for which `f` returns `false` is encountered.

Call signature of `f` is: `f(t::AbstractGraph)`.
"""
allat(f::Function, t::AbstractGraph, N::Int) = allat(f, t, N, 1)

"""
    allall(f::Function, t::AbstractGraph)

Determine whether predicate `f` returns `true` for all elements at all levels,
returning `false` as soon as the first item in `t` for which `f` returns `false`
is encountered.

Call signature of `f` is: `f(t::AbstractGraph)`.
"""
function allall(f::Function, t::AbstractGraph)
    f(t) || return false
    for p ∈ t
        allall(f, p.second) || return false
    end
    return true
end

"""
    allfrom(f::Function, t::AbstractGraph, N::Int, C::Int)

Determine whether predicate `f` returns `true` for all elements at all
levels after and including `N`. That is, proceed to `N - 1`th level,
begin recursion which acts on all levels, the first of which is `N`.
Caller is responsible for `C + 1` ≤ `N`.

Call signature of `f` is: `f(t::AbstractGraph)`.
"""
function allfrom(f::Function, t::AbstractGraph, N::Int, C::Int)
    C̃ = C + 1
    if C̃ < N
        for p ∈ t
            allfrom(f, p.second, N, C̃) || return false
        end
    else#if C̃ == N
        for p ∈ t
            allall(f, p.second) || return false
        end
    end
    return true
end

allfrom(f::Function, t::AbstractGraph, N::Int) = allfrom(f, t, N, 1)

"""
    allthrough(f::Function, t::AbstractGraph, N::Int, C::Int)

Determine whether predicate `f` returns `true` for all elements on all levels
up to and including `N`, returning `false` as soon as the first item in `t`
for which `f` returns `false` is encountered. Caller is responsible for `C + 1` ≤ `N`.

Call signature of `f` is: `f(t::AbstractGraph)`.
"""
function allthrough(f::Function, t::AbstractGraph, N::Int, C::Int)
    C̃ = C + 1
    f(t) || return false
    if C̃ ≤ N
        for p ∈ t
            allthrough(f, p.second, N, C̃) || return false
        end
    end
    return true
end

"""
    allthrough(f::Function, t::AbstractGraph, N::Int)

Determine whether predicate `f` returns `true` for all elements on all levels
up to and including `N`, returning `false` as soon as the first item in `t`
for which `f` returns `false` is encountered.

Call signature of `f` is: `f(t::AbstractGraph)`.
"""
allthrough(f::Function, t::AbstractGraph, N::Int) = allthrough(f, t, N, 1)

"""
    allupto(f::Function, t::AbstractGraph, N::Int, C::Int)

Determine whether predicate `f` returns `true` for all elements on all levels
up to _but not_ including `N` (i.e. up to and including `N - 1`), returning
`false` as soon as the first item in `t` for which `f` returns `false` is encountered.
Caller is responsible for `C + 1` ≤ `N`.

Call signature of `f` is: `f(t::AbstractGraph)`.
"""
allupto(f::Function, t::AbstractGraph, N::Int, C::Int) = allthrough(f, t, N - 1, C)

"""
    allupto(f::Function, t::AbstractGraph, N::Int)

Determine whether predicate `f` returns `true` for all elements on all levels
up to _but not_ including `N` (i.e. up to and including `N - 1`), returning
`false` as soon as the first item in `t` for which `f` returns `false` is encountered.

Call signature of `f` is: `f(t::AbstractGraph)`.
"""
allupto(f::Function, t::AbstractGraph, N::Int) = allthrough(f, t, N - 1)

# sg() = SimpleDiGraph()
# a = sg()
# get!.(sg, Ref(a), [1,2,3,4,5])
# push!.(getproperty.(getindex.(Ref(a), [1,2,3,4,5]), :data), 1)
# g = i -> i.data[1] == 1
# @code_warntype allat(g, a, 2, 1)
# a[5].data[1] = 1
# allat(g, a, 2)
# #
# push!(a.data, 1)
# a.data[1] = 1
# u = i -> i.data[1] == 1
# @benchmark allthrough(u, a, 2, 1)
# allupto(u, a, 2)
# #
# h = i -> i.data[1] == 2
# anyat(h, a, 2, 1)
# @benchmark anyat(h, a, 2, 1)

"""
    anyat(f::Function, t::AbstractGraph, N::Int, C::Int)

Determine whether predicate `f` returns `true` for any element at the given level `N`,
starting at the level `C`, returning `true` as soon as the first item in `t`
for which `f` returns `true` is encountered. Caller is responsible for `C + 1` ≤ `N`.

Call signature of `f` is: `f(t::AbstractGraph)`.
"""
function anyat(f::Function, t::AbstractGraph, N::Int, C::Int)
    isempty(t) && return false
    C̃ = C + 1
    if C̃ < N
        for p ∈ t
            anyat(f, p.second, N, C̃) && return true
        end
    else#if C̃ == N
        for p ∈ t
            f(p.second) && return true
        end
    end
    return false
end

"""
    anyat(f::Function, t::AbstractGraph, N::Int)

Determine whether predicate `f` returns `true` for any element at the given level `N`,
starting at the level `C`, returning `true` as soon as the first item in `t`
for which `f` returns `true` is encountered.

Call signature of `f` is: `f(t::AbstractGraph)`.
"""
anyat(f::Function, t::AbstractGraph, N::Int) = anyat(f, t, N, 1)

"""
    anyall(f::Function, t::AbstractGraph)

Determine whether predicate `f` returns `true` for any element at any level,
returning `true` as soon as the first item in `t` for which `f` returns `true`
is encountered.

Call signature of `f` is: `f(t::AbstractGraph)`.
"""
function anyall(f::Function, t::AbstractGraph)
    f(t) && return true
    for p ∈ t
        anyall(f, p.second) || return true
    end
    return false
end

"""
    anyfrom(f::Function, t::AbstractGraph, N::Int, C::Int)

Determine whether predicate `f` returns `true` for any element at any
level after and including `N`. That is, proceed to `N - 1`th level,
begin recursion which acts on all levels, the first of which is `N`.
Caller is responsible for `C + 1` ≤ `N`.

Call signature of `f` is: `f(t::AbstractGraph)`.
"""
function anyfrom(f::Function, t::AbstractGraph, N::Int, C::Int)
    C̃ = C + 1
    if C̃ < N
        for p ∈ t
            anyfrom(f, p.second, N, C̃) || return false
        end
    else#if C̃ == N
        for p ∈ t
            anyall(f, p.second) || return false
        end
    end
    return true
end

anyfrom(f::Function, t::AbstractGraph, N::Int) = anyfrom(f, t, N, 1)

"""
    anythrough(f::Function, t::AbstractGraph, N::Int, C::Int)

Determine whether predicate `f` returns `true` for any elements on any levels
up to and including `N`, returning `true` as soon as the first item in `t`
for which `f` returns `true` is encountered. Canyer is responsible for `C + 1` ≤ `N`.

Cany signature of `f` is: `f(t::AbstractGraph)`.
"""
function anythrough(f::Function, t::AbstractGraph, N::Int, C::Int)
    C̃ = C + 1
    f(t) && return true
    if C̃ ≤ N
        for p ∈ t
            anythrough(f, p.second, N, C̃) && return true
        end
    end
    return false
end

"""
    anythrough(f::Function, t::AbstractGraph, N::Int)

Determine whether predicate `f` returns `true` for any elements on any levels
up to and including `N`, returning `true` as soon as the first item in `t`
for which `f` returns `true` is encountered.

Cany signature of `f` is: `f(t::AbstractGraph)`.
"""
anythrough(f::Function, t::AbstractGraph, N::Int) = anythrough(f, t, N, 1)

"""
    anyupto(f::Function, t::AbstractGraph, N::Int, C::Int)

Determine whether predicate `f` returns `true` for any elements on any levels
up to _but not_ including `N` (i.e. up to and including `N - 1`), returning
`true` as soon as the first item in `t` for which `f` returns `true` is encountered.
Canyer is responsible for `C + 1` ≤ `N`.

Cany signature of `f` is: `f(t::AbstractGraph)`.
"""
anyupto(f::Function, t::AbstractGraph, N::Int, C::Int) = anythrough(f, t, N - 1, C)

"""
    anyupto(f::Function, t::AbstractGraph, N::Int)

Determine whether predicate `f` returns `true` for any elements on any levels
up to _but not_ including `N` (i.e. up to and including `N - 1`), returning
`true` as soon as the first item in `t` for which `f` returns `true` is encountered.

Cany signature of `f` is: `f(t::AbstractGraph)`.
"""
anyupto(f::Function, t::AbstractGraph, N::Int) = anythrough(f, t, N - 1)

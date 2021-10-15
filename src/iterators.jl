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

Apply `f` to all nodes at a given level, `N`.
"""
function foreachat(f::Function, t::AbstractNode, N::Int, C::Int)
    # if C < N - 1
    #     if !isempty(t)
    #         for p in t
    #             foreachat!(f, p.second, N, C + 1)
    #         end
    #     end
    # elseif C == N - 1
    #     if !isempty(t)
    #         for p in t
    #             f(p)
    #         end
    #     end
    # end
    # nothing
    # A terse definition
    isempty(t) && return nothing
    C̃ = C + 1
    if C̃ < N
        for p in t
            foreachat(f, p.second, N, C̃)
        end
    elseif C̃ == N
        # for p in t
        #     f(p)
        # end
        # Alternative that may lead to more inlining:
        foreach(f, t)
    end
    nothing
end

"""
    foreachat(f::Function, t::AbstractNode, N::Int)

Apply `f` to all nodes at a given level, `N`.
"""
function foreachat(f::Function, t::AbstractNode, N::Int)
    foreachat(f, t, N, 1)
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

# Not exactly the most value, since foreachat covers this, but, perhaps, there
# is minor value in providing a convenience iterator that applies f directly
# to values rather than to pairs of (K => V)
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

function foreachfilterat(f::Function, fs::Vector{Function}, t::AbstractNode, N::Int, C::Int)
    isempty(t) && return nothing
    C̃ = C + 1
    g = fs[C]
    if C̃ < N
        for p in t
            g(p) && foreachfilterat(f, fs, p.second, N, C̃)
        end
    elseif C̃ == N
        for p in t
            g(p) && f(p)
        end
    end
    nothing
end

function foreachfilterat(f::Function, fs::Vector{Function}, t::AbstractNode, N::Int)
    foreachfilterat(f, fs, t, N, 1)
end

function foreachfilterat(f::Function, fs::Vector{Function}, t::AbstractNode)
    foreachfilterat(f, fs, t, length(fs), 1)
end

####
# modified 2021-10-03; see p. 480-481
function findallpathsat!(f::Function, A::Vector{<:Vector}, ks::Vector, t::AbstractNode,
                         N::Int, C::Int)
    # if C < N - 1
    #     if !isempty(t)
    #         for p in t
    #             findallat!(f, A, setindex!(ks, p.first, C), p.second, N, C + 1)
    #         end
    #     end
    # elseif C == N - 1
    #     if !isempty(t)
    #         for p in t
    #             setindex!(ks, p.first, C)
    #             f(p) && push!(A, copyto!(similar(ks), ks))
    #         end
    #     end
    # end
    # return A
    # Terse form
    isempty(t) && return A
    C̃ = C + 1 #tmp = N - 1
    if C̃ < N # equivalent to C < N - 1
        for p in t
            @inbounds setindex!(ks, p.first, C)
            findallat!(f, A, ks, p.second, N, C̃)
        end
    elseif C̃ == N
        for p in t
            # @inbounds setindex!(ks, p.first, C)
            # f(p) && push!(A, copyto!(similar(ks), ks))
            f(p) && (setindex!(ks, p.first, C); push!(A, copyto!(similar(ks), ks)))
        end
    end
    return A
end

function findallpathsat(f::Function, t::AbstractNode, N::Int)
    ks = Vector{Any}(undef, N - 1)
    A = Vector{Vector{Any}}()
    C = 1
    findallpathsat!(f, A, ks, t, N, C)
end

####
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
#     countat!(f, A, t, N, 1)
#     A[1]
# end

# _rettrue(p) = true
# @benchmark countat(_rettrue, t, 10)
# countat(x -> first(x) % 3 == 0, t, 10)

# Technically, faster by ≈ 25%, but loses convenience of countat!
# function _countat(f::Function, tmp::Int, t::AbstractNode, N::Int, C::Int)
#     isempty(t) && return tmp
#     C̃ = C + 1
#     if C̃ < N
#         tmp2 = 0
#         for p in t
#             tmp2 += _countat(f, tmp, p.second, N, C̃)
#             # tmp += _countat(f, 0, p.second, N, C̃)
#         end
#         return tmp2
#     elseif C̃ == N
#         tmp3 = 0
#         for p in t
#             f(p) && (tmp3 += 1)
#             # f(p) && (tmp += 1)
#         end
#         return tmp3
#     end
#     tmp
# end
# _countat(f::Function, t::AbstractNode, N::Int) = _countat(f, 0, t, N, 1)

# @benchmark _countat(_rettrue, t, 10)
# # _countat(x -> first(x) % 3 == 0, t, 10)

function countat(f::Function, t::AbstractNode, N::Int, C::Int)
    tmp = 0
    isempty(t) && return tmp
    C̃ = C + 1
    if C̃ < N
        for p in t
            tmp += countat(f, p.second, N, C̃)
        end
    elseif C̃ == N
        for p in t
            f(p) && (tmp += 1)
        end
    end
    return tmp
end

countat(f::Function, t::AbstractNode, N::Int) = countat(f, t, N, 1)

####
function forall_depthfirst!(f::Function, t::AbstractNode)
    if !isempty(t)
        for p in t
            forall_depthfirst!(f, p.second)
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

function forall_breadthfirst!(f::Function, t::AbstractNode)
    f(t)
    if !isempty(t)
        for p in t
            forall_breadthfirst!(f, p.second)
        end
    end
    return t
    # Terse form
    # f(t)
    # isempty(t) && return nothing
    # for p in t
    #     forall_breadthfirst!(f, p.second)
    # end
    # return nothing
end

function forallfrom!(f::Function, t::AbstractNode, N::Int, C::Int)
    # if C < N - 1
    #     if !isempty(t)
    #         for p in t
    #             forallfrom!(f, p.second, N, C + 1)
    #         end
    #     end
    # elseif C == N - 1
    #     if !isempty(t)
    #         for p in t
    #             forall_breadthfirst!(f, p.second)
    #         end
    #     end
    # end
    # nothing
    # Terse form
    isempty(t) && return nothing
    C̃ = C + 1
    if C̃ < N
        for p in t
            forallfrom!(f, p.second, N, C̃)
        end
    elseif C̃ == N
        for p in t
            forall_breadthfirst!(f, p.second)
        end
    end
    nothing
end

function forallfrom!(f::Function, t::AbstractNode, N::Int)
    forallfrom!(f, t, N, 1)
end

function forallthrough!(f::Function, t::AbstractNode, N::Int, C::Int)
    f(t)
    isempty(t) && return nothing
    C̃ = C + 1
    if C̃ < N
        for p in t
            forallthrough!(f, p.second, N, C̃)
        end
    elseif C̃ == N
        for p in t
            f(p.second)
        end
    end
    nothing
    # Works, but enters function for each node at last level. Might be less efficient.
    # f(t)
    # isempty(t) && return nothing
    # C̃ = C + 1
    # if C̃ ≤ N
    #     for p in t
    #         forallthrough!(f, p.second, N, C̃)
    #     end
    # end
    # return nothing
end

function forallthrough(f::Function, t::AbstractNode, N::Int)
    forallthrough!(f, t, N, 1)
end

####
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

# Much more efficient -- 16 bytes vs. 112
function countall(f::Function, t::AbstractNode)
    tmp = 0
    f(t) && (tmp +=1)
    isempty(t) && return tmp
    # isempty(t) && return f(t) ? 1 : 0
    for p in t
        tmp += countall(f, p.second)
    end
    return tmp
end



function countallfrom!(f::Function, A::Vector{Int}, t::AbstractNode, N::Int, C::Int)
    # if C < N - 1
    #     if !isempty(t)
    #         for p in t
    #             countallfrom!(f, A, p.second, N, C + 1)
    #         end
    #     end
    # elseif C == N - 1
    #     if !isempty(t)
    #         for p in t
    #             countall!(f, A, p.second)
    #         end
    #     end
    # end
    # return A
    # Terse form
    # A[1] += f(t) # Not necessary as countall! will perform this, even if empty
    isempty(t) && return A
    C̃ = C + 1
    if C̃ < N
        for p in t
            countallfrom!(f, A, p.second, N, C̃)
        end
    elseif C̃ == N
        for p in t
            countall!(f, A, p.second)
        end
    end
    return A
end

function countallfrom(f::Function, t::AbstractNode, N::Int)
    A = Int[0]
    C = 1
    countallfrom!(f, A, t, N, C)
    A[1]
end

function countallthrough!(f::Function, A::Vector{Int}, t::AbstractNode, N::Int, C::Int)
    A[1] += f(t)
    isempty(t) && return A
    C̃ = C + 1
    if C̃ < N
        for p in t
            countallthrough!(f, A, p.second, N, C̃)
        end
    elseif C̃ == N
        for p in t
            A[1] += f(p.second)
        end
    end
    return A
    # Works, but enters function for each node at last level. Might be less efficient.
    # A[1] += f(t)
    # isempty(t) && return A
    # C̃ = C + 1
    # if C̃ ≤ N
    #     for p in t
    #         countallthrough!(f, A, p.second, N, C̃)
    #     end
    # end
    # return A
end

function countallthrough(f::Function, t::AbstractNode, N::Int)
    A = Int[0]
    C = 1
    countallthrough!(f, A, t, N, C)
    A[1]
end

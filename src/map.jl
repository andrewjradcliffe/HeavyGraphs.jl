#
# Date created: 2021-10-12
# Author: aradclif
#
#
############################################################################################
#### analyzeallupto!, analyzeat! using subset functions, p.500-507, 2021-10-12
# Alas, it seems to be slightly faster to use Arrays for access,
# but for subsequent operations, e.g. a .+ b, Tuple{Vararg{Array}} is much faster,
# or, quantitatively, about 12x faster.
# Note: Tuple{Vararg{Array{T} where T}} ≡ Tuple{Vararg{Array{T}} where T}

#### A draft of analyzeat!, p. 483, 2021-10-03; p. 500-507, 2021-10-12
"""
    mapat!(f::Function, dest::Tuple{Vararg{Array{T}} where T}, t::AbstractNode, N::Int, C::Int)

Transform the graph/tree `t` by applying `f` to each element at a given level, `N`.
Results are stored in `dest`; it is assumed that the number and dimensions of arrays
are sufficient to store the results of calling `f`.

Call signature of `f` is: `f(dest, p::Pair)`.

See also: [`mapat`](@ref), [`mapfilterat`](@ref), [`mapfilterat!`](@ref)

"""
function mapat!(f::Function, dest::Tuple{Vararg{Array{T}} where T}, t::AbstractNode, N::Int, C::Int)
    isempty(t) && return dest
    C̃ = C + 1
    if C̃ < N
        for p in t
            mapat!(f, dest, p.second, N, C̃)
        end
    elseif C̃ == N
        for p in t
            f(dest, p)
        end
    end
    dest
end

"""
    mapat(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S}, t::AbstractNode, N::Int)

Transform the graph/tree `t` by applying `f` to each element at a given level, `N`.
Results are stored in a tuple of zero-initialized arrays, the number and dimension of which
are specified by `dims`.

Call signature of `f` is: `f(dest, p::Pair)`.

See also: [`mapat!`](@ref), [`mapfilterat`](@ref), [`mapfilterat!`](@ref)

"""
function mapat(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S}, t::AbstractNode, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    mapat!(f, dest, t, N, 1)
end

#### filter
"""
    mapfilterat!(f::Function, fs::Vector{Function}, dest::Tuple{Vararg{Array{T}} where T},
                 t::AbstractNode, N::Int, C::Int)

Performs a filtered traversal in which a subset is formed at each level from
the corresponding element (note: linear indexed) of the filter functions `fs`.

Call signature of `f` is: `f(dest, p::Pair)`.
Call signature of `fs[C]` is: `fs[C](p::Pair)`.

See also: [`mapfilterat`](@ref), [`mapat`](@ref), [`mapat!`](@ref)
"""
function mapfilterat!(f::Function, fs::Vector{Function}, dest::Tuple{Vararg{Array{T}} where T},
                      t::AbstractNode, N::Int, C::Int)
    isempty(t) && return dest
    C̃ = C + 1
    g = fs[C]
    if C̃ < N
        for p in t
            g(p) && mapfilterat!(f, fs, dest, p.second, N, C̃)
        end
    elseif C̃ == N
        for p in t
            g(p) && f(dest, p)
        end
    end
    dest
end

"""
    mapfilterat(f::Function, fs::Vector{Function}, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                t::AbstractNode, N::Int)

Performs a filtered traversal in which a subset is formed at each level from
the corresponding element (note: linear indexed) of the filter functions `fs`.

Call signature of `f` is: `f(dest, p::Pair)`.
Call signature of `fs[C]` is: `fs[C](p::Pair)`.

See also: [`mapfilterat!`](@ref), [`mapat`](@ref), [`mapat!`](@ref)
"""
function mapfilterat(f::Function, fs::Vector{Function}, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                     t::AbstractNode, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    mapfilterat!(f, fs, dest, t, N, 1)
end

################ reduction of vector of into single dest
function mapat!(f::Function, dest::Tuple{Vararg{Array{T}} where T},
                ts::Vector{<:AbstractNode}, N::Int, C::Int)
    for t in ts
        mapat!(f, dest, t, N, C)
    end
    dest
end

function mapat(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S},
               ts::Vector{<:AbstractNode}, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    mapat!(f, dest, ts, N, 1)
end

function tmapat(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                ts::Vector{<:AbstractNode}, L::Int)
    N = length(ts)
    M = Threads.threads()
    ranges = equalranges(N, M)
    A = Vector{Tuple{(Array{Int, length(d)} for d in dims)...}}(undef, M)
    Threads.@threads for m = 1:M
        A[m] = mapat(f, dims, ts[ranges[m]], L)
    end
    return A
end

#### filter
function mapfilterat!(f::Function, fs::Vector{Function}, dest::Tuple{Vararg{Array{T}} where T},
                      ts::Vector{<:AbstractNode}, N::Int, C::Int)
    for t in ts
        mapfilterat!(f, fs, dest, t, N, C)
    end
    dest
end

function mapfilterat(f::Function, fs::Vector{Function}, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                     ts::Vector{<:AbstractNode}, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    mapfilterat!(f, fs, dest, ts, N, 1)
end

function tmapfilterat(f::Function, fs::Vector{Function}, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                      ts::Vector{<:AbstractNode}, L::Int)
    N = length(ts)
    M = Threads.threads()
    ranges = equalranges(N, M)
    A = Vector{Tuple{(Array{Int, length(d)} for d in dims)...}}(undef, M)
    Threads.@threads for m = 1:M
        A[m] = mapfilterat(f, fs, dims, ts[ranges[m]], L)
    end
    return A
end

################################################################

"""
    mapupto!(f::Function, dest::Tuple{Vararg{Array{T}} where T}, t::AbstractNode, N::Int)

Transform the graph/tree `t` by applying `f` to each element on each level up to a
given level `N`. Results are stored in `dest`; it is assumed that the number and
dimensions of arrays are sufficient to store the results of calling `f`.

Call signature of `f` is: `f(dest, t::AbstractNode)`.

See also: [`mapupto`](@ref), [`mapfilterupto`](@ref), [`mapfilterupto!`](@ref)
"""
function mapupto!(f::Function, dest::Tuple{Vararg{Array{T}} where T}, t::AbstractNode, N::Int, C::Int)
    C̃ = C + 1
    C̃ < N || return dest
    f(dest, t)
    isempty(t) && return dest
    for p in t
        mapupto!(f, dest, p.second, N, C̃)
    end
    dest
end

"""
    mapupto(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S}, t::AbstractNode, N::Int)

Transform the graph/tree `t` by applying `f` to each element on each level up to a
given level `N`. Results are stored in a tuple of zero-initialized arrays, the number
and dimension of which are specified by `dims`.

Call signature of `f` is: `f(dest, t::AbstractNode)`.

See also: [`mapupto!`](@ref), [`mapfilterupto`](@ref), [`mapfilterupto!`](@ref)
"""
function mapupto(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S}, t::AbstractNode, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    mapupto!(f, dest, t, N, 1)
end

#### filter
"""
    mapfilterupto!(f::Function, fs::Vector{Function}, dest::Tuple{Variant{Array{T}} where T},
                   t::AbstractNode, N::Int, C::Int)

Performs a filtered traversal in which a subset is formed at each level from
the corresponding element (note: linear indexed) of the filter functions `fs`.

Call signature of `f` is: `f(dest, t::AbstractNode)`.
Call signature of `fs[C]` is: `fs[C](p::Pair)`.

See also: [`mapfilterupto`](@ref), [`mapupto`](@ref), [`mapupto!`](@ref)
"""
function mapfilterupto!(f::Function, fs::Vector{Function}, dest::Tuple{Vararg{Array{T}} where T},
                        t::AbstractNode, N::Int, C::Int)
    C̃ = C + 1
    C̃ < N || return dest
    f(dest, t)
    isempty(t) && return dest
    g = fs[C]
    for p in t
        g(p) && mapfilterupto!(f, fs, dest, p.second, N, C̃)
    end
    dest
end

"""
    mapfilterupto(f::Function, fs::Vector{Function}, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                  t::AbstractNode, N::Int)

Performs a filtered traversal in which a subset is formed at each level from
the corresponding element (note: linear indexed) of the filter functions `fs`.

Call signature of `f` is: `f(dest, t::AbstractNode)`.
Call signature of `fs[C]` is: `fs[C](p::Pair)`.

See also: [`mapfilterupto!`](@ref), [`mapupto`](@ref), [`mapupto!`](@ref)
"""
function mapfilterupto(f::Function, fs::Vector{Function},
                       dims::Tuple{Vararg{NTuple{S, Int}} where S}, t::AbstractNode, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    mapfilterupto!(f, fs, dest, t, N, 1)
end

#### levs_ks
"""
    mapupto!(f::Function, dest::Tuple{Vararg{Array{T}} where T},
             t::AbstractNode, N::Int, C::Int, levs_ks::Vector{Vector{Any}})

Provides the given level `N`, the current level `C`, and the vector of
level-respective index sets, `levs_ks` as arguments to `f`.

Call signature of `f` is: `f(dest, t::AbstractNode, N, C, levs_ks)`.
"""
function mapupto!(f::Function, dest::Tuple{Vararg{Array{T}} where T},
                  t::AbstractNode, N::Int, C::Int, levs_ks::Vector{Vector{Any}})
    C̃ = C + 1
    C̃ < N || return dest
    f(dest, t, N, C, levs_ks)
    isempty(t) && return dest
    for p in t
        mapupto!(f, dest, p.second, N, C̃, levs_ks)
    end
    dest
end

"""
    mapupto(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S},
            t::AbstractNode, N::Int, levs_ks::Vector{Vector{Any}})

Provides the given level `N`, the current level `C`, and the vector of
level-respective index sets, `levs_ks` as arguments to `f`.

Call signature of `f` is: `f(dest, t::AbstractNode, N, C, levs_ks)`.
"""
function mapupto(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                 t::AbstractNode, N::Int, levs_ks::Vector{Vector{Any}})
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    mapupto!(f, dest, t, N, 1, levs_ks)
end

#### filter and levs_ks
"""
    mapfilterupto!(f::Function, fs::Vector{Function}, dest::Tuple{Vararg{Array{T}} where T},
                   t::AbstractNode, N::Int, C::Int, levs_ks::Vector{Vector{Any}})

Performs filtered traversal; provides the given level `N`,
the current level `C`, and the vector of level-respective index sets,
`levs_ks` as arguments to `f`.

Call signature of `f` is: `f(dest, t::AbstractNode, N, C, levs_ks)`.
Call signature of `fs[C]` is: `fs[C](p::Pair)`.
"""
function mapfilterupto!(f::Function, fs::Vector{Function}, dest::Tuple{Vararg{Array{T}} where T},
                        t::AbstractNode, N::Int, C::Int, levs_ks::Vector{Vector{Any}})
    C̃ = C + 1
    C̃ < N || return dest
    f(dest, t, N, C, levs_ks)
    isempty(t) && return dest
    g = fs[C]
    for p in t
        g(p) && mapfilterupto!(f, fs, dest, p.second, N, C̃, levs_ks)
    end
    dest
end

"""
    mapfilterupto(f::Function, fs::Vector{Function}, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                  t::AbstractNode, N::Int, levs_ks::Vector{Vector{Any}})

Performs filtered traversal; provides the given level `N`,
the current level `C`, and the vector of level-respective index sets,
`levs_ks` as arguments to `f`.

Call signature of `f` is: `f(dest, t::AbstractNode, N, C, levs_ks)`.
Call signature of `fs[C]` is: `fs[C](p::Pair)`.
"""
function mapfilterupto(f::Function, fs::Vector{Function},
                       dims::Tuple{Vararg{NTuple{S, Int}} where S}, t::AbstractNode,
                       N::Int, levs_ks::Vector{Vector{Any}})
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    mapfilterupto!(f, fs, dest, t, N, 1, levs_ks)
end

################ reduction of vector of into single dest
function mapupto!(f::Function, dest::Tuple{Vararg{Array{T}} where T},
                  ts::Vector{<:AbstractNode}, N::Int, C::Int)
    for t in ts
        mapupto!(f, dest, t, N, C)
    end
    dest
end

function mapupto(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                 ts::Vector{<:AbstractNode}, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    mapupto!(f, dest, ts, N, 1)
end

function tmapupto(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                  ts::Vector{<:AbstractNode}, L::Int)
    N = length(ts)
    M = Threads.threads()
    ranges = equalranges(N, M)
    A = Vector{Tuple{(Array{Int, length(d)} for d in dims)...}}(undef, M)
    Threads.@threads for m = 1:M
        A[m] = mapupto(f, dims, ts[ranges[m]], L)
    end
    return A
end

#### filter
function mapfilterupto!(f::Function, fs::Vector{Function}, dest::Tuple{Vararg{Array{T}} where T},
                        ts::Vector{<:AbstractNode}, N::Int, C::Int)
    for t in ts
        mapfilterupto!(f, fs, dest, t, N, C)
    end
    dest
end

function mapfilterupto(f::Function, fs::Vector{Function},
                       dims::Tuple{Vararg{NTuple{S, Int}} where S},
                       ts::Vector{<:AbstractNode}, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    mapfilterupto!(f, fs, dest, ts, N, 1)
end

function tmapfilterupto(f::Function, fs::Vector{Function},
                        dims::Tuple{Vararg{NTuple{S, Int}} where S},
                        ts::Vector{<:AbstractNode}, L::Int)
    N = length(ts)
    M = Threads.threads()
    ranges = equalranges(N, M)
    A = Vector{Tuple{(Array{Int, length(d)} for d in dims)...}}(undef, M)
    Threads.@threads for m = 1:M
        A[m] = mapfilterupto(f, fs, dims, ts[ranges[m]], L)
    end
    return A
end

#### levs_ks

function mapupto!(f::Function, dest::Tuple{Vararg{Array{T}} where T},
                  ts::Vector{<:AbstractNode}, N::Int, C::Int, levs_kss::Vector{Vector{Vector{Any}}})
    for i in eachindex(ts)
        mapupto!(f, dest, ts[i], N, C, levs_kss[i])
    end
    dest
end

function mapupto(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                 ts::Vector{<:AbstractNode}, N::Int, levs_kss::Vector{Vector{Vector{Any}}})
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    mapupto!(f, dest, ts, N, 1, levs_kss)
end

function tmapupto(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                  ts::Vector{<:AbstractNode}, L::Int, levs_kss::Vector{Vector{Vector{Any}}})
    N = length(ts)
    M = Threads.threads()
    ranges = equalranges(N, M)
    A = Vector{Tuple{(Array{Int, length(d)} for d in dims)...}}(undef, M)
    Threads.@threads for m = 1:M
        A[m] = mapupto(f, dims, ts[ranges[m]], L, levs_kss[ranges[m]])
    end
    return A
end

#### filter and levs_ks
function mapfilterupto!(f::Function, fs::Vector{Function}, dest::Tuple{Vararg{Array{T}} where T},
                        ts::Vector{<:AbstractNode}, N::Int, C::Int,
                        levs_kss::Vector{Vector{Vector{Any}}})
    for i in eachindex(ts)
        mapfilterupto!(f, fs, dest, ts[i], N, C, levs_kss[i])
    end
    dest
end

function mapfilterupto(f::Function, fs::Vector{Function},
                       dims::Tuple{Vararg{NTuple{S, Int}} where S},
                       ts::Vector{<:AbstractNode}, N::Int, levs_kss::Vector{Vector{Vector{Any}}})
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    mapfilterupto!(f, fs, dest, ts, N, 1, levs_kss)
end

function tmapfilterupto(f::Function, fs::Vector{Function},
                        dims::Tuple{Vararg{NTuple{S, Int}} where S},
                        ts::Vector{<:AbstractNode}, L::Int, levs_kss::Vector{Vector{Vector{Any}}})
    N = length(ts)
    M = Threads.threads()
    ranges = equalranges(N, M)
    A = Vector{Tuple{(Array{Int, length(d)} for d in dims)...}}(undef, M)
    Threads.@threads for m = 1:M
        A[m] = mapfilterupto(f, fs, dims, ts[ranges[m]], L, levs_kss[ranges[m]])
    end
    return A
end

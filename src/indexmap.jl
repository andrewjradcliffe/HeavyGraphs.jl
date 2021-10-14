#
# Date created: 2021-10-13
# Author: aradclif
#
#
############################################################################################
"""
    indexmapat!(f::Function, dest::Tuple{Vararg{Array{T}} where T}, ks::Vector,
                t::AbstractNode, N::Int, C::Int)

Analogy to `mapat!` which tracks Cartesian key-index of iteration, passing
it as an additional argument to `f`.

Call signature of `f` is: `f(dest, ks, p::Pair)`.

See also: [`mapat!`](@ref), [`indexmapat`](@ref), [`indexmapfilterat`](@ref), [`indexmapfilterat!`](@ref)

"""
function indexmapat!(f::Function, dest::Tuple{Vararg{Array{T}} where T}, ks::Vector,
                     t::AbstractNode, N::Int, C::Int)
    isempty(t) && return dest
    C̃ = C + 1
    if C̃ < N
        for p in t
            setindex!(ks, p.first, C)
            indexmapat!(f, dest, ks, p.second, N, C̃)
        end
    elseif C̃ == N
        for p in t
            setindex!(ks, p.first, C)
            f(dest, ks, p)
        end
    end
    dest
end

"""
    indexmapat(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S}, t::AbstractNode, N::Int)

Analogy to `mapat` which tracks Cartesian key-index of iteration, passing
it as an additional argument to `f`.

Call signature of `f` is: `f(dest, ks, p::Pair)`.

See also: [`mapat`](@ref), [`indexmapat`](@ref), [`indexmapfilterat`](@ref), [`indexmapfilterat!`](@ref)

"""
function indexmapat(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                    t::AbstractNode, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    ks = Vector{Any}(undef, N - 1)
    indexmapat!(f, dest, ks, t, N, 1)
end

#### filter
"""
    indexmapfilterat!(f::Function, fs::Vector{Function}, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                      ks::Vector, t::AbstractNode, N::Int, C::Int)

Indexing analogy to `mapfilterat!` with filtered traversal.

Call signature of `f` is: `f(dest, ks, p::Pair)`.
Call signature of `fs[C]` is: fs[C](p::Pair).

See also: [`mapfilterat!`](@ref), [`indexmapat`](@ref), [`indexmapat!`](@ref)
"""
function indexmapfilterat!(f::Function, fs::Vector{Function}, dest::Tuple{Vararg{Array{T}} where T},
                           ks::Vector, t::AbstractNode, N::Int, C::Int)
    isempty(t) && return dest
    C̃ = C + 1
    g = fs[C]
    if C̃ < N
        for p in t
            g(p) && (setindex!(ks, p.first, C); indexmapfilterat!(f, fs, dest, p.second, N, C̃))
        end
    elseif C̃ == N
        for p in t
            g(p) && (setindex!(ks, p.first, C); f(dest, ks, p))
        end
    end
    dest
end

"""
    indexmapfilterat(f::Function, fs::Vector{Function}, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                     t::AbstractNode, N::Int)

Indexing analogy to `mapfilterat` with filtered traversal.

Call signature of `f` is: `f(dest, ks, p::Pair)`.
Call signature of `fs[C]` is: fs[C](p::Pair).

See also: [`mapfilterat`](@ref), [`indexmapat`](@ref), [`indexmapat!`](@ref)
"""
function indexmapfilterat(f::Function, fs::Vector{Function},
                          dims::Tuple{Vararg{NTuple{S, Int}} where S}, t::AbstractNode, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    ks = Vector{Any}(undef, N - 1)
    indexmapfilterat!(f, fs, dest, ks, t, N, 1)
end

################ reduction of vector of into single dest
function indexmapat!(f::Function, dest::Tuple{Vararg{Array{T}} where T}, ks::Vector,
                     ts::Vector{<:AbstractNode}, N::Int, C::Int)
    for t in ts
        indexmapat!(f, dest, ks, t, N, C)
    end
    dest
end

function indexmapat(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                    ts::Vector{<:AbstractNode}, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    ks = Vector{Any}(undef, N - 1)
    indexmapat!(f, dest, ks, ts, N, 1)
end

function tindexmapat(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                     ts::Vector{<:AbstractNode}, L::Int)
    N = length(ts)
    M = Threads.threads()
    ranges = equalranges(N, M)
    A = Vector{Tuple{(Array{Int, length(d)} for d in dims)...}}(undef, M)
    Threads.@threads for m = 1:M
        A[m] = indexmapat(f, dims, ts[ranges[m]], L)
    end
    return A
end

#### filter
function indexmapfilterat!(f::Function, fs::Vector{Function}, dest::Tuple{Vararg{Array{T}} where T},
                           ks::Vector, ts::Vector{<:AbstractNode}, N::Int, C::Int)
    for t in ts
        indexmapfilterat!(f, fs, dest, ks, t, N, C)
    end
    dest
end

function indexmapfilterat(f::Function, fs::Vector{Function},
                          dims::Tuple{Vararg{NTuple{S, Int}} where S},
                          ts::Vector{<:AbstractNode}, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    ks = Vector{Any}(undef, N - 1)
    indexmapfilterat!(f, fs, dest, ks, ts, N, 1)
end

function tindexmapfilterat(f::Function, fs::Vector{Function},
                           dims::Tuple{Vararg{NTuple{S, Int}} where S},
                           ts::Vector{<:AbstractNode}, L::Int)
    N = length(ts)
    M = Threads.threads()
    ranges = equalranges(N, M)
    A = Vector{Tuple{(Array{Int, length(d)} for d in dims)...}}(undef, M)
    Threads.@threads for m = 1:M
        A[m] = indexmapfilterat(f, fs, dims, ts[ranges[m]], L)
    end
    return A
end

################################################################
"""
    indexmapupto!(f::Function, dest::Tuple{Vararg{Array{T}} where T}, ks::Vector
                  t::AbstractNode, N::Int, C::Int)

Analogy to `mapupto!` which tracks Cartesian key-index of iteration, passing
it as an additional argument to `f`.

Call signature of `f` is: `f(dest, ks, t::AbstractNode, N, C)`.

See also: [`mapupto!`](@ref), [`indexmapfilterupto`](@ref), [`indexmapfilterupto!`](@ref)

"""
function indexmapupto!(f::Function, dest::Tuple{Vararg{Array{T}} where T}, ks::Vector,
                       t::AbstractNode, N::Int, C::Int)
    C̃ = C + 1
    C̃ < N || return dest
    f(dest, ks, t, N, C)
    isempty(t) && return dest
    for p in t
        setindex!(ks, p.first, C)
        indexmapupto!(f, dest, ks, p.second, N, C̃)
    end
    dest
end

"""
    indexmapupto(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S}, t::AbstractNode, N::Int)

Analogy to `mapupto!` which tracks Cartesian key-index of iteration, passing
it as an additional argument to `f`.

Call signature of `f` is: `f(dest, ks, t::AbstractNode, N, C)`.

See also: [`mapupto`](@ref), [`indexmapfilterupto`](@ref), [`indexmapfilterupto!`](@ref)

"""
function indexmapupto(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                      t::AbstractNode, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    ks = Vector{Any}(undef, N - 1)
    indexmapupto!(f, dest, ks, t, N, 1)
end

#### filter
"""
    indexmapfilterupto!(f::Function, fs::Vector{Function}, dest::Tuple{Vararg{Array{T}} where T},
                        ks::Vector, t::AbstractNode, N::Int, C::Int)

Filtered traversal analogy to `mapfilterupto!`.

Call signature of `f` is: `f(dest, ks, t::AbstractNode, N, C)`.
Call signature of `fs[C]` is: fs[C](p::Pair).

See also: [`mapfilterupto!`](@ref), [`indexmapupto`](@ref), [`indexmapupto!`](@ref)
"""
function indexmapfilterupto!(f::Function, fs::Vector{Function},
                             dest::Tuple{Vararg{Array{T}} where T},
                             ks::Vector, t::AbstractNode, N::Int, C::Int)
    C̃ = C + 1
    C̃ < N || return dest
    f(dest, ks, t, N, C)
    isempty(t) && return dest
    g = fs[C]
    for p in t
        g(p) && (setindex!(ks, p.first, C); indexmapfilterupto!(f, fs, dest, ks, p.second, N, C̃))
    end
    dest
end

"""
    indexmapfilterupto(f::Function, fs::Vector{Function},
                       dims::Tuple{Vararg{NTuple{S, Int}} where S},
                       t::AbstractNode, N::Int)

Filtered traversal analogy to `mapfilterupto`.

Call signature of `f` is: `f(dest, ks, t::AbstractNode, N, C)`.
Call signature of `fs[C]` is: fs[C](p::Pair).

See also: [`mapfilterupto`](@ref), [`indexmapupto`](@ref), [`indexmapupto!`](@ref)

"""
function indexmapfilterupto(f::Function, fs::Vector{Function},
                            dims::Tuple{Vararg{NTuple{S, Int}} where S},
                            t::AbstractNode, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    ks = Vector{Any}(undef, N - 1)
    indexmapfilterupto!(f, fs, dest, ks, t, N, 1)
end

#### levs_ks

"""
    indexmapupto!(f::Function, dest::Tuple{Vararg{Array{T}} where T}, ks::Vector,
                  t::AbstractNode, N::Int, C::Int, levs_ks::Vector{Vector{Any}})

Level index set-respective analogy to `mapupto!`.

Call signature of `f` is: `f(dest, ks, t::AbstractNode, N, C, levs_ks)`.
"""
function indexmapupto!(f::Function, dest::Tuple{Vararg{Array{T}} where T}, ks::Vector,
                       N::Int, C::Int, levs_ks::Vector{Vector{Any}})
    C̃ = C + 1
    C̃ < N || return dest
    f(dest, ks, t, N, C, levs_ks)
    isempty(t) && return dest
    for p in t
        setindex!(ks, p.first, C)
        indexmapupto!(f, dest, ks, p.second, N, C̃, levs_ks)
    end
    dest
end

"""
    indexmapupto(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                 t::AbstractNode, N::Int, levs_ks::Vector{Vector{Any}})

Level index set-respective analogy to `mapupto`.

Call signature of `f` is: `f(dest, ks, t::AbstractNode, N, C, levs_ks)`.
"""
function indexmapupto(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                      t::AbstractNode, N::Int, levs_ks::Vector{Vector{Any}})
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    ks = Vector{Any}(undef, N - 1)
    indexmapupto!(f, dest, ks, t, N, 1, levs_ks)
end

#### filter and levs_ks
"""
    indexmapfilterupto!(f::Function, fs::Vector{Function}, dest::Tuple{Vararg{Array{T}} where T},
                        ks::Vector, t::AbstractNode, N::Int, C::Int, levs_ks::Vector{Vector{Any}})

Filtered traversal, with level-respective index sets analogy to `mapfilterupto!`.

Call signature of `f` is: `f(dest, ks, t::AbstractNode, N, C, levs_ks)`.
Call signature of `fs[C]` is: fs[C](p::Pair).
"""
function indexmapfilterupto!(f::Function, fs::Vector{Function},
                             dest::Tuple{Vararg{Array{T}} where T}, ks::Vector,
                             t::AbstractNode, N::Int, C::Int, levs_ks::Vector{Vector{Any}})
    C̃ = C + 1
    C̃ < N || return dest
    f(dest, ks, t, N, C, levs_ks)
    isempty(t) && return dest
    g = fs[C]
    for p in t
        g(p) && (setindex!(ks, p.first, C); indexmapfilterupto!(f, fs, dest, ks, p.second, N, C̃, levs_ks))
    end
    dest
end

"""
    indexmapfilterupto(f::Function, fs::Vector{Function}, dest::Tuple{Vararg{Array{T}} where T},
                       t::AbstractNode, N::Int, levs_ks::Vector{Vector{Any}})

Filtered traversal, with level-respective index sets analogy to `mapfilterupto`.

Call signature of `f` is: `f(dest, ks, t::AbstractNode, N, C, levs_ks)`.
Call signature of `fs[C]` is: fs[C](p::Pair).
"""
function indexmapfilterupto(f::Function, fs::Vector{Function},
                            dims::Tuple{Vararg{NTuple{S, Int}} where S},
                            t::AbstractNode, N::Int, levs_ks::Vector{Vector{Any}})
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    ks = Vector{Any}(undef, N - 1)
    indexmapfilterupto!(f, fs, dest, ks, t, N, 1, levs_ks)
end

################ reduction of vector of into single dest
function indexmapupto!(f::Function, dest::Tuple{Vararg{Array{T}} where T}, ks::Vector,
                       ts::Vector{<:AbstractNode}, N::Int, C::Int)
    for t in ts
        indexmapupto!(f, dest, ks, t, N, C)
    end
    dest
end

function indexmapupto(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                      ts::Vector{<:AbstractNode}, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    ks = Vector{Any}(undef, N - 1)
    indexmapupto!(f, dest, ks, ts, N, 1)
end

function tindexmapupto(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                       ts::Vector{<:AbstractNode}, L::Int)
    N = length(ts)
    M = Threads.threads()
    ranges = equalranges(N, M)
    A = Vector{Tuple{(Array{Int, length(d)} for d in dims)...}}(undef, M)
    Threads.@threads for m = 1:M
        A[m] = indexmapupto(f, dims, ts[ranges[m]], L)
    end
    return A
end

#### filter
function indexmapfilterupto!(f::Function, fs::Vector{Function},
                             dest::Tuple{Vararg{Array{T}} where T},
                             ks::Vector, ts::Vector{<:AbstractNode}, N::Int, C::Int)
    for t in ts
        indexmapfilterupto!(f, fs, dest, ks, t, N, C)
    end
    dest
end

function indexmapfilterupto(f::Function, fs::Vector{Function},
                            dims::Tuple{Vararg{NTuple{S, Int}} where S},
                            ts::Vector{<:AbstractNode}, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    ks = Vector{Any}(undef, N - 1)
    indexmapfilterupto!(f, fs, dest, ks, ts, N, 1)
end

function tindexmapfilterupto(f::Function, fs::Vector{Function},
                             dims::Tuple{Vararg{NTuple{S, Int}} where S}, ts::Vector{<:AbstractNode},
                             L::Int)
    N = length(ts)
    M = Threads.threads()
    ranges = equalranges(N, M)
    A = Vector{Tuple{(Array{Int, length(d)} for d in dims)...}}(undef, M)
    Threads.@threads for m = 1:M
        A[m] = indexmapfilterupto(f, fs, dims, ts[ranges[m]], L)
    end
    return A
end

#### levs_ks

function indexmapupto!(f::Function, dest::Tuple{Vararg{Array{T}} where T}, ks::Vector,
                       ts::Vector{<:AbstractNode}, N::Int, C::Int,
                       levs_kss::Vector{Vector{Vector{Any}}})
    for i in eachindex(ts)
        indexmapupto!(f, dest, ks, ts[i], N, C, levs_kss[i])
    end
    dest
end

function indexmapupto(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                      ts::Vector{<:AbstractNode}, N::Int, levs_kss::Vector{Vector{Vector{Any}}})
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    ks = Vector{Any}(undef, N - 1)
    indexmapupto!(f, dest, ks, ts, N, 1, levs_kss)
end

function tindexmapupto(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                       ts::Vector{<:AbstractNode}, L::Int, levs_kss::Vector{Vector{Vector{Any}}})
    N = length(ts)
    M = Threads.threads()
    ranges = equalranges(N, M)
    A = Vector{Tuple{(Array{Int, length(d)} for d in dims)...}}(undef, M)
    Threads.@threads for m = 1:M
        A[m] = indexmapupto(f, dims, ts[ranges[m]], L, levs_kss[ranges[m]])
    end
    return A
end

#### filter and levs_ks
function indexmapfilterupto!(f::Function, fs::Vector{Function},
                             dest::Tuple{Vararg{Array{T}} where T},
                             ks::Vector, ts::Vector{<:AbstractNode}, N::Int, C::Int,
                             levs_kss::Vector{Vector{Vector{Any}}})
    for i in eachindex(ts)
        indexmapfilterupto!(f, fs, dest, ks, ts[i], N, C, levs_kss[i])
    end
    dest
end

function indexmapfilterupto(f::Function, fs::Vector{Function},
                            dims::Tuple{Vararg{NTuple{S, Int}} where S},
                            ts::Vector{<:AbstractNode}, N::Int,
                            levs_kss::Vector{Vector{Vector{Any}}})
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    ks = Vector{Any}(undef, N - 1)
    indexmapfilterupto!(f, fs, dest, ks, ts, N, 1, levs_kss)
end

function tindexmapfilterupto(f::Function, fs::Vector{Function},
                             dims::Tuple{Vararg{NTuple{S, Int}} where S},
                             ts::Vector{<:AbstractNode}, L::Int,
                             levs_kss::Vector{Vector{Vector{Any}}})
    N = length(ts)
    M = Threads.threads()
    ranges = equalranges(N, M)
    A = Vector{Tuple{(Array{Int, length(d)} for d in dims)...}}(undef, M)
    Threads.@threads for m = 1:M
        A[m] = indexmapfilterupto(f, fs, dims, ts[ranges[m]], L, levs_kss[ranges[m]])
    end
    return A
end

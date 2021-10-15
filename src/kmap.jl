#
# Date created: 2021-10-13
# Author: aradclif
#
#
############################################################################################
"""
    kmapat!(f::Function, dest::Tuple{Vararg{Array{T}} where T}, ks::Vector,
            t::AbstractNode, N::Int, C::Int)

Analogy to `mapat!` which tracks Cartesian key-index of iteration, passing
it as an additional argument to `f`.

Call signature of `f` is: `f(dest, ks, p::Pair)`.

See also: [`mapat!`](@ref), [`kmapat`](@ref), [`kmapfilterat`](@ref), [`kmapfilterat!`](@ref)

"""
function kmapat!(f::Function, dest::Tuple{Vararg{Array{T}} where T}, ks::Vector,
                 t::AbstractNode, N::Int, C::Int)
    isempty(t) && return dest
    C̃ = C + 1
    if C̃ < N
        for p in t
            setindex!(ks, p.first, C)
            kmapat!(f, dest, ks, p.second, N, C̃)
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
    kmapat(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S}, t::AbstractNode, N::Int)

Analogy to `mapat` which tracks Cartesian key-index of iteration, passing
it as an additional argument to `f`.

Call signature of `f` is: `f(dest, ks, p::Pair)`.

See also: [`mapat`](@ref), [`kmapat`](@ref), [`kmapfilterat`](@ref), [`kmapfilterat!`](@ref)

"""
function kmapat(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                t::AbstractNode, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    ks = Vector{Any}(undef, N - 1)
    kmapat!(f, dest, ks, t, N, 1)
end

#### filter
"""
    kmapfilterat!(f::Function, fs::Vector{Function}, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                  ks::Vector, t::AbstractNode, N::Int, C::Int)

Indexing analogy to `mapfilterat!` with filtered traversal.

Call signature of `f` is: `f(dest, ks, p::Pair)`.
Call signature of `fs[C]` is: `fs[C](p::Pair)`.

See also: [`mapfilterat!`](@ref), [`kmapat`](@ref), [`kmapat!`](@ref)
"""
function kmapfilterat!(f::Function, fs::Vector{Function}, dest::Tuple{Vararg{Array{T}} where T},
                       ks::Vector, t::AbstractNode, N::Int, C::Int)
    isempty(t) && return dest
    C̃ = C + 1
    g = fs[C]
    if C̃ < N
        for p in t
            g(p) && (setindex!(ks, p.first, C); kmapfilterat!(f, fs, dest, p.second, N, C̃))
        end
    elseif C̃ == N
        for p in t
            g(p) && (setindex!(ks, p.first, C); f(dest, ks, p))
        end
    end
    dest
end

"""
    kmapfilterat(f::Function, fs::Vector{Function}, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                 t::AbstractNode, N::Int)

Indexing analogy to `mapfilterat` with filtered traversal.

Call signature of `f` is: `f(dest, ks, p::Pair)`.
Call signature of `fs[C]` is: `fs[C](p::Pair)`.

See also: [`mapfilterat`](@ref), [`kmapat`](@ref), [`kmapat!`](@ref)
"""
function kmapfilterat(f::Function, fs::Vector{Function},
                      dims::Tuple{Vararg{NTuple{S, Int}} where S}, t::AbstractNode, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    ks = Vector{Any}(undef, N - 1)
    kmapfilterat!(f, fs, dest, ks, t, N, 1)
end

################ reduction of vector of into single dest
function kmapat!(f::Function, dest::Tuple{Vararg{Array{T}} where T}, ks::Vector,
                 ts::Vector{<:AbstractNode}, N::Int, C::Int)
    for t in ts
        kmapat!(f, dest, ks, t, N, C)
    end
    dest
end

function kmapat(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                ts::Vector{<:AbstractNode}, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    ks = Vector{Any}(undef, N - 1)
    kmapat!(f, dest, ks, ts, N, 1)
end

function tkmapat(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                 ts::Vector{<:AbstractNode}, L::Int)
    N = length(ts)
    M = Threads.threads()
    ranges = equalranges(N, M)
    A = Vector{Tuple{(Array{Int, length(d)} for d in dims)...}}(undef, M)
    Threads.@threads for m = 1:M
        A[m] = kmapat(f, dims, ts[ranges[m]], L)
    end
    return A
end

#### filter
function kmapfilterat!(f::Function, fs::Vector{Function}, dest::Tuple{Vararg{Array{T}} where T},
                       ks::Vector, ts::Vector{<:AbstractNode}, N::Int, C::Int)
    for t in ts
        kmapfilterat!(f, fs, dest, ks, t, N, C)
    end
    dest
end

function kmapfilterat(f::Function, fs::Vector{Function}, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                      ts::Vector{<:AbstractNode}, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    ks = Vector{Any}(undef, N - 1)
    kmapfilterat!(f, fs, dest, ks, ts, N, 1)
end

function tkmapfilterat(f::Function, fs::Vector{Function},
                       dims::Tuple{Vararg{NTuple{S, Int}} where S},
                       ts::Vector{<:AbstractNode}, L::Int)
    N = length(ts)
    M = Threads.threads()
    ranges = equalranges(N, M)
    A = Vector{Tuple{(Array{Int, length(d)} for d in dims)...}}(undef, M)
    Threads.@threads for m = 1:M
        A[m] = kmapfilterat(f, fs, dims, ts[ranges[m]], L)
    end
    return A
end

################################################################
"""
    kmapupto!(f::Function, dest::Tuple{Vararg{Array{T}} where T}, ks::Vector
              t::AbstractNode, N::Int, C::Int)

Analogy to `mapupto!` which tracks Cartesian key-index of iteration, passing
it as an additional argument to `f`.

Call signature of `f` is: `f(dest, ks, t::AbstractNode, N, C)`.

See also: [`mapupto!`](@ref), [`kmapfilterupto`](@ref), [`kmapfilterupto!`](@ref)

"""
function kmapupto!(f::Function, dest::Tuple{Vararg{Array{T}} where T}, ks::Vector,
                   t::AbstractNode, N::Int, C::Int)
    C̃ = C + 1
    C̃ < N || return dest
    f(dest, ks, t, N, C)
    isempty(t) && return dest
    for p in t
        setindex!(ks, p.first, C)
        kmapupto!(f, dest, ks, p.second, N, C̃)
    end
    dest
end

"""
    kmapupto(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S}, t::AbstractNode, N::Int)

Analogy to `mapupto!` which tracks Cartesian key-index of iteration, passing
it as an additional argument to `f`.

Call signature of `f` is: `f(dest, ks, t::AbstractNode, N, C)`.

See also: [`mapupto`](@ref), [`kmapfilterupto`](@ref), [`kmapfilterupto!`](@ref)

"""
function kmapupto(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S}, t::AbstractNode, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    ks = Vector{Any}(undef, N - 1)
    kmapupto!(f, dest, ks, t, N, 1)
end

#### filter
"""
    kmapfilterupto!(f::Function, fs::Vector{Function}, dest::Tuple{Vararg{Array{T}} where T},
                    ks::Vector, t::AbstractNode, N::Int, C::Int)

Filtered traversal analogy to `mapfilterupto!`.

Call signature of `f` is: `f(dest, ks, t::AbstractNode, N, C)`.
Call signature of `fs[C]` is: `fs[C](p::Pair)`.

See also: [`mapfilterupto!`](@ref), [`kmapupto`](@ref), [`kmapupto!`](@ref)
"""
function kmapfilterupto!(f::Function, fs::Vector{Function}, dest::Tuple{Vararg{Array{T}} where T},
                         ks::Vector, t::AbstractNode, N::Int, C::Int)
    C̃ = C + 1
    C̃ < N || return dest
    f(dest, ks, t, N, C)
    isempty(t) && return dest
    g = fs[C]
    for p in t
        g(p) && (setindex!(ks, p.first, C); kmapfilterupto!(f, fs, dest, ks, p.second, N, C̃))
    end
    dest
end

"""
    kmapfilterupto(f::Function, fs::Vector{Function}, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                   t::AbstractNode, N::Int)

Filtered traversal analogy to `mapfilterupto`.

Call signature of `f` is: `f(dest, ks, t::AbstractNode, N, C)`.
Call signature of `fs[C]` is: `fs[C](p::Pair)`.

See also: [`mapfilterupto`](@ref), [`kmapupto`](@ref), [`kmapupto!`](@ref)

"""
function kmapfilterupto(f::Function, fs::Vector{Function},
                        dims::Tuple{Vararg{NTuple{S, Int}} where S}, t::AbstractNode, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    ks = Vector{Any}(undef, N - 1)
    kmapfilterupto!(f, fs, dest, ks, t, N, 1)
end

#### levs_ks

"""
    kmapupto!(f::Function, dest::Tuple{Vararg{Array{T}} where T}, ks::Vector,
              t::AbstractNode, N::Int, C::Int, levs_ks::Vector{Vector{Any}})

Level index set-respective analogy to `mapupto!`.

Call signature of `f` is: `f(dest, ks, t::AbstractNode, N, C, levs_ks)`.
"""
function kmapupto!(f::Function, dest::Tuple{Vararg{Array{T}} where T}, ks::Vector,
                   N::Int, C::Int, levs_ks::Vector{Vector{Any}})
    C̃ = C + 1
    C̃ < N || return dest
    f(dest, ks, t, N, C, levs_ks)
    isempty(t) && return dest
    for p in t
        setindex!(ks, p.first, C)
        kmapupto!(f, dest, ks, p.second, N, C̃, levs_ks)
    end
    dest
end

"""
    kmapupto(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S},
             t::AbstractNode, N::Int, levs_ks::Vector{Vector{Any}})

Level index set-respective analogy to `mapupto`.

Call signature of `f` is: `f(dest, ks, t::AbstractNode, N, C, levs_ks)`.
"""
function kmapupto(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                  t::AbstractNode, N::Int, levs_ks::Vector{Vector{Any}})
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    ks = Vector{Any}(undef, N - 1)
    kmapupto!(f, dest, ks, t, N, 1, levs_ks)
end

#### filter and levs_ks
"""
    kmapfilterupto!(f::Function, fs::Vector{Function}, dest::Tuple{Vararg{Array{T}} where T},
                    ks::Vector, t::AbstractNode, N::Int, C::Int, levs_ks::Vector{Vector{Any}})

Filtered traversal, with level-respective index sets analogy to `mapfilterupto!`.

Call signature of `f` is: `f(dest, ks, t::AbstractNode, N, C, levs_ks)`.
Call signature of `fs[C]` is: `fs[C](p::Pair)`.
"""
function kmapfilterupto!(f::Function, fs::Vector{Function}, dest::Tuple{Vararg{Array{T}} where T},
                         ks::Vector, t::AbstractNode, N::Int, C::Int, levs_ks::Vector{Vector{Any}})
    C̃ = C + 1
    C̃ < N || return dest
    f(dest, ks, t, N, C, levs_ks)
    isempty(t) && return dest
    g = fs[C]
    for p in t
        g(p) && (setindex!(ks, p.first, C); kmapfilterupto!(f, fs, dest, ks, p.second, N, C̃, levs_ks))
    end
    dest
end

"""
    kmapfilterupto(f::Function, fs::Vector{Function}, dest::Tuple{Vararg{Array{T}} where T},
                   t::AbstractNode, N::Int, levs_ks::Vector{Vector{Any}})

Filtered traversal, with level-respective index sets analogy to `mapfilterupto`.

Call signature of `f` is: `f(dest, ks, t::AbstractNode, N, C, levs_ks)`.
Call signature of `fs[C]` is: `fs[C](p::Pair)`.
"""
function kmapfilterupto(f::Function, fs::Vector{Function},
                        dims::Tuple{Vararg{NTuple{S, Int}} where S},
                        t::AbstractNode, N::Int, levs_ks::Vector{Vector{Any}})
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    ks = Vector{Any}(undef, N - 1)
    kmapfilterupto!(f, fs, dest, ks, t, N, 1, levs_ks)
end

################ reduction of vector of into single dest
function kmapupto!(f::Function, dest::Tuple{Vararg{Array{T}} where T}, ks::Vector,
                   ts::Vector{<:AbstractNode}, N::Int, C::Int)
    for t in ts
        kmapupto!(f, dest, ks, t, N, C)
    end
    dest
end

function kmapupto(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                  ts::Vector{<:AbstractNode}, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    ks = Vector{Any}(undef, N - 1)
    kmapupto!(f, dest, ks, ts, N, 1)
end

function tkmapupto(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                   ts::Vector{<:AbstractNode}, L::Int)
    N = length(ts)
    M = Threads.threads()
    ranges = equalranges(N, M)
    A = Vector{Tuple{(Array{Int, length(d)} for d in dims)...}}(undef, M)
    Threads.@threads for m = 1:M
        A[m] = kmapupto(f, dims, ts[ranges[m]], L)
    end
    return A
end

#### filter
function kmapfilterupto!(f::Function, fs::Vector{Function}, dest::Tuple{Vararg{Array{T}} where T},
                         ks::Vector, ts::Vector{<:AbstractNode}, N::Int, C::Int)
    for t in ts
        kmapfilterupto!(f, fs, dest, ks, t, N, C)
    end
    dest
end

function kmapfilterupto(f::Function, fs::Vector{Function},
                        dims::Tuple{Vararg{NTuple{S, Int}} where S},
                        ts::Vector{<:AbstractNode}, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    ks = Vector{Any}(undef, N - 1)
    kmapfilterupto!(f, fs, dest, ks, ts, N, 1)
end

function tkmapfilterupto(f::Function, fs::Vector{Function},
                         dims::Tuple{Vararg{NTuple{S, Int}} where S}, ts::Vector{<:AbstractNode},
                         L::Int)
    N = length(ts)
    M = Threads.threads()
    ranges = equalranges(N, M)
    A = Vector{Tuple{(Array{Int, length(d)} for d in dims)...}}(undef, M)
    Threads.@threads for m = 1:M
        A[m] = kmapfilterupto(f, fs, dims, ts[ranges[m]], L)
    end
    return A
end

#### levs_ks

function kmapupto!(f::Function, dest::Tuple{Vararg{Array{T}} where T}, ks::Vector,
                   ts::Vector{<:AbstractNode}, N::Int, C::Int, levs_kss::Vector{Vector{Vector{Any}}})
    for i in eachindex(ts)
        kmapupto!(f, dest, ks, ts[i], N, C, levs_kss[i])
    end
    dest
end

function kmapupto(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                  ts::Vector{<:AbstractNode}, N::Int, levs_kss::Vector{Vector{Vector{Any}}})
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    ks = Vector{Any}(undef, N - 1)
    kmapupto!(f, dest, ks, ts, N, 1, levs_kss)
end

function tkmapupto(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                   ts::Vector{<:AbstractNode}, L::Int, levs_kss::Vector{Vector{Vector{Any}}})
    N = length(ts)
    M = Threads.threads()
    ranges = equalranges(N, M)
    A = Vector{Tuple{(Array{Int, length(d)} for d in dims)...}}(undef, M)
    Threads.@threads for m = 1:M
        A[m] = kmapupto(f, dims, ts[ranges[m]], L, levs_kss[ranges[m]])
    end
    return A
end

#### filter and levs_ks
function kmapfilterupto!(f::Function, fs::Vector{Function}, dest::Tuple{Vararg{Array{T}} where T},
                         ks::Vector, ts::Vector{<:AbstractNode}, N::Int, C::Int,
                         levs_kss::Vector{Vector{Vector{Any}}})
    for i in eachindex(ts)
        kmapfilterupto!(f, fs, dest, ks, ts[i], N, C, levs_kss[i])
    end
    dest
end

function kmapfilterupto(f::Function, fs::Vector{Function},
                        dims::Tuple{Vararg{NTuple{S, Int}} where S},
                        ts::Vector{<:AbstractNode}, N::Int, levs_kss::Vector{Vector{Vector{Any}}})
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    ks = Vector{Any}(undef, N - 1)
    kmapfilterupto!(f, fs, dest, ks, ts, N, 1, levs_kss)
end

function tkmapfilterupto(f::Function, fs::Vector{Function},
                         dims::Tuple{Vararg{NTuple{S, Int}} where S},
                         ts::Vector{<:AbstractNode}, L::Int, levs_kss::Vector{Vector{Vector{Any}}})
    N = length(ts)
    M = Threads.threads()
    ranges = equalranges(N, M)
    A = Vector{Tuple{(Array{Int, length(d)} for d in dims)...}}(undef, M)
    Threads.@threads for m = 1:M
        A[m] = kmapfilterupto(f, fs, dims, ts[ranges[m]], L, levs_kss[ranges[m]])
    end
    return A
end

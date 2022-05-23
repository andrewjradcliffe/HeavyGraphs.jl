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
# Note: Tuple{Vararg{Array{T} where T}} ≡ Tuple{Vararg{AbstractArray}}

#### A draft of analyzeat!, p. 483, 2021-10-03; p. 500-507, 2021-10-12
"""
    mapat!(f::Function, dest::Tuple{Vararg{AbstractArray}}, t::AbstractGraph, N::Int, C::Int)

Transform the graph/tree `t` by applying `f` to each element at a given level, `N`.
Results are stored in `dest`; it is assumed that the number and dimensions of arrays
are sufficient to store the results of calling `f`.

Call signature of `f` is: `f(dest, t::AbstractGraph)`.

See also: [`mapat`](@ref), [`mapfilterat`](@ref), [`mapfilterat!`](@ref)

"""
function mapat!(f::Function, dest::Tuple{Vararg{AbstractArray}}, t::AbstractGraph, N::Int, C::Int)
    isempty(t) && return dest
    C̃ = C + 1
    if C̃ < N
        for p ∈ t
            mapat!(f, dest, p.second, N, C̃)
        end
    elseif C̃ == N
        for p ∈ t
            f(dest, p.second)
        end
    end
    dest
end

"""
    mapat(f::Function, dims::Tuple{Vararg{Tuple{Vararg{Int}}}}, t::AbstractGraph, N::Int)

Transform the graph/tree `t` by applying `f` to each element at a given level, `N`.
Results are stored in a tuple of zero-initialized arrays, the number and dimension of which
are specified by `dims`.

Call signature of `f` is: `f(dest, t::AbstractGraph)`.

See also: [`mapat!`](@ref), [`mapfilterat`](@ref), [`mapfilterat!`](@ref)

"""
function mapat(f::Function, dims::Tuple{Vararg{Tuple{Vararg{Int}}}}, t::AbstractGraph, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    mapat!(f, dest, t, N, 1)
end

#### filter
"""
    mapfilterat!(f::Function, fs::Vector{Function}, dest::Tuple{Vararg{AbstractArray}},
                 t::AbstractGraph, N::Int, C::Int)

Performs a filtered traversal in which a subset is formed at each level from
the corresponding element (note: linear indexed) of the filter functions `fs`.

Call signature of `f` is: `f(dest, t::AbstractGraph)`.
Call signature of `fs[C]` is: `fs[C](p::Pair)`.

See also: [`mapfilterat`](@ref), [`mapat`](@ref), [`mapat!`](@ref)
"""
function mapfilterat!(f::Function, fs::Vector{Function}, dest::Tuple{Vararg{AbstractArray}},
                      t::AbstractGraph, N::Int, C::Int)
    isempty(t) && return dest
    C̃ = C + 1
    g = fs[C]
    if C̃ < N
        for p ∈ t
            g(p) && mapfilterat!(f, fs, dest, p.second, N, C̃)
        end
    elseif C̃ == N
        for p ∈ t
            g(p) && f(dest, p.second)
        end
    end
    dest
end

"""
    mapfilterat(f::Function, fs::Vector{Function}, dims::Tuple{Vararg{Tuple{Vararg{Int}}}},
                t::AbstractGraph, N::Int)

Performs a filtered traversal in which a subset is formed at each level from
the corresponding element (note: linear indexed) of the filter functions `fs`.

Call signature of `f` is: `f(dest, t::AbstractGraph)`.
Call signature of `fs[C]` is: `fs[C](p::Pair)`.

See also: [`mapfilterat!`](@ref), [`mapat`](@ref), [`mapat!`](@ref)
"""
function mapfilterat(f::Function, fs::Vector{Function}, dims::Tuple{Vararg{Tuple{Vararg{Int}}}},
                     t::AbstractGraph, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    mapfilterat!(f, fs, dest, t, N, 1)
end

################ reduction of vector of into single dest
function mapat!(f::Function, dest::Tuple{Vararg{AbstractArray}},
                ts::Vector{<:AbstractGraph}, N::Int, C::Int)
    for t ∈ ts
        mapat!(f, dest, t, N, C)
    end
    dest
end

function mapat(f::Function, dims::Tuple{Vararg{Tuple{Vararg{Int}}}},
               ts::Vector{<:AbstractGraph}, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    mapat!(f, dest, ts, N, 1)
end

function tmapat(f::Function, dims::Tuple{Vararg{Tuple{Vararg{Int}}}},
                ts::Vector{<:AbstractGraph}, L::Int, M::Int=Threads.nthreads())
    N = length(ts)
    # M = Threads.threads()
    ranges = equalranges(N, M)
    A = Vector{Tuple{(Array{Int, length(d)} for d ∈ dims)...}}(undef, M)
    Threads.@threads for m = 1:M
        A[m] = mapat(f, dims, ts[ranges[m]], L)
    end
    # B = ntuple(i -> zeros(Int, dims[i]), length(dims))
    # for m ∈ eachindex(A)
    #     for i ∈ eachindex(B)
    #         sum!(B[i], A[m][i])
    #     end
    # end
    return reduce(.+, A)
end

#### filter
function mapfilterat!(f::Function, fs::Vector{Function}, dest::Tuple{Vararg{AbstractArray}},
                      ts::Vector{<:AbstractGraph}, N::Int, C::Int)
    for t ∈ ts
        mapfilterat!(f, fs, dest, t, N, C)
    end
    dest
end

function mapfilterat(f::Function, fs::Vector{Function}, dims::Tuple{Vararg{Tuple{Vararg{Int}}}},
                     ts::Vector{<:AbstractGraph}, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    mapfilterat!(f, fs, dest, ts, N, 1)
end

function tmapfilterat(f::Function, fs::Vector{Function}, dims::Tuple{Vararg{Tuple{Vararg{Int}}}},
                      ts::Vector{<:AbstractGraph}, L::Int, M::Int=Threads.nthreads())
    N = length(ts)
    # M = Threads.threads()
    ranges = equalranges(N, M)
    A = Vector{Tuple{(Array{Int, length(d)} for d ∈ dims)...}}(undef, M)
    Threads.@threads for m = 1:M
        A[m] = mapfilterat(f, fs, dims, ts[ranges[m]], L)
    end
    return reduce(.+, A)
end

################################################################
"""
    mapat_pairs!(f::Function, dest::Tuple{Vararg{AbstractArray}}, t::AbstractGraph,
                 N::Int, C::Int)

Analogous to `mapat!`, but call signature of `f` is: `f(dest, p::Pair)`.

See also: [`mapat!`](@ref), [`mapfilterat_pairs!`]

"""
function mapat_pairs!(f::Function, dest::Tuple{Vararg{AbstractArray}}, t::AbstractGraph,
                      N::Int, C::Int)
    isempty(t) && return dest
    C̃ = C + 1
    if C̃ < N
        for p ∈ t
            mapat_pairs!(f, dest, p.second, N, C̃)
        end
    elseif C̃ == N
        for p ∈ t
            f(dest, p)
        end
    end
    dest
end

"""
    mapat_pairs(f::Function, dims::Tuple{Vararg{Tuple{Vararg{Int}}}}, t::AbstractGraph, N::Int)

Analogous to `mapat`, but call signature of `f` is: `f(dest, p::Pair)`.

See also: [`mapat`](@ref), [`mapfilterat_pairs`](@ref)

"""
function mapat_pairs(f::Function, dims::Tuple{Vararg{Tuple{Vararg{Int}}}}, t::AbstractGraph,
                     N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    mapat_pairs!(f, dest, t, N, 1)
end

#### filter
"""
    mapfilterat_pairs!(f::Function, fs::Vector{Function}, dest::Tuple{Vararg{AbstractArray}},
                       t::AbstractGraph, N::Int, C::Int)

Analogous to `mapfilterat!`, but call signature of `f` is: `f(dest, p::Pair)`.
Call signature of `fs[C]` remains unchanged: `fs[C](p::Pair)`.

See also: [`mapfilterat!`](@ref), [`mapfilterat_pairs`](@ref)
"""
function mapfilterat_pairs!(f::Function, fs::Vector{Function}, dest::Tuple{Vararg{AbstractArray}},
                            t::AbstractGraph, N::Int, C::Int)
    isempty(t) && return dest
    C̃ = C + 1
    g = fs[C]
    if C̃ < N
        for p ∈ t
            g(p) && mapfilterat_pairs!(f, fs, dest, p.second, N, C̃)
        end
    elseif C̃ == N
        for p ∈ t
            g(p) && f(dest, p)
        end
    end
    dest
end

"""
    mapfilterat_pairs(f::Function, fs::Vector{Function},
                      dims::Tuple{Vararg{Tuple{Vararg{Int}}}}, t::AbstractGraph, N::Int)

Analogous to `mapfilterat`, but call signature of `f` is: `f(dest, p::Pair)`.
Call signature of `fs[C]` remains unchanged: `fs[C](p::Pair)`.

See also: [`mapfilterat`](@ref), [`mapfilterat_pairs!`](@ref)
"""
function mapfilterat_pairs(f::Function, fs::Vector{Function},
                           dims::Tuple{Vararg{Tuple{Vararg{Int}}}}, t::AbstractGraph, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    mapfilterat_pairs!(f, fs, dest, t, N, 1)
end

################ reduction of vector of into single dest
function mapat_pairs!(f::Function, dest::Tuple{Vararg{AbstractArray}},
                      ts::Vector{<:AbstractGraph}, N::Int, C::Int)
    for t ∈ ts
        mapat_pairs!(f, dest, t, N, C)
    end
    dest
end

function mapat_pairs(f::Function, dims::Tuple{Vararg{Tuple{Vararg{Int}}}},
                     ts::Vector{<:AbstractGraph}, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    mapat_pairs!(f, dest, ts, N, 1)
end

function tmapat_pairs(f::Function, dims::Tuple{Vararg{Tuple{Vararg{Int}}}},
                      ts::Vector{<:AbstractGraph}, L::Int, M::Int=Threads.nthreads())
    N = length(ts)
    # M = Threads.threads()
    ranges = equalranges(N, M)
    A = Vector{Tuple{(Array{Int, length(d)} for d ∈ dims)...}}(undef, M)
    Threads.@threads for m = 1:M
        A[m] = mapat_pairs(f, dims, ts[ranges[m]], L)
    end
    # B = ntuple(i -> zeros(Int, dims[i]), length(dims))
    # for m ∈ eachindex(A)
    #     for i ∈ eachindex(B)
    #         sum!(B[i], A[m][i])
    #     end
    # end
    return reduce(.+, A)
end

#### filter
function mapfilterat_pairs!(f::Function, fs::Vector{Function}, dest::Tuple{Vararg{AbstractArray}},
                            ts::Vector{<:AbstractGraph}, N::Int, C::Int)
    for t ∈ ts
        mapfilterat_pairs!(f, fs, dest, t, N, C)
    end
    dest
end

function mapfilterat_pairs(f::Function, fs::Vector{Function},
                           dims::Tuple{Vararg{Tuple{Vararg{Int}}}},
                           ts::Vector{<:AbstractGraph}, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    mapfilterat_pairs!(f, fs, dest, ts, N, 1)
end

function tmapfilterat_pairs(f::Function, fs::Vector{Function},
                            dims::Tuple{Vararg{Tuple{Vararg{Int}}}},
                            ts::Vector{<:AbstractGraph}, L::Int, M::Int=Threads.nthreads())
    N = length(ts)
    # M = Threads.threads()
    ranges = equalranges(N, M)
    A = Vector{Tuple{(Array{Int, length(d)} for d ∈ dims)...}}(undef, M)
    Threads.@threads for m = 1:M
        A[m] = mapfilterat_pairs(f, fs, dims, ts[ranges[m]], L)
    end
    return reduce(.+, A)
end

################################################################

"""
    mapupto!(f::Function, dest::Tuple{Vararg{AbstractArray}}, t::AbstractGraph, N::Int)

Transform the graph/tree `t` by applying `f` to each element on each level up to a
given level `N`. Results are stored in `dest`; it is assumed that the number and
dimensions of arrays are sufficient to store the results of calling `f`.

Call signature of `f` is: `f(dest, t::AbstractGraph)`.

See also: [`mapupto`](@ref), [`mapfilterupto`](@ref), [`mapfilterupto!`](@ref)
"""
function mapupto!(f::Function, dest::Tuple{Vararg{AbstractArray}}, t::AbstractGraph, N::Int, C::Int)
    # C̃ = C + 1
    # C̃ < N || return dest
    # f(dest, t)
    # isempty(t) && return dest
    # C < N || return dest # necessary to ensure safety, but technically optional
    f(dest, t)
    isempty(t) && return dest
    C̃ = C + 1
    if C̃ < N #|| return dest
        for p ∈ t
            mapupto!(f, dest, p.second, N, C̃)
        end
    end
    return dest
end

"""
    mapupto(f::Function, dims::Tuple{Vararg{Tuple{Vararg{Int}}}}, t::AbstractGraph, N::Int)

Transform the graph/tree `t` by applying `f` to each element on each level up to a
given level `N`. Results are stored in a tuple of zero-initialized arrays, the number
and dimension of which are specified by `dims`.

Call signature of `f` is: `f(dest, t::AbstractGraph)`.

See also: [`mapupto!`](@ref), [`mapfilterupto`](@ref), [`mapfilterupto!`](@ref)
"""
function mapupto(f::Function, dims::Tuple{Vararg{Tuple{Vararg{Int}}}}, t::AbstractGraph, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    mapupto!(f, dest, t, N, 1)
end

#### filter
"""
    mapfilterupto!(f::Function, fs::Vector{Function}, dest::Tuple{Variant{Array{T}} where T},
                   t::AbstractGraph, N::Int, C::Int)

Performs a filtered traversal in which a subset is formed at each level from
the corresponding element (note: linear indexed) of the filter functions `fs`.

Call signature of `f` is: `f(dest, t::AbstractGraph)`.
Call signature of `fs[C]` is: `fs[C](p::Pair)`.

See also: [`mapfilterupto`](@ref), [`mapupto`](@ref), [`mapupto!`](@ref)
"""
function mapfilterupto!(f::Function, fs::Vector{Function}, dest::Tuple{Vararg{AbstractArray}},
                        t::AbstractGraph, N::Int, C::Int)
    # C̃ = C + 1
    # f(dest, t)
    # C̃ < N || return dest
    # isempty(t) && return dest
    # g = fs[C]
    # for p ∈ t
    #     g(p) && mapfilterupto!(f, fs, dest, p.second, N, C̃)
    # end
    # dest
    # C < N || return dest # necessary to ensure safety, but technically optional
    f(dest, t)
    isempty(t) && return dest
    C̃ = C + 1
    g = fs[C]
    if C̃ < N #|| return dest
        for p ∈ t
            g(p) && mapfilterupto!(f, fs, dest, p.second, N, C̃)
        end
    end
    return dest
end

"""
    mapfilterupto(f::Function, fs::Vector{Function}, dims::Tuple{Vararg{Tuple{Vararg{Int}}}},
                  t::AbstractGraph, N::Int)

Performs a filtered traversal in which a subset is formed at each level from
the corresponding element (note: linear indexed) of the filter functions `fs`.

Call signature of `f` is: `f(dest, t::AbstractGraph)`.
Call signature of `fs[C]` is: `fs[C](p::Pair)`.

See also: [`mapfilterupto!`](@ref), [`mapupto`](@ref), [`mapupto!`](@ref)
"""
function mapfilterupto(f::Function, fs::Vector{Function},
                       dims::Tuple{Vararg{Tuple{Vararg{Int}}}}, t::AbstractGraph, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    mapfilterupto!(f, fs, dest, t, N, 1)
end

#### levs_ks
"""
    mapupto!(f::Function, dest::Tuple{Vararg{AbstractArray}},
             t::AbstractGraph, N::Int, C::Int, levs_ks::Vector{Vector{Any}})

Provides the given level `N`, the current level `C`, and the vector of
level-respective index sets, `levs_ks` as arguments to `f`.

Call signature of `f` is: `f(dest, t::AbstractGraph, N, C, levs_ks)`.
"""
function mapupto!(f::Function, dest::Tuple{Vararg{AbstractArray}},
                  t::AbstractGraph{T}, N::Int, C::Int, levs_ks::Vector{Vector{T}}) where {T}
    # C̃ = C + 1
    # f(dest, t, N, C, levs_ks)
    # C̃ < N || return dest
    # isempty(t) && return dest
    # for p ∈ t
    #     mapupto!(f, dest, p.second, N, C̃, levs_ks)
    # end
    # dest
    # C < N || return dest # necessary to ensure safety, but technically optional
    f(dest, t, N, C, levs_ks)
    isempty(t) && return dest
    C̃ = C + 1
    if C̃ < N #|| return dest
        for p ∈ t
            mapupto!(f, dest, p.second, N, C̃, levs_ks)
        end
    end
    return dest
end

"""
    mapupto(f::Function, dims::Tuple{Vararg{Tuple{Vararg{Int}}}},
            t::AbstractGraph, N::Int, levs_ks::Vector{Vector{Any}})

Provides the given level `N`, the current level `C`, and the vector of
level-respective index sets, `levs_ks` as arguments to `f`.

Call signature of `f` is: `f(dest, t::AbstractGraph, N, C, levs_ks)`.
"""
function mapupto(f::Function, dims::Tuple{Vararg{Tuple{Vararg{Int}}}},
                 t::AbstractGraph{T}, N::Int, levs_ks::Vector{Vector{T}}) where {T}
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    mapupto!(f, dest, t, N, 1, levs_ks)
end

#### filter and levs_ks
"""
    mapfilterupto!(f::Function, fs::Vector{Function}, dest::Tuple{Vararg{AbstractArray}},
                   t::AbstractGraph, N::Int, C::Int, levs_ks::Vector{Vector{Any}})

Performs filtered traversal; provides the given level `N`,
the current level `C`, and the vector of level-respective index sets,
`levs_ks` as arguments to `f`.

Call signature of `f` is: `f(dest, t::AbstractGraph, N, C, levs_ks)`.
Call signature of `fs[C]` is: `fs[C](p::Pair)`.
"""
function mapfilterupto!(f::Function, fs::Vector{Function}, dest::Tuple{Vararg{AbstractArray}},
                        t::AbstractGraph{T}, N::Int, C::Int, levs_ks::Vector{Vector{T}}) where {T}
    # C̃ = C + 1
    # f(dest, t, N, C, levs_ks)
    # C̃ < N || return dest
    # isempty(t) && return dest
    # g = fs[C]
    # for p ∈ t
    #     g(p) && mapfilterupto!(f, fs, dest, p.second, N, C̃, levs_ks)
    # end
    # dest
    # C < N || return dest # necessary to ensure safety, but technically optional
    f(dest, t, N, C, levs_ks)
    isempty(t) && return dest
    C̃ = C + 1
    g = fs[C]
    if C̃ < N #|| return dest
        for p ∈ t
            g(p) && mapfilterupto!(f, fs, dest, p.second, N, C̃, levs_ks)
        end
    end
    return dest
end

"""
    mapfilterupto(f::Function, fs::Vector{Function}, dims::Tuple{Vararg{Tuple{Vararg{Int}}}},
                  t::AbstractGraph, N::Int, levs_ks::Vector{Vector{Any}})

Performs filtered traversal; provides the given level `N`,
the current level `C`, and the vector of level-respective index sets,
`levs_ks` as arguments to `f`.

Call signature of `f` is: `f(dest, t::AbstractGraph, N, C, levs_ks)`.
Call signature of `fs[C]` is: `fs[C](p::Pair)`.
"""
function mapfilterupto(f::Function, fs::Vector{Function},
                       dims::Tuple{Vararg{Tuple{Vararg{Int}}}}, t::AbstractGraph{T},
                       N::Int, levs_ks::Vector{Vector{T}}) where {T}
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    mapfilterupto!(f, fs, dest, t, N, 1, levs_ks)
end

################ reduction of vector of into single dest
function mapupto!(f::Function, dest::Tuple{Vararg{AbstractArray}},
                  ts::Vector{<:AbstractGraph}, N::Int, C::Int)
    for t ∈ ts
        mapupto!(f, dest, t, N, C)
    end
    dest
end

function mapupto(f::Function, dims::Tuple{Vararg{Tuple{Vararg{Int}}}},
                 ts::Vector{<:AbstractGraph}, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    mapupto!(f, dest, ts, N, 1)
end

function tmapupto(f::Function, dims::Tuple{Vararg{Tuple{Vararg{Int}}}},
                  ts::Vector{<:AbstractGraph}, L::Int, M::Int=Threads.nthreads())
    N = length(ts)
    # M = Threads.threads()
    ranges = equalranges(N, M)
    A = Vector{Tuple{(Array{Int, length(d)} for d ∈ dims)...}}(undef, M)
    Threads.@threads for m = 1:M
        A[m] = mapupto(f, dims, ts[ranges[m]], L)
    end
    return reduce(.+, A)
end

#### filter
function mapfilterupto!(f::Function, fs::Vector{Function}, dest::Tuple{Vararg{AbstractArray}},
                        ts::Vector{<:AbstractGraph}, N::Int, C::Int)
    for t ∈ ts
        mapfilterupto!(f, fs, dest, t, N, C)
    end
    dest
end

function mapfilterupto(f::Function, fs::Vector{Function},
                       dims::Tuple{Vararg{Tuple{Vararg{Int}}}},
                       ts::Vector{<:AbstractGraph}, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    mapfilterupto!(f, fs, dest, ts, N, 1)
end

function tmapfilterupto(f::Function, fs::Vector{Function},
                        dims::Tuple{Vararg{Tuple{Vararg{Int}}}},
                        ts::Vector{<:AbstractGraph}, L::Int, M::Int=Threads.nthreads())
    N = length(ts)
    # M = Threads.threads()
    ranges = equalranges(N, M)
    A = Vector{Tuple{(Array{Int, length(d)} for d ∈ dims)...}}(undef, M)
    Threads.@threads for m = 1:M
        A[m] = mapfilterupto(f, fs, dims, ts[ranges[m]], L)
    end
    return reduce(.+, A)
end

#### levs_ks

function mapupto!(f::Function, dest::Tuple{Vararg{AbstractArray}},
                  ts::Vector{T}, N::Int, C::Int, levs_kss::Vector{Vector{Vector{U}}}) where {U, T<:AbstractGraph{U}}
    for i ∈ eachindex(ts)
        mapupto!(f, dest, ts[i], N, C, levs_kss[i])
    end
    dest
end

function mapupto(f::Function, dims::Tuple{Vararg{Tuple{Vararg{Int}}}},
                 ts::Vector{T}, N::Int, levs_kss::Vector{Vector{Vector{U}}}) where {U, T<:AbstractGraph{U}}
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    mapupto!(f, dest, ts, N, 1, levs_kss)
end

function tmapupto(f::Function, dims::Tuple{Vararg{Tuple{Vararg{Int}}}},
                  ts::Vector{T}, L::Int, levs_kss::Vector{Vector{Vector{U}}},
                  M::Int=Threads.nthreads()) where {U, T<:AbstractGraph{U}}
    N = length(ts)
    # M = Threads.threads()
    ranges = equalranges(N, M)
    A = Vector{Tuple{(Array{Int, length(d)} for d ∈ dims)...}}(undef, M)
    Threads.@threads for m = 1:M
        A[m] = mapupto(f, dims, ts[ranges[m]], L, levs_kss[ranges[m]])
    end
    return reduce(.+, A)
end

#### filter and levs_ks
function mapfilterupto!(f::Function, fs::Vector{Function}, dest::Tuple{Vararg{AbstractArray}},
                        ts::Vector{T}, N::Int, C::Int,
                        levs_kss::Vector{Vector{Vector{U}}}) where {U, T<:AbstractGraph{U}}
    for i ∈ eachindex(ts)
        mapfilterupto!(f, fs, dest, ts[i], N, C, levs_kss[i])
    end
    dest
end

function mapfilterupto(f::Function, fs::Vector{Function},
                       dims::Tuple{Vararg{Tuple{Vararg{Int}}}},
                       ts::Vector{T}, N::Int, levs_kss::Vector{Vector{Vector{U}}}) where {U, T<:AbstractGraph{U}}
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    mapfilterupto!(f, fs, dest, ts, N, 1, levs_kss)
end

function tmapfilterupto(f::Function, fs::Vector{Function},
                        dims::Tuple{Vararg{Tuple{Vararg{Int}}}},
                        ts::Vector{T}, L::Int, levs_kss::Vector{Vector{Vector{U}}},
                        M::Int=Threads.nthreads()) where {U, T<:AbstractGraph{U}}
    N = length(ts)
    # M = Threads.threads()
    ranges = equalranges(N, M)
    A = Vector{Tuple{(Array{Int, length(d)} for d ∈ dims)...}}(undef, M)
    Threads.@threads for m = 1:M
        A[m] = mapfilterupto(f, fs, dims, ts[ranges[m]], L, levs_kss[ranges[m]])
    end
    return reduce(.+, A)
end

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

function mapat(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S}, t::AbstractNode, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    mapat!(f, dest, t, N, 1)
end

# Variant with filter
function mapat!(f::Function, fs::Vector{Function}, dest::Tuple{Vararg{Array{T}} where T},
                t::AbstractNode, N::Int, C::Int)
    isempty(t) && return dest
    C̃ = C + 1
    g = fs[C]
    if C̃ < N
        for p in t
            g(p) && mapat!(f, fs, dest, p.second, N, C̃)
        end
    elseif C̃ == N
        for p in t
            g(p) && f(dest, p)
        end
    end
    dest
end

function mapat(f::Function, fs::Vector{Function}, dims::Tuple{Vararg{NTuple{S, Int}} where S},
               t::AbstractNode, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    mapat!(f, fs, dest, t, N, 1)
end

# Variant with ks
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
            f(dest, p, ks)
        end
    end
    dest
end

function indexmapat(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                    t::AbstractNode, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    ks = Vector{Any}(undef, N - 1)
    indexmapat!(f, dest, ks, t, N, 1)
end

# Variant with ks and filter
function indexmapat!(f::Function, fs::Vector{Function}, dest::Tuple{Vararg{Array{T}} where T},
                     ks::Vector, t::AbstractNode, N::Int, C::Int)
    isempty(t) && return dest
    C̃ = C + 1
    g = fs[C]
    if C̃ < N
        for p in t
            g(p) && (setindex!(ks, p.first, C); indexmapat!(f, fs, dest, p.second, N, C̃))
        end
    elseif C̃ == N
        for p in t
            g(p) && (setindex!(ks, p.first, C); f(dest, p, ks))
        end
    end
    dest
end

function indexmapat(f::Function, fs::Vector{Function}, dims::Tuple{Vararg{NTuple{S, Int}} where S},
               t::AbstractNode, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    ks = Vector{Any}(undef, N - 1)
    indexmapat!(f, fs, dest, ks, t, N, 1)
end

####

# p. 501
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

function mapupto(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S}, t::AbstractNode, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    mapupto!(f, dest, t, N, 1)
end

# filter variant
function mapupto!(f::Function, fs::Vector{Function}, dest::Tuple{Vararg{Array{T}} where T},
                  t::AbstractNode, N::Int, C::Int)
    C̃ = C + 1
    C̃ < N || return dest
    f(dest, t)
    isempty(t) && return dest
    g = fs[C]
    for p in t
        g(p) && mapupto!(f, fs, dest, p.second, N, C̃)
    end
    dest
end

function mapupto(f::Function, fs::Vector{Function}, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                 t::AbstractNode, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    mapupto!(f, fs, dest, t, N, 1)
end

# levs_ks variant
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

function mapupto(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                 t::AbstractNode, N::Int, levs_ks::Vector{Vector{Any}})
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    mapupto!(f, dest, t, N, 1, levs_ks)
end

# filter and levs_ks
function mapupto!(f::Function, fs::Vector{Function}, dest::Tuple{Vararg{Array{T}} where T},
                  t::AbstractNode, N::Int, C::Int, levs_ks::Vector{Vector{Any}})
    C̃ = C + 1
    C̃ < N || return dest
    f(fs, dest, t, N, C, levs_ks)
    isempty(t) && return dest
    g = fs[C]
    for p in t
        g(p) && mapupto!(f, fs, dest, p.second, N, C̃, levs_ks)
    end
    dest
end

function mapupto(f::Function, fs::Vector{Function}, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                 t::AbstractNode, N::Int, levs_ks::Vector{Vector{Any}})
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    mapupto!(f, fs, dest, t, N, 1, levs_ks)
end

# p. 500
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

function indexmapupto(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                      t::AbstractNode, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    ks = Vector{Any}(undef, N - 1)
    indexmapupto!(f, dest, ks, t, N, 1)
end

# Variant with levs_ks p. 502 -- dispatch on lev_aks
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

function indexmapupto(f::Function, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                      t::AbstractNode, N::Int, levs_ks::Vector{Vector{Any}})
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    ks = Vector{Any}(undef, N - 1)
    indexmapupto!(f, dest, ks, t, N, 1, levs_ks)
end

#### filter variants -- p. 505-507

# Variant with filter
function indexmapupto!(f::Function, fs::Vector{Function}, dest::Tuple{Vararg{Array{T}} where T},
                       ks::Vector, t::AbstractNode, N::Int, C::Int)
    C̃ = C + 1
    C̃ < N || return dest
    f(fs, dest, ks, t, N, C)
    isempty(t) && return dest
    g = fs[C]
    for p in t
        g(p) && (setindex!(ks, p.first, C); indexmapupto!(f, fs, dest, ks, p.second, N, C̃))
    end
    dest
end

function indexmapupto(f::Function, fs::Vector{Function}, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                      t::AbstractNode, N::Int)
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    ks = Vector{Any}(undef, N - 1);
    indexmapupto!(f, fs, dest, ks, t, N, 1)
end

# Variant with filter and levs_ks
function indexmapupto!(f::Function, fs::Vector{Function}, dest::Tuple{Vararg{Array{T}} where T},
                       ks::Vector, t::AbstractNode, N::Int, C::Int, levs_ks::Vector{Vector{Any}})
    C̃ = C + 1
    C̃ < N || return dest
    f(fs, dest, ks, t, N, C, levs_ks)
    isempty(t) && return dest
    g = fs[C]
    for p in t
        g(p) && (setindex!(ks, p.first, C); indexmapupto!(f, fs, dest, ks, p.second, N, C̃, levs_ks))
    end
    dest
end

function indexmapupto(f::Function, fs::Vector{Function}, dims::Tuple{Vararg{NTuple{S, Int}} where S},
                      t::AbstractNode, N::Int, levs_ks::Vector{Vector{Any}})
    dest = ntuple(i -> zeros(Int, dims[i]), length(dims))
    ks = Vector{Any}(undef, N - 1)
    indexmapupto!(f, fs, dest, ks, t, N, 1, levs_ks)
end

#
# Date created: 2021-10-14
# Author: aradclif
#
#
############################################################################################
#### p. 532-533, 2021-10-20
# countabsent! as a closure:
# Given a f::Function, ca = (dest, N, C, levs_ks) -> countabsent!(f, dest, N, C, levs_ks)
# Given a f::Function and fs::Vector{Function},
# fca = (dest, N, C, levs_ks) -> countabsent!(f, dest, N, C, levs_ks)
##
# Then, it is used as:
# mapupto(ca, dims, t, N)                OR    mapupto(ca, dims, t, N, levs_ks)
# mapfilterupto(fca, mfs, dims, t, N)    OR    mapfilterupto(fca, mfs, dims, t, N, levs_ks)
##
# The f itself is a closure around some increment!(g, dest, ν, ks) function, e.g.
# f = (dest, ν, ks) -> incrementrows!(g, dest, ν, ks)
# Or, given that we know dest to be a Tuple{Vararg{Array{T}} where T}, it might be:
# f = (dest, ν, ks) -> incrementrows!(g, dest[2], ν, ks)
################
function countabsent!(f::Function, dest::Tuple{Vararg{Array{T}} where T},
                      t::AbstractGraph, N::Int, C::Int, levs_ks::Vector{Vector{Any}})
    Ñ = N - 1
    if C == Ñ
        lev_ks = setdiff(levs_ks[C], keys(t))
        ν = 1
        isempty(lev_ks) || f(dest, ν, lev_ks)
    elseif C < Ñ
        ν = multabsentupto(t, N, C, levs_ks)
        iszero(ν) || f(dest, ν, levs_ks[Ñ])
    end
    dest
end

function countabsent!(f::Function, fs::Vector{Function}, dest::Tuple{Vararg{Array{T}} where T},
                      t::AbstractGraph, N::Int, C::Int, levs_ks::Vector{Vector{Any}})
    Ñ = N - 1
    if C == Ñ
        lev_ks = filter!(fs[C], setdiff(levs_ks[C], keys(t)))
        ν = 1
        isempty(lev_ks) || f(dest, ν, lev_ks)
    elseif C < Ñ
        ν = multabsentupto(fs, t, N, C, levs_ks)
        iszero(ν) || f(dest, ν, filter(fs[Ñ], levs_ks[Ñ]))
    end
    dest
end
################################################################
#### p. 540-545, 2021-10-25
# levs_ks::Vector{Vector{Any}} ∈ 𝔻²
# dest::Tuple{Vararg{Array{T}} where T} ∈ 𝔻², with dest[2] ∈ ℕᴵ, dest[1] ∈ ℕᴵˣᴶ
# f::Function : closure around some increment!(g, dest, ν, ks)
## e.g. of countabsent!
# db::Dict{String,Int}
# g = let db = db
#     k -> db[k]
# end
# f = (dest, ν, ks) -> incrementrows!(g, dest[2], ν, ks)
## e.g. of countstatus!
# ds::Dict{String, Int}
# h = let ds = ds
#     k -> ds[k]
# end
# f = (dest, ν, rk, ck) -> incrementrowcol!(g, h, dest[1], ν, rk, ck)
################
function countstatus!(f::Function, dest::Tuple{Vararg{Array{T}} where T}, k, t::AbstractGraph)
    status = t.data[1]
    f(dest, 1, k, status)
    dest
end
function countstatus!(f::Function, dest::Tuple{Vararg{Array{T}} where T}, p::Pair)
    countstatus!(f, dest, p.first, p.second)
end
################
## e.g. of kcountstatus!
# dd::Dict{NTuple{2, Int}, Int}
# d = let dd = dd
#     (k) -> dd[k]
# end
# db::Dict{String, Int}
# g = let db = db
#     (k) -> db[k]
# end
# ds::Dict{NTuple{3, String}, Int}
# h = let ds = ds
#     (k) -> ds[k]
# end
# f = (dest, ν, ks...) -> incrementnd!([d, g, h], dest[1], ν, ks...)
function kcountstatus!(f::Function, dest::Tuple{Vararg{Array{T}} where T}, ks::Vector,
                       g::AbstractGraph)
    status = g.data[1]
    f(dest, 1, ks..., status)
    dest
end
#### Usage example
# cs = (dest, ks, g) -> kcountstatus!(f, dest, ks, g)
# mapat(cs, ((282, 47, 7),), g, 3)
############################################################################################
#### 2021-11-04: p. 564-571
# Increment count(s) for the absent vertices, indexing each level to the respective
# dimension in a multidimensional array. Only the Cᵗʰ level edges are used for iteration,
# with the higher dimensions treated in a generic fashion. Lower dimensions
# are ignored by creating a view which includes the current level and all higher dimensions.
##
# This version counts on all dimensions.
#### Usage example
# fs = [g, h]
# ca = (dest, ks, x, N, C, levs_ks) -> kcountabsent!(fs, dest[1], x, ks, N, C, levs_ks)
# mapupto(ca, ((440, 47),), g, 3)
function kcountabsent!(fs::Vector{Function}, A::AbstractArray, x::AbstractGraph, ks::Vector{Any},
                       N::Int, C::Int, levs_ks::Vector{Vector{Any}})
    # Option 1: single view
    mks = setdiff(levs_ks[C], keys(x))
    isempty(mks) && return A
    idxs = ntuple(i -> fs[i](ks[i]), C - 1)
    colons = ntuple(i -> :, N - C - 1)
    for k ∈ mks
        idx = fs[C](k)
        Ã = view(A, idxs..., idx, colons...)
        for i ∈ eachindex(Ã)
            Ã[i] += 1
        end
    end
    return A
    # Option 2: multiple view
    # mks = setdiff(levs_ks[C], keys(x))
    # isempty(mks) && return A
    # idxs = ntuple(i -> fs[i](ks[i]), C - 1)
    # colons = ntuple(i -> :, N - C)
    # Ã = view(A, idxs..., colons...)
    # colons₂ = ntuple(i -> :, N - C - 1)
    # for k ∈ mks
    #     idx = fs[C](k)
    #     Ã̃ = view(Ã, idx, colons2...)
    #     for i ∈ eachindex(Ã̃)
    #         Ã̃[i] += 1
    #     end
    # end
    # return A
    # Option 3: dispatching to increment function for type-stability
    # _vinc!(tf[C], A, mks, idxs, colons, 1)
end

# A = reshape([1:24;], 3, 4, 2)
# idxs = (1,)
# colons = (:, :)
# Ã = view(A, idxs..., colons...)
# colons2 = (:,)
# idx = 3
# Ã̃ = view(Ã, idx, colons2...)
# for i ∈ eachindex(Ã̃)
#     Ã̃[i] += 1
# end
# Ã₂ = view(A, idxs..., idx, colons2...)

################ The concept
# The idea is to specify the indices preceding the dimension corresponding to the
# current level, then iterate through the current-level-dimension, handling
# all higher dimensions implicitly (via a view).
# A ∈ ℝᴵˣᴶˣᴷˣ⋯ → ndims(A) = N - 1
# indices ∈ 𝔻ᴺ⁻¹
# ρ = preceding indices    ∈ 𝔻ᶜ⁻¹
# ξ = current index        ∈ 𝔻¹
# η = higher indices       ∈ 𝔻ᴺ⁻ᶜ⁻¹, noting: N - 1 = (C - 1) + 1 + H ⇒ H = N - C - 1
# With E ∈ 𝔻ᴺ⁻¹, under full indexing of the result, the ρ, ξ, η can always be constructed.
## What if the number of dimensions of A and dimensionality of E do not match?
# e.g E = edge1-edge2-edge3-edge4-edge5 ∈ 𝔻⁵⁼ᴺ⁻¹, but A ∈ ℝᴵˣᴶˣᴷ,
# with I => edge1, J => edge4, K => edge5.
# This is in fact a commonly desired result: indexing of counts to only a subset of
# the indices. One can deal with this in generic manner by building on the convention
# outlined above, assigning unit dimensions to each of the collapsed dimensions.
# In the motivating example, A would become A ∈ ℝᴵˣ¹ˣ¹ˣᴶˣᴷ. This enables the
# ρ, ξ, η convention to work, given fs that return 1 for the dimensions that are unity.
# Thus, the ρ can always be formed, and so can ξ and η.
# Now the increment factor needs to be dynamically computed, as the collapsed
# dimensions must still be included. In essence, a multiplier should
# cover the unit dimensions between C and the next non-unit dimension.
# How to know the next non-unit dimension? It can be deduced from the
# dimensions of A. Once the index of the next-non-unit (NNU) dimension is known,
# computing the multiplier which encompasses elements between C and NNU
# is straightforward, given the levs_ks.
#### Extensions
# Apply this concept to handle all multi-indexing cases: kcountstatus!, and so on.
################

function nextnonunit(dims::NTuple{M, Int}, C::Int) where {M}
    # dims[C + 1] != 1 ? C + 1 : nextnonunit2(dims, C + 1)
    C̃ = C + 1
    dims[C̃] != 1 ? C̃ : nextnonunit(dims, C̃)
end

function dimsmultiplier(dims::NTuple{M, Int}, N::Int, C::Int, levs_ks::Vector{Vector{Any}}) where {M}
    # Could optimize via rearrangement
    # nnu = C ≥ N - 1 ? N - 1 : nextnonunit(dims, C)
    # C̃ = C + 1
    # ν = 1
    # # Slightly slower
    # C̃ = C + 1
    # ν = 1
    # C̃ ≥ N && return ν
    # Slightly faster (≈ 10%)
    C ≥ N - 1 && return 1
    C̃ = C + 1
    ν = 1
    # Slightly faster, but less than above (≈ 5%)
    # C̃ = C + 1
    # C̃ ≥ N && return 1
    # ν = 1
    nnu = nextnonunit(dims, C)
    while C̃ < nnu
        ν *= length(levs_ks[C̃])
        C̃ += 1
    end
    return ν
end

# diml = (440, 1, 1, 32, 32)
# levs_ks5 = Vector{Any}[[1:440;], [1:10;], [1:10;], [1:32;], [1:32;]];
# @benchmark dimsmultiplier(diml, 5, 2, levs_ks5)

#### 2021-11-05: p. 572-575
function kcountabsent!(tfs::Vector{Function}, dims::NTuple{M, Int}, A::Array{T, M},
                       x::AbstractGraph, ks::Vector{Any},
                       N::Int, C::Int, levs_ks::Vector{Vector{Any}}) where {M} where {T<:Number}
    mks = setdiff(levs_ks[C], keys(x))
    isempty(mks) && return A
    idxs = ntuple(i -> tfs[i](ks[i]), C - 1)
    colons = ntuple(i -> :, N - C - 1)
    ν = dimsmultiplier(dims, N, C, levs_ks)
    for k ∈ mks
        idx = tfs[C](k)
        Ã = view(A, idxs..., idx, colons...)
        for i ∈ eachindex(Ã)
            Ã[i] += ν
        end
    end
    return A
    # _vinc!(tf[C], A, mks, idxs, colons, 1)
end
#### Usage example
# fs = [g, h]
# ca = (dest, ks, x, N, C, levs_ks) -> kcountabsent!(fs, (440, 47), dest[1], x, ks, N, C, levs_ks)
# mapupto(ca, ((440, 47),), g, 3)
@code_warntype kcountabsent!(fsₐ₂, (282, 47), As2₂[1], gti, ks_e, 3, 1, levs_ks)

################################################################
#### 2021-11-05: metaprogramming experiments
@generated function _idxs(tfs::NTuple{N, Function}, ks::Vector{Any}) where {N}
    quote
        Base.Cartesian.@ntuple $N i -> tfs[i](ks[i])
    end
end

ks_e = Any[kk₁..., chainu]
tf = (uf, nufₗᵢₙ, h)

_idxs(tf, ks_e)
################

@generated function _nidxs(tfs::NTuple{N, Function}, ks::Vector{Any}, C::Val{M})::NTuple{M, Int} where {N} where {M}
    quote
        Base.Cartesian.@ntuple $M i -> tfs[i](ks[i])
    end
end

@benchmark _nidxs(tf, ks_e, Val(2))

# # Not type-stable.
# @generated function _nidxs(tfs::Vector{Function}, ks::Vector{Any}, C::Val{M})::NTuple{M, Int} where {N} where {M}
#     quote
#         Base.Cartesian.@ntuple $M i -> tfs[i](ks[i])
#     end
# end
# _nidxs(fsᵢ, ks_e, Val(2))

function _nidxsplain(tfs::NTuple{N, Function}, ks::Vector{Any}, C::Int) where {N}
    ntuple(i -> tfs[i](ks[i]), C)
end
# function _nidxsplain(tfs::Vector{Function}, ks::Vector{Any}, C::Int) where {N}
#     ntuple(i -> tfs[i](ks[i]), C)
# end
@benchmark _nidxsplain(tf, ks_e, 2)

function kcountabsent!(tfs::NTuple{M, Function}, dims::NTuple{M, Int}, A::Array{T, M},
                       x::AbstractGraph, ks::Vector{Any},
                       N::Int, C::Int, levs_ks::Vector{Vector{Any}}) where {M} where {T<:Number}
    mks = setdiff(levs_ks[C], keys(x))
    isempty(mks) && return A
    idxs = _nidxs(tfs, ks, Val(C - 1))
    colons = ntuple(i -> :, N - C - 1)
    ν = dimsmultiplier(dims, N, C, levs_ks)
    # tf = tfs[C]
    # for k ∈ mks
    #     # idx = tfs[C](k)
    #     idx::Int = tf(k)
    #     Ã = view(A, idxs..., idx, colons...)
    #     for i ∈ eachindex(Ã)
    #         Ã[i] += ν
    #     end
    # end
    # return A
    _vinc!(tf[C], A, mks, idxs, colons, ν)
end

function _vinc!(tf::Function, A::Array{T, N}, mks::Vector{S}, idxs::NTuple{M, Int}, colons::NTuple{H, Colon}, ν::Int) where {T, N} where {S} where {M} where {H}
    for k ∈ mks
        idx = tf(k)
        Ã = view(A, idxs..., idx, colons...)
        for i ∈ eachindex(Ã)
            Ã[i] += ν
        end
    end
    return A
end
@code_warntype _vinc!(tf₂[1], Aex, levs_ks[1], (), (:,), 1)
@benchmark _vinc!(tf₂[1], Aex, levs_ks[1], (), (:,), 1)
@timev _vinc!(tf₂[1], Aex, levs_ks[1], (), (:,), 1);



tf₂ = (uf, nufₗᵢₙ)
Aex = As2₂[1]
@code_warntype kcountabsent!(tf₂, (282, 91), Aex, gti, ks_e, 3, 1, levs_ks)

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
cs = (dest, ks, g) -> kcountstatus!(f, dest, ks, g)
mapat(cs, ((282, 47, 7),), g, 3)
############################################################################################
#### 2021-11-04: p. 564-571
# Increment count(s) for the absent vertices, indexing each level to the respective
# dimension in a multidimensional array. Only the Cᵗʰ level edges are used for iteration,
# with the higher dimensions treated in a generic fashion. Lower dimensions
# are ignored by creating a view which includes the current level and all higher dimensions.
##
# This version counts on all dimensions.
#### Usage example
fs = [g, h]
ca = (dest, ks, x, N, C, levs_ks) -> kcountabsent!(fs, dest[1], x, ks, N, C, levs_ks)
mapupto(ca, ((440, 47),), g, 3)
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
end

A = reshape([1:24;], 3, 4, 2)
idxs = (1,)
colons = (:, :)
Ã = view(A, idxs..., colons...)
colons2 = (:,)
idx = 3
Ã̃ = view(Ã, idx, colons2...)
for i ∈ eachindex(Ã̃)
    Ã̃[i] += 1
end
Ã₂ = view(A, idxs..., idx, colons2...)

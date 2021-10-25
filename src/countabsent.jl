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
                      t::AbstractNode, N::Int, C::Int, levs_ks::Vector{Vector{Any}})
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
                      t::AbstractNode, N::Int, C::Int, levs_ks::Vector{Vector{Any}})
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
function countstatus!(f::Function, dest::Tuple{Vararg{Array{T}} where T}, k, t::AbstractNode)
    status = t.val[1]
    f(dest, 1, k, status)
    dest
end
function countstatus!(f::Function, dest::Tuple{Vararg{Array{T}} where T}, p::Pair)
    countstatus!(f, dest, p.first, p.second)
end

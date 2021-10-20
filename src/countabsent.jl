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
# The f itself is a closure around some increment!(g, dest, ks, ν) function, e.g.
# f = (dest, ks, ν) -> incrementrows!(g, dest, ks, ν)
# Or, given that we know dest to be a Tuple{Vararg{Array{T}} where T}, it might be:
# f = (dest, ks, ν) -> incrementrows!(g, dest[2], ks, ν)
################
function countabsent!(f::Function, dest::Tuple{Vararg{Array{T}} T},
                      t::AbstractNode, N::Int, C::Int, levs_ks::Vector{Vector{Any}})
    Ñ = N - 1
    if C == Ñ
        lev_ks = setdiff(levs_ks[C], keys(t))
        ν = 1
        isempty(lev_ks) || f(dest, lev_ks, ν)
    elseif C < Ñ
        ν = multabsentupto(t, N, C, levs_ks)
        iszero(ν) || f(dest, levs_ks[Ñ], ν)
    end
    dest
end

function countabsent!(f::Function, fs::Vector{Function}, dest::Tuple{Vararg{Array{T}} T},
                      t::AbstractNode, N::Int, C::Int, levs_ks::Vector{Vector{Any}})
    Ñ = N - 1
    if C == Ñ
        lev_ks = filter!(fs[C], setdiff(levs_ks[C], keys(t)))
        ν = 1
        isempty(lev_ks) || f(dest, lev_ks, ν)
    elseif C < Ñ
        ν = multabsentupto(fs, t, N, C, levs_ks)
        iszero(ν) || f(dest, filter(fs[Ñ], levs_ks[Ñ]), ν)
    end
    dest
end

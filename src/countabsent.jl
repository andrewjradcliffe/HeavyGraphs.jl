#
# Date created: 2021-10-14
# Author: aradclif
#
#
############################################################################################

function countabsent!(f::Function, dest::Tuple{Vararg{Array{T}} T},
                      t::AbstractNode, N::Int, C::Int, levs_ks::Vector{Vector{Any}})
    Ñ = N - 1
    if C == Ñ
        lev_ks = setdiff(levs_ks[C], keys(t))
        ν = 1
        isempty(ks) || f(dest, lev_ks, ν)
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
        isempty(ks) || f(dest, lev_ks, ν)
    elseif C < Ñ
        ν = multabsentupto(fs, t, N, C, levs_ks)
        iszero(ν) || f(dest, filter(fs[Ñ], levs_ks[Ñ]), ν)
    end
    dest
end

#
# Date created: 2021-10-13
# Author: aradclif
#
#
############################################################################################
#### p. 517 - 522, 2021-10-13

"""
    multabsentupto(t::AbstractNode, N::Int, C::Int, levs_ks::Vector{Vector{Any}})

Return the number of `N-1`th level nodes which would be present if the current
level `C`'s absent branches were instantiated at full density. The fully dense
graph/tree would be the Cartesian product of the sets `levs_ks`; the linear
index of each component of `levs_ks` corresponds to the level index, with
total levels counted as `length(levs_ks) + 1 ≡ N`.

See also: [`multabsent`](@ref)

"""
function multabsentupto(t::AbstractNode, N::Int, C::Int, levs_ks::Vector{Vector{Any}})
    Ñ = N - 1
    C ≥ Ñ && return 0
    # ν = length(setdiff(levs_ks[C], keys(t)))
    # much more efficient
    ν = length(levs_ks[C]) - length(t) + count(∉(levs_ks[C]), keys(t))
    C̃ = C + 1
    while C̃ < Ñ
        ν *= length(levs_ks[C̃])
        C̃ += 1
    end
    ν
end

"""
    multabsentupto(fs::Vector{Function}, t::AbstractNode, N::Int, C::Int,
                   levs_ks::Vector{Vector{Any}})

Return the number of `N-1`th level nodes which would be present if the current
level `C`'s absent branches were instantiated using the subsets produced by applying
the level-indexed filtering functions `fs`.

"""
function multabsentupto(fs::Vector{Function}, t::AbstractNode, N::Int, C::Int,
                        levs_ks::Vector{Vector{Any}})
    Ñ = N - 1
    C ≥ Ñ && return 0
    ν = count(fs[C], setdiff(levs_ks[C], keys(t)))
    C̃ = C + 1
    while C̃ < Ñ
        ν *= count(fs[C̃], levs_ks[C̃])
        C̃ += 1
    end
    ν
end

multabsentupto(t::AbstractNode, levs_ks::Vector{Vector{Any}}) =
    multabsentupto(t, length(levs_ks), 1, levs_ks)

multabsentupto(fs::Vector{Function}, t::AbstractNode, levs_ks::Vector{Vector{Any}}) =
    multabsentupto(fs, t, length(levs_ks), 1, levs_ks)

# Optimized option for non-filtering version
function multabsentupto(t::AbstractNode, N::Int, C::Int, levs_ks::Vector{Int})
    Ñ = N - 1
    C ≥ Ñ && return 0
    ν = levs_ks[C] - length(t)
    C̃ = C + 1
    while C̃ < Ñ
        ν *= levs_ks[C̃]
        C̃ += 1
    end
    ν
end
####

"""
    multabsent(t::AbstractNode, N::Int, C::Int, levs_ks::Vector{Vector{Any}})

Return the number of `N`th level nodes which would be present if the current
level `C`'s absent branches were instantiated at full density. The fully dense
graph/tree would be the Cartesian product of the sets `levs_ks`; the linear
index of each component of `levs_ks` corresponds to the level index, with
total levels counted as `length(levs_ks) + 1 ≡ N`.

See also: [`multabsentupto`](@ref)

"""
function multabsent(t::AbstractNode, N::Int, C::Int, levs_ks::Vector{Vector{Any}})
    Ñ = N - 1
    C > Ñ && return 0
    # ν = length(setdiff(levs_ks[C], keys(t)))
    # much more efficient
    ν = length(levs_ks[C]) - length(t) + count(∉(levs_ks[C]), keys(t))
    C̃ = C + 1
    while C̃ ≤ Ñ
        ν *= length(levs_ks[C̃])
        C̃ += 1
    end
    ν
end

"""
    multabsent(fs::Vector{Function}, t::AbstractNode, N::Int, C::Int,
               levs_ks::Vector{Vector{Any}})

Return the number of `N`th level nodes which would be present if the current
level `C`'s absent branches were instantiated using the subsets produced by applying
the level-indexed filtering functions `fs`.

"""
function multabsent(fs::Vector{Function}, t::AbstractNode, N::Int, C::Int,
                    levs_ks::Vector{Vector{Any}})
    Ñ = N - 1
    C > Ñ && return 0
    ν = count(fs[C], setdiff(levs_ks[C], keys(t)))
    C̃ = C + 1
    while C̃ ≤ Ñ
        ν *= count(fs[C̃], levs_ks[C̃])
        C̃ += 1
    end
    ν
end

multabsent(t::AbstractNode, levs_ks::Vector{Vector{Any}}) =
    multabsent(t, length(levs_ks), 1, levs_ks)

multabsent(fs::Vector{Function}, t::AbstractNode, levs_ks::Vector{Vector{Any}}) =
    multabsent(fs, t, length(levs_ks), 1, levs_ks)

# Optimized option for non-filtering version
function multabsent(t::AbstractNode, N::Int, C::Int, levs_ks::Vector{Int})
    Ñ = N - 1
    C > Ñ && return 0
    ν = levs_ks[C] - length(t)
    C̃ = C + 1
    while C̃ ≤ Ñ
        ν *= levs_ks[C̃]
        C̃ += 1
    end
    ν
end
####

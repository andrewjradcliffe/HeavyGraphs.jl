#
# Date created: 2021-10-13
# Author: aradclif
#
#
############################################################################################
#### p. 517 - 522, 2021-10-13; see also: p. 530-533, 2021-10-20
# The index of levs_ks corresponds to the true level depth. The set held at a given
# index represents the set of all keys which a node at that level may hold.
# If the terminal node were at N = 5, then the multiplier up to it would be:
# |level_1| × |level_2| × |level_3|
# If one desired the multiplier for all levels, it would be:
# |level_1| × |level_2| × |level_3| × |level_4|
# Why? Because the Nᵗʰ level is fully circumscribed by the 1,…,N-1 levels preceding it.
# If one were at the N=4ᵗʰ level, then the multiplier up to 5 is automatically 0
# as one is sitting on the last level before 5. Correspondingly, the full multiplier
# at N=4ᵗʰ would be |level_4 \ edges(node)|.
##
# The use of countabsent! always involves the multiplier up to the Nᵗʰ level as the Nᵗʰ
# level is the first unique level. Hence, one needs to determine the identities of the missing
# nodes at the Nᵗʰ level in order to take appropriate actions. -- The multiplier up to the Nᵗʰ
# is useful in that it can be combined with a set of identities. The multiplier up to the Nᵗʰ
# is the number of time each element in that set would be instantiated.
#     multiplier up to Nᵗʰ × {Nᵗʰ level identities} ≡ νₙ₋₁ × {Nᵗʰ level identities}
#     multiplier through Nᵗʰ = νₙ₋₁ × |{Nᵗʰ level identities}| ≡ νₙ
# Thus, the multiplier through Nᵗʰ would destroy the unique identities of the Nᵗʰ level.
##
# The above simplifies the concept and code.
# multabsentupto(t, N, C, levs_ks) = multabsent(t, N - 1, C, levs_ks)
#                                  ↕
#     multabsent(t, N, C, levs_ks) = multabsentupto(t, N + 1, C, levs_ks)
##
# Concept: use in countabsent!
# Determine if already at unique level, i.e. C == Ñ. If so, then compute
# {level_N-1} \ {edges{node}}. Otherwise, find the multiplier up to Nᵗʰ level, νₙ₋₁;
# automatically, one knows that the set of identities is {level_N-1} which is the set of
# identities stored at the N-1ᵗʰ index of levs_ks.
################
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
    # much more efficient -- or: setdiffcard(levs_ks[C], keys(t))
    # ν = length(levs_ks[C]) - length(t) + count(∉(levs_ks[C]), keys(t))
    ν = length(levs_ks[C]) - count(∈(levs_ks[C]), keys(t)) # |A| - |A ∩ B|
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
multabsent(t::AbstractNode, N::Int, C::Int, levs_ks::Vector{Vector{Any}}) =
    multabsentupto(t, N + 1, C, levs_ks)
# function multabsent(t::AbstractNode, N::Int, C::Int, levs_ks::Vector{Vector{Any}})
#     Ñ = N - 1
#     C > Ñ && return 0
#     # ν = length(setdiff(levs_ks[C], keys(t)))
#     # much more efficient -- or: setdiffcard(levs_ks[C], keys(t))
#     # ν = length(levs_ks[C]) - length(t) + count(∉(levs_ks[C]), keys(t))
#     ν = length(levs_ks[C]) - count(∈(levs_ks[C]), keys(t)) # |A| - |A ∩ B|
#     C̃ = C + 1
#     while C̃ ≤ Ñ
#         ν *= length(levs_ks[C̃])
#         C̃ += 1
#     end
#     ν
# end

"""
    multabsent(fs::Vector{Function}, t::AbstractNode, N::Int, C::Int,
               levs_ks::Vector{Vector{Any}})

Return the number of `N`th level nodes which would be present if the current
level `C`'s absent branches were instantiated using the subsets produced by applying
the level-indexed filtering functions `fs`.

"""
multabsent(fs::Vector{Function}, t::AbstractNode, N::Int, C::Int, levs_ks::Vector{Vector{Any}}) =
    multabsentupto(fs, t, N + 1, C, levs_ks)
# function multabsent(fs::Vector{Function}, t::AbstractNode, N::Int, C::Int,
#                     levs_ks::Vector{Vector{Any}})
#     Ñ = N - 1
#     C > Ñ && return 0
#     ν = count(fs[C], setdiff(levs_ks[C], keys(t)))
#     C̃ = C + 1
#     while C̃ ≤ Ñ
#         ν *= count(fs[C̃], levs_ks[C̃])
#         C̃ += 1
#     end
#     ν
# end

multabsent(t::AbstractNode, levs_ks::Vector{Vector{Any}}) =
    multabsent(t, length(levs_ks), 1, levs_ks)

multabsent(fs::Vector{Function}, t::AbstractNode, levs_ks::Vector{Vector{Any}}) =
    multabsent(fs, t, length(levs_ks), 1, levs_ks)

# Optimized option for non-filtering version
multabsent(t::AbstractNode, N::Int, C::Int, levs_ks::Vector{Int}) =
    multabsentupto(t, N + 1, C, levs_ks)
# function multabsent(t::AbstractNode, N::Int, C::Int, levs_ks::Vector{Int})
#     Ñ = N - 1
#     C > Ñ && return 0
#     ν = levs_ks[C] - length(t)
#     C̃ = C + 1
#     while C̃ ≤ Ñ
#         ν *= levs_ks[C̃]
#         C̃ += 1
#     end
#     ν
# end
####

#
# Date created: 2021-10-20
# Author: aradclif
#
#
############################################################################################
#### Methods of grow_(:field)! : see p. 443-446, 451-455, 2021-09-14/15
# grow! -> get! only has 2 signatures:
# (f::Function, t::AbstractNode, p::AbstractEdges)
# (f::Function, t::AbstractNode, p::AbstractEdges, itr)
# grow_val! -> get_pushval! only has 2 signatures:
# (f::Function, t::AbstractNode, v::AbstractEdges, p::AbstractEdges)
# (f::Function, t::AbstractNode, v::AbstractEdges, p::AbstractEdges, itr)
# grow_spec! -> get_pushspec! has 2 signatures, same as grow_val!
# grow_sval! -> get_pushsval! has 2 signatures, same as grow_val!

############################################################################################
#### Methods of grow_(:field)! : see p. 443-446, 451-455, 2021-09-14/15
function grow!(f::Function, g::AbstractGraph, p::AbstractEdges)
    for x ∈ p
        get!(f, g, x...)
    end
    return g
end

# Replaces splatting method
@generated function grow!(f::Function, g::AbstractGraph, p::AbstractEdges{U, N}, itr) where {U, N}
    quote
        # for item ∈ itr
        #     tmp = g
        #     Base.Cartesian.@nexprs $N i -> tmp = get!(f, tmp, p[i](item))
        # end
        Base.Cartesian.@nextract $N e i -> p[i]
        for item ∈ itr
            tmp = g
            Base.Cartesian.@nexprs $N i -> tmp = get!(f, tmp, e_i(item))
        end
        return g
    end
end

# Convenience wrappers, but useful nonetheless
grow(f::Function, p::AbstractEdges) = grow!(f, f(), p)
grow(f::Function, p::AbstractEdges, itr) = grow!(f, f(), p, itr)

function datagrow!(f::Function, g::AbstractGraph, v::T, p::AbstractEdges) where {T<:Union{AbstractLabel, AbstractLabels}}
    for x ∈ p
        get_datapush!(f, g, v(x), x...)
    end
    return g
end

@generated function datagrow!(f::Function, g::AbstractGraph, v::T,
                              p::AbstractEdges{U, N}, itr) where {U, N} where {T<:Union{AbstractLabel, AbstractLabels}}
    quote
        # for item ∈ itr
        #     tmp = g
        #     Base.Cartesian.@nexprs $N i -> tmp = get!(f, tmp, p[i](item))
        #     push!(tmp.data, v(item))
        # end
        Base.Cartesian.@nextract $N e i -> p[i]
        for item ∈ itr
            tmp = g
            Base.Cartesian.@nexprs $N i -> tmp = get!(f, tmp, e_i(item))
            push!(tmp.data, v(item))
        end
        return g
    end
end


# Convenience wrappers, but useful nonetheless
datagrow(f::Function, v, p) = datagrow!(f, f(), v, p)
datagrow(f::Function, v, p, itr) = datagrow!(f, f(), v, p, itr)

# 2022-01-10: multi-step growth
function datagrow!(f::Function, x::AbstractGraph, vs::Vector{<:AbstractEdge},
                   ps::Vector{<:AbstractEdges}, itrs::Vector)
    for i ∈ eachindex(vs, ps, itrs)
        datagrow!(f, x, vs[i], ps[i], itrs[i])
    end
    return x
end
datagrow(f::Function, vs::Vector{<:AbstractEdge}, ps::Vector{<:AbstractEdges}, itrs::Vector) =
    datagrow!(f, f(), vs, ps, itrs)

# 2022-01-11: alternative meta
# @generated function datagrow!(f::Function, x::AbstractGraph,
#                               vs::Tuple{Vararg{T, N} where {T<:Union{AbstractLabel, AbstractLabels}}},
#                               ps::Tuple{Vararg{S, N} where {S<:AbstractEdges}}, itrs) where {N}
#     quote
#         Base.Cartesian.@nexprs $N i -> datagrow!(f, x, vs[i], ps[i], itrs[i])
#     end
# end
@generated function datagrow!(f::Function, x::AbstractGraph,
                              vs::Tuple{Vararg{Union{AbstractLabel, AbstractLabels}, N}},
                              ps::Tuple{Vararg{AbstractEdges, N}}, itrs) where {N}
    quote
        Base.Cartesian.@nexprs $N i -> datagrow!(f, x, vs[i], ps[i], itrs[i])
    end
end
# Convenience wrapper
# function datagrow(f::Function,
#                   vs::Tuple{Vararg{T, N} where {T<:Union{AbstractLabel, AbstractLabels}}},
#                   ps::Tuple{Vararg{S, N} where {S<:AbstractEdges}}, itrs) where {N}
#     datagrow!(f, f(), vs, ps, itrs)
# end
function datagrow(f::Function,
                  vs::Tuple{Vararg{Union{AbstractLabel, AbstractLabels}, N}},
                  ps::Tuple{Vararg{AbstractEdges, N}}, itrs) where {N}
    datagrow!(f, f(), vs, ps, itrs)
end


################################
# Growth from non-flat sources - p. 475, 2021-09-22
function datagrow!(f::Function, g::AbstractGraph, v::T, p::AbstractEdges, vitr, pitr) where {T<:Union{AbstractLabel, AbstractLabels}}
    x = Vector{Any}(undef, p.N)
    for (a, b) ∈ zip(vitr, pitr)
        p(x, b)
        get_datapush!(f, g, v(a), x...)
    end
    return g
end

# Possible parallel growth method
function tdatagrow!(f::Function, g::AbstractGraph, v::T, p::AbstractEdges,
                    itrsource::AbstractDict) where {T<:Union{AbstractLabel, AbstractLabels}}
    @sync for pₜ ∈ t
        Threads.@spawn datagrow!(f, pₜ.second, v, p, eachcol(itrsource[pₜ.first]))
        # Alternative:
        # let itr = eachcol(itrsource[p.first])
        #     Threads.@spawn datagrow!(f, p.second, v, p, itr)
        # end
    end
    return g
end

################
# 2021-11-10: map constructors for vectors of iterables
function mapgrow!(f::Function, gs::Vector{A}, p::AbstractEdges, itrs::Vector{T}) where {T} where {A<:AbstractGraph}
    @inbounds for i ∈ eachindex(gs, itrs)
        grow!(f, gs[i], p, itrs[i])
    end
    return gs
end
function mapgrow(f::Function, p::AbstractEdges, itrs::Vector{T}) where {T}
    gs = [f() for _ = 1:length(itrs)]
    mapgrow!(f, gs, p, itrs)
end

function tmapgrow!(f::Function, gs::Vector{A}, p::AbstractEdges, itrs::Vector{T}) where {T} where {A<:AbstractGraph}
    @sync @inbounds for i ∈ eachindex(gs, itrs)
        Threads.@spawn grow!(f, gs[i], p, itrs[i])
    end
    return gs
end

function tmapgrow(f::Function, p::AbstractEdges, itrs::Vector{T}) where {T}
    gs = [f() for _ = 1:length(itrs)]
    tmapgrow!(f, gs, p, itrs)
end

function tʳmapgrow!(f::Function, gs::Vector{A}, p::AbstractEdges, itrs::Vector{T}, M::Int=Threads.nthreads()) where {T} where {A<:AbstractGraph}
    ranges = equalranges(length(itrs), M)
    @sync for r ∈ ranges
        Threads.@spawn mapgrow!(f, gs[r], p, itrs[r])
    end
    return gs
end

function tʳmapgrow(f::Function, p::AbstractEdges, itrs::Vector{T}, M::Int=Threads.nthreads()) where {T}
    gs = [f() for _ = 1:length(itrs)]
    tʳmapgrow!(f, gs, p, itrs, M)
end

####
function mapdatagrow!(f::Function, gs::Vector{A}, v::U, p::AbstractEdges, itrs::Vector{T}) where {T} where {A<:AbstractGraph} where {U<:Union{AbstractLabel, AbstractLabels}}
    @inbounds for i ∈ eachindex(gs, itrs)
        datagrow!(f, gs[i], v, p, itrs[i])
    end
    return gs
end
function mapdatagrow(f::Function, v::U, p::AbstractEdges, itrs::Vector{T}) where {T} where {U<:Union{AbstractLabel, AbstractLabels}}
    gs = [f() for _ = 1:length(itrs)]
    mapdatagrow!(f, gs, v, p, itrs)
end

function tmapdatagrow!(f::Function, gs::Vector{A}, v::U, p::AbstractEdges, itrs::Vector{T}) where {T} where {A<:AbstractGraph} where {U<:Union{AbstractLabel, AbstractLabels}}
    @sync @inbounds for i ∈ eachindex(gs, itrs)
        Threads.@spawn datagrow!(f, gs[i], v, p, itrs[i])
    end
    return gs
end

function tmapdatagrow(f::Function, v::U, p::AbstractEdges, itrs::Vector{T}) where {T} where {U<:Union{AbstractLabel, AbstractLabels}}
    gs = [f() for _ = 1:length(itrs)]
    tmapdatagrow!(f, gs, v, p, itrs)
end

function tʳmapdatagrow!(f::Function, gs::Vector{A}, v::U, p::AbstractEdges, itrs::Vector{T}, M::Int=Threads.nthreads()) where {T} where {A<:AbstractGraph} where {U<:Union{AbstractLabel, AbstractLabels}}
    ranges = equalranges(length(itrs), M)
    @sync for r ∈ ranges
        Threads.@spawn mapdatagrow!(f, gs[r], v, p, itrs[r])
    end
    return gs
end

function tʳmapdatagrow(f::Function, v::U, p::AbstractEdges, itrs::Vector{T}, M::Int=Threads.nthreads()) where {T} where {U<:Union{AbstractLabel, AbstractLabels}}
    gs = [f() for _ = 1:length(itrs)]
    tʳmapdatagrow!(f, gs, v, p, itrs, M)
end


# function mapdatagrow!(f::Function, gs::Vector{A}, v::AbstractEdge, p::AbstractEdges, itrs::Vector{T}) where {T} where {A<:AbstractGraph}
#     for i ∈ eachindex(gs, itrs)
#         datagrow!(f, gs[i], v, p, itrs[i])
#     end
#     return gs
# end
# function mapdatagrow(f::Function, v::AbstractEdge, p::AbstractEdges, itrs::Vector{T}) where {T}
#     gs = Vector{typeof(f())}(undef, length(itrs))
#     # @inbounds for n ∈ eachindex(itrs)
#     #     gs[n] = datagrow(f, v, p, eachcol(itrs[n]))
#     # end
#     mapdatagrow!(f, gs, v, p, itrs)
#     return gs
# end

# function tmapdatagrow!(f::Function, gs::Vector{A}, v::AbstractEdge, p::AbstractEdges, itrs::Vector{T}, M::Int=Threads.nthreads()) where {T} where {A<:AbstractGraph}
#     N = length(itrs)
#     ranges = equalranges(N, M)
#     @inbounds @sync for m ∈ eachindex(ranges)
#         Threads.@spawn mapdatagrow!(f, gs[ranges[m]], v, p, itrs[ranges[m]])
#     end
#     return gs
# end

# function tmapdatagrow(f::Function, v::AbstractEdge, p::AbstractEdges, itrs::Vector{T},
#                       M::Int=Threads.nthreads()) where {T}
#     N = length(itrs)
#     ranges = equalranges(N, M)
#     gs = Vector{Vector{typeof(f())}}(undef, M)
#     # Threads.@threads for m ∈ 1:M #eachindex(ranges)
#     @inbounds @sync for m ∈ eachindex(ranges)
#         Threads.@spawn gs[m] = mapdatagrow(f, v, p, itrs[ranges[m]])
#     end
#     return vcat(gs...)
# end

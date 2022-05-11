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

function grow!(f::Function, t::AbstractNode, p::AbstractEdges)
    for x ∈ p
        get!(f, t, x...)
    end
    return t
end
function grow!(f::Function, t::AbstractNode, p::AbstractEdges, itr)
    x = Vector{Any}(undef, p.N)
    for item ∈ itr
        p(x, item)
        get!(f, t, x...)
    end
    return t
end

# Alternately, just remove the type on v

# function _valgrow!(f::Function, t::AbstractNode, v, p::AbstractEdges)
#     for x ∈ p
#         get_valpush!(f, t, v(x), x...)
#     end
#     return t
# end
# function _valgrow!(f::Function, t::AbstractNode, v, p::AbstractEdges, itr)
#     x = Vector{Any}(undef, p.N)
#     for item ∈ itr
#         p(x, item)
#         get_valpush!(f, t, v(item), x...)
#     end
#     return t
# end
# valgrow!(f::Function, t::AbstractNode, v::AbstractEdges, p::AbstractEdges) =
#     _valgrow!(f, t, v, p)
# valgrow!(f::Function, t::AbstractNode, v::AbstractEdge, p::AbstractEdges) =
#     _valgrow!(f, t, v, p)
# valgrow!(f::Function, t::AbstractNode, v::AbstractEdges, p::AbstractEdges, itr) =
#     _valgrow!(f, t, v, p, itr)
# valgrow!(f::Function, t::AbstractNode, v::AbstractEdge, p::AbstractEdges, itr) =
#     _valgrow!(f, t, v, p, itr)

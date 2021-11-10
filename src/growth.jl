#
# Date created: 2021-10-20
# Author: aradclif
#
#
############################################################################################
#### Methods of grow_(:field)! : see p. 443-446, 451-455, 2021-09-14/15
# grow! -> get! only has 2 signatures:
# (f::Function, t::AbstractNode, p::AbstractPathKeys)
# (f::Function, t::AbstractNode, p::AbstractPathKeys, itr)
# grow_val! -> get_pushval! only has 2 signatures:
# (f::Function, t::AbstractNode, v::AbstractPathKeys, p::AbstractPathKeys)
# (f::Function, t::AbstractNode, v::AbstractPathKeys, p::AbstractPathKeys, itr)
# grow_spec! -> get_pushspec! has 2 signatures, same as grow_val!
# grow_sval! -> get_pushsval! has 2 signatures, same as grow_val!

function grow!(f::Function, t::AbstractNode, p::AbstractPathKeys)
    for x ∈ p
        get!(f, t, x...)
    end
    return t
end
function grow!(f::Function, t::AbstractNode, p::AbstractPathKeys, itr)
    x = Vector{Any}(undef, p.N)
    for item ∈ itr
        p(x, item)
        get!(f, t, x...)
    end
    return t
end

# Alternately, just remove the type on v

# function _valgrow!(f::Function, t::AbstractNode, v, p::AbstractPathKeys)
#     for x ∈ p
#         get_valpush!(f, t, v(x), x...)
#     end
#     return t
# end
# function _valgrow!(f::Function, t::AbstractNode, v, p::AbstractPathKeys, itr)
#     x = Vector{Any}(undef, p.N)
#     for item ∈ itr
#         p(x, item)
#         get_valpush!(f, t, v(item), x...)
#     end
#     return t
# end
# valgrow!(f::Function, t::AbstractNode, v::AbstractPathKeys, p::AbstractPathKeys) =
#     _valgrow!(f, t, v, p)
# valgrow!(f::Function, t::AbstractNode, v::AbstractPathKey, p::AbstractPathKeys) =
#     _valgrow!(f, t, v, p)
# valgrow!(f::Function, t::AbstractNode, v::AbstractPathKeys, p::AbstractPathKeys, itr) =
#     _valgrow!(f, t, v, p, itr)
# valgrow!(f::Function, t::AbstractNode, v::AbstractPathKey, p::AbstractPathKeys, itr) =
#     _valgrow!(f, t, v, p, itr)

#
# Date created: 2021-09-08
# Author: aradclif
#
#
############################################################################################
#### 2021-09-08: Experiments with Types, Interfaces, iterators for Node
# - Experiment with Interfaces
################################################################
#### Experiment with Interfaces -- see p. 404-405, 2021-09-08
abstract type AbstractNode{T, S} end
struct NodeWrap{T<:AbstractDict, S<:AbstractVector} <: AbstractNode{T, S}
    link::T
    val::S
end
NodeWrap(link::T) where {T<:AbstractDict} = NodeWrap(link, Any[])
NodeWrap() = NodeWrap(Dict{Any,Any}())
Base.iterate(x::A where {A<:NodeWrap{T,S}}) where {T,S} = iterate(x.link)
Base.iterate(x::A where {A<:NodeWrap{T,S}}, i) where {T,S} = iterate(x.link, i)
Base.eltype(x::A where {A<:NodeWrap{T,S}}) where {T,S} = Pair{Any, NodeWrap{T,S}}
Base.length(x::A where {A<:NodeWrap{T, S}}) where {T,S} = length(x.link)
Base.isequal(x::A where {A<:NodeWrap{T, S}}, y::B where {B<:NodeWrap{T, S}}) where {T,S} =
    x.link == y.link && x.val == y.val
Base.:(==)(x::A where {A<:NodeWrap{T, S}}, y::B where {B<:NodeWrap{T, S}}) where {T,S} =
    Base.isequal(x, y)
Base.keys(x::A where {A<:NodeWrap{T, S}}) where {T,S} = keys(x.link)
Base.values(x::A where {A<:NodeWrap{T, S}}) where {T,S} = values(x.link)
Base.pairs(x::A where {A<:NodeWrap{T, S}}) where {T,S} = pairs(x.link)
# Experimental ports
Base.get(x::A where {A<:NodeWrap{T, S}}, key, default) where {T,S} =
    get(x.link, key, default)
# get2(x::A where {A<:NodeWrap{T, S}}, key, default) where {T,S} =
#     get(x.link, key, default)::NodeWrap{T,S}
Base.get!(x::A where {A<:NodeWrap{T, S}}, key, default) where {T,S} =
    get!(x.link, key, default)
Base.get!(f::Function, x::A where {A<:NodeWrap{T, S}}, key) where {T,S} =
    get!(f, x.link, key)
#
a = NodeWrap()
b = NodeWrap(Dict{String,Int}(), Float64[])
c = NodeWrap(Dict{String,Int}("$(2i)" => i for i=1:7), Float64[])
# test of simple iterator -- works.
next = iterate(c)
while next !== nothing
    (item, state) = next
    println("Item is: ", item, " State is: ", state)
    next = iterate(c, state)
end
# test of container iterator -- works
for p in c
    println("K is: ", p.first, " V is: ", p.second)
end
# test of stateful iterator -- cannot access state during iteration (likely for safety)
# for (k, v, s) in c
#     println("Item is: ", k, " State is: ", s)
# end
####
d = NodeWrap(Dict{Any,Any}("$i" => NodeWrap() for i = 1:3))
eltype(d)
dp = ("1" => NodeWrap())
dp in d
d1 = first(iterate(d))
@code_warntype dp == dp
for p in d
    println("K is: ", p.first, " V is: ", p.second)
    println(typeof(p))
end
fcol(d) = [p.second for p in d]
@code_warntype fcol(d)
# allocation test -- same time and allocation for iteration
function checklength(x)
    v::Int = 0
    for p in x
        v += length(p.first)
    end
    return v
end
@timev checklength(d)
@code_warntype checklength(d)
dd = d.link; @timev checklength(dd)
# Type-stability during iteration
function iterprint(iter)
    for p in iter
        println("K is: ", p.first, " V is: ", p.second)
    end
    nothing
end
@code_warntype iterprint(d)
# test of enumerate
for (i, x) in enumerate(d)
    println("index is: ", i, " x is: ", x)
end
################ Experimental test of buildnode!
function buildnode!(x::A where {A<:NodeWrap{T, S}}, ks, v) where {T,S}
    # local node::NodeWrap{T, S}
    # node = x
    node::NodeWrap{T,S} = x
    for k in ks
        node = get!(node, k, NodeWrap())
    end
    push!(node.val, v)
    # nothing
    return x
end
dd = deepcopy(d)
ks = ["1", "2", "3"]
v = [4, 5, 6]
@code_warntype buildnode!(dd, ks, v)
@benchmark buildnode!(dd, ks, v)
l1 = get(dd, "1", "stuff")
l2 = get(l1, "2", "stuff")
l3 = get(l2, "3", "stuff")
#### As a recursion -- maybe slightly worse performance, but could just be machine noise.
# Certainly, it involves fewer lines of LLVM code; perhaps the 5 fewer lines are
# worthwhile, if it makes the difference between inlining.
function rbuildnode!(x::A where {A<:NodeWrap{T, S}}, ks, v, N::Int, C::Int) where {T,S}
    local node::NodeWrap{T,S}
    if C < N
        node = get!(x, ks[C], NodeWrap())
        rbuildnode!(node, ks, v, N, C + 1)
    else # C == N
        node = get!(x, ks[C], NodeWrap())
        push!(node.val, v)
    end
    # nothing
    return x
end
rbuildnode!(x::A where {A<:NodeWrap{T, S}}, ks, v) where {T,S} = rbuildnode!(x, ks, v, length(ks), 1)
ddd = deepcopy(d)
@code_warntype rbuildnode!(ddd, ks, v, 3, 1)
@benchmark rbuildnode!(ddd, ks, v, 3, 1)
@code_warntype rbuildnode!(ddd, ks, v)
@benchmark rbuildnode!(ddd, ks, v)
#### An even more basic recursion -- worse performance of the 3
# function rbuildnode2!(x::A where {A<:NodeWrap{T, S}}, ks, v, N::Int, C::Int) where {T,S}
#     local node::NodeWrap{T,S}
#     if C < N
#         rbuildnode2!(get!(x, ks[C], NodeWrap()), ks, v, N, C + 1)
#     else # C == N
#         push!(get!(x, ks[C], NodeWrap()).val, v)
#     end
#     nothing
# end
# rbuildnode2!(x::A where {A<:NodeWrap{T, S}}, ks, v) where {T,S} =
#     rbuildnode2!(x, ks, v, length(ks), 1)
# ddd = deepcopy(d)
# @code_warntype rbuildnode2!(ddd, ks, v, 3, 1)
# @benchmark rbuildnode2!(ddd, ks, v, 3, 1)
# @code_warntype rbuildnode2!(ddd, ks, v)
# @benchmark rbuildnode2!(ddd, ks, v)
################################################################
#### Further experimentation: abstract types
abstract type AbstractNode{T, S} end
struct NodeWrap{T<:AbstractDict, S<:AbstractVector} <: AbstractNode{T, S}
    link::T
    val::S
end
NodeWrap(link::T) where {T<:AbstractDict} = NodeWrap(link, Any[])
NodeWrap() = NodeWrap(Dict{Any,Any}())
struct NodeWrap2{T<:AbstractDict, S<:AbstractVector, U<:AbstractDict} <: AbstractNode{T,S}
    link::T
    val::S
    spec::U
end
NodeWrap2(link::T, val::S) where {T<:AbstractDict, S<:AbstractVector} =
    NodeWrap2(link, val, Dict{String,String}())
NodeWrap2(link::T) where {T<:AbstractDict} = NodeWrap2(link, Any[])
NodeWrap2() = NodeWrap2(Dict{Any,Any}())
NodeWrap2(x::NodeWrap{T, S}) where {T,S} = NodeWrap2(x.link, x.val) # constructor
struct NotNode{T<:AbstractDict}
    link::T
end
NotNode() = NotNode(Dict{Int,Int}())
####
eltype2(x::A where {A<:AbstractNode{T,S}}) where {T,S} = println(T, S)
eltype3(x::A where {A<:AbstractNode{T,S}}) where {T,S} = println(A{T,S}) # Does not work
eltype4(x::A) where {A<:AbstractNode{T,S}} where {T,S} = A
eltype8(::A) where {A<:AbstractNode} = A # works on Any subtype of AbstractNode
eltype9(::A) where {A<:AbstractNode{T,S}} where {T,S} = (T, S)
eltype10(::A) where {A<:AbstractNode{T,S}} where {T,S} = A
####
a = NodeWrap()
a2 = NodeWrap2()
an = NotNode()
#
eltype4(a), eltype8(a), eltype9(a)
eltype4(a2), eltype8(a2), eltype9(a2)
####
Base.promote_rule(::Type{NodeWrap2}, ::Type{NodeWrap}) = NodeWrap2
Base.promote_rule(::Type{NodeWrap2{A,B,C}}, ::Type{NodeWrap{T,S}}) where {A,B,C,T,S} =
    NodeWrap2{promote_type(A, T), promote_type(B, S), C}
promote_type(NodeWrap, NodeWrap2)
Base.convert(::Type{NodeWrap2{A,B,C}}, x::NodeWrap{T,S}) where {A,B,C, T,S} =
    NodeWrap2(x.link, x.val)
Base.convert(::Type{NodeWrap2}, x::NodeWrap{T,S}) where {T,S} = NodeWrap2(x)

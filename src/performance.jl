#
# Date created: 2021-10-25
# Author: aradclif
#
#
############################################################################################
#### Separate performance tests
# default iterables
mat = reshape([1:200000;], (4, 50000));
pmat = reshape([1:200000;], (20, 10000));
qmat3 = reshape([1:200000;], (50, 4000));
qmat4 = reshape([1:200000;], (100, 2000));
qmat5 = reshape([1:200000;], (200, 1000));
qmat6 = reshape([1:200000;], (10, 20000));

# default Edges
pni = Edges([IndexedEdge(i) for i = 1:4]);
pm = Edges([[IndexedEdge(i) for i = 1:3]; IndexedEdge(string, 4)]);
pni2 = Edges([IndexedEdge(i) for i = 1:20]);
pm2 = Edges([[IndexedEdge(i) for i = 1:10]; IndexedEdge(string, 11); [IndexedEdge(i) for i = 12:20]]);
pni3 = Edges([IndexedEdge(i) for i = 1:50]);
pm3 = Edges([[IndexedEdge(i) for i = 1:10]; IndexedEdge(string, 11); [IndexedEdge(i) for i = 12:50]]);
pni4 = Edges([IndexedEdge(i) for i = 1:100]);
pm4 = Edges([[IndexedEdge(i) for i = 1:10]; IndexedEdge(string, 11); [IndexedEdge(i) for i = 12:100]]);
pni5 = Edges([IndexedEdge(i) for i = 1:200]);
pm5 = Edges([[IndexedEdge(i) for i = 1:10]; IndexedEdge(string, 11); [IndexedEdge(i) for i = 12:200]]);
pni6 = Edges([IndexedEdge(i) for i = 1:10]);
pm6 = Edges([[IndexedEdge(i) for i = 1:9]; IndexedEdge(string, 10)]);

# gf() = SimpleNode()
sdg() = SimpleDiGraph()
sg() = SimpleGraph()
#### SimpleNode vs. AryNode
# Result: growth. SimpleNode seems to be moderately faster
@benchmark grow!(gf, SimpleNode(), pni3, eachcol(qmat3)) seconds=30
@benchmark grow!(gf, SimpleNode(), pm3, eachcol(qmat3)) seconds=30
@benchmark grow!(af, AryNode(), pni3, eachcol(qmat3)) seconds=30
@benchmark grow!(af, AryNode(), pm3, eachcol(qmat3)) seconds=30
@benchmark grow!(sdg, SimpleDiGraph(), pni3, eachcol(qmat3)) seconds=30
@benchmark grow!(sdg, SimpleDiGraph(), pm3, eachcol(qmat3)) seconds=30
@benchmark grow(sdg, pni3, eachcol(qmat3)) seconds=30

# t0 = grow!(gf, SimpleNode(), pni, eachcol(mat));
# a0 = grow!(af, AryNode(), pni, eachcol(mat));
g0 = grow!(sdg, SimpleDiGraph(), pni, eachcol(mat));
# t3 = grow!(gf, SimpleNode(), pni3, eachcol(qmat3));
# a3 = grow!(af, AryNode(), pni3, eachcol(qmat3));
g3 = grow!(sdg, SimpleDiGraph(), pni3, eachcol(qmat3));
g4 = grow!(sdg, SimpleDiGraph(), pni3, eachcol(qmat3));
# t6 = grow!(gf, SimpleNode(), pni6, eachcol(qmat6));
# a6 = grow!(af, AryNode(), pni6, eachcol(qmat6));
g6 = grow!(sdg, SimpleDiGraph(), pni6, eachcol(qmat6));

# Result: traversal tests.
# 4x50000: AryNode is faster overall
# 50x4000: Other than getindex, SimpleNode is faster
# 10x20000: Other than getindex, SimpleNode is faster
# @benchmark maxdepth(t3)
@benchmark depth(g3)
# @benchmark maxbreadth(t3)
@benchmark maxbreadth(g3)
# @benchmark size(t3)
@benchmark size(g3)
# @benchmark size(t3, 10)
@benchmark size(g3, 10)
# @benchmark size(t3, 50)
@benchmark size(g3, 50)
# @benchmark rlength(t3)
@benchmark rlength(g3)
z = [1:10;];
# @benchmark t3[z...]
@benchmark g3[z...]

# Result: map, reduce, etc.
@benchmark countall(isempty, t3)
@benchmark countall(isempty, g3)
fl!(dest, x) = (dest[1][1] += length(x); dest)
fl_v2!(dest, x) = (dest[1][1] = dest[1][1] + length(x))
fl2!(dest, x) = (dest[1][1] += length(x.val[1]); dest)
fl3!(dest, x) = (dest[1][1] += length(unique(y -> y.val[1], values(x))))
fl4!(dest, x) = (dest[1][1] = max(dest[1][1], length(x)))
fl3_v2!(dest, x) = (dest[1][1] += length(unique(p -> p.second.val[1], x)))
fl3_v3!(dest, x) = (dest[1][1] += length(unique!([p.second.val[1] for p in x])))
fl3_v3g!(dest, x) = (dest[1][1] += length(unique!([p.second.data[1] for p in x])))
@benchmark mapat(fl!, ((5,),), t3, 10)
@benchmark mapat(fl!, ((5,),), a3, 10)
@benchmark mapat(fl!, ((5,),), g3, 10)
@benchmark mapupto(fl!, ((5,),), g3, 40)

################
prefix = "/nfs/site/home/aradclif/my_cit_scratch1/diagnosis/1276/x76se/gt_g1m/2021-09-15";
pni7 = Edges([IndexedEdge(i) for i = 1:10]);
suffix8 = "D117E4C0_36.arrow" # 218MiB
file8 = joinpath(prefix, suffix8);
dmat8 = permutedims(hcat([convert(Vector{Any}, copy(col)) for col in Arrow.Table(file8)]...));
####
#### Comparison given optimized Arrow files
cellinst(x) = @inbounds (x[1], pathprefix(x[2]))
prefix10 = "/nfs/site/home/aradclif/my_cit_scratch1/diagnosis/1276/x76se/gt_g1m/2021-09-24";
suffix10 = "D117E4C0_36.arrow";
file10 = joinpath(prefix10, suffix10);
dmat10 = permutedims(hcat([convert(Vector{Any}, copy(col)) for col in Arrow.Table(file10)]...));
pni10 = Edges([IndexedEdge(x -> Tuple(x), [1, 2]); [IndexedEdge(i) for i = 3:8]])
vi10 = IndexedEdge(cellinst, [10, 9])
pni10_4 = Edges4((IndexedEdge(x -> Tuple(x), [1, 2]), (IndexedEdge(i) for i = 3:8)...))
vi10_4 = IndexedEdge(cellinst, [10, 9])
####

# Result: growth on heterogeneous keys. SimpleNode is faster
@benchmark grow!(gf, SimpleNode(), pni7, eachcol(dmat8)) seconds=240
@benchmark grow!(af, AryNode(), pni7, eachcol(dmat8)) seconds=240
@benchmark grow!(sdg, SimpleDiGraph(), pni7, eachcol(dmat8)) seconds=240
@benchmark valgrow!(gf, SimpleNode(), vi10, pni10, eachcol(dmat10)) seconds=240
@benchmark valgrow!(af, AryNode(), vi10, pni10, eachcol(dmat10)) seconds=240
@benchmark datagrow!(sdg, SimpleDiGraph(), vi10, pni10, eachcol(dmat10)) seconds=240

t7 = grow!(gf, SimpleNode(), pni7, eachcol(dmat8));
a7 = grow!(af, AryNode(), pni7, eachcol(dmat8));
g7 = grow!(sdg, SimpleDiGraph(), pni7, eachcol(dmat8));
t10 = valgrow!(gf, SimpleNode(), vi10, pni10, eachcol(dmat10));
a10 = valgrow!(af, AryNode(), vi10, pni10, eachcol(dmat10));
@timev datagrow!(sdg, SimpleDiGraph(), vi10, pni10, eachcol(dmat10));
@timev datagrow!(sdg, SimpleDiGraph(), vi10_4, pni10_4, eachcol(dmat10));
@timev datagrow(sdg, vi10, pni10, eachcol(dmat10));
@timev datagrow(sdg, vi10_4, pni10_4, eachcol(dmat10));

# Result: traversal tests
# t7: 4387231 nodes. AryNode is moderately faster
# t10: 1673512 nodes. AryNode is faster: (290 - 170) / 290, generally by ≈40%
# @benchmark maxdepth(t10)
@benchmark depth(g10)
# @benchmark maxbreadth(t10)
@benchmark maxbreadth(g10)
# @benchmark size(t10)
@benchmark size(g10)
# @benchmark size(t10, 8)
@benchmark size(g10, 8)
# @benchmark rlength(t10)
@benchmark rlength(g10)
# z10 = pni10(first(eachcol(dmat8)));
z10 = pni10(first(eachcol(dmat10)));
# @benchmark t10[z10...]
@benchmark g10[z10...]

# Result: map, reduce, etc.
# @benchmark countall(isempty, t10)
@benchmark countall(isempty, g10)
# @benchmark mapat(fl!, ((5,),), t7, 10)
@benchmark mapat(fl!, ((5,),), g7, 10)
# @benchmark mapupto(fl!, ((5,),), t7, 10)
@benchmark mapupto(fl!, ((5,),), g7, 10)

fs = [x -> true, x -> true, x -> true, x -> x[1] ∈ (1,2,3)]
# @benchmark mapat(fl3_v3!, ((5,),), t10, 7)
@benchmark mapat(fl3_v3g!, ((5,),), g10, 7)
# @benchmark mapfilterat(fl!, fs, ((5,),), t10, 5)
@benchmark mapfilterat(fl!, fs, ((5,),), g10, 5)
####
n7 = sum(size(t7))
n10 = sum(size(t10))
time7 = [753, 508]
time10 = [123, 57]
n7 ./ time7
n10 ./ time10
mem7 = [48.4 * 2^20, 160]
mem10 = [6.99 * 2^20, 160]
mpn7 = mem7 ./ n7
mpn10 = mem10 ./ n10
bpernode_max = maximum(maximum, [mpn7, mpn10])
bpernode_min = minimum(minimum, [mpn7, mpn10])
filemem = 218 * 2^20
npm = n10 / filemem
Nₘ = 26 * 2^30
Nₜ = Nₘ * npm
Nₜ * bpernode_max
Nₜ * bpernode_min
4000 * 160
# eval of unique function on suspects: 875ms, 167.85MiB
time = 875
Tₘ = 167.85 * 2^20
su = 1215165
sps = su / .875
mps = Tₘ / su
filemps = filemem / su
Tₛᵤ = Nₘ / filemps
allocunique = Tₛᵤ * mps
timeunique = Tₛᵤ / sps
############################################################################################
#### 2021-11-10: checking get, get!
retnot() = nothing
sdg() = SimpleDiGraph()
a1 = sdg()
get!(sdg, a1, 1,2,3,4,5)
@benchmark get(retnot, a1, 1,2,3,4,5)
sg() = SimpleGraph()
a2 = sg()
bget!(sg, a2, 1,2,3,4,5)
# @benchmark get(retnot, a2, 1,2,3,4,5)
@benchmark bget(retnot, a2, 1,2,3,4,5)

############################################################################################
#### 2022-01-11: Add performance tests of meta-unrolling
p4 = Edges(Tuple(IndexedEdge(i) for i = 1:9))
v4 = IndexedEdge(10)

# @timev datagrow(sdg, v4, p4, eachcol(m_dg));


p1 = Edges([IndexedEdge(i) for i = 1:9])
v1 = IndexedEdge(10)
@timev datagrow(sdg, v1, p1, eachcol(m_dg));
@timev datagrow3(sdg, v1, p1, eachcol(m_dg));

L = 100
rmat = rand(Int, L, 2000);

p4_l = Edges(Tuple(IndexedEdge(i) for i = 1:L-1));
v4_l = IndexedEdge(L);

@timev datagrow(sdg, v4_l, p4_l, eachcol(rmat));
@code_warntype datagrow!(sdg, sdg(), v4_l, p4_l, eachcol(rmat), Val(100))

p1_l = Edges([IndexedEdge(i) for i = 1:L-1]);
v1_l = IndexedEdge(L);
@timev datagrow(sdg, v1_l, p1_l, eachcol(rmat));
@timev datagrow3(sdg, v1_l, p1_l, eachcol(rmat));

L = 10
rmat = rand(Int, L, 2000);
p4_l = Edges(Tuple(IndexedEdge(i) for i = 1:L-1));
v4_l = IndexedEdge(L);

@benchmark datagrow(sdg, v4_l, p4_l, eachcol(rmat))

p2_l = Edges2(Tuple(IndexedEdge(i) for i = 1:L-1));

@benchmark datagrow2(sdg, v4_l, p2_l, eachcol(rmat))

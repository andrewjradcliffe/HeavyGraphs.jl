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

# default PathKeys
pni = PathKeys([IndexedPathKey(i) for i = 1:4]);
pm = PathKeys([[IndexedPathKey(i) for i = 1:3]; IndexedPathKey(string, 4)]);
pni2 = PathKeys([IndexedPathKey(i) for i = 1:20]);
pm2 = PathKeys([[IndexedPathKey(i) for i = 1:10]; IndexedPathKey(string, 11); [IndexedPathKey(i) for i = 12:20]]);
pni3 = PathKeys([IndexedPathKey(i) for i = 1:50]);
pm3 = PathKeys([[IndexedPathKey(i) for i = 1:10]; IndexedPathKey(string, 11); [IndexedPathKey(i) for i = 12:50]]);
pni4 = PathKeys([IndexedPathKey(i) for i = 1:100]);
pm4 = PathKeys([[IndexedPathKey(i) for i = 1:10]; IndexedPathKey(string, 11); [IndexedPathKey(i) for i = 12:100]]);
pni5 = PathKeys([IndexedPathKey(i) for i = 1:200]);
pm5 = PathKeys([[IndexedPathKey(i) for i = 1:10]; IndexedPathKey(string, 11); [IndexedPathKey(i) for i = 12:200]]);
pni6 = PathKeys([IndexedPathKey(i) for i = 1:10]);
pm6 = PathKeys([[IndexedPathKey(i) for i = 1:9]; IndexedPathKey(string, 10)]);

#### SimpleNode vs. AryNode
# Result: growth. SimpleNode seems to be moderately faster
@benchmark grow!(gf, SimpleNode(), pni3, eachcol(qmat3)) seconds=30
@benchmark grow!(gf, SimpleNode(), pm3, eachcol(qmat3)) seconds=30
@benchmark grow!(af, AryNode(), pni3, eachcol(qmat3)) seconds=30
@benchmark grow!(af, AryNode(), pm3, eachcol(qmat3)) seconds=30

t0 = grow!(gf, SimpleNode(), pni, eachcol(mat));
a0 = grow!(af, AryNode(), pni, eachcol(mat));
t3 = grow!(gf, SimpleNode(), pni3, eachcol(qmat3));
a3 = grow!(af, AryNode(), pni3, eachcol(qmat3));
t6 = grow!(gf, SimpleNode(), pni6, eachcol(qmat6));
a6 = grow!(af, AryNode(), pni6, eachcol(qmat6));

# Result: traversal tests.
# 4x50000: AryNode is faster overall
# 50x4000: Other than getindex, SimpleNode is faster
# 10x20000: Other than getindex, SimpleNode is faster
@benchmark maxdepth(t3)
@benchmark maxdepth(a3)
@benchmark maxbreadth(t3)
@benchmark maxbreadth(a3)
@benchmark size(t3)
@benchmark size(a3)
@benchmark size(t3, 10)
@benchmark size(a3, 10)
@benchmark size(t3, 50)
@benchmark size(a3, 50)
z = [1:10;];
@benchmark t3[z...]
@benchmark a3[z...]

# Result: map, reduce, etc.
@benchmark countall(isempty, t3)
@benchmark countall(isempty, a3)

# Result: growth on heterogeneous keys. SimpleNode is faster
@benchmark grow!(gf, SimpleNode(), pni7, eachcol(dmat8)) seconds=240
@benchmark grow!(af, AryNode(), pni7, eachcol(dmat8)) seconds=240
@benchmark valgrow!(gf, SimpleNode(), vi10, pni10, eachcol(dmat10)) seconds=240
@benchmark valgrow!(af, AryNode(), vi10, pni10, eachcol(dmat10)) seconds=240

t7 = grow!(gf, SimpleNode(), pni7, eachcol(dmat8));
a7 = grow!(af, AryNode(), pni7, eachcol(dmat8));
t10 = valgrow!(gf, SimpleNode(), vi10, pni10, eachcol(dmat10));
a10 = valgrow!(af, AryNode(), vi10, pni10, eachcol(dmat10));

# Result: traversal tests
# t7: 4300000 nodes. AryNode is moderately faster
@benchmark maxdepth(t7)
@benchmark maxdepth(a7)
@benchmark maxbreadth(t7)
@benchmark maxbreadth(a7)
@benchmark size(t7)
@benchmark size(a7)
@benchmark size(t7, 10)
@benchmark size(a7, 10)
z7 = pni7(first(eachcol(dmat8)));
z10 = pni10(first(eachcol(dmat10)));
@benchmark t7[z...]
@benchmark a7[z...]

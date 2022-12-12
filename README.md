# HeavyGraphs

## Installation
```julia
using Pkg
Pkg.add("HeavyGraphs")    # the appropriate registry must be added for this to work
```

## Description
Associative array-backed (hence, "heavy") graph types. Provides
parameterizable constructors, growth, traversals; map, mapreduce:
indexed, non-indexed, and all the derived higher-order functions (any,
all, count, find, filter, etc.). The interface supports iteration,
Cartesian indexing, set operations, etc. -- everything you'd expect of
a fully-featured (recursive algebraic) data type. Notably, the
higher-order functions make use of the idea of depth relative to
entry point into the graph in order to provide consistent notions of
what something like "map" means at any given point during a graph
traversal.

This package provides the `SimpleDiGraph` and `SimpleGraph` types,
which correspond to direct and undirected graphs. Cycles are
inherently supported,, but whether a graph exhibits cycles is left to
the user -- presence or lack thereof is not checked for or enforced by
this package. Thus, if one wishes to construct directed acyclic
graphs, one would utilize the `SimpleDiGraph` type and simply not
introduce cycles. It is worth mentioning that many of the higher-order
functions derived from map and mapreduce are likely to infinite loop
if cycles are present, unless said cycles are handled by the caller.

## Usage

Given a flat representation of a graph, the connectivity of which is
specified, for example, using a walk.
```julia
julia> x,y = [1,2,3], [1,2,4];

julia> g = grow(() -> SimpleDiGraph{Int}(), Edges(i for i = 1:3), (x, y))
SimpleDiGraph{Int64}(Dict{Int64, SimpleDiGraph{Int64}}(1 => SimpleDiGraph{Int64}(Dict{Int64, SimpleDiGraph{Int64}}(2 => SimpleDiGraph{Int64}(Dict{Int64, SimpleDiGraph{Int64}}(4 => SimpleDiGraph{Int64}(Dict{Int64, SimpleDiGraph{Int64}}(), Any[]), 3 => SimpleDiGraph{Int64}(Dict{Int64, SimpleDiGraph{Int64}}(), Any[])), Any[])), Any[])), Any[])
```

Demonstrates the use of novel design patterns which utilize multiple
dispatch, metaprogramming, recursion, Varargs functions and functional
constructors to provide an interface which permits the use of
Cartesian indexing for access, assignment, adjacency, growth, pruning,
etc. at arbitrary depth (depth relative to the node one happens to be
on, or, equivalently: depth relative to the entry point of the
caller).

This incorporates functionality from `AbstractArray` and
`AbstractDict`, enabling flexible creation of graphs with arbitrary
labels on edges (and metadata on nodes, if desired).

See the source code functions `get!`, `getindex` to see the pattern
implemented in a manner which optimizes over the Varargs types,
leading to a recursion which dispatches to a fast, type-stable loop at
the first occurrence of an NTuple, otherwise reverting to a recursion
which consumes two elements of the Varargs tuple before repeating the
cycle.

# HeavyGraphs
Associative array-backed (hence, "heavy") graphs. Demonstrates the use of
novel design patterns which utilize multiple dispatch, recursion, Varargs
functions and functional constructors to provide an interface which
permits the use of Cartesian indexing for access, assignment, adjacency,
growth, pruning, etc. at arbitrary depth (depth relative to the node one
happens to be on).
See the source code functions `get!`, `getindex` to see the pattern implemented
in a manner which optimizes over the Varargs types, leading to a
recursion which dispatches to a fast, type-stable loop at the first occurrence
of an NTuple, otherwise reverting to a recursion which consumes two elements
of the Varargs tuple before repeating the cycle.

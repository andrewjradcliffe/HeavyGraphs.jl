module HeavyGraphs

using EqualRanges
# using SetDiffCard

# abstractedges.jl
export AbstractEdge, AbstractIndexedEdge, IndexedEdge,
    AbstractEdges, Edges

# abstractlabels.jl
export AbstractLabel, Label,
    AbstractIndexedLabel, IndexedLabel,
    AbstractLabels, Labels

# abstractgraph.jl
export AbstractGraph, AbstractSimpleDiGraph, AbstractSimpleGraph,
    SimpleDiGraph, SimpleGraph, sdg, sg, usdg, usg

export haspath, depth, maxbreadth, findmaxbreadth, rlength, rlength2, maxsize

export isbidirectional, bget!, bget, hasbipath

export get_datapush!

# growth.jl
export grow!, grow, datagrow!, datagrow, tdatagrow!

export mapgrow!, mapgrow, tmapgrow!, tmapgrow, tʳmapgrow!, tʳmapgrow,
    mapdatagrow!, mapdatagrow, tmapdatagrow!, tmapdatagrow, tʳmapdatagrow!, tʳmapdatagrow

# ndops.jl
export ndadd!, ndadd1!, ndpush!

# map.jl
export mapat, mapat!, tmapat,
    mapfilterat, mapfilterat!, tmapfilterat,
    mapat_pairs, mapat_pairs!, tmapat_pairs,
    mapfilterat_pairs, mapfilterat_pairs!, tmapfilterat_pairs,
    mapupto, mapupto!, tmapupto,
    mapfilterupto, mapfilterupto!, tmapfilterupto

# kmap.jl
export kmapat, kmapat!, tkmapat,
    kmapfilterat, kmapfilterat!, tkmapfilterat,
    kmapupto, kmapupto!, tkmapupto,
    kmapfilterupto, kmapfilterupto!, tkmapfilterupto

# iterators.jl
export foreachat, foreachfilterat, foreach_depthfirst, foreach_breadthfirst,
    foreachfrom, foreachthrough, foreachupto,
    countat, countall, countfrom, countthrough, countupto,
    findpathsat, findpathsat!, findpathsall, findpathsall!,
    findpathsfrom, findpathsfrom!, findpathsthrough, findpathsthrough!, findpathsupto, findpathsupto!,
    allat, allall, allfrom, allthrough,
    anyat, anyany, anyfrom, anythrough,
    mapall, mapall!, mapfrom, mapfrom!, mapthrough, mapthrough!

# multipliers.jl
export multabsentupto, multabsent

# countabsent.jl
export countabsent!, countstatus!, kcountstatus!, kcountabsent!,
    nextnonunit, dimsmultiplier

include("abstractedges.jl")
include("abstractlabels.jl")
include("abstractgraph.jl")
include("growth.jl")

include("ndops.jl")

include("map.jl")
include("kmap.jl")
include("iterators.jl")

include("multipliers.jl")
include("countabsent.jl")

end

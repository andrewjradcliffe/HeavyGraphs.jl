module HeavyGraphs

using EqualRanges
using SetDiffCard

# abstractpathkeys.jl
export AbstractPathKey, AbstractIndexedPathKey, IndexedPathKey,
    AbstractPathKeys, PathKeys

# abstractgraph.jl
export AbstractGraph, AbstractSimpleDiGraph, AbstractSimpleGraph,
    SimpleDiGraph, SimpleGraph

export haspath, depth, maxbreadth, findmaxbreadth, rlength, rlength2, maxsize

export isbidirectional, bget!, bget, hasbipath

export get_datapush!

export grow!, grow, datagrow!, datagrow, tdatagrow!

# map.jl
export mapat, mapat!, tmapat,
    mapfilterat, mapfilterat!, tmapfilterat,
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
    findpathsfrom, findpathsfrom!, findpathsthrough, findpathsthrough!, findpathsupto, findpathsupto!

# multipliers.jl
export multabsentupto, multabsent

# countabsent.jl
export countabsent!, countstatus!, kcountstatus!, kcountabsent!,
    nextnonunit, dimsmultiplier, ndadd!

include("abstractpathkeys.jl")
include("abstractgraph.jl")

include("map.jl")
include("kmap.jl")
include("iterators.jl")

include("multipliers.jl")
include("countabsent.jl")

end

#! /usr/bin/env Rscript
'Generate a scale-free, undirected, bipartite graph
Usage:
    bigraph_r.r --beta=<float> --nodes=<int> --degrees=<int> [--output=<./bipiartite_sf.csv>] [--edges=<int>] [--seed=<int>]
    
Options:
    -h --help  Show this screen.
    --beta=<float>  fitness of node, gamma = 1 / beta + 1, gamma is scale-free exponent 
    --nodes=<int>  2*n is the number of nodes in graph. Restrictions: number of nodes in group 1 = number of nodes in group 2. 
    --degrees=<int>  average degree of nodes 
    --output=<str> output filename with path [default: ./bipartite_sf.csv]
    --edges=<int>  number of (directed) edges
    --seed=<int>  set a seed for random graph generation 
' -> doc

withCallingHandlers({
    suppressMessages({
        library(StabEco, quietly = T)
        library(igraph, quietly = T)
        library(docopt, quietly = T)
})})


opts = docopt(doc)
for (opt in c('nodes', 'edges', 'degrees', 'beta', 'seed')){
    val = opts[[opt]]
    if (!is.null(val)){
        opts[[opt]] = as.numeric(val)
        }
}

bg_sf<-function(n1, beta, k, fn, m, seed){
    set.seed(seed)
    B = BiGraph$new(n1=n1, beta=beta, k=k, m=m, 
                    type = 'bipartite_sf', directed = F, is_adj = T)
    G = B$get_graph()
    write.csv(G, file = fn)
}

bg_sf(n1=opts[['nodes']], beta=opts[['beta']], k = opts[['degrees']], fn = opts[['output']], m = opts[['edges']], 
     seed = opts[['seed']])

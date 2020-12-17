
suppressPackageStartupMessages({
    library('SymSim')
    library(StabEco)
    library(igraph)
})

BiGraph()

bg_sf<-function(...){
    B = BiGraph$new(type = 'bipartite_sf', ...)
    G = B$get_graph()
    return(G)
}




beta = 2, n1 = 100, k = 90

B$





ngenes<-500
ncells<-300

# sample from a vector of gene lengths
gene_len <- sample(SymSim::gene_len_pool, ngenes, replace = FALSE)

# true transcript counts (intrinsic and extrinsic variation)
phyla <- Phyla5() # create tree structure for 5 separate populations (can create your own, see vignett)
tcr <- SimulateTrueCounts(ncells_total=ncells, min_popsize=50, i_minpop=2, ngenes=ngenes, nevf=10, 
                                      evf_type="discrete", n_de_evf=9, vary="s", Sigma=0.4, phyla=phyla, randseed=0)

# add technical variation
observed_counts <- True2ObservedCounts(true_counts=true_counts_res[[1]], meta_cell=true_counts_res[[3]], protocol="nonUMI", alpha_mean=0.1, alpha_sd=0.05, gene_len=gene_len, depth_mean=1e5, depth_sd=3e3)

# add technical variation
observed_counts <- True2ObservedCounts(true_counts=true_counts_res[[1]], meta_cell=true_counts_res[[3]], protocol="UMI", alpha_mean=0.05, alpha_sd=0.02, gene_len=gene_len, depth_mean=5e4, depth_sd=3e3)



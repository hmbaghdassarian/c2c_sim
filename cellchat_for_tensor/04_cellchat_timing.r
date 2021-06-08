'Time how long it takes to run cellchat given some samples
Usage:
    Rscript 04_cellchat_timing.r --number=<int> --group=<bool> [--seed=<int>]
    
Options:
    -h --help  Show this screen.
    --number=<int>  number of samples to include
    --group=<bool>  whether to group by context (1) or keep separated by samples (0)
    --seed=<int>  set a seed for random graph generation 
' -> doc

suppressPackageStartupMessages({
    library(CellChat, quietly = T)
    library(patchwork, quietly = T)
    library(RhpcBLASctl, quietly = T)
    library(Matrix, quietly = T)
    library(docopt, quietly = T)
})
options(stringsAsFactors = FALSE)
expression_data_path = '/data2/hratch/immune_CCI/covid/covid_atlas/interim/umi_for_timing/'
input_path = '/data2/hratch/immune_CCI/covid/covid_atlas/interim/timing_inputs/'

RhpcBLASctl::blas_set_num_threads(1) # no multithreading

cell_grouper<-'majorType'

opts = docopt(doc)
for (opt in c('number', 'group', 'seed')){
    val = opts[[opt]]
    if (!is.null(val)){
        opts[[opt]] = as.numeric(val)
        }

set.seed(seed)

samples<-read.csv(paste0(input_path, 'samples_for_timing.csv'))
sample.names<-unlist(strsplit(samples[samples$n_samples == number, 'sample_names'], '; '))

# load LR pairs
# filter for the LR pairs used by Tensor cell2cell
# lr_pairs<-read.csv(paste0(input_data_path,'Tensor-cell2cell-LRpairs.csv'))
# lr_pairs<-lr_pairs$interaction_name
humandb<-CellChatDB.human
# humandb$interaction<-CellChatDB.human$interaction[CellChatDB.human$interaction$interaction_name %in% lr_pairs, ] 
# saveRDS(humandb, paste0(output_results_path, 'humandb.rds'))

# load metadata
md.cell <- read.csv(paste0(input_path,'metadata_for_timing.csv'), row.names = 1)
sample_context_map<-readRDS(paste0(input_path, 'sample_context_map.rds'))
contexts<-unique(sample_context_map)

# load the UMI counts
if (!group){
    counts<-lapply(setNames(sample.names, sample.names), function(sn) t(read.csv(paste0(expression_data_path, sn, '.csv'), row.names = 1)))
}else{ # group by context
    by.context<-lapply(setNames(contexts, contexts), function(context) names(sample_context_map[sample.names][sample_context_map[sample.names] == context]))
    
    group_by_context<-function(context){
        sns<-by.context[[context]]       
        counts<-lapply(sns, function(sn) t(read.csv(paste0(expression_data_path, sn, '.csv'), row.names = 1)))    
        counts<-do.call(cbind, counts)
        return (counts)
    }
    counts<-lapply(setNames(contexts, contexts), function(context) group_by_context(context))                   
                    
}

# create cellchat object for each sample or sample.name
covid.list<-list()
for (sample.name in names(counts)){
    # loop through each sample.name and create a cell type future
    expr<-CellChat::normalizeData(counts[[sample.name]])
    cellchat<-createCellChat(object = as(expr, "dgCMatrix"), meta = md.cell[colnames(expr)], 
                                   group.by = cell_grouper)
    cellchat@DB <- humandb # human organism

    cellchat <- subsetData(cellchat) # subset the expression data of signaling genes, assign to @data.signalling 
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat) # generate @ LR slot used by computeCommunProb
    cellchat <- projectData(cellchat, PPI.human) # shallow sequencing depth
    
    cellchat <- computeCommunProb(cellchat, raw.use = F, type = 'triMean', trim = NULL, seed.use = seed, 
                                 population.size = F) 
    
    # The functional similarity analysis requires the same cell population composition between two datasets.
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    cellchat <- computeCommunProbPathway(cellchat)
    covid.list[[sample.name]]<-cellchat
}

# merge and analyze
cellchat <- mergeCellChat(covid.list, add.names = names(covid.list))
cellchat <- computeNetSimilarityPairwise(cellchat, type = 'structural')
cellchat <- netEmbedding(cellchat, type = 'structural')
cellchat <- netClustering(cellchat, type = 'structural',  do.parallel = F, do.plot = F)
# Manifold learning of the signaling networks for datasets 
pairwise_comparisons<-sapply(as.data.frame(combn(seq(1:length(covid.list)),2)), 
                         function(x) as.numeric(x), simplify = F) # pairwise combination of elements

for (pc in names(pairwise_comparisons)){
    path.dist<-rankSimilarity_(cellchat, type = 'structural', comparison1 = 1:length(covid.list),
                               comparison2 = pairwise_comparisons[[pc]]) 
}                      

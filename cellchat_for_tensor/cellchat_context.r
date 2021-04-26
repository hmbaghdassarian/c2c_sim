
suppressPackageStartupMessages({
    library(CellChat)
    library(patchwork)
    library(RhpcBLASctl)
    library(Matrix)
})
options(stringsAsFactors = FALSE)
input_data_path = '/data2/hratch/immune_CCI/covid_data/'
output_results_path = '/data2/hratch/immune_CCI/covid_results/'

#' Rank the similarity of the shared signaling pathways based on their joint manifold learning
#'
#' @param object CellChat object
#' @param slot.name the slot name of object that is used to compute centrality measures of signaling networks
#' @param type "functional","structural"
#' @param comparison1 a numerical vector giving the datasets for comparison. This should be the same as `comparison` in `computeNetSimilarityPairwise`
#' @param comparison2 a numerical vector with two elements giving the datasets for comparison.
#'
#' If there are more than 2 datasets defined in `comparison1`, `comparison2` can be defined to indicate which two datasets used for computing the distance.
#' e.g., comparison2 = c(1,3) indicates the first and third datasets defined in `comparison1` will be used for comparison.
#' @import ggplot2
#' @importFrom methods slot
#' @return
#' @export
#'
#' @examples
rankSimilarity_ <- function(object, slot.name = "netP", type = c("functional","structural"), comparison1 = NULL,  
                           comparison2 = c(1,2)) {
  type <- match.arg(type)

  if (is.null(comparison1)) {
    comparison1 <- 1:length(unique(object@meta$datasets))
  }
  comparison.name <- paste(comparison1, collapse = "-")
  cat("Compute the distance of signaling networks between datasets", as.character(comparison1[comparison2]), '\n')
  comparison2.name <- names(methods::slot(object, slot.name))[comparison1[comparison2]]

  Y <- methods::slot(object, slot.name)$similarity[[type]]$dr[[comparison.name]]
  group <- sub(".*--", "", rownames(Y))
  data1 <- Y[group %in% comparison2.name[1], ]
  data2 <- Y[group %in% comparison2.name[2], ]
  rownames(data1) <- sub("--.*", "", rownames(data1))
  rownames(data2) <- sub("--.*", "", rownames(data2))

  pathway.show = as.character(intersect(rownames(data1), rownames(data2)))
  data1 <- data1[pathway.show, ]
  data2 <- data2[pathway.show, ]
  euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
  dist <- NULL
  for(i in 1:nrow(data1)) dist[i] <- euc.dist(data1[i,],data2[i,])
  df <- data.frame(name = pathway.show, dist = dist, row.names = pathway.show)
  df <- df[order(df$dist), , drop = F]
  df$name <- factor(df$name, levels = as.character(df$name))

  return(df)
}

#' Run the CellChat pipeline on a single sample from Covid datasets
#'
#' @param data.input <- normalized data matrix (genes x cells) for one sample; dgCMatrix
#' @param meta <- associated metadata for data.input cells; one column must be "celltype"
#' @param organism <- string; default 'human' (only current option)
#' @param n.cores <- number of cores for parallelization (int) default NULL
#' @param score.type <- way to calculate average cell expression, see CellChat::computeCommunProb for details, str, default 'trimean'
#' @param trim <- float [0,1] accompanying score.type, see CellChat::computeCommunProb for details, default NULL
#' @param seed <- int, default NULL
#' @param comparison.type <- string, either functional or 'structural', see CellChat::computeNetSimilarityPairwise, default 'structural'

get_cellchat_sample<-function(data.input, meta, organism = 'human', n.cores = NULL, score.type = 'triMean', 
                              trim = NULL, seed = NULL){
    
    # loop through each context and create a cell type future
    cellchat <- createCellChat(object = data.input, meta = meta, group.by = 'celltype')
    
    if (organism == 'human'){
        cellchat@DB <- CellChatDB.human # human organism
    }else{stop('Only human organism considered currently')}
    # parallelization
    if (!(is.null(n.cores)) | (!n.cores<=1)){
        RhpcBLASctl::blas_set_num_threads(n.cores)
        future::plan("multiprocess", workers = n.cores) 
    }else{
        RhpcBLASctl::blas_set_num_threads(1)
    }
    cellchat <- subsetData(cellchat) # subset the expression data of signaling genes, assign to @data.signalling 


    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat) # generate @ LR slot used by computeCommunProb
    cellchat <- projectData(cellchat, PPI.human) # shallow sequencing depth
    
    # raw.use = F assumes shallow sequencing depth (10x used on COVID dataset)
    # population.size = T considers effect of cell proportion; since a CDR correction is applied, we set to F
    cellchat <- computeCommunProb(cellchat, raw.use = F, type = score.type, trim = trim, seed.use = seed, 
                                 population.size = F) 
    
    # The functional similarity analysis requires the same cell population composition between two datasets.
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    cellchat <- computeCommunProbPathway(cellchat)
    return(cellchat)
    }

#' Run the CellChat pipeline on multiple samples (contexts across Covid datasets)
#'
#' @param object.list <- list of cellchat objects, each for a different sample, output from get_cellchat_sample
#' @param n.cores <- number of cores for parallelization (int) default NULL
#' @param comparison.type <- string, either functional or 'structural', see CellChat::computeNetSimilarityPairwise, default 'structural'

get_cellchat_context<-function(object.list, n.cores = NULL, comparison.type = 'structural'){
    cellchat <- mergeCellChat(object.list, add.names = names(object.list))
    cellchat <- computeNetSimilarityPairwise(cellchat, type = comparison.type)
    cellchat <- netEmbedding(cellchat, type = comparison.type)
    # Manifold learning of the signaling networks for datasets 
    if (!(is.null(n.cores)) | (!n.cores<=1)){
        cellchat <- netClustering(cellchat, type = comparison.type, nCores = n.cores, do.plot = F)
    }else{
        cellchat <- netClustering(cellchat, type = comparison.type, do.parallel = F, do.plot = F)
    }
#     path.dist<-rankSimilarity_(cellchat, type = comparison.type) #CHANGE to pairwise comparison
    pairwise_comparisons<-sapply(as.data.frame(combn(seq(1:length(object.list)),2)), 
                             function(x) as.numeric(x), simplify = F) # pairwise combination of elements
    path_dis_res = list()
    for (pc in names(pairwise_comparisons)){
        print(pc)
        path.dist<-rankSimilarity_(cellchat, type = comparison.type, comparison1 = 1:length(object.list),
                                   comparison2 = pairwise_comparisons[[pc]]) 
        path_dis_res[[paste(names(object.list)[pairwise_comparisons[[pc]]], collapse = '_')]]<-path.dist
    }                            
#     ifl<-rankNet(cellchat, mode = 'comparison', do.stat = T, return.data = T) 
                                 
    if_res = list()
    for (pc in names(pairwise_comparisons)){
        ifl<-rankNet(cellchat, mode = 'comparison', do.stat = T, return.data = T, 
                    comparison = pairwise_comparisons[[pc]]) 
        if_res[[paste(names(object.list)[pairwise_comparisons[[pc]]], collapse = '_')]]<-ifl$signaling.contribution
    }
    res = list()
    res[['pathway_distance']] = path_dis_res
    res[['information.flow']] = if_res
    return(res)
}

organism = 'human'
n.cores = 17
score.type = 'triMean'
trim = NULL
comparison.type = 'structural'
# if by sample, will create a separate cellchat object for each sample, otherwise, will first merge by 
# the condition, and create a separate cellchat object for each condition
by.sample<-FALSE 
if (by.sample){seed<-18}else{seed<-20}



# compare across contexts
if (by.sample){
    covid.list<-readRDS(paste0(output_results_path, 'cellchat_list_bysample.rds')) 
    res<-get_cellchat_context(object.list = covid.list, n.cores = length(covid.list), comparison.type = 'structural')
    saveRDS(res, paste0(output_results_path, 'cellchat_context_results_bysample.rds')) # checkpoint
}else{
    covid.list<-readRDS(paste0(output_results_path, 'cellchat_list_bycontext.rds')) 
    res<-get_cellchat_context(object.list = covid.list, n.cores = length(covid.list), comparison.type = 'structural')
    saveRDS(res, paste0(output_results_path, 'cellchat_context_results_bycontext.rds')) # checkpoint
}

# res<-readRDS(paste0(output_results_path, 'cellchat_context_results.rds')) 



suppressPackageStartupMessages({
    library(CellChat)
    library(patchwork)
    library(RhpcBLASctl)
    library(Matrix)
    library(fgsea)
    library(psych)
})
options(stringsAsFactors = FALSE)
input_data_path = '/data2/hratch/immune_CCI/covid_data/'
output_results_path = '/data2/hratch/immune_CCI/covid_results/timing'

by.sample<-TRUE 
fns<-readRDS(paste0(output_results_path, 'fns.rds'))


# create cellchat object for each sample or context
if (!by.sample){ # by context
    # merge the metadata files
    print('Read in metadata files')
    meta_all <- read.csv(paste0(input_data_path, fns[[1]]$Meta))
    for (sample.name in names(fns)[2:length(fns)]){
        meta_ <- read.csv(paste0(input_data_path, fns[[sample.name]]$Meta))
        meta_all<-rbind(meta_all,meta_)
    }
    
    if (length(unique(meta_all$ID)) != dim(meta_all)[[1]]){
        stop('Overlapping cell barcodes across samples')
    }
    
    # map samples to context
    context_map<-meta_all$sample
    names(context_map)<-meta_all$group
    context_map<-context_map[!duplicated(context_map)]
    context_map<-setNames(names(context_map), context_map)
    contexts = unique(context_map)
    
    # separate metadata by context
    meta_map<-list()
    for (context in contexts){
        meta_map[[context]]<-meta_all[meta_all$group==context,]
    }
    
    # create context-merged expression matrices
    data_input_map<-rep(list(NULL), each = 3)
    names(data_input_map)<-contexts
    counter<-1
    print('Merge expression matrices by context: ')
    for (sample.name in names(fns)){
        print(paste0(counter, ' of ', length(fns)))
        context = context_map[[sample.name]]
        sample.data<-read.csv(paste0(input_data_path, fns[[sample.name]]$DGE))
        rownames(sample.data)<-sample.data$Gene
        sample.data<-sample.data[ , !(colnames(sample.data) %in% c('Gene'))]
        if (is.null(data_input_map[[context]])){ 
            data_input_map[[context]]<-sample.data
        }
        else{
            if (dim(sample.data)[[1]] != dim(data_input_map[[context]])[[1]]){ # genes are same in all dfs so don't have to worry about this
                stop('Not the same genes')
            }else{sample.data<-sample.data[rownames(data_input_map[[context]]),]}
            data_input_map[[context]]<-cbind(data_input_map[[context]], sample.data)
        }
        counter<-counter + 1
        print('------------------------------')    
    }
    
    # create context-specific cellchat objects
    print('Format for inputs to cellchat function')
    counter<-1
    for (context in contexts){
        print(paste0(counter, ' of ', length(contexts)))
        data.input<-data_input_map[[context]]
        meta<-meta_map[[context]]

        cell.names<-colnames(data.input)
        gene.names<-rownames(data.input)
        print('Convert to sparse dgcmatrix')
        data.input<-Reduce(cbind2, lapply(data.input, Matrix, sparse = TRUE))# convert to sparse dgcmatrix, slow
        colnames(data.input) = cell.names
        rownames(data.input) = gene.names

        rownames(meta) = meta$ID
        meta = meta[colnames(data.input), ] # order
        
        data_input_map[[context]]<-data.input
        meta_map[[context]]<-meta
        counter<-counter + 1
        print('------------------------------')
    }
    saveRDS(data_input_map, paste0(output_results_path, 'data_input_map_by_context.rds'))
    saveRDS(meta_map, paste0(output_results_path, 'meta_map_by_context.rds'))
}else{ # by sample
    data_input_map<-list()
    meta_map<-list()
    counter<-1
    for (sample.name in names(fns)){
        print(paste0(counter, ' of ', length(fns)))
        print('Read in data')
        data.input<-read.csv(paste0(input_data_path, fns[[sample.name]]$DGE))
        rownames(data.input)<-data.input$Gene
        data.input<-data.input[ , !(colnames(data.input) %in% c('Gene'))]
     
        cell.names<-colnames(data.input)
        gene.names<-rownames(data.input)
        data.input<-Reduce(cbind2, lapply(data.input, Matrix, sparse = TRUE))# convert to sparse dgcmatrix, slow
        colnames(data.input) = cell.names
        rownames(data.input) = gene.names

        meta = read.csv(paste0(input_data_path, fns[[sample.name]]$Meta))
        rownames(meta) = meta$ID
        meta = meta[colnames(data.input), ] # order
        data_input_map[[sample.name]]<-data.input
        meta_map[[sample.name]]<-meta

        counter<-counter + 1
        print('------------------------------')
        }
    saveRDS(data_input_map, paste0(output_results_path, 'data_input_map_by_sample.rds'))
    saveRDS(meta_map, paste0(output_results_path, 'meta_map_by_sample.rds'))
}

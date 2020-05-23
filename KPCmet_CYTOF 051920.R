# CYTOF DATA ANALYSIS FOR KPC MET MODEL
# Adapting CyTOF pipeline from https://f1000research.com/articles/6-748/v2 to Won's KPC data

# NOTES
# We used mean instead of median for differential analysis and for diagnostic plots
# The input data were spillover corrected using CATALYST
# There is a mild batch effect that we couldn't corrected, so we used batches as a random-effect in differential analysis

rm(list = ls())

library(XML)
library(readxl)
library(flowCore)
library(matrixStats)
library(stringi)
library(reshape2)
library(dplyr)
library(limma)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)
library(FlowSOM)
library(ConsensusClusterPlus)
library(ggridges)
library(ComplexHeatmap)
library(umap)
library(cowplot)
library(lme4)
library(multcomp)

# Specify input files and directories
projectDir = 'C:/Users/who10/Desktop/Research/Rfiles/CYTOF/KPC Mets' # working directory
metaDataFile = paste0(projectDir,'/TME_metadata_cp.xlsx') # meat data
panel_filename <- paste0(projectDir,"/TME_panel_LD.xlsx") # panel with metals, isotopes, and markers
dataFileFolder = paste0(projectDir,"/Corrected_data") # input fcs file (corrected for spillover)
clusterMergeFile = paste0(projectDir,"/TME_cluster_merging.xlsx") # cluster merging
resultsDir = paste0(projectDir,'/Analysis') # folder with output plots and results

# reduced the threshold from 0.05 to 0.01 to make it more stringent
FDR_cutoff <- 0.01

## Define cluster colors (here there are 30 colors)
color_clusters <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
                    "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
                    "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
                    "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
                    "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
                    "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")


# check if the resultsDir exists, create if not
if(!file.exists(resultsDir)) dir.create(resultsDir)


# read in metadata
md <- read_excel(metaDataFile)

## Make sure the condition variables are factors with the right levels
md$condition <- factor(md$condition, levels = c('NL', 'KPC'))

## Make sure the tissue variables are factors
md$tissue <- factor(tolower(md$tissue))

# add batch
md$batch = md$sample_id
for(i in 1:5)
	md$batch[which(md$batch %in% paste0(c("KPCliver", "KPClung", "NLliver", "NLlung"),i))] = paste0('batch',i)
md$batch <- factor(md$batch)

rownames(md) = md$sample_id

## Define colors for the conditions
color_conditions <- c('#d11141', '#00aedb')
names(color_conditions) <- levels(md$condition)

## Define colors for the tissues
color_tissues <- c("#899DA4", "#C93312")
names(color_tissues) <- levels(md$tissue)
shape_tissues <- c(16,17)
names(shape_tissues) <- levels(md$tissue)

## Load content of .fcs files.
fcs_raw <- read.flowSet(md$file_name, path = dataFileFolder, transformation = FALSE, truncate_max_range = FALSE)
## Generate sample IDs corresponding to each cell in the 'expr' matrix
sample_ids <- rep(md$sample_id, fsApply(fcs_raw, nrow))

## Read panel info file
panel <- read_excel(panel_filename)
head(data.frame(panel))

## Replace problematic characters
panel$Metal <- gsub('-', '_', panel$Metal)
panel_fcs <- pData(parameters(fcs_raw[[1]]))
panel_fcs$desc <- gsub('-', '_', panel_fcs$desc)
panel_fcs$desc
rownames(panel_fcs) = panel_fcs$name
# use panel$Antigen to fix description in panel_fcs
# use metal+isotope as mapping between panel from xlsx and panel from the fcs files
panel_fcs[paste0(panel$Metal,panel$Isotope,'Di'),2] = panel$Antigen
#======================================
# TODO try to extract isotope from panel_fcs$desc to use it for mapping
#======================================

## Replace paramater data in flowSet
pData(parameters(fcs_raw[[1]])) <- panel_fcs

## Define variables indicating marker types
subtype_markers <- panel$Antigen[panel$Subtype == 1]
functional_markers <- panel$Antigen[panel$Functional == 1]

## Spot checks
all(subtype_markers %in% panel_fcs$desc)
all(functional_markers %in% panel_fcs$desc)
## arcsinh transformation and column subsetting
fcs <- fsApply(fcs_raw, function(x, cofactor = 5){
  colnames(x) <- panel_fcs$desc
  expr <- exprs(x)
  expr <- asinh(expr[, union(subtype_markers, functional_markers)] / cofactor)
  exprs(x) <- expr
  x
})
fcs

## Extract expression
expr <- fsApply(fcs, exprs)
dim(expr)

# Save matrix for CoGAPS run
#save(expr, file = 'kpcMets_expr.rda')

## Scale expression of all markers to values between 0 and 1
rng <- colQuantiles(expr, probs = c(0.01, 0.99))
expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
expr01[expr01 < 0] <- 0
expr01[expr01 > 1] <- 1

#pdf('Analysis/boxplots.pdf')
#boxplot(expr ~ sample_ids)
#dev.off()
#===============================================
# violin plot markers by batches to study batch effect
ggdf <- data.frame(sample_id = sample_ids, expr[, functional_markers])
mm <- match(ggdf$sample_id, md$sample_id)
ggdf$batch <- factor(md$batch[mm])

ggdf <- melt(ggdf, id.var = c('sample_id','batch'), value.name = 'expression', 
             variable.name = 'antigen')
mm <- match(ggdf$sample_id, md$sample_id)
ggdf$condition <- md$condition[mm]

pdf(paste0(resultsDir,'/KPC_diagnostic_by_func_marker_batch.pdf'), height = 15, width = 10)
print(	ggplot(ggdf, aes(x = batch, y = expression, color = condition, fill = condition)) +
	 geom_violin() +
	 stat_summary(fun.y=mean, geom="point", shape=16, size=2) +
	 facet_wrap(~ antigen, nrow = 7, scales = 'free') 
  )
dev.off()


## Plot per-sample marker expression distributions
    #save all diagnostic plots into a file
makeDiagnosticPlots = function(exprData, fileName = 'diagnostics.pdf', tit = '', fun = mean)
{
	pdf(file = fileName)
	
	# plot 1
		ggdf <- data.frame(sample_id = sample_ids, exprData)
		ggdf <- melt(ggdf, id.var = 'sample_id', value.name = 'expression', 
					 variable.name = 'antigen')
		mm <- match(ggdf$sample_id, md$sample_id)
		ggdf$condition <- md$condition[mm]
		print(ggplot(ggdf, aes(x = expression, color = condition, group = sample_id)) + 
		  geom_density() +
		  facet_wrap(~ antigen, nrow = 4, scales = 'free') + theme_bw() +
		  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
				strip.text = element_text(size = 7),
				axis.text = element_text(size = 5)) + 
		  scale_color_manual(values = color_conditions) )
	# plot 2

		  ## Spot check - number of cells per sample
		cell_table <- table(sample_ids)
		ggdf <- data.frame(sample_id = names(cell_table), 
						   cell_counts = as.numeric(cell_table))
		mm <- match(ggdf$sample_id, md$sample_id)
		ggdf$condition <- md$condition[mm]
		print(ggplot(ggdf, aes(x = sample_id, y = cell_counts, fill = condition)) + 
		  geom_bar(stat = 'identity') + 
		  geom_text(aes(label = cell_counts), hjust = 0.5, vjust = -0.5, size = 2.5) + 
		  theme_bw() +
		  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +  
		  scale_fill_manual(values = color_conditions, drop = FALSE) + 
		  scale_x_discrete(drop = FALSE))

		## Multi-dimensional scaling plot to show similarities between samples
	# plot 3

		## Get the mean marker expression per sample##############################################
		expr_mean_sample_tbl <- data.frame(sample_id = sample_ids, exprData) %>%
		  group_by(sample_id) %>%  summarize_all(funs(fun))
		expr_mean_sample <- t(expr_mean_sample_tbl[, -1])
		colnames(expr_mean_sample) <- expr_mean_sample_tbl$sample_id
		mds <- plotMDS(expr_mean_sample, plot = FALSE)
		ggdf <- data.frame(MDS1 = mds$x, MDS2 = mds$y,
						   sample_id = colnames(expr_mean_sample))
		mm <- match(ggdf$sample_id, md$sample_id)
		ggdf$condition <- md$condition[mm]
		ggdf$tissue <- md$tissue[mm]
		print(ggplot(ggdf, aes(x = MDS1, y = MDS2, color = condition, shape = tissue)) +
		  geom_point(size = 2, alpha = 0.8) +
		  geom_label_repel(aes(label = sample_id)) +
		  theme_bw() +
		  scale_color_manual(values = color_conditions) +
		  scale_shape_manual(values = shape_tissues) +
		  coord_fixed())

	# plot 4
		## Can see differences between tissues as well as conditions
		## Column annotation for the heatmap
		mm <- match(colnames(expr_mean_sample), md$sample_id)
		annotation_col <- data.frame(condition = md$condition[mm], 
									 row.names = colnames(expr_mean_sample))
		annotation_colors <- list(condition = color_conditions)

		## Colors for the heatmap
		color <- colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(100)
		pheatmap(expr_mean_sample, color = color, display_numbers = TRUE,
				 number_color = "black", fontsize_number = 5, 
				 annotation_col = annotation_col, main = tit,
				 annotation_colors = annotation_colors, clustering_method = "average")

	# plot 5

				 ## Define a function that calculates the Non-Redundancy Score per sample
		NRS <- function(x, ncomp = 3){
		  pr <- prcomp(x, center = TRUE, scale. = FALSE)
		  score <- rowSums(outer(rep(1, ncol(x)), pr$sdev[1:ncomp]^2) * 
							 abs(pr$rotation[,1:ncomp]))
		  return(score)
		}

		## Calculate the score
		## May want to do the same with other markers
		nrs_sample <- fsApply(fcs[, subtype_markers], NRS, use.exprs = TRUE)
		rownames(nrs_sample) <- md$sample_id
		nrs <- colMeans(nrs_sample, na.rm = TRUE)

		## Plot the NRS for ordered markers
		## May be helpful to look at tissue instead of condition
		subtype_markers_ord <- names(sort(nrs, decreasing = TRUE))
		nrs_sample <- data.frame(nrs_sample)
		nrs_sample$sample_id <- rownames(nrs_sample)
		ggdf <- melt(nrs_sample, id.var = "sample_id",
					 value.name = "nrs", variable.name = "antigen")
		ggdf$antigen <- factor(ggdf$antigen, levels = subtype_markers_ord)
		mm <- match(ggdf$sample_id, md$sample_id)
		ggdf$condition <- md$condition[mm]
		print(ggplot(ggdf, aes(x = antigen, y = nrs)) +
		  geom_point(aes(color = condition), alpha = 0.9,
					 position = position_jitter(width = 0.3, height = 0)) +
		  geom_boxplot(outlier.color = NA, fill = NA) +
		  stat_summary(fun.y = "mean", geom = "point", shape = 21, fill = "white") +
		  theme_bw() + ggtitle(tit)+ 
		  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) )
		  # scale_color_manual(values = color_conditions)

	dev.off()
}

makeDiagnosticPlots(expr, fileName = paste0(resultsDir,'/KPC_diagnostics_expr_mean.pdf'), 'unscaled')
makeDiagnosticPlots(expr01, fileName = paste0(resultsDir,'/KPC_diagnostics_expr01_mean.pdf'), 'scaled (0-1)')
#makeDiagnosticPlots(expr, fileName = paste0(resultsDir,'/CyTOF_KPC_diagnostics_expr_median.pdf', 'unscaled', fun = median)


#======================BEGIN ANALYSIS HERE=============================#                  
# Clustering

## Cell population identification with FlowSOM and ConsensusClusterPlus
fsom <- ReadInput(fcs, transform = FALSE, scale = FALSE)
set.seed(1234)
som <- BuildSOM(fsom, colsToUse = subtype_markers)

## Get the cell clustering into 100 SOM codes
cell_clustering_som <- som$map$mapping[,1]

## Metaclustering into 20 clusters with ConsensusClusterPlus
codes <- som$map$codes
plot_outdir <- "consensus_plots"
nmc <- 20
mc <- ConsensusClusterPlus(t(codes), maxK = nmc, reps = 100,
                           pItem = 0.9, pFeature = 1, title = plot_outdir, 
                           plot = "png", clusterAlg = "hc", 
                           innerLinkage = "average", finalLinkage = "average",
                           distance = "euclidean", seed = 1234)

## Get cluster ids for each cell
code_clustering1 <- mc[[nmc]]$consensusClass
cell_clustering1 <- code_clustering1[cell_clustering_som]

## Investigate characteristics of identified clusters

####====wrapper function for the heatmap
## Create wrapper function for heatmap to reuse later
plot_clustering_heatmap_wrapper <- function(expr, expr01, cell_clustering, 
                                            color_clusters, cluster_merging = NULL, fileName = 'clustering.pdf'){
  
  ## Calculate the mean expression##################################################
 pdf(fileName) 
  expr_mean <- data.frame(expr, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  expr01_mean <- data.frame(expr01, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  
  ## Calculate cluster frequencies
  
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  
  ## Sort the cell clusters with hierarchical clustering
  
  d <- dist(expr_mean[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_mean[, colnames(expr01)])
  rownames(expr_heat) <- expr01_mean$cell_clustering
  
  ## Colors for the heatmap
  
  color_heat <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100)
  legend_breaks = seq(from = 0, to = 1, by = 0.2)
  labels_row <- paste0(expr01_mean$cell_clustering, " (", clustering_prop ,
                       "%)")
  
  ## Annotation for the original clusters
  
  annotation_row <- data.frame(Cluster = factor(expr01_mean$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  color_clusters1 <- color_clusters[1:nlevels(annotation_row$Cluster)]
  names(color_clusters1) <- levels(annotation_row$Cluster)
  annotation_colors <- list(Cluster = color_clusters1)
  
  ## Annotation for the merged clusters
  
  if(!is.null(cluster_merging)){
    cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
    annotation_row$Merged <- cluster_merging$new_cluster
    color_clusters2 <- color_clusters[1:nlevels(cluster_merging$new_cluster)]
    names(color_clusters2) <- levels(cluster_merging$new_cluster)
    annotation_colors$Merged <- color_clusters2
  }
  
  pheatmap(expr_heat, color = colorRampPalette(brewer.pal(n = 9, name = "YlOrRd"))(100), 
           cluster_cols = FALSE,
           cluster_rows = cluster_rows, 
           labels_row = labels_row,
           #scale="column",
           display_numbers = FALSE, number_color = "black",
           fontsize = 9, fontsize_number = 6,  
           #legend_breaks = legend_breaks,
           annotation_row = annotation_row, 
           annotation_colors = annotation_colors,
           cellwidth = 8,
           cellheight = 8
  )
dev.off()  
}

plot_clustering_heatmap_wrapper(expr = expr[, subtype_markers],
                                expr01 = expr01[, subtype_markers],
                                cell_clustering = cell_clustering1, 
                                color_clusters = color_clusters, fileName = paste0(resultsDir,'/KPC_clustering_20.pdf'))

####=======manual cluster merging based on initial metaclustering results from the wrapper function for heatmap
cluster_merging <- read_excel(clusterMergeFile)
data.frame(cluster_merging)

plot_clustering_heatmap_wrapper(expr = expr[, subtype_markers],
                                expr01 = expr01[, subtype_markers],
                                cell_clustering = cell_clustering1, 
                                color_clusters = color_clusters, 
								cluster_merging = cluster_merging,
								fileName = paste0(resultsDir,'/KPC_clustering_mergedREVISED.pdf'))


## New clustering1m
mm <- match(cell_clustering1, cluster_merging$original_cluster)
cell_clustering1m <- cluster_merging$new_cluster[mm]

mm <- match(code_clustering1, cluster_merging$original_cluster)
code_clustering1m <- cluster_merging$new_cluster[mm]

plot_clustering_heatmap_wrapper(expr = expr[, subtype_markers],
                                expr01 = expr01[, subtype_markers],
                                cell_clustering = cell_clustering1, 
                                color_clusters = color_clusters)

## Create wrapper function for density plot
plot_clustering_distr_wrapper <- function(expr, cell_clustering){
  ## Calculate the mean expression
  cell_clustering <- factor(cell_clustering)
  expr_mean <- data.frame(expr, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  
  ## Sort the cell clusters with hierarchical clustering
  
  d <- dist(expr_mean[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  ## Calculate cluster frequencies
  
  freq_clust <- table(cell_clustering)
  freq_clust <- round(as.numeric(freq_clust)/sum(freq_clust)*100, 2)
  cell_clustering <- factor(cell_clustering,
                            labels = paste0(levels(cell_clustering),
                                            "  (", freq_clust, "%)"))
  
  ## Data organized per cluster
  
  ggd <- melt(data.frame(cluster = cell_clustering, expr),
              id.vars = "cluster", value.name = "expression",
              variable.name = "antigen")
  ggd$antigen <- factor(ggd$antigen, levels = colnames(expr))
  ggd$reference <- "no"
  
  ## The reference data
  
  ggd_bg <- ggd
  ggd_bg$cluster <- "reference"
  ggd_bg$reference <- "yes"
  
  ggd_plot <- rbind(ggd, ggd_bg)
  ggd_plot$cluster <- factor(ggd_plot$cluster,
                             levels = c(levels(cell_clustering)
                                        [rev(cluster_rows$order)], "reference"))
  
  ggplot() +
    geom_density_ridges(data = ggd_plot, aes(x = expression, y = cluster,
                                             color = reference, 
                                             fill = reference), alpha = 0.3) +
    facet_wrap( ~ antigen, scales = "free_x", nrow = 2) +
    theme_ridges() +
    theme(axis.text = element_text(size = 7),  
          strip.text = element_text(size = 7), legend.position = "none")
  
}

#plot_clustering_distr_wrapper(expr = expr[, subtype_markers], 
#                              cell_clustering = cell_clustering1m)

## Create wrapper function for heatmap
plot_clustering_heatmap_wrapper2 <- function(expr, expr01, subtype_markers, 
                                             functional_markers = NULL, sample_ids = NULL, cell_clustering,
                                             color_clusters, cluster_merging = NULL, 
                                             plot_cluster_annotation = TRUE){
  
  ## Calculate the mean expression of lineage markers
  
  expr_mean <- data.frame(expr[, subtype_markers],
                          cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  expr01_mean <- data.frame(expr01[, subtype_markers],
                            cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  
  ## Calculate cluster frequencies
  
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  
  ## Sort the cell clusters with hierarchical clustering
  
  d <- dist(expr_mean[, subtype_markers], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_mean[, subtype_markers])
  
  ## mean expression of activation markers in each sample per cluster
  
  expr_mean_sample_cluster_tbl <- data.frame(expr01[, functional_markers, 
                                                    drop = FALSE], 
                                             sample_id = sample_ids, 
                                             cluster = cell_clustering) %>%
    group_by(sample_id, cluster) %>% summarize_all(funs(mean))
  
  ## Colors for the heatmap
  
  color_heat <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100)
  legend_breaks = seq(from = 0, to = 1, by = 0.2)
  labels_row <- paste0(expr01_mean$cell_clustering, 
                       " (", clustering_prop , "%)") 
  
  ## Annotation for the original clusters
  
  annotation_row1 <- data.frame(Cluster = 
                                  factor(expr01_mean$cell_clustering))
  color_clusters1 <- color_clusters[1:nlevels(annotation_row1$Cluster)]
  names(color_clusters1) <- levels(annotation_row1$Cluster)
  
  ## Annotation for the merged clusters
  
  if(!is.null(cluster_merging)){
    mm <- match(annotation_row1$Cluster, cluster_merging$original_cluster)
    annotated_row2 <- data.frame(Cluster_merging =
                                   factor(cluster_merging$new_cluster[mm]))
    color_clusters2 <-color_clusters[1:nlevels(annotated_row2$Cluster_merging)]
    names(color_clusters2) <- levels(annotated_row2$Cluster_merging)
  }
  
  
  ## Heatmap annotation for the original clusters
  
  ha1 <- Heatmap(annotation_row1, name = "Cluster",
                 col = color_clusters1, cluster_columns = FALSE,
                 cluster_rows = cluster_rows, row_dend_reorder = FALSE,
                 show_row_names = FALSE, width = unit(0.5, "cm"),
                 rect_gp = gpar(col = "grey"))
  
  ## Heatmap annotation for the merged clusters
  
  if(!is.null(cluster_merging)){
    ha2 <- Heatmap(annotated_row2, name = "Cluster \nmerging",
                   col = color_clusters2, cluster_columns = FALSE,
                   cluster_rows = cluster_rows, row_dend_reorder = FALSE,
                   show_row_names = FALSE, width = unit(0.5, "cm"),
                   rect_gp = gpar(col = "grey"))
  }
  
  ## Cluster names and sizes - text
  
  ha_text <- rowAnnotation(text = row_anno_text(labels_row,
                                                gp = gpar(fontsize = 6)), 
                           width = max_text_width(labels_row))
  ## Cluster sizes - barplot
  
  ha_bar <- rowAnnotation("Freq (%)" = row_anno_barplot(
    x = clustering_prop, border = FALSE, axis = TRUE,
    axis_gp = gpar(fontsize = 5), gp = gpar(fill = "#696969", col = "#696969"),
    bar_width = 0.9), width = unit(0.7, "cm"), show_annotation_name = TRUE,
    annotation_name_rot = 0, annotation_name_offset = unit(5, "mm"),
    annotation_name_gp = gpar(fontsize = 7))
  
  ## Heatmap for the lineage markers
  
  ht1 <- Heatmap(expr_heat, name = "Expr", column_title = "Lineage markers",
                 col = color_heat, cluster_columns = FALSE, 
                 cluster_rows = cluster_rows, row_dend_reorder = FALSE, 
                 heatmap_legend_param = list(at = legend_breaks,
                                             labels = legend_breaks, 
                                             color_bar = "continuous"),
                 show_row_names = FALSE, row_dend_width = unit(2, "cm"),
                 rect_gp = gpar(col = "grey"), 
                 column_names_gp = gpar(fontsize = 8))
  
  if(plot_cluster_annotation){
    draw_out <- ha1
  }else{
    draw_out <- NULL
  }
  if(!is.null(cluster_merging)){
    draw_out <- draw_out + ha2 + ht1 + ha_bar + ha_text
  }else{
    draw_out <- draw_out + ht1 + ha_bar + ha_text
  }
  
  ## Heatmaps for the functional markers
  
  if(!is.null(functional_markers)){
    for(i in 1:length(functional_markers)){
      
      ## Rearrange so the rows represent clusters
      
      expr_heat_fun <- as.matrix(dcast(expr_mean_sample_cluster_tbl[,
                                                                    c("sample_id", "cluster", 
                                                                      functional_markers[i])],
                                       cluster ~ sample_id, 
                                       value.var = functional_markers[i])[, -1])
      draw_out <- draw_out + Heatmap(expr_heat_fun,
                                     column_title = functional_markers[i], 
                                     col = color_heat,
                                     cluster_columns = FALSE, 
                                     cluster_rows = cluster_rows,
                                     row_dend_reorder = FALSE, 
                                     show_heatmap_legend = FALSE,
                                     show_row_names = FALSE, 
                                     rect_gp = gpar(col = "grey"),
                                     column_names_gp = gpar(fontsize = 8))
    }
  }
  draw(draw_out, row_dend_side = "left")
}



pdf(paste0(resultsDir,"/1_KPC_heatmaps_per_funMarker.pdf"), width=10, height=4)
for(i in functional_markers)
plot_clustering_heatmap_wrapper2(expr = expr, expr01 = expr01,
                                 subtype_markers = subtype_markers, 
                                 functional_markers = i,
                                 sample_ids = sample_ids, 
                                 cell_clustering = cell_clustering1m,
                                 color_clusters = color_clusters, 
                                 cluster_merging = NULL)

dev.off() 

#===========================================================
## UMAP
#===========================================================

## Find and skip duplicates
dups <- which(!duplicated(expr[, union(subtype_markers,functional_markers)]))

## Data subsampling: create indices by sample
inds <- split(1:length(sample_ids), sample_ids)

## How many cells to downsample per-sample
umap_ncells <- pmin(table(sample_ids), 2000)
umap_ncells5k <- pmin(table(sample_ids), 5000)

## Get subsampled indices

getUmapInds = function(seed = 1234, inds, umap_ncells, dups)
{
	set.seed(seed)
	umap_inds <- lapply(names(inds), function(i){
	  s <- sample(inds[[i]], umap_ncells[i], replace = FALSE)
	  intersect(s, dups)
	})
	umap_inds <- unlist(umap_inds)
	return(umap_inds)
}
umap_inds = getUmapInds(seed = 1234, inds, umap_ncells, dups)
umap_inds5k = getUmapInds(seed = 1234, inds, umap_ncells5k, dups)

## Run umap

custom.settings = umap.defaults
custom.settings$seed = 1234
custom.settings$n.neighbors = 10

umap_out <- umap(expr[umap_inds, union(subtype_markers,functional_markers)], config = custom.settings, method = 'naive')

umapRes = data.frame(umap1 = umap_out$layout[, 1], umap2 = umap_out$layout[, 2],
					expr[umap_inds, union(subtype_markers,functional_markers)],
					sample_id = sample_ids[umap_inds], cell_clustering = factor(cell_clustering1m[umap_inds]))
mm <- match(umapRes$sample_id, md$sample_id)
umapRes$condition <- md$condition[mm]
umapRes$tissue  <- md$tissue[mm]
#adding a batch column 
umapRes$batch <- md$batch[mm]

#=============================
# create dataset with 5k randomly sampled from each sample
# run UMAP and CoGAPS on the same dataset
#umap_out5k <- umap(expr[umap_inds5k, union(subtype_markers,functional_markers)], config = custom.settings, method = 'naive')
#umapRes5k = data.frame(umap1 = umap_out5k$layout[, 1], umap2 = umap_out5k$layout[, 2],
#					expr01[umap_inds5k, union(subtype_markers,functional_markers)],
#					sample_id = sample_ids[umap_inds5k], cell_clustering = factor(cell_clustering1m[umap_inds5k]))
#mm <- match(umapRes5k$sample_id, md$sample_id)
#umapRes5k$condition <- md$condition[mm]
#umapRes5k$tissue = md$tissue[mm]
#umapRes5k$batch<- md$batch[mm]
#==================
# save downsampled data for CoGAPS
#expr5k = expr[umap_inds5k, union(subtype_markers,functional_markers)]
#save(expr5k, file = 'kpcMets_expr_5k.rda')
#rm(expr5k)

## Plot umap colored by clusters

plotUmap = function(umapRes, fileName = 'umap.pdf')
{
	pdf(fileName)
	
	ggp <- ggplot(umapRes,  aes(x = umap1, y = umap2, color = cell_clustering)) +
	  geom_point(size = 0.25) +
	  theme_bw() +
	  scale_color_manual(values = color_clusters) +
	  guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))

	print(ggp)
	#other options
	#ggp + facet_wrap(~ sample_id)
	#ggp + facet_wrap(~ condition)
	#ggp + facet_wrap(~ tissue)
	print(ggp + facet_wrap(~ condition ~tissue))
	print(ggp + facet_wrap(~ batch))
 dev.off()
}
plotUmap(umapRes,paste0(resultsDir,'/3_cluster_umap.pdf'))


##visualizing umap just for condition == KPC
umapResKPC <- umapRes[umapRes$condition=="KPC",]
ggp <- ggplot(umapResKPC,  aes(x = umap1, y = umap2, color = cell_clustering)) +
  geom_point(size = 0.25) +
  theme_bw() +
  scale_color_manual(values = color_clusters) +
  labs(color = "CLUSTERS") +
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 1))
pdf("umap kpc only.pdf", width=5.5, height=4)
ggp + facet_wrap(~ condition ~tissue)
dev.off()

##visualizing umap just for tissue == LUNG
umapResLg <- umapRes[umapRes$tissue=="lung",]
ggp <- ggplot(umapResLg,  aes(x = umap1, y = umap2, color = cell_clustering)) +
  geom_point(size = 0.25) +
  theme_bw() +
  scale_color_manual(values = color_clusters) +
  labs(color = "CLUSTERS") +
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 1))
pdf("umap lung only.pdf", width=5.5, height=4)
ggp + facet_wrap(~ condition ~tissue)
dev.off()

##visualizing umap just for tissue == LUNG
umapResLr <- umapRes[umapRes$tissue=="liver",]
ggp <- ggplot(umapResLr,  aes(x = umap1, y = umap2, color = cell_clustering)) +
  geom_point(size = 0.25) +
  theme_bw() +
  scale_color_manual(values = color_clusters) +
  labs(color = "CLUSTERS") +
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 1))
pdf("umap liver only.pdf", width=5.5, height=4)
ggp + facet_wrap(~ condition ~tissue)
dev.off()


#plotUmap(umapRes5k,paste0(resultsDir,'/umap_unscaled_5k.pdf'))

############################################################################
## Plot umap colored by marker expression

pdf(paste0(resultsDir,"/2_marker_umaps.pdf"), width=9, height=7)
for(i in functional_markers)
{
	ggp <- ggplot(umapRes,  aes(x = umap1, y = umap2, color = umapRes[,i])) +
	  geom_point(size = 0.25) +
	  theme_bw() +
	  scale_color_gradientn(i, colours = 
							  colorRampPalette(rev(brewer.pal(n = 11, 
															  name = "Spectral")))(50))
	print(ggp + facet_wrap(~ condition ~tissue))
}
dev.off()


ggp <- ggplot(umapRes,  aes(x = umap1, y = umap2, color = Lag3)) +
  geom_point(size = 0.25) +
  theme_bw() +
  scale_color_gradientn("Lag3", 
                        values = rev(c(1.2,1,0.8,0.6,0.4,0.3,0.25,0.2,0.15,0.1,0.05,0)),
                        colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50)
                        ,limits = c(0,1.3)
  )
ggp + facet_wrap(~ condition ~tissue)





#===========================================================================
## Plot clusters
#===========================================================================

counts_table <- table(cell_clustering1m, sample_ids)
props_table <- t(t(counts_table) / colSums(counts_table)) * 100
counts <- as.data.frame.matrix(counts_table)
props <- as.data.frame.matrix(props_table)
###### relative abudance of each cell population
ggdf <- melt(data.frame(cluster = rownames(props), props),
             id.vars = "cluster", value.name = "proportion", variable.name = "sample_id")

## Add condition and tissue info
mm <- match(ggdf$sample_id, md$sample_id)
ggdf$condition <- factor(md$condition[mm])
ggdf$tissue <- factor(md$tissue[mm])

barplot <- ggplot(ggdf, aes(x = condition, y = proportion, fill = cluster)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ tissue, scales = "free_x") +
  theme_bw() +
  ylab("% of CD45")+
  xlab("")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_blank()) +
  scale_fill_manual(values = color_clusters)

ggdfkpc <- ggdf[ggdf$condition=="KPC",]
pdf("KPClivervslung.pdf",width=3.5,height=4)
ggplot(ggdfkpc, aes(x = sample_id, y = proportion, fill = cluster)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ tissue, scales = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = color_clusters, name="") +
  xlab("") +
  ylab("Proportion")
dev.off()

ggdf2 <- ggdf[ggdf$condition=="KPC",]
ggdf <- ggdf[ggdf$cluster!="CD62L+Ly6C+ cells",]

boxplot <- ggplot(ggdf, aes(x=tissue, y=proportion, fill=condition))+
  geom_boxplot(outlier.size=0, lwd=0.25, color="black")+
  facet_wrap(~cluster, scales="free", ncol=6)+
  ylab("% of CD45")+
  scale_fill_manual(values=color_tissues)+
  theme(axis.text.x = element_text(size=8, angle=30, vjust=1, hjust=1),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.2),
        axis.line.y = element_line(size=0.2),
        axis.title.y = element_text(size=10),
        axis.ticks.y = element_line(size=0.2),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=7),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(1,'lines'),
        legend.text = element_text(size=5),
        legend.key = element_rect(fill="white")
  )
pdf(paste0(resultsDir,'/KPC_cellprops_Boxplot.pdf'),width=7, height=4);boxplot;dev.off()
pdf(paste0(resultsDir,'/KPC_cellprops_Barplot.pdf'),width=2.5, height=4);barplot;dev.off()


#===========================================================================
## Differential Analyses
#===========================================================================

## Create contrasts

contrast_names <- c("NormalvsKPC")
k1 <- c(0, 1)
K <- matrix(k1, nrow = 1, byrow = TRUE, dimnames = list(contrast_names))

contrast_names2 <- c("LivervsLung")

#=================================================================
## Differential Cell Population Abundance
#=================================================================


# add random-effect for batch
formula_glmer_binomial1 <- y/total ~ condition + (1|sample_id)
formula_glmer_binomial2 <- y/total ~ tissue + (1|sample_id)

## Wrapper function, input is data frame with cell counts 
## (Row is a population, column is a sample), the metadata table, and formula
## Performs diff analysis for each population separately (with contrast K)
## Returns a table with non-adjusted and adjusted p-values

differential_abundance_wrapper <- function(counts, md, formula, K, 
                                           contrast_names){
  ## Fit the GLMM for each cluster separately
  
  ntot <- colSums(counts)
  fit_binomial <- lapply(1:nrow(counts), function(i){
    
    data_tmp <- data.frame(y = as.numeric(counts[i, md$sample_id]),
                           total = ntot[md$sample_id], md)
    
    fit_tmp <- glmer(formula, weights = total, family = binomial,
                     data = data_tmp)
    
    ## Fit contrasts one by one
    
    out <- apply(K, 1, function(k){
      contr_tmp <- glht(fit_tmp, linfct = matrix(k, 1))
      summ_tmp <- summary(contr_tmp)
      pval <- summ_tmp$test$pvalues
      return(pval)
    })
    return(out)
  })
  pvals <- do.call(rbind, fit_binomial)
  colnames(pvals) <- paste0("pval_", contrast_names)
  rownames(pvals) <- rownames(counts)
  
  ## Adjust the p-values
  
  adjp <- apply(pvals, 2, p.adjust, method = "BH")
  colnames(adjp) <- paste0("adjp_", contrast_names)
  return(list(pvals = pvals, adjp = adjp))
}

## First GLMM, looking at differences between conditions (subsetting by tissue)
## Includes observation-level random effects
k1 <- c(0, 1)
K <- matrix(k1, nrow = 1, byrow = TRUE, dimnames = list(contrast_names))

da_out_liver <- differential_abundance_wrapper(counts, 
                                               md = subset(md, tissue == "liver"),
                                               formula = formula_glmer_binomial1,
                                               K = K, 
                                               contrast_names = contrast_names)
apply(da_out_liver$adjp < FDR_cutoff, 2, table)

da_out_lung <- differential_abundance_wrapper(counts, 
                                              md = subset(md, tissue == "lung"),
                                              formula = formula_glmer_binomial1, 
                                              K = K, 
                                              contrast_names = contrast_names)
apply(da_out_lung$adjp < FDR_cutoff, 2, table)

## Second GLMM, looking at variance between tissues (subsetting by condition)
## Includes observation-level random effects

da_out_normal <- differential_abundance_wrapper(counts, 
                                                md = subset(md, condition == "NL"),
                                                formula = formula_glmer_binomial2,
                                                K = K,
                                                contrast_names = contrast_names2)
apply(da_out_normal$adjp < FDR_cutoff, 2, table)

da_out_KPC <- differential_abundance_wrapper(counts, 
                                             md = subset(md, condition == "KPC"),
                                             formula = formula_glmer_binomial2, 
                                             K = K, 
                                             contrast_names = contrast_names2)
apply(da_out_KPC$adjp < FDR_cutoff, 2, table)

## Output tables containing observed cell population proportions

da_output_liver <- data.frame(cluster = rownames(props), 
                              subset(props, select = c(KPCliver1:KPCliver5, NLliver1:NLliver5)),
                              da_out_liver$pvals, da_out_liver$adjp, row.names = NULL)
print(head(da_output_liver), digits = 2)
write.csv(da_output_liver,paste0(resultsDir,"/DA_res_liver_clusters.csv"))

da_output_lung <- data.frame(cluster = rownames(props), 
                             subset(props, select = c(KPClung1:KPClung5, NLlung1:NLlung5)),
                             da_out_lung$pvals, da_out_lung$adjp, row.names = NULL)
print(head(da_output_lung,15),digits = 2)
write.csv(da_output_lung,paste0(resultsDir,"/DA_res_lung_clusters.csv"))

da_output_normal <- data.frame(cluster = rownames(props), 
                               subset(props, select = c(NLliver1:NLlung5)),
                               da_out_normal$pvals, da_out_normal$adjp, 
                               row.names = NULL)
print(head(da_output_normal), digits = 2)
write.csv(da_output_normal,paste0(resultsDir,"/DA_res_normal_clusters.csv"))

da_output_KPC <- data.frame(cluster = rownames(props), 
                            subset(props, select = c(KPCliver1:KPClung5)),
                            da_out_KPC$pvals, da_out_KPC$adjp,
                            row.names = NULL)
print(head(da_output_KPC), digits = 2)
write.csv(da_output_KPC,paste0(resultsDir,"/DA_res_kpc_clusters.csv"))

#heatmap wrapper of p values

#normalize first
normalization_wrapper <- function(expr, th = 2.5){
  expr_norm <- apply(expr, 1, function(x){
    sdx <- sd(x, na.rm = TRUE)
    if(sdx == 0){
      x <- (x - mean(x, na.rm = TRUE))
    }else{
      x <- (x - mean(x, na.rm = TRUE)) / sdx
    }
    x[x > th] <- th
    x[x < -th] <- -th
    return(x)
  })
  expr_norm <- t(expr_norm)
}

#differential heatmap wrapper

plot_differential_heatmap_wrapper <- function(expr_norm, sign_adjp,
                                              condition, color_conditions, th = 2.5){
  ## Order samples by condition
  oo <- order(condition)
  condition <- condition[oo]
  expr_norm <- expr_norm[, oo, drop = FALSE]
  
  ## Create the row labels with adj p-values and other objects for pheatmap
  labels_row <- paste0(rownames(expr_norm), " (", sprintf( "%.02e", sign_adjp), ")")
  labels_col <- colnames(expr_norm)
  annotation_col <- data.frame(condition = factor(condition))
  rownames(annotation_col) <- colnames(expr_norm) 
  annotation_colors <- list(condition = color_conditions)
  
  color <- colorRampPalette(rev(brewer.pal(n=11,name="PuOr")))(100)
  
  breaks = seq(from = -th, to = th, length.out = 101)
  legend_breaks = seq(from = -round(th), to = round(th), by = 1)
  gaps_col <- as.numeric(table(annotation_col$condition))
  
  pheatmap(expr_norm, color = color, breaks = breaks,
           legend_breaks = legend_breaks, cluster_cols = FALSE, cluster_rows = FALSE,
           labels_col = labels_col, labels_row = labels_row, gaps_col = gaps_col,
           annotation_col = annotation_col, annotation_colors = annotation_colors,
           fontsize = 8)
}

## Apply the arcsine-square-root transformation
asin_table <- asin(sqrt((t(t(counts_table) / colSums(counts_table))))) 
asin <- as.data.frame.matrix(asin_table)

#create pdf for p value heatmaps for cluster comparisons
pdf(paste0(resultsDir,"/4 cluster p value heatmaps.pdf"), width=9, height=7)

	## Heatmap for (1) liver
	## Keep significant clusters and sort them by significance
	sign_clusters_liver <- names(which(sort(da_out_liver$adjp[, "adjp_NormalvsKPC"]) < FDR_cutoff))
	## Get the adjusted p-values
	sign_adjp_liver <- da_out_liver$adjp[sign_clusters_liver , "adjp_NormalvsKPC", drop=FALSE]
	## Keep samples for the tissue type and normalize to mean = 0 and sd = 1
	asin_norm <- normalization_wrapper(asin[sign_clusters_liver, ])
	subsettemp <- subset(asin_norm, select = c(KPCliver1:KPCliver5, NLliver1:NLliver5)) #subset is different for each comparison
	mm <- match(colnames(subsettemp), md$sample_id)

	plot_differential_heatmap_wrapper(expr_norm = subsettemp, 
									  sign_adjp = sign_adjp_liver,
									  condition = md$condition[mm], #make sure to check this stratification
									  color_conditions = color_conditions)


	## Heatmap for (2) lung
	## Keep significant clusters and sort them by significance
	sign_clusters_lung <- names(which(sort(da_out_lung$adjp[, "adjp_NormalvsKPC"]) < FDR_cutoff))
	## Get the adjusted p-values
	sign_adjp_lung <- da_out_lung$adjp[sign_clusters_lung , "adjp_NormalvsKPC", drop=FALSE]
	## Keep samples for the tissue type and normalize to mean = 0 and sd = 1
	asin_norm <- normalization_wrapper(asin[sign_clusters_lung, ])
	subsettemp <- subset(asin_norm, select = c(KPClung1:KPClung5, NLlung1:NLlung5)) #subset is different for each comparison
	mm <- match(colnames(subsettemp), md$sample_id)

	plot_differential_heatmap_wrapper(expr_norm = subsettemp, 
									  sign_adjp = sign_adjp_lung,
									  condition = md$condition[mm], #make sure to check this stratification
									  color_conditions = color_conditions)



	## Heatmap for (3) normal
	## Keep significant clusters and sort them by significance
	sign_clusters_normal <- names(which(sort(da_out_normal$adjp[, "adjp_LivervsLung"]) < FDR_cutoff))
	## Get the adjusted p-values
	sign_adjp_normal <- da_out_normal$adjp[sign_clusters_normal , "adjp_LivervsLung", drop=FALSE]
	## Keep samples for the tissue type and normalize to mean = 0 and sd = 1
	asin_norm <- normalization_wrapper(asin[sign_clusters_normal, ])
	subsettemp <- subset(asin_norm, select = c(NLliver1:NLlung5)) #subset is different for each comparison
	mm <- match(colnames(subsettemp), md$sample_id)

	plot_differential_heatmap_wrapper(expr_norm = subsettemp, 
									  sign_adjp = sign_adjp_normal,
									  condition = md$tissue[mm], #make sure to check this stratification
									  color_conditions = color_tissues)


	## Heatmap for (4) KPC
	## Keep significant clusters and sort them by significance
	sign_clusters_KPC <- names(which(sort(da_out_KPC$adjp[, "adjp_LivervsLung"]) < FDR_cutoff))
	## Get the adjusted p-values
	sign_adjp_KPC <- da_out_KPC$adjp[sign_clusters_KPC , "adjp_LivervsLung", drop=FALSE]
	## Keep samples for the tissue type and normalize to mean = 0 and sd = 1
	asin_norm <- normalization_wrapper(asin[sign_clusters_KPC, ])
	subsettemp <- subset(asin_norm, select = c(KPCliver1:KPClung5)) #subset is different for each comparison
	mm <- match(colnames(subsettemp), md$sample_id)

	plot_differential_heatmap_wrapper(expr_norm = subsettemp, 
									  sign_adjp = sign_adjp_KPC,
									  condition = md$tissue[mm], #make sure to check this stratification
									  color_conditions = color_tissues)

dev.off()

#=============================================================================
#########Differential analysis of marker expression stratified by cell population (using random-effect for batch)
#============================================================================

## Get mean marker expression per sample and cluster
expr_mean_sample_cluster_tbl <- data.frame(expr[, functional_markers],
                                           sample_id = sample_ids, cluster = cell_clustering1m) %>%
  group_by(sample_id, cluster) %>%
  summarize_all(funs(mean))

expr_mean_sample_cluster_tbl$condition <- factor(c(rep("KPC",130),rep("NL",130)),levels=c("NL","KPC"))
expr_mean_sample_cluster_tbl$tissue <- factor(rep(c(rep("Liver",65),rep("Lung",65)),2),levels=c("Liver","Lung"))
expr_mean_sample_cluster_tbl$combined <- paste(expr_mean_sample_cluster_tbl$condition, expr_mean_sample_cluster_tbl$tissue, sep=" ")
expr_mean_sample_cluster_tbl$combined <- factor(expr_mean_sample_cluster_tbl$combined, levels = c("NL Liver", "KPC Liver", "NL Lung", "KPC Lung"))

## Melt
expr_mean_sample_cluster_melt <- melt(expr_mean_sample_cluster_tbl,
                                      id.vars = c("sample_id", "cluster","condition","tissue", "combined"), value.name = "mean_expression",
                                      variable.name = "antigen")
ggdffunc <- expr_mean_sample_cluster_melt
ggdffunc2<- ggdffunc[ggdffunc$condition=="KPC",]
ggdffunc2<- ggdffunc2[ggdffunc2$cluster!="CD62L+Ly6C+ cells",]

ggp <- ggplot(ggdffunc2, aes(x=cluster, y=mean_expression, fill=tissue))+
  facet_wrap(~antigen, scales='free', ncol=2)+
  scale_fill_manual(values=c("#d11141","#00aedb"))+
  geom_boxplot(outlier.size=0, lwd=0.25, color="black")+
  theme(axis.text.x = element_text(size=8, angle = 45, hjust=1, vjust=1),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.5),
        axis.line.y = element_line(size=0.5),
        axis.title.y = element_blank(),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=10),
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(3,'lines'),
        legend.text = element_text(size=10),
        legend.key = element_rect(fill="white")
  )
pdf("funcmarkersplots.pdf", width=10, height=20);ggp;dev.off()

## Rearange so the rows represent clusters and markers
expr_mean_sample_cluster <- dcast(expr_mean_sample_cluster_melt,
                                  cluster + antigen ~ sample_id,  value.var = "mean_expression")
rownames(expr_mean_sample_cluster) <- paste0(expr_mean_sample_cluster$cluster,
                                             "_", expr_mean_sample_cluster$antigen)


## Eliminate clusters with low frequency
clusters_keep <- names(which((rowSums(counts < 5) == 0)))
keepLF <- expr_mean_sample_cluster$cluster %in% clusters_keep
expr_mean_sample_cluster <- expr_mean_sample_cluster[keepLF, ]
## Eliminate cases with zero expression in all samples
keep0 <- rowSums(expr_mean_sample_cluster[, md$sample_id]) > 0
expr_mean_sample_cluster <- expr_mean_sample_cluster[keep0, ]





## WRAPPER FUNCTION 
differential_expression_wrapper <- function(expr_mean, md, model = "lmer", formula, K, contrast_names){
  ## Fit LMM or LM for each marker separately
  fit_gaussian <- lapply(1:nrow(expr_mean), function(i){
    data_tmp <- data.frame(y = as.numeric(expr_mean[i, md$sample_id]), md)
    switch(model,
           lmer = {
             fit_tmp <- lmer(formula, data = data_tmp)
           },
           lm = {
             fit_tmp <- lm(formula, data = data_tmp)
           })
    ## Fit contrasts one by one
    out <- apply(K, 1, function(k){
      contr_tmp <- glht(fit_tmp, linfct = matrix(k, 1))
      summ_tmp <- summary(contr_tmp)
      pval <- summ_tmp$test$pvalues
      return(pval)
    })
    return(out)
  })
  pvals <- do.call(rbind, fit_gaussian)
  colnames(pvals) <- paste0("pval_", contrast_names)
  rownames(pvals) <- rownames(expr_mean)
  ## Adjust the p-values
  adjp <- apply(pvals, 2, p.adjust, method = "BH")
  colnames(adjp) <- paste0("adjp_", contrast_names)
  return(list(pvals = pvals, adjp = adjp))
}

# formula without random effect
formula_lm <- y ~ condition 
# and with random effect of batches
formula_lmer <- y ~ condition + (1|batch)

K <- matrix(k1, nrow = 1, byrow = TRUE, dimnames = list(contrast_names))

										 
runDE = function(mean_expr, tis, nBatch = 5, fileName = 'deRes.csv',FDR_cutoff = 0.05, ...)
{
	md_tissue <- subset(md, tissue==tis)
	print(md_tissue)
	de_out1 <- differential_expression_wrapper(
	  expr_mean = mean_expr, #subset here
	  md = md_tissue, #subset here
	  ...)

	apply(de_out1$adjp < FDR_cutoff, 2, table)

	de_output1 <- data.frame(mean_expr,
							 de_out1$pvals, de_out1$adjp, row.names = NULL)
	write.csv(data.frame(mean_expr[,1:2], 
						mean_normal = apply(mean_expr[,paste0(rep(paste0("NL",tis), nBatch), 1:nBatch)],1,mean),
						mean_KPC = apply(mean_expr[,paste0(rep(paste0("KPC",tis), nBatch), 1:nBatch)],1,mean),
						de_out1$pvals, de_out1$adjp, row.names = NULL), file = fileName)
	print(head(de_output1), digits = 2)
	return(de_out1)
}


#create a file of all mean expression
write.csv(data.frame(expr_mean_sample_cluster),file="Analysis/DEmarkers_all.csv")


# fit the model with the batch as random-effect
### (1) liver normal vs. kpc
de_out_lv_r = runDE(expr_mean_sample_cluster, tis = 'liver', fileName = paste0(resultsDir,'/KPC_DE_markers_liver_normalVsKPC_r.csv'), 
	model = 'lmer', formula = formula_lmer, K = K, contrast_names = contrast_names)
apply(de_out_lv_r$adjp < FDR_cutoff, 2, table)

# no random effect
#de_out_lv = runDE(expr_mean_sample_cluster, tis = 'liver', fileName = 'Analysis/KPC_DE_markers_liver_normalVsKPC.csv', 
#	model = 'lm', formula = formula_lm, K = K, contrast_names = contrast_names)
#apply(de_out_lv$adjp < FDR_cutoff, 2, table)
#head(data.frame(de_out_lv))

### (2) ling normal vs. kpc
de_out_lg_r = runDE(expr_mean_sample_cluster, tis = 'lung', fileName = paste0(resultsDir,'/KPC_DE_markers_lung_normalVsKPC_r.csv'), 
	model = 'lmer', formula = formula_lmer, K = K, contrast_names = contrast_names)
apply(de_out_lg_r$adjp < FDR_cutoff, 2, table)

plot_differential_marker_heatmap_wrapper = function(de_out, # output from runDE
	adjpCol = "adjp_NormalvsKPC", # adjusted p-value column in de_out
	expr_mean, # mean expression matrix
	expr_mean_melt, # melt matrix
	md, # metadata
	samp, # samples to plot
	colBy = color_conditions, # colors to use for conditions
	plotBy = 'condition', # split by condtions
	fileName = 'heatmap.pdf', # file name
	h = 7, # file height
	w = 7) # file width
{
	## Keep the significant markers, sort them by significance and group by cluster
	sign_clusters_markers <- names(which(de_out$adjp[, adjpCol] < FDR_cutoff)) #specify the correct comparison
	oo <- order(expr_mean[sign_clusters_markers, "cluster"],
				de_out$adjp[sign_clusters_markers, adjpCol])
	sign_clusters_markers <- sign_clusters_markers[oo]

	## Get the significant adjusted p-values
	sign_adjp <- de_out$adjp[sign_clusters_markers , adjpCol]

	## Normalize expression to mean = 0 and sd = 1
	expr_s <- expr_mean[sign_clusters_markers,md$sample_id]
	expr_s_tissue <- subset(expr_s, select = c(samp)) #subset here
	expr_mean_sample_cluster_norm <- normalization_wrapper(expr_s_tissue) #make sure to specify correct expr_s

	mm <- match(colnames(expr_mean_sample_cluster_norm), md$sample_id)

	pdf(fileName, height = h, width = w)
	##PLOT HEATMAP
	plot_differential_heatmap_wrapper(expr_norm = expr_mean_sample_cluster_norm,
									  sign_adjp = sign_adjp, condition = unlist(md[mm,plotBy]),
									  color_conditions = colBy)

	## PLOT BY MARKER
		# create sample names to be plotted using tissue and nBatch
	expr_mean_sample_cluster_tissue <- subset(expr_mean_melt, subset = 
											expr_mean_melt$sample_id %in%samp)
	ggdf <- expr_mean_sample_cluster_tissue[expr_mean_sample_cluster_tissue$cluster  #make sure to specify correct expr_s
										 %in% unique(expr_mean$cluster), ]
	mm <- match(ggdf$sample_id, md$sample_id)
	ggdf$condition <- factor(unlist(md[mm,plotBy]))
	ggdf$sample_id <- factor(md$sample_id[mm])

	print(	ggplot(ggdf) +
		  geom_boxplot(aes(x = antigen, y = mean_expression,
						   color = condition),
					   position = position_dodge(), alpha = 0.5, outlier.color = NA) +
		  geom_point(aes(x = antigen, y = mean_expression, color = condition,
						 shape = sample_id), alpha = 0.8, position = position_jitterdodge(),
					 size = 0.7) +
		  facet_wrap(~ cluster, scales = "free_y", ncol=2) +
		  theme_bw() +
		  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
		  ylab("mean expression") +
		  labs(shape="samples",color="group") +
		  scale_color_manual(values = colBy) +
		  scale_shape_manual(values = c(1:20))
	)
	dev.off()

}

plot_differential_marker_heatmap_wrapper(de_out_lv_r,adjpCol = "adjp_NormalvsKPC", expr_mean_sample_cluster, expr_mean_sample_cluster_melt, md, 
	samp = paste0(rep(paste0(c("KPC","NL"),'liver'), each = 5), 1:5),
	fileName = paste0(resultsDir,'/KPC_DE_markers_liver_normalVsKPC_heatmap.pdf'), h = 15, w = 10)
plot_differential_marker_heatmap_wrapper(de_out_lv_r,adjpCol = "adjp_NormalvsKPC", expr_mean_sample_cluster, expr_mean_sample_cluster_melt, md, 
	samp = paste0(rep(paste0(c("KPC","NL"),'lung'), each = 5), 1:5),
	fileName = paste0(resultsDir,'/KPC_DE_markers_lung_normalVsKPC_heatmap.pdf'), h = 15, w = 10)

										
### (3) Normal liver vs. lung

K <- matrix(k1, nrow = 1, byrow = TRUE, dimnames = list(contrast_names2))

#Subset the relevant groups only
samp = paste0(rep(c("NLliver","NLlung"), each = 5), 1:5)
expr_mean_sample_cluster_nl <- subset(expr_mean_sample_cluster_melt, subset = 
                                    expr_mean_sample_cluster_melt$sample_id %in%samp)
md_nl <- subset(md, condition=="NL")
expr_mean_sample_cluster_nl_1 <- dcast(expr_mean_sample_cluster_nl,
                                  cluster + antigen ~ sample_id,  value.var = "mean_expression")
rownames(expr_mean_sample_cluster_nl_1) <- paste0(expr_mean_sample_cluster_nl_1$cluster,
                                             "_", expr_mean_sample_cluster_nl_1$antigen)

de_out1 <- differential_expression_wrapper(
  expr_mean = expr_mean_sample_cluster_nl_1,
  md = md_nl, #subset here
  model = "lmer", formula = y ~ tissue + (1|batch), K = K, contrast_names = contrast_names2)

de_output1 <- data.frame(expr_mean_sample_cluster_nl_1,
                         de_out1$pvals, de_out1$adjp, row.names = NULL)
print(head(de_output1), digits = 2)
write.csv(de_output1,file="Analysis/KPC_DE_markers_normal_lungvsliver_r.csv")

plot_differential_marker_heatmap_wrapper(de_out1,adjpCol = "adjp_LivervsLung", expr_mean_sample_cluster, expr_mean_sample_cluster_melt, md, 
	samp = paste0(rep(paste0("NL",c('liver','lung')), each = 5), 1:5),
	colBy = color_tissues, plotBy = 'tissue',
	fileName = paste0(resultsDir,'/KPC_DE_markers_normal_liverVsLung_heatmap.pdf'), h = 15, w = 10)

### (4) KPC liver vs. lung

#Subset the relevant groups only
samp = paste0(rep(c("KPCliver","KPClung"), each = 5), 1:5)
expr_mean_sample_cluster_kpc <- subset(expr_mean_sample_cluster_melt, subset = 
                                        expr_mean_sample_cluster_melt$sample_id %in%samp)

md_kpc <- subset(md, condition=="KPC")
expr_mean_sample_cluster_kpc_1 <- dcast(expr_mean_sample_cluster_kpc,
                                       cluster + antigen ~ sample_id,  value.var = "mean_expression")
rownames(expr_mean_sample_cluster_kpc_1) <- paste0(expr_mean_sample_cluster_kpc_1$cluster,
                                                  "_", expr_mean_sample_cluster_kpc_1$antigen)
de_out1 <- differential_expression_wrapper(
  expr_mean = expr_mean_sample_cluster_kpc_1,
  md = md_kpc, #subset here
  model = "lmer", formula = y ~ tissue + (1|batch), K = K, contrast_names = contrast_names2)

de_output1 <- data.frame(expr_mean_sample_cluster_kpc_1,
                         de_out1$pvals, de_out1$adjp, row.names = NULL)
print(head(de_output1), digits = 2)
write.csv(de_output1,file="Analysis/KPC_DE_markers_kpc_lungvsliver_r.csv")

plot_differential_marker_heatmap_wrapper(de_out1,adjpCol = "adjp_LivervsLung", expr_mean_sample_cluster, expr_mean_sample_cluster_melt, md, 
                                         samp = paste0(rep(paste0("KPC",c('liver','lung')), each = 5), 1:5),
                                         colBy = color_tissues, plotBy = 'tissue',
                                         fileName = paste0(resultsDir,'/KPC_DE_markers_KPC_liverVsLung_heatmap.pdf'), h = 15, w = 10)



#==================================
# all normals vs all KPC
K <- matrix(k1, nrow = 1, byrow = TRUE, dimnames = list(contrast_names))

de_out1 <- differential_expression_wrapper(
  expr_mean = expr_mean_sample_cluster, 
  md = md,
  model = "lmer", formula = y ~ tissue + (1|batch), K = K, contrast_names = contrast_names)

apply(de_out1$adjp < FDR_cutoff, 2, table)

de_output1 <- data.frame(expr_mean_sample_cluster,
                         de_out1$pvals, de_out1$adjp, row.names = NULL)
write.csv(data.frame(expr_mean_sample_cluster[,1:2],
                         de_out1$pvals, de_out1$adjp, row.names = NULL), file = paste0(resultsDir,'/KPC_DE_markers_normalVsKPC.csv'))

print(head(de_output1), digits = 2)
samp = paste0(rep(as.vector(sapply(c("KPC","NL"),function(x) paste0(x,c('liver','lung')))), each = 5), 1:5)
plot_differential_marker_heatmap_wrapper(de_out1,adjpCol = "adjp_NormalvsKPC", expr_mean_sample_cluster, expr_mean_sample_cluster_melt, md, samp,
	fileName = paste0(resultsDir,'/KPC_DE_markers_normalVsKPC_heatmap.pdf'), h = 15, w = 10)




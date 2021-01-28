

library(monocle)
library(devtools)

#########################################################################################################################################
#########################################################################################################################################

# Extract data, phenotype data, and feature data from the SeuratObject
data <- as(as.matrix(immune.combined@assays$integrated@data), 'sparseMatrix')

pd <- new('AnnotatedDataFrame', data = immune.combined@meta.data)

fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))

fd <- new('AnnotatedDataFrame', data = fData)

# Construct monocle cds
HSMM <- newCellDataSet(data,
                       phenoData = pd,
                       featureData = fd,
                       #lowerDetectionLimit = 0.5,
                       expressionFamily = uninormal())# since I have already normalized, thresholded and scalled in Suerat v3.0.0.9150

# View data
pData(HSMM)
fData(HSMM)


# Run ordering algorithm

var_genes <- seurat_object[["integrated"]]@var.features
ordering_genes <- var_genes

HSMM <- setOrderingFilter(HSMM, ordering_genes)
print(dim(exprs(HSMM)))

# reduce dimension - do not normalize or include pseudo count. Use monocle scaling

HSMM <- reduceDimension(HSMM,norm_method="none", 
                        reduction_method="DDRTree",
                        max_components=2,
                        scaling=TRUE,
                        verbose=TRUE,
                        pseudo_expr=0)

# First decide what you want to color your cells by

print(head(pData(HSMM)))

# order cells

HSMM <- orderCells(HSMM,reverse = F)

# Plot trajectory

plot_cell_trajectory(HSMM, 
                     color_by = "Pseudotime",
                     theta = -15,
                     show_branch_points = FALSE,
                     show_tree = TRUE,
                     cell_size = 1) + theme(legend.position = "right") + facet_wrap(~time, ncol = 3)


saveRDS(HSMM, file = "C:/Users/Utente/Desktop/MONOCLE.rds")

library(knitr)
library(markdown)
library(rmarkdown)

cluster_order <- list(Cardiomyocytes=c(0,1,3,6,7,8,14),
                      Endothelial_cells=c(5),
                      Fibroblasts=c(4,11,12,16),
                      Myeloid=c(9),
                      Lymphoid=c(13,18),
                      Adipocytes=c(17),
                      SMC=c(10),
                      PC=c(2),
                      Neuronal=c(15)) #combined some cluster to match cell types from hca.

path <- "/media/Helios_scStorage/SingleCell/Human/WholeHeart/Heart_cell_atlas_aortic_stenosis_harmonized/Harmonized_whole_dataset_only_as_only_septum_only_nuc/"

for (cell_cluster in 1:length(cluster_order)){
  cell_type_name <- names(cluster_order)[cell_cluster]
  rmarkdown::render(paste0(path, "template.Rmd"),
                    output_file =  paste0("Human_Heart_snRNAseq_Adult-Aortic-Stenosis-HF_", cell_type_name, ".html"), 
                    output_dir = paste0(path, cell_type_name))
}

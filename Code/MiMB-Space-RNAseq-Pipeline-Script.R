###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
#################################################### MiMB Space RNAseq Pipeline Script ####################################################
###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
#
# DESCRIPTION:  This script is linked to the "Transcriptomics in space: a basic R pipeline for analysing bulk RNA-sequencing data" chapter 
#               as part of the "Biomedicine and Biotechnology in Space - Methods and Protocols in Systems Biology" book in the "Methods in 
#               Molecular Biology" series by Springer-Nature.
#
# SECTIONS:     3.2. Source RNA-seq reads
#               3.3. Quality control check of RNA-seq reads
#               3.4. Quantification of gene expression from RNA-seq reads
#               3.5. Downstream analysis of RNA-seq counts
#
# NOTE:         Section/subsection numbers in the script match the associated section/subsection numbers in the chapter text.
#
###########################################################################################################################################
######################################################## 3.2. Source RNA-seq reads ########################################################
###########################################################################################################################################

###########################################################################################################################################
###################################################### 3.2.1. Create metadata table #######################################################

### Initialise table with sample run IDs ###
meta.data = data.frame('Run' = c('SRR26332897', 'SRR26332898', 'SRR26332902', 'SRR26332903', 'SRR26332906', 'SRR26332907'))

### Add condition column ###
meta.data$Condition = ifelse(meta.data$Run %in% c('SRR26332906', 'SRR26332902', 'SRR26332897'), 'GC', 'SF')

### Update row names as run IDs ###
rownames(meta.data) = meta.data$Run

###########################################################################################################################################
####################################################### 3.2.2. Obtain RNA-seq reads #######################################################

### Install required R packages ###

# Install remotes R package (v2.5.0) #
install.packages(pkgs = 'remotes')
remotes::install_version(package = 'remotes', version = '2.5.0', force = T)

# Install BiocManager R package (v1.30.25) #
remotes::install_version(package = 'BiocManager', version = '1.30.25', force = T)

# Install GEOfastq R package (v1.14.0) #
BiocManager::install(pkgs = 'GEOfastq', version = '3.20', force = T, ask = F)

### Create directory to store fastq files in ###
dir.create(path = './Data/Reads', recursive = T)

### Get urls for reads ###

# Define part of url for general EBI-ENA database #
url_pre = 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq'

# Add url for each read mate to meta data table #
meta.data$url_R1 = sapply(meta.data$Run, function(x) paste0(url_pre, '/', GEOfastq::get_dldir(x), '/', x, '_1.fastq.gz'))
meta.data$url_R2 = sapply(meta.data$Run, function(x) paste0(url_pre, '/', GEOfastq::get_dldir(x), '/', x, '_2.fastq.gz'))

### Download reads ###
read_urls = c(meta.data$url_R1, meta.data$url_R2)
read_destinations = sapply(read_urls, function(x) paste0('./Data/Reads/', tail(unlist(strsplit(x, split = '/')), n = 1)))
Map(function(u, d) download.file(url = u, destfile = d, method = 'curl'), read_urls, read_destinations)

###########################################################################################################################################
############################################### 3.3. Quality control check of RNA-seq reads ###############################################
###########################################################################################################################################

### Install fastqcr R package (v0.1.3) ###
remotes::install_version(package = 'fastqcr', version = '0.1.3', force = T)

### Install FastQC (v0.12.1) ###

# Create directory to store FastQC software (and the others) #
dir.create(path = './Software', recursive = T) 

# Define download url #
fastqc_url = 'https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip'

# Download and unzip #
fastqcr::fastqc_install(url = fastqc_url, dest.dir = "./Software")

### Run FastQC via the fastqcr package ###

# Create directory to save outputs #
dir.create(path = './Outputs/FastQC', recursive = T)

# Run fastQC #
fastqcr::fastqc(fq.dir = './Data/Reads', qc.dir = './Outputs/FastQC', fastqc.path = './Software/FastQC/fastqc')

# Generate single qc report across all samples #
fastqcr::qc_report(qc.path = './Outputs/FastQC', result.file = './Outputs/FastQC/multi-qc-report')

###########################################################################################################################################
######################################## 3.4. Quantification of gene expression from RNA-seq reads ########################################
###########################################################################################################################################

###########################################################################################################################################
######################################## 3.4.1. Pseudo-alignment to quantify transcript abundances ########################################

### Install kallisto ###

# Define download url #
kallisto_url = 'https://github.com/pachterlab/kallisto/releases/download/v0.50.0/kallisto_mac_m1-v0.50.0.tar.gz'

# Download and unzip #
temp = tempfile()
download.file(url = kallisto_url, destfile = temp, method = 'curl', extra = '-L')
untar(tarfile = temp, exdir = './Software')
unlink(temp)

### Obtain reference transcriptome ###

# Create directory to store Ensembl cdna file #
dir.create(path = './Data/Ensembl_cdna', recursive = T)

# Define url for Ensembl cdna file #
cdna_url = 'https://ftp.ensembl.org/pub/release-113/fasta/caenorhabditis_elegans/cdna/Caenorhabditis_elegans.WBcel235.cdna.all.fa.gz'

# Download #
download.file(url = cdna_url, destfile = './Data/Ensembl_cdna/Caenorhabditis_elegans.WBcel235.cdna.all.fa.gz', method = 'curl')

### Build transcriptome index ###

# Create directory to store kallisto index #
dir.create(path = './Outputs/kallisto/Index', recursive = T)

# Run kallisto index via system integration #
kallisto_path = './Software/kallisto/kallisto'
cdna_input_path = './Data/Ensembl_cdna/Caenorhabditis_elegans.WBcel235.cdna.all.fa.gz'
index_out_path = './Outputs/kallisto/Index/WBcel235.cdna.all.transcripts.idx'
system(paste(kallisto_path, 'index', '-i', index_out_path, cdna_input_path))

### Run quantification algorithm ###

# Create directory to store outputs #
dir.create(path = './Outputs/kallisto/Quants', recursive = T)

# Run kallisto quant via system integration #
runs = meta.data$Run
R1s = paste0('./Data/Reads', '/', meta.data$Run, '_1.fastq.gz')
R2s = paste0('./Data/Reads', '/', meta.data$Run, '_2.fastq.gz')
out_dirs = paste0('./Outputs/kallisto/Quants', '/', meta.data$Run)
opt_args = c('--rf-stranded')
Map(function(run, R1, R2, out) system(paste(kallisto_path, 'quant', '-i', index_out_path, '-o', out, opt_args, R1, R2)), runs, R1s, R2s, out_dirs)

###########################################################################################################################################
######################################## 3.4.2. Inferring gene expression from transcript abundances ######################################

### Obtain transcript-to-gene ID mappings ###

# Install biomaRt R package (v2.62.1) #
BiocManager::install(pkgs = 'biomaRt', version = '3.20', force = T, ask = F)

# Define gene annotations database for C. elegans from Biomart #
biomart_db_ce = biomaRt::useEnsembl(biomart = 'genes', dataset = 'celegans_gene_ensembl', version = 113)

# Pull the gene annotations #
gene_annot = biomaRt::getBM(mart =  biomart_db_ce, attributes = c('ensembl_transcript_id', 'ensembl_gene_id', 'external_gene_name'), uniqueRows = T)

### Infer gene expression from transcript abundances ###

# define core transcript-to-gene mapping data frame #
t2g = data.frame('transcript' = gene_annot$ensembl_transcript_id, 'gene' = gene_annot$ensembl_gene_id)

# Install tximport R package (v1.34.0) #
BiocManager::install(pkgs = 'tximport', version = '3.20', force = T, ask = F)

# Import transcript-level estimates and summarise as gene counts #
txi_kallisto = tximport::tximport(files = paste0('./Outputs/kallisto/Quants', '/', meta.data$Run, '/', 'abundance.tsv'), type = 'kallisto', tx2gene = t2g, dropInfReps = T)

###########################################################################################################################################
################################################ 3.5. Downstream analysis of RNA-seq counts ###############################################
###########################################################################################################################################

###########################################################################################################################################
######################################################## 3.5.1. Intermediatory steps ######################################################

### Initialise DESeq Dataset (DDS) object ###

# Install DESeq2 R package (v1.46.0) #
BiocManager::install(pkgs = 'DESeq2', version = '3.20', force = T, ask = F)

# Create DDS object from tximport object #
dds = DESeq2::DESeqDataSetFromTximport(txi = txi_kallisto, colData = meta.data, design = ~ Condition)

# Add gene annotation data to DDS object #
GenomicRanges::mcols(dds)[, c('ensembl_gene_id', 'external_gene_name')] = gene_annot[match(names(dds), gene_annot$ensembl_gene_id), c('ensembl_gene_id', 'external_gene_name')]

### Gene pre-filtering ###

# Remove lowly expressed genes #
dds = dds[rowSums(DESeq2::counts(dds) >= 10) == ncol(DESeq2::counts(dds)), ]

###########################################################################################################################################
#################################################### 3.5.2. Principal component analysis ##################################################

### Install PCAtools R package (v2.18.0) ###
BiocManager::install(pkgs = 'PCAtools', version = '3.20', force = T, ask = F)

### Extract principal component analysis data ###
pca_dat = PCAtools::pca(mat = SummarizedExperiment::assay(DESeq2::vst(dds)), metadata = SummarizedExperiment::colData(dds), removeVar = 0.75)

### Construct principal component analysis plot ###
pca_plot = PCAtools::biplot(pcaobj = pca_dat, colby = 'Condition', pointSize = 6, legendPosition = 'right', labSize = 4.0)

### Save principal component analysis plot ###

# Create directory to save plot in #
dir.create(path = './Outputs/Analysis/Plots', recursive = T)

# Save plot as PDF file #
ggplot2::ggsave(plot = pca_plot, filename = 'PCA_Plot.pdf',path = './Outputs/Analysis/Plots', device = "pdf", height = 8, width = 8, units = "in")

###########################################################################################################################################
############################################### 3.5.3. Differential gene expression analysis ##############################################

### Calculate differential expression ###
dds = DESeq2::DESeq(dds)

### Extract initial results ###
de_res = DESeq2::results(object = dds, contrast = c('Condition', 'SF', 'GC'), cooksCutoff = F, independentFiltering = F)

### Perform log2-fold change shrinkage ###

# Install ashr R package (v2.2-63) #
remotes::install_version(package = 'ashr', version = '2.2-63', force = T)

# Apply ashr log2-fold change shrinkage method #
de_res = DESeq2::lfcShrink(dds = dds, res = de_res, type = 'ashr')

### Save differential expression results as a text file ###

# Create new directory to store text file results #
dir.create(path = './Outputs/Analysis/Text_Files', recursive = T)

# Create results dataframe including gene IDs #
de_res_df = cbind(GenomicRanges::mcols(dds)[, c('ensembl_gene_id', 'external_gene_name')], de_res)

# Write to text file #
write.table(x = de_res_df, file = './Outputs/Analysis/Text_Files/SF_vs_GC_DESeq2_Results.csv', sep = ',', row.names = F)

### Volcano plot of  differential gene expression results ###

# Install EnhancedVolcano R package (v1.24.0) #
BiocManager::install(pkgs = 'EnhancedVolcano', version = '3.20', force = T, ask = F)

# Construct volcano plot #
volcano_plot = EnhancedVolcano::EnhancedVolcano(toptable = de_res_df, 
                                                x = 'log2FoldChange', 
                                                y = 'padj',
                                                lab = paste0("italic('", de_res_df$external_gene_name, "')"),
                                                parseLabels = T,
                                                xlab = 'log2-FC',
                                                ylab = '-log10(padj)',
                                                pCutoffCol = 'padj',
                                                pCutoff = 	0.05,
                                                FCcutoff = 3,
                                                pointSize = 2.5,
                                                labSize = 4.0,
                                                boxedLabels = T,
                                                drawConnectors = T,
                                                colAlpha = 0.75,
                                                legendLabels = c('NS', 'log2-FC', 'padj', 'padj and log2-FC'),
                                                col = c("grey", "grey", "lightgreen", "forestgreen"),
                                                title = 'SF versus GC',
                                                subtitle = 'Volcano plot')

# Save volcano plot as PDF file #
ggplot2::ggsave(plot = volcano_plot, filename = 'Volcano_Plot.pdf',path = './Outputs/Analysis/Plots', device = "pdf", height = 8, width = 8, units = "in")

###########################################################################################################################################
################################################# 3.5.4. Gene Ontology enrichment analysis ################################################

### Install required R packages ###

# Install clusterProfiler R package (v4.14.6) #
BiocManager::install(pkgs = 'clusterProfiler', version = '3.20', force = T, ask = F)

# Install org.Ce.eg.db R package (v3.20.0) #
BiocManager::install(pkgs = 'org.Ce.eg.db', version = '3.20', force = T, ask = F)

### Define input gene lists ###

# differentially expressed gene lists #
ur_genes = de_res_df[de_res_df$padj < 0.05 & de_res_df$log2FoldChange > 0, 'ensembl_gene_id']
dr_genes = de_res_df[de_res_df$padj < 0.05 & de_res_df$log2FoldChange < 0, 'ensembl_gene_id']

# Background gene list #
bg_genes = de_res_df$ensembl_gene_id

### Run Gene Ontology enrichment via over-representation analysis ###
ur_go_bp = clusterProfiler::enrichGO(gene = ur_genes, universe = bg_genes, OrgDb = org.Ce.eg.db::org.Ce.eg.db, keyType = 'ENSEMBL', ont = 'BP', readable = T) 
dr_go_bp = clusterProfiler::enrichGO(gene = dr_genes, universe = bg_genes, OrgDb = org.Ce.eg.db::org.Ce.eg.db, keyType = 'ENSEMBL', ont = 'BP', readable = T) 

### Save results as text files ###
write.table(x = ur_go_bp@result, file = './Outputs/Analysis/Text_Files/SF_vs_GC_Upregulated_GO_BP_Results.csv', sep = ',', row.names = F)
write.table(x = dr_go_bp@result, file = './Outputs/Analysis/Text_Files/SF_vs_GC_Downregulated_GO_BP_Results.csv', sep = ',', row.names = F)

### Visualise Gene Ontology enrichment results as enrichment maps ###

# Install enrichplot R package (v1.26.6) #
BiocManager::install(pkgs = 'enrichplot', version = '3.20', force = T, ask = F)

# Create enrichment map plots #
go_bp_plot_ur = enrichplot::emapplot(enrichplot::pairwise_termsim(ur_go_bp)) 
go_bp_plot_dr = enrichplot::emapplot(enrichplot::pairwise_termsim(dr_go_bp)) 

# Save enrichment map plots as PDF files #
ggplot2::ggsave(plot = go_bp_plot_ur, filename = 'EMAP_Plot_GO_BP_Upregulated.pdf',path = './Outputs/Analysis/Plots', device = "pdf", height = 10, width = 10, units = "in")
ggplot2::ggsave(plot = go_bp_plot_dr, filename = 'EMAP_Plot_GO_BP_downregulated.pdf',path = './Outputs/Analysis/Plots', device = "pdf", height = 10, width = 10, units = "in")

###########################################################################################################################################
######################################################### 3.5.5. Network analysis #########################################################

### Install STRINGdb R package (v2.18.0) ###
BiocManager::install(pkgs = 'STRINGdb', version = '3.20', force = T, ask = F)

### Create a new STRING database object ###
string_db = STRINGdb::STRINGdb$new(version="12.0", species = 6239, score_threshold = 900)

### Map differential expression results to string ID database ###
de_res_df_mapped = string_db$map(my_data_frame = de_res_df, my_data_frame_id_col_name = 'ensembl_gene_id', removeUnmappedRows = TRUE)

### Extract String IDs significantly upregulated and significantly downregulated genes ###
ur_genes_str =  de_res_df_mapped[de_res_df_mapped$padj < 0.05 & de_res_df_mapped$log2FoldChange > 0, 'STRING_id']
dr_genes_str =  de_res_df_mapped[de_res_df_mapped$padj < 0.05 & de_res_df_mapped$log2FoldChange < 0, 'STRING_id']

### Run protein-protein interaction network clustering ###

# Perform clustering #
ur_genes_str_clusts = string_db$get_clusters(ur_genes_str)
dr_genes_str_clusts = string_db$get_clusters(dr_genes_str)

# Name the clusters #
names(ur_genes_str_clusts) = paste0('UR_Clust_', 1:length(ur_genes_str_clusts))
names(dr_genes_str_clusts) = paste0('DR_Clust_', 1:length(dr_genes_str_clusts))

### Save cluster data as text file ###

# Extract cluster data as data frame #
str_clust_df = stack(c(ur_genes_str_clusts, dr_genes_str_clusts))

# Map cluster assignments per gene #
de_res_df_mapped$Cluster = str_clust_df[match(de_res_df_mapped$STRING_id, str_clust_df$values), 'ind']

# Save cluster assignments as text file #
write.table(x = de_res_df_mapped[!(is.na(de_res_df_mapped$Cluster)), ], file = './Outputs/Analysis/Text_Files/SF_vs_GC_String_Clusters.csv', sep = ',', row.names = F)

### Plot larger clusters ###

# Extract clusters with ≥ 10 features #
clusts_to_plot = c(ur_genes_str_clusts[lengths(ur_genes_str_clusts) >= 10], dr_genes_str_clusts[lengths(dr_genes_str_clusts) >= 10])

# Plot the clusters and save as PDF files #
for(clust in names(clusts_to_plot)){pdf(file = paste0('./Outputs/Analysis/Plots', '/', clust, '_Plot.pdf'), width = 10, height = 10)
                                    string_db$plot_network(clusts_to_plot[[clust]])
                                    dev.off()}

###########################################################################################################################################
###########################################################################################################################################

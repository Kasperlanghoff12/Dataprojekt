load_GOI <- function(ctrl_dir, sample_dir, annot_gr, ctrl, sample, GOI, verbose=TRUE, verbose2=TRUE){
  
  ##########################################################
  #### (1) Indlæs gen-specifikt data
  ########################################################## 
  
  # Output the index of the gene of interest (GOI) and its name if verbose or verbose2 is TRUE
  if (verbose | verbose2){
    cat( GOI, sep='\n')
  }
  
  # Output a tab or newline character based on the value of verbose
  if (verbose){
    cat(ifelse(verbose, '\t', '\n'))
  }
  
  # Create an empty list to store data for each replicate
  rt_list = list()

  strand_options = c("+"="plus", "-"="minus")
  
  # Subsetting annotation on gene name
  gene_annot = subset(annot_gr, gene_name == GOI)
  chrom_no = as.character(seqnames(gene_annot)@values)
  strand_sign = as.character(strand(gene_annot)@values)
  chrom_no_loop = paste0("chr", chrom_no) # TEMP FIX
  
  # Coordinates
  start_coord = min(start(gene_annot))
  end_coord = max(end(gene_annot))
  
  # Find coord of next gene
  next_gene_annot = sort(subset(annot_gr, seqnames %in% c(chrom_no, chrom_no_loop) & strand == strand_sign))
  if (strand_sign == "+") {
    s = start(next_gene_annot)
    e = end(next_gene_annot)
    next_gene_coord = s[which(s > end_coord)[1]]
    prev_gene_coord = e[rev(which(e < start_coord))[1]]
  } else {
    s = end(next_gene_annot)
    e = start(next_gene_annot)
    next_gene_coord = s[rev(which(s < start_coord))[1]]
    prev_gene_coord = e[which(e > end_coord)[1]]
  }
  
  start_coord_ext =  ifelse(strand_sign == "+", min(start(gene_annot)) - usExt, next_gene_coord + 1)
  end_coord_ext = ifelse(strand_sign == "+", next_gene_coord - 1, max(end(gene_annot)) + usExt)
  
  
  # Empty DF with correct length
  indexes <- seq(start_coord_ext, end_coord_ext)
  GOI_data <- data.frame(index = indexes)
  
  # Loop for loading replicates and saving to DF
  for (fname in ctrl) {
    ctrl_fname = paste0(ctrl_dir, fname, "_", strand_options[strand_sign], ".bw")
    ctrl_bw = import(ctrl_fname, 
                     which = GenomicRanges::GRanges(seqnames = chrom_no_loop, 
                                                    ranges = IRanges::IRanges(start = start_coord_ext, end = end_coord_ext)), 
                     as = "NumericList")[[1]]
    GOI_data[[fname]] = ctrl_bw
  }
  
  # Loop for loading replicates and saving to DF
  for (fname in sample) {
    sample_fname = paste0(sample_dir, fname, "_", strand_options[strand_sign], ".bw")
    sample_bw = import(sample_fname, 
                       which = GenomicRanges::GRanges(seqnames = chrom_no_loop, 
                                                      ranges = IRanges::IRanges(start = start_coord_ext, end = end_coord_ext)), 
                       as = "NumericList")[[1]]
    GOI_data[[fname]] = sample_bw
  }
  
  # Setting the DF indices to the positions of the reads and deleting index column
  if (strand_sign == "-"){
    GOI_data = GOI_data[order(GOI_data$index, decreasing=TRUE), ]
  }
  rownames(GOI_data) = GOI_data$index
  GOI_data$index = NULL

  return (GOI_data)
}

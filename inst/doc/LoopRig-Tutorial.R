## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  out.width = "100%",
  tidy.opts = list(width.cutoff = 80),
  tidy = TRUE
)

## ----LoopsToRanges-------------------------------------------------------
library(LoopRig)

# Load example files for chromatin loop data. These are BEDPE files with no extra columns.
ovary_loops <- system.file("extdata/loops", "ovary_hg19.bedpe", package = "LoopRig", mustWork = TRUE)
pancreas_loops <- system.file("extdata/loops", "pancreas_hg19.bedpe", package = "LoopRig", mustWork = TRUE)
spleen_loops <- system.file("extdata/loops", "spleen_hg19.bedpe", package = "LoopRig", mustWork = TRUE)

# Call LoopsToRanges() on all files at once. Since there are custom columns, we indicate this by custom_cols = 0.
loops <- LoopsToRanges(ovary_loops, pancreas_loops, spleen_loops, custom_cols = 0, loop_names = c("ovary", "pancreas", "spleen"))

# View LoopRanges object for first loop dataset (ovary). 
head(loops, 1)

# We can see that each element in the LoopRanges list contains the loops from a given dataset, in this case the ovarian loops. 

# We can confirm the class of our object by using the class() function.
class(loops)

## ----ElementsToRanges----------------------------------------------------
# Load example files for genomic element data. These are BED4 files, indicating they have one extra column for metadata.
enhancers <- system.file("extdata/elements", "enhancers.bed", package = "LoopRig", mustWork = TRUE)
promoters <- system.file("extdata/elements", "promoters.bed", package = "LoopRig", mustWork = TRUE)

# Call ElementsToRanges() on all files at once. We must specify that there is one extra column by custom_cols = 1. We want to use this column as metadata for our object, so we indicate this by custom_mcols = 4, meaning that the metadata columns for our ElementRanges object will be from column 4 of the original BED files.
element_ranges <- ElementsToRanges(enhancers, promoters, element_names = c("enhancers", "promoters"), custom_cols = 1, custom_mcols = 4)

# View ElementRanges object for first element type (enhancers).  
head(element_ranges, 1)

# As we can see, we have metadata for the GRanges object that stores element information.  

# Confirm object class.
class(element_ranges)

## ----ConsensusLoops------------------------------------------------------
# Call the ConsensusLoops() function. Stringency indicates how many datasets a loop must be present in to be considered a consensus loop, and the overlap threshold (nucleotides) defines a *hit* across datasets. In this case, a loop must be present in at least 2 datasets based on overlap with another of 1 nucleotide.
consensus_loops <- ConsensusLoops(loop_ranges = loops, stringency = 2, overlap_threshold = 1)

# View consensus_loops. 
head(consensus_loops, 1)

# As we can see, we still have a LoopRanges object, but it only comprises of one GRangesList object containing the consensus loop data. We can double-check this using the length() function.
length(loops)
length(consensus_loops)

# We can further check class to ensure LoopRanges object is intact.
class(consensus_loops)

## ----DropLoops-----------------------------------------------------------
# Call DropLoops() and indicate the 'type' parameter to subset loops by loop size and use the 'size' parameter to specify lengths to keep (0-10 Kb).
consensus_loops_dropped <- DropLoops(loop_ranges = consensus_loops, type = "loop_size", size = c(100000, 200000))

# We can view the change by determining how many loops are in each object. The first index must be specified because we are still using LoopRanges list objects, and we must index an anchor to determine a count.
length(consensus_loops[[1]][[1]])
length(consensus_loops_dropped[[1]][[1]])

# Our total number of loops reduced from 68 to 41, and these 41 loops all have total loop sizes (anchor end-to-end) between 100-200 Kb.  

# Recheck class of subset object.
class(consensus_loops_dropped)

## ----LinkedElements------------------------------------------------------
# Call LinkedElements() for the promoter and enhancer ranges in element_ranges and use the default output parameters (dataframe). We must subset the ElementRanges object we created initially, as it has both the enhancer (index 1) and promoter (index 2) data in one object. 
linked_elements <- LinkedElements(loop_ranges = consensus_loops, element_ranges_x = element_ranges[[1]], element_ranges_y = element_ranges[[2]])

linked_elements

# Using our consensus set of loops, enhancer data, and promoter data, we find that the only linked elements are Enhancer 3762 and the promoter corresponding to the KSR1 gene. 

# We can do the exact same analysis, but this time output the promoters subset by this linkage instead of the dataframe indicating the links (i.e. the range corresponding to the KSR1 promoter).
linked_elements_promoters <- LinkedElements(loop_ranges = consensus_loops, element_ranges_x = element_ranges[[1]], element_ranges_y = element_ranges[[2]], range_out_y = TRUE)

linked_elements_promoters

class(linked_elements_promoters)

# As we can see, the output this time is an ElementRanges object containing the promoter data. 

## ----StackedElements-----------------------------------------------------
# Call StackedElements() for the promoter and enhancer data by subsetting the element_ranges object and using the consensus_loops object for the looping data. The default output will be a dataframe object.
stacked_elements <- StackedElements(loop_ranges = consensus_loops, element_ranges_x = element_ranges[[1]], element_ranges_y = element_ranges[[2]])

stacked_elements

# From here, we can see that we have three instances of stacked promoters and enhancers on our consensus loop anchors. 

# We can output which enhancers correspond to these stacks by using the range_out_x parameter. This will output these enhancers as an ElementRanges object.
stacked_elements_enhancers <- StackedElements(loop_ranges = consensus_loops, element_ranges_x = element_ranges[[1]], element_ranges_y = element_ranges[[2]], range_out_x = TRUE)

stacked_elements_enhancers

class(stacked_elements_enhancers)

## ----ScaffoldElemenets---------------------------------------------------
# Call ScaffoldElements() for the promoter and enhancer data using indexing and use the consensus_loops object for the loop data. Use default options so the output is a dataframe object.
scaffold_elements <- ScaffoldElements(loop_ranges = consensus_loops, element_ranges_x = element_ranges[[1]], element_ranges_y = element_ranges[[2]])

scaffold_elements

# Similar to LinkedElements() and StackedElements(), we can output either the linked enhancers or promoters as ElementRanges objects. 

# Let's output the promoters using the range_out_y parameter.
scaffold_elements_promoters <- ScaffoldElements(loop_ranges = consensus_loops, element_ranges_x = element_ranges[[1]], element_ranges_y = element_ranges[[2]], range_out_y = TRUE)

scaffold_elements_promoters

class(scaffold_elements_promoters)

## ----ExportBED, eval = FALSE---------------------------------------------
#  # Call ExportBED() for the linked_elements_promoters objects. Specify element index, file name/location, and indicate mcol=TRUE to save object metadata with the BED file.
#  
#  # NOT RUN
#  ExportBED(obj = linked_elements_promoters, index = 1, mcol = TRUE, file_name = "promoter_linked_enhancers.bed")


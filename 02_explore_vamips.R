# Updated on 2023-05-01 # Additional batches considered
# Started on 2023-04-11 by Rodrigo García-López for Prof. Carlos Arias Laboratory at the Instituto de Biotecnología of the UNAM
# R version 4.2.3 (2023-03-15) -- "Shortstop Beagle"
# Under GNU GPLv3 license
# This script contains the processing of the VAMIP table created from the ivar tables from the whole CoViGen collection as deposited in the IBT cluster and CoViGen server.
# As of 2023-05-01, this includes the main 80 batches from IMSS (spanning February, 2021 to January, 2023) as well as these extra items: CIENI and Cardio (genomes processed by [INER] from the Centro de Investigación en Enfermedades Infecciosas and INC, Instituto Nacional de Cardiología), CIAD items from Sinaloa and Nayarit mostly, and the P1000 project, a previous project with samples from 2020 (the first year of the pandemic). This latter and the CIENI set are the only sets containing any items from the first year.
# The input table has the following fields (columns) Genome Position In Out  Alt_Freq DP Ref_rev Alt_rev and each mutation is independently registered for each genome. DP is for depth, Ref is the same as Wuhan's reference (no mutation) and Alt is any alternative substitution. Rev is for reverse (how many are in the complementary strand, which are ideally not zero).

# ### Some functions ###
find_items <- function(queryStrings, dataf, colname){ # IMPORTANT. This is used to avoid partial index matching. Use it to hash a table and use rownames to match a vector. colname may be numeric
	searchVect <- dataf[,colname]
	names(searchVect) <- rownames(dataf)
	out <- searchVect[unlist(queryStrings)]
	return(out)
}
strip <- function(string){ # Takes a string object and strips rare characters
	string <- iconv(string,from="UTF-8",to="ASCII//TRANSLIT") # Remove accents and rare chars
# 	string <- tolower(string) # change to lower
	string <- gsub("(?<=[\\s])\\s*|^\\s+|\\s+$", "", string, perl=TRUE) # Remove multiple spaces
	return(string)
}

# ### Preprocess data ###
setwd("/home/rod/Documents/01_Projects/SARS/VAMIP/") # This is a local folder containing the input table
df <- read.table("06_VAMIP_tables/All_mut_NoSingl_filtered.tsv", header=F, sep ='\t', stringsAsFactors = FALSE,fill=T) # input table (mutation x ivar specs)
dim(df)
# [1] 2717614       8
df[!complete.cases(df),] # Some line may have NAs, these should be ignored (this is not ideal but it is like adding an N, basically)
#                            V1    V2 V3 V4 V5       V6 V7 V8
# 762142 L030_202101206551_S199 21973  T  + NA 0.333333  3  3
# 762143                      0    NA       NA       NA NA NA
df <- df[complete.cases(df),] # remove these faulty items
dim(df)
# [1] 2717612       8
names(df) <- c("Genome", "Position", "In", "Out", "Alt_Freq", "DP", "Ref_rev", "Alt_rev") # Rename ivar columns
# We will now remove the _SXXX part of the name (this was kept until now to avoid repeated items but there should be removed at this point if any are still present)
temp_names <- unique(df[,"Genome"]) # Now, derreplicate genome names as they are (with the library number)
length(temp_names)
# [1] 30517 # These are the total items that are considered for the analysis, several will be removed and only those matching metadata may be used later
dup_for_removal <- names(which(table(temp_names)>1)) # Flag repeated items if analysis
if(length(dup_for_removal)!=0){ # in case there are any repeated items
	remove_list <- temp_names[as.vector(unlist(sapply(dup_for_removal,function(x){grep(x,temp_names)})))] # Go through the list and get all repeated items
	row_for_removal <- as.vector(unlist(sapply(remove_list,function(x){grep(x,df[,"Genome"])}))) # Get their row numbers
	dim(df[row_for_removal,]) # Total items to be removed (all mutations per sample)
	write.table(df[row_for_removal,], "07_metadata/RepeatedGenomes_vamips_removed.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)
	df <- df[-row_for_removal,]
	dim(df)
}

# Now that there are no genomes that have the same name (not repeated but faulty called the same), we can replace the _SXX suffix
df[,"Position"] <- sprintf("%05d", df[,"Position"]) # add 0s to fix position sorting
df[,"DP_Alt"] <- round(df[,"DP"]*df[,"Alt_Freq"]) # add the depth of the alternative vars
df[,"DP_Ref"] <- df[,"DP"]-df[,"DP_Alt"] # and the depth of reference (original) vars
df[,"Mutation"] <- paste0(df[,"Position"],":",df[,"In"],"->",df[,"Out"]) # This will hold mutations in a new nomenclature such as 00021:C->T with position:ref->alt
df <- df[order(df[,"Position"]),] # now, sort observations by position

# LOAD METADATA
meta <- read.table("07_metadata/2023-05-01_All_Metadata_Fastq_CoViGen_variantInfo.tsv", header=T, sep ='\t',stringsAsFactors = FALSE, check.names=F, fill=T, quote="") # Load metadata with id folio for xref. This was the result of matching all items with their GISAID-deposited metadata. However, the current (2023-04-28) version has no info for several fields, including CT, and such. I'm currently using a more complete metadata version which has extra info on variants and has some issues already fixed (such as having repeated items
dim(meta)
# [1] 29717    59
# Fix items in the main df that do not have the right names
temp <- grep("Cardio", meta[,"Deprecated_name"]) # First, Cardio genomes
dict <- meta[temp,"File basename"]
names(dict) <- meta[temp,"Deprecated_name"]
temp <- dict[df[,"Genome"]] # match the items
df[!is.na(temp),"Genome"] <- temp[!is.na(temp)] # Only replace matches
# Similarly, we have to fix some items having UIBMZ instead of UIBZ as part of their names
df[,"Genome"] <- sub("UIBMZ","UIBZ",df[,"Genome"])

meta[meta[,"VOCs"]=="B.1.1.519","VOCs"] <- "519" # We'll simplify variant name for 519 and 222
meta[meta[,"VOCs"]=="B.1.1.222","VOCs"] <- "222"
rownames(meta) <- meta[,"File basename"] # and use the sample name as row identifier
meta[,"Rename"] <- paste0(meta[,"Fecha de recoleccion"],"|", meta[,"Accession ID"],"|",meta[,"New_Lineage"],"(",meta[,"Clado Nexstrain"],"_",meta[,"VOCs"],")","|",meta[,"Edad"],"|",meta[,"Genero"],"|",meta[,"Estado"]) # append a new column with a date|name|PANGO(Clade)|age|gender|State
temp <- meta[,"Rename"]; names(temp) <- rownames(meta) # IMPORTANT: create an auxiliary vector to avoid partial matches (df has partial string matching)
# Now, match both sets

df[,"Names"] <- temp[df[,"Genome"]] # add names now
bad <- df[is.na(df[,"Names"]),] # Get those that have no metadata available
nrow(bad) # Mutations with no matching genome
# [1] 83382
nonmatched_Genomes <- unique(bad[,"Genome"])
length(nonmatched_Genomes) # Non-matching genomes (these have no metadata)
# [1] 842
write.table(bad, "07_metadata/Mutations_withNO_metadata.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA) # The complete collection of mutations from those items failing
write.table(sort(nonmatched_Genomes), "07_metadata/Genomes_withNO_metadata.tsv", sep="\t", quote=FALSE, row.names=F, col.names=F) # Print a list of items failing
write.table(table(df[!is.na(df[,"Names"]),"Genome"]),"07_metadata/test_genomes_inVAMIP_wMetadata.tsv", sep="\t", quote=FALSE, row.names=F, col.names=F) # Test for good items.
df <- df[!is.na(df[,"Names"]),] # Now, filter those that could not be identified
dim(df)
# [1] 2634230      12
meta <- meta[names(table(df[,"Genome"])),] # IMPORTANT: Create a new metadata table, keeping only those items that were actually there.
dim(meta)
# [1] 29675    60
write.table(meta, "07_metadata/2023-05-01_Metadata_all_in_VAMIP_set_that_matched.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)
dim(df)
# [1] 2634230      12

# ### Filters ###
df[grep("\\+",df[,"Out"]),"Alt_rev"] <- 1 # This offset is a small fix to avoid loosing all indels (since they always have 0 in the Alt_rev column (only shown in the Fwd strand)
df[grep("\\-",df[,"Out"]),"Alt_rev"] <- 1
# temp <- df[,"Ref_rev"]/df[,"DP_Ref"]; temp[is.na(temp)] <- 0
temp <- df[,"Alt_rev"]/df[,"DP_Alt"]; temp[is.na(temp)] <- 0 # Calculate the frequency of items in the reverse strand out of all alternative observations (those having a mutation). These should not be only in one strand to be valid.
temp <- as.logical((temp > 0) * (temp < 1)) # Create a mask to detect those having any observations in both strands
sum(!temp) # How many have either all observations on the F or R strand?
# [1] 289174
sum(df[,"Alt_Freq"] < 0.05) # How many are not in at least 0.05 of observations?
# [1] 475182
df[!temp,"Alt_Freq"] <- 0.000001 # Those with no reads in both strands are marked as 0 frequency + an offset
df[df[,"Alt_Freq"] < 0.05,"Alt_Freq"] <- 0.000001 # Those will have there values set to 0 + an offset
sum(df[,"DP"] < 20) # How many have less than 20 reads for depth
# [1] 193644
df[df[,"DP"] < 20,"Alt_Freq"] <- 0.000001 # Same for those with low DP
sum(df[,"Alt_Freq"]==0.000001) # How many were marked as bad?
# [1] 716951
sum(df[,"Alt_Freq"]>0.000001) # How many passed all filters
# [1] 1917279
# ### Reduce table size ###
# UPDATE 2023-05-02: Due to the large size of the table, I changed some filters to this section, as to avoid unncessary memory usage
df <- df[df[,"Alt_Freq"]>=0.05,] # First, remove items not passing filters from the full table
dim(df)
# [1] 1917279      12
save.image("temp.Rdata")
# Then, count how many mutations are observed in less than 5 genomes (this will be included to track
test <- table(df[,"Mutation"])<5
length(test) # How many different mutations are observed
# [1] 50911
sum(test) # How many mutations are only seen in 4 or less genomes?
# [1] 35938
# We will only ignore only the IPMAVs in those with less than this cutoff as the ones that are already fixated are actually confirmed
sum((test[df[,"Mutation"]])*(df[,"Alt_Freq"]<0.5)) # How many items would be removed?
# [1] 50526
sum(!(test[df[,"Mutation"]])*(df[,"Alt_Freq"]<0.5)) # How many would pass?
# [1] 1866753
df <- df[!as.logical((test[df[,"Mutation"]])*(df[,"Alt_Freq"]<0.5)),] # create a mask for those items in rare (<5 observations) AND not having fixated mutations (freq<0.5)
dim(df)
# [1] 1866753      12

save.image("first_full_filtered_Alt_whole_table.Rdata")
write.table(df, "Alt_Whole_table-Filtered.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)

# ### Extract tables ###
# First, frequency
load("first_full_filtered_Alt_whole_table.Rdata")
genomes <- names(table(df[,"Names"])) # Extract the names
# genomes <- genomes[order(sub(".*\\|","",genomes))] # Sort them by lineage
genomes <- sort(genomes) # sort by date
mutations <- names(table(df[,"Mutation"])) # and extract the mutations
out <- as.data.frame(matrix("", ncol=length(genomes), nrow=length(mutations))); colnames(out) <- genomes; rownames(out) <- mutations # Create a void matrix container
save.image("first_full_alt_variant_freq_table.Rdata")
# load("first_full_alt_variant_freq_table.Rdata")
for(i in 1:nrow(df)){ # and survey those that are present (the rest will be NAs)
	out[df[i,"Mutation"],df[i,"Names"]] <- df[i,"Alt_Freq"]
}
dim(out)
# [1] 80087 29661
save.image("first_full_alt_variant_freq_table.Rdata")
# load("first_full_alt_variant_freq_table.Rdata")
dim(out) # How many remain?
# [1]  7340 17927
write.table(out, "Alt_variant_Freq.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)


# Next, depth
load("first_full_filtered_Alt_whole_table.Rdata")
genomes2 <- names(table(df[,"Names"])) # Extract the names
# genomes2 <- genomes2[order(sub(".*\\|","",genomes2))] # Sort them by lineage
genomes2 <- sort(genomes2) # Sort them by date
mutations2 <- names(table(df[,"Mutation"])) # and extract the mutations2
out2 <- as.data.frame(matrix("", ncol=length(genomes2), nrow=length(mutations2))); colnames(out2) <- genomes2; rownames(out2) <- mutations2 # Create a void matrix container
for(i in 1:nrow(df)){ # and survey those that are present (the rest will be NAs)
	out2[df[i,"Mutation"],df[i,"Names"]] <- df[i,"DP"]
}
dim(out2)
# [1] 46877 17927
save.image("first_full_depth_table.Rdata")
write.table(out2, "Alt_variant_TotalDepth.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)

# And now, with alt allelic depth
load("first_full_filtered_Alt_whole_table.Rdata")
genomes2 <- names(table(df[,"Names"])) # Extract the names based on Selene's table
genomes2 <- genomes2[order(sub(".*\\|","",genomes2))] # Sort them by lineage
mutations2 <- names(table(df[,"Mutation"])) # and extract the mutations2
out2 <- as.data.frame(matrix("", ncol=length(genomes2), nrow=length(mutations2))); colnames(out2) <- genomes2; rownames(out2) <- mutations2 # Create a void matrix container
for(i in 1:nrow(df)){ # and survey those that are present (the rest will be NAs)
	out2[df[i,"Mutation"],df[i,"Names"]] <- df[i,"DP_Alt"]
}
dim(out2)
# [1] 46877 17927
save.image("first_full_alt_variant_AltDepth_table.Rdata")
write.table(out2, "Alt_variant_AltAllelicDepth.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)

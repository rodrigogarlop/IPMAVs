# Updated 2023-05-03 Calculate main mutations per genome
# Updated 2022-10-06 Aded recombinants
# Started 2022-09-29 by Rodrigo Garcia-Lopez for GAL, iBT, UNAM: We will get percentage of VAMIPS vs non vamips
# by Rodrigo García-López for Carlos Arias's Virology Group at IBt, UNAM, Cuernavaca, Mexico as part of the CoViGen-Mex SARS-CoV-2 survillance in Mexico.
# This script was created to analyze all mutations (IPMAVs or fixated mutations per sample instead of separately as script 03_Analyze_vamip_contingency_Bar_n_boxplots). It is used to construct haplotypes for variants and VOCs in the set
# R version 4.2.3 (2023-03-15) -- "Shortstop Beagle"
# Under GNU GPLv3 license
# Disclamer: This script was written for use with specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.
# ### LOAD LIBRARIES & INPUT ###
library("pheatmap")
library("stringr") # Methods have been written for replacing this, to avoid version issues
df <- read.table("Alt_variant_Freq.tsv", sep="\t",header=T, skip=0, comment.char='',quote="",fill=F,row.names=1, check.names=FALSE) # Load input table
# Each column was named as follows: "2023-01-27|EPI_ISL_16983502|XBB.1.5(20B_Others)|065|F|Zacatecas" and they have been sorted before
df<- df[,order(colnames(df))] # In case they were not sorted, we do it here
sum(!is.na(df))
# [1] 1866675
dim(df)
# [1] 25053 29661
save.image("Genomes-VAMIPsVsGenomes_chkpt1.Rdata")

# UPDATE 2023-05-03: We want to get which are the main mutations per genome
# ### Functions ###
strip <- function(string){ # Takes a string object and strips rare characters
	string <- iconv(string,from="UTF-8",to="ASCII//TRANSLIT") # Remove accents and rare chars
# 	string <- tolower(string) # change to lower
	string <- gsub("(?<=[\\s])\\s*|^\\s+|\\s+$", "", string, perl=TRUE) # Remove multiple spaces
	string <- gsub("\\+","",string) # remove + characters
	return(string)
}
find_items <- function(queryStrings, dataf, colname){ # IMPORTANT. This is used to avoid partial index matching. Use it to hash a table and use rownames to match a vector. colname may be numeric
	searchVect <- dataf[,colname]
	names(searchVect) <- rownames(dataf)
	out <- searchVect[unlist(queryStrings)]
	return(out)
}
group_map <- function(vect) { # Get a mutation table (columns = genomes; rows = mutations), and a vector of the same rowsize. Then use the vector to define how groups should be splitted and build a list of df with the table observations per group (a map of positions)
	group_names <- levels(as.factor(vect))
	group_positions <- sapply(group_names,function(x){grep(paste0("^",strip(x),"$"),strip(vect))})
	return(group_positions)
}
group_prevalence <- function(gmap,inmat) { # Get the group map from the corresponding function and get the prevalence per row item. This is calculated based on the total NAs.
	nam <- names(gmap) # Save names for indices
	out_raw <- matrix(NA, nrow=nrow(inmat), ncol=length(nam)); rownames(out_raw) <- rownames(inmat); colnames(out_raw) <- nam # Create a void matrix container
	out_rel <- out_raw # make them two
	for(i in nam){ # now, for each table
# 		print(i)
		single_grp <- inmat[,gmap[[i]]] # First, subset the input table using the group map for a single group
		out_raw[,i] <- apply(single_grp,1,function(x){length(x[!is.na(x)])})  # and calculate the total items
		out_rel[,i] <- out_raw[,i]*100/ncol(single_grp)
	}
	out <- list(out_raw,out_rel)
	names(out) <- c("raw","rel")
	return(out)
}
variant_completeness <- function(freqs_vect,pos_vect) { # extract only the items in pos_vect positions of vector freqs_vect, then calculate which have not NAs in those positions, the output is a single %
	test_vect <- freqs_vect[pos_vect]
	MutObs <- sum(!is.na(test_vect))
	MutExp <- length(test_vect)
	ObsPerc <- round(MutObs/MutExp*100,2)
	out <- list(MutObs, MutExp, ObsPerc)
	names(out) <- c("MutObs", "MutExp", "ObsPerc")
	return(out)
}
# ### MAIN ###


for(mat in names(multiMatrix)){
	out_mat <- t(apply(df, 2, function(genome) {apply(multiMatrix[[mat]], 2, function(var){variant_completeness(genome,var)$ObsPerc})})) # test all variants on each genome and create a table bearing the completeness
	write.table(out_mat, paste0("10_Analysis_per_genome/genome_variant_completeness/genome_variant_completeness_", mat, "_confidence.tsv"), sep="\t", quote=FALSE, row.names=T, col.names=NA)
}

sum(!is.na(df)) # we test the total number of items
# [1] 1866675
dim(df)
# [1]  25053 29661
# Now, load the metadata to determine lineages
# PREPARE METADATA
meta <- read.table("07_metadata/2023-05-01_Metadata_all_in_VAMIP_set_that_matched-derep_sortByTime.tsv", header=T, sep ='\t',stringsAsFactors = FALSE, check.names=F, fill=T, quote="") # Load metadata with id folio for xref. IDs should be unique
rownames(meta) <- meta[,"Accession ID"]
dim(meta)
# [1] 29661    60
meta[,"VOCs_three"] <- sub("XBB\\..*","XBB",meta[,"VOCs_three"]) # Update 2023-08-05: Not enough variation was seen among XBB variants that justify analyzing them separately.
write.table(meta, "07_metadata/2023-05-01_Metadata_all_in_VAMIP_set_that_matched-derep_sortByTime_xbbFix.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)
# Process the name to extract basic info and the Folio ID
# names are processed into a table: Date|Folio ID|Pango(clade)|age|gender|State
# Example: "2021-03-04|L004_202101083013|B.1.1.222(20B_Others)|046|M|Puebla"
# minimetatab <- t(as.data.frame(str_split(colnames(df), "\\|")))
minimetatab <- t(as.data.frame(sapply(colnames(df),function(x){strsplit(x,"\\|")[[1]]}))) #UPDATE 2023-05-02: After updating to R.4.2.3, the stringr, the str_split function stopped working so I changed it back to R base
rownames(minimetatab) <- NULL # row names are removed
# This may be a larger table but only matching items will be kept
meta <- meta[unique(minimetatab[,2]),]
dim(meta) # This should not change if all were found
# [1] 29661    60

all_lineages <- find_items(minimetatab[,2],meta, "Linaje Pangolin") # Load the lineage metadata (this is drawn from the "meta" object since the one in the colname may be oudated
lineage_map <- group_map(all_lineages) # Create a map of which columns belong to which lineages
lineage_totals <- sort(unlist(lapply(lineage_map, length)), decreasing=TRUE) # Now, count them
n=20 # Define the least number of genomes for each variant to be considered for the analysis
pass_yes <- names(which(lineage_totals>=n)) # Get which items have at least n genomes
pass_not <- names(which(lineage_totals<n)) # And which ones don't
lineage_map_filt <- lineage_map[pass_yes] # Start by creating a list with those items that pass the filter
temp <- unlist(lineage_map[pass_not[grep("B\\.1\\.631",pass_not)]]);names(temp) <- NULL; lineage_map_filt[["B.1.631"]] <- temp # and append additional items for the variant from which XBB originates (along with B.1.634)
compare_variant_fil <- group_prevalence(lineage_map_filt,df) # Now, get the prevalence of each mutation per variant
write.table(compare_variant_fil[[1]], "10_Analysis_per_genome/genome_variant_completeness/2023-05-03_Main_mutations_per_variant-filtered.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)
compare_variant_fil_rel <- round(compare_variant_fil[[2]],2)/100
write.table(compare_variant_fil_rel, "10_Analysis_per_genome/genome_variant_completeness/2023-05-03_Main_mutations_per_variant-filtered-rel.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)
# We'll be using a frequency cutoff of 0.6 for all items to consider them key mutations

multiMatrix <- list() # create a new list, each item will be a matrix with varying numbers of representative mutations per variant (filtered at different permissive levels [percentages] of the genomes that have them per variant)
for (i in rev(seq(.50,.95,.05))) { # We will just be using those in most (>50% samples)
	print(i)
	multiMatrix[[paste0("m",i)]] <- compare_variant_fil_rel>=i # each output in the list will be a logical matrix with various levels of permissiveness
};rm(i)

# We will now do the variant analysis (main mutations) in a more general way, using for example Delta, Omicron, etc.
# UPDATE 2023-05-03: This now includes recombinants
meta[meta[,"VOCs_three"]=="B.1.1.222","VOCs_three"] <- "222" # simplify 222's name
meta[meta[,"VOCs_three"]=="B.1.1.519","VOCs_three"] <- "519" # and 519's
all_VOCs <- find_items(minimetatab[,2],meta, "VOCs_three") # Load the VOC metadata (this is drawn from the "meta" object since the one in the colname may be oudated
VOC_map <- group_map(all_VOCs) # Create a map of which columns belong to which VOCs
VOC_totals <- sort(unlist(lapply(VOC_map, length)), decreasing=TRUE) # Now, count them
n=16 # Define the least number of genomes for each variant to be considered for the analysis
pass_yes <- names(which(VOC_totals>=n)) # Get which items have at least n genomes
pass_not <- names(which(VOC_totals<n)) # And which ones don't
VOC_map_filt <- VOC_map[pass_yes] # Start by creating a list with those items that pass the filter
compare_VOC_fil <- group_prevalence(VOC_map_filt,df) # Now, get the prevalence of each mutation per variant
write.table(compare_VOC_fil[[1]], "10_Analysis_per_genome/genome_VOC_completeness/2023-05-03_Main_mutations_per_VOC-filtered.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)
compare_VOC_fil <- compare_VOC_fil[,c("222", "519", "Alpha", "Gamma", "Delta", "Omicron_BA.1.x", "Omicron_BA.2.x", "Omicron_BA.4.x", "Omicron_BA.5.x", "XB", "XBB", "XAS", "MultLin", "Others")]
compare_VOC_fil_rel <- round(compare_VOC_fil[[2]],2)/100
compare_VOC_fil_rel <- compare_VOC_fil_rel[,c("222", "519", "Alpha", "Gamma", "Delta", "Omicron_BA.1.x", "Omicron_BA.2.x", "Omicron_BA.4.x", "Omicron_BA.5.x", "XB", "XBB", "XAS", "MultLin", "Others")]
write.table(compare_VOC_fil_rel, "10_Analysis_per_genome/genome_VOC_completeness/2023-05-03_Main_mutations_per_VOC-filtered-rel.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)

multiMatrix <- list() # create a new list, each item will be a matrix with varying numbers of representative mutations per variant (filtered at different permissive levels [percentages] of the genomes that have them per variant)
for (i in rev(seq(.50,.95,.05))) { # We will just be using those in most (>50% samples)
	print(i)
	multiMatrix[[paste0("m",i)]] <- compare_VOC_fil_rel>=i # each output in the list will be a logical matrix with various levels of permissiveness
};rm(i)

for(mat in names(multiMatrix)){
	out_mat <- t(apply(df, 2, function(genome) {apply(multiMatrix[[mat]], 2, function(var){variant_completeness(genome,var)$ObsPerc})})) # test all variants on each genome and create a table bearing the completeness
	write.table(out_mat, paste0("10_Analysis_per_genome/genome_VOC_completeness/genome_VOC_completeness_", mat, "_confidence.tsv"), sep="\t", quote=FALSE, row.names=T, col.names=NA)
}

#UPDATE 2022-10-06: We now want to evaluate genome by genome to determine the % of VAMIPs
load("Genomes-VAMIPsVsGenomes_chkpt1.Rdata")
dir.create("10_Analysis_per_genome/genomePos_vs_SNPnVAMIP_contents", showWarnings=FALSE)
df[is.na(df)] <- 0 # We need to remove NAs before anything else
VAMIPs <- SNPs <- df
SNPs[SNPs<0.5] <- 0 # First, get two matrices, one with only SNPs
VAMIPs[VAMIPs>=0.5] <- 0 # and the other one with only VAMIPs
# Now, we need to collate each table where mutations are summed by position. We only want to keep substitutions and indels separated
mut_names <- rownames(df) # save the mutations' names
rm(df)
save.image("Genomes-VAMIPsVsGenomes_chkpt2a.Rdata") # We need to flush the memory, which can be done by saving and reloading the image
# load("Genomes-VAMIPsVsGenomes_chkpt2.Rdata")
fill_missing_items <- function(mat,vector){ # Gets a table and a vector with all rows that should be present, outputs an expanded row collection with missing dates for complete calendar in that range. rownames should have date format as %Y-%m-%d
	xtable <- as.data.frame(matrix(0,nrow=length(vector),ncol=ncol(mat)), stringsAsFactors = FALSE) # Create empty vessel for output
	rownames(xtable) <- vector # use the vector as rownames
	colnames(xtable) <- colnames(mat) # and inherit the names
	invisible(sapply(rownames(mat),function(x) {xtable[x,] <<- mat[x,]}))	# append the original values in the corresponding places (write to higher env variable
	return(xtable)
}
strip <- function(string){ # Takes a string object and strips rare characters
	string <- iconv(string,from="UTF-8",to="ASCII//TRANSLIT") # Remove accents and rare chars
# 	string <- tolower(string) # change to lower
	string <- gsub("(?<=[\\s])\\s*|^\\s+|\\s+$", "", string, perl=TRUE) # Remove multiple spaces
	string <- gsub("\\+","",string) # remove + characters
	string <- gsub("\\|", "_-_",string)
	string <- gsub(" ", "_",string)
	string <- gsub("\\(", "_",string)
	string <- gsub("\\)", "_",string)
	string <- gsub("\\/", "_",string)
	string <- gsub("__", "_",string)
	return(string)
}
all_mut <- as.data.frame(cbind("id"=1:length(mut_names), "Pos"=sub(":.*","",mut_names), "MutType"=rep("_Sub", length(mut_names))));rownames(all_mut) <- mut_names # we start with the full collection of mutation names, and the type of mutations, by default, these are marked by a _Sub suffix
# now, fill the indels accordingly to identify the type of mutation
all_mut[grep(">\\+",mut_names),3] <- "-Ins"
all_mut[grep(">-",mut_names),3] <- "-Del"
SNPs <- rowsum(SNPs, group=paste0(all_mut[,2],all_mut[,3])) # Collate by the postition and type of mutation # This will substiute the SNPs table
dim(SNPs)
# [1]  19189 29661
VAMIPs <- rowsum(VAMIPs, group=paste0(all_mut[,2],all_mut[,3])) # Collate by the postition and type of mutation. This will substitute the VAMIP table
dim(VAMIPs)
# [1]  19189 29661
SNPs[SNPs>1] <- 1; VAMIPs[VAMIPs>1] <- 1 # Due to rounding, some items may have values slightly larger than 1 (<1e-5), thus, we will truncate those to 1s

id_Sub <- grep("Sub",rownames(SNPs)) # Locate where each type of mutation is seen
id_Del <- grep("Del",rownames(SNPs))
id_Ins <- grep("Ins",rownames(SNPs))
length(c(id_Sub,id_Del,id_Ins)) # This checkpoint should be the same as the number of rows
# [1] 19189

save.image("Genomes-VAMIPsVsGenomes_chkpt2.Rdata") # Update the image
# load("Genomes-VAMIPsVsGenomes_chkpt2.Rdata") # Reload to flush memory
all_pos <- sprintf("%05d", 1:29903) # Add a vector for all available positions
# for(i in 1:ncol(SNPs)){ # For each genome
# for(i in 2270:15000){ # For each genome ### SUSPENDED 3048 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# for(i in 15207:ncol(SNPs)){ # For each genome ### SUSPENDED 15958
# for(i in 1:2){ # For each genome
mut_position_per_genome <- function(numbers){
	for(i in numbers){
		print(i)
		name <- colnames(SNPs)[i] # This holds the name of the genome
		all_mut <- cbind("SNPs"=SNPs[,i], "VAMIPs"=VAMIPs[,i]); rownames(all_mut) <- rownames(SNPs) # this holds the SNPs and VAMIPs for the genome
		sub <- all_mut[id_Sub,]; rownames(sub) <- sub("_Sub","",rownames(sub)) # The table is split in 3 subtables
		del <- all_mut[id_Del,]; rownames(del) <- sub("-Del","",rownames(del))
		ins <- all_mut[id_Ins,]; rownames(ins) <- sub("-Ins","",rownames(ins))
		sub <- fill_missing_items(sub,all_pos); colnames(sub) <- paste("Sub",colnames(sub)) # create a full table containing 0s for missing items, then rename the columns
		del <- fill_missing_items(del,all_pos); colnames(del) <- paste("Del",colnames(del))
		ins <- fill_missing_items(ins,all_pos); colnames(ins) <- paste("Ins",colnames(ins))
		data <- cbind(sub, del, ins)
#	 write.table(data[rowSums(data)>0,], paste0("genomePos_vs_SNPnVAMIP_contents/", str_split(name, "\\|")[[1]][2], ".tsv"), sep="\t", quote=FALSE, row.names=T, col.names=NA)
		write.table(data[rowSums(data)>0,], paste0("10_Analysis_per_genome/genomePos_vs_SNPnVAMIP_contents/", gsub("/", "_", gsub("\\|", "_-_",name)), ".tsv"), sep="\t", quote=FALSE, row.names=T, col.names=NA)
		pdf(paste0("10_Analysis_per_genome/genomePos_vs_SNPnVAMIP_contents/", gsub("/", "_", gsub("\\|", "_-_",name)), ".pdf"),width=18, height=10)
		cols <- c("turquoise","blue","pink","violetred","chartreuse3","forestgreen")
		ltys <- c(1,1,2,2,4,4)
		par(oma = c(1.5, 0, 1, 6)) # This is just a creative fix to plot the legend outside
		par(mar=c(5.1, 5.1, 4.1, 2.1))
		matplot(las=1, data[,ncol(data):1],type='l', col=rev(cols),lty=rev(ltys),lwd=, ylim=c(0,1), xlim=c(1,29903), ylab="Mutation frequency per position", main=name, xaxt='n', yaxt='n',frame.plot=FALSE, xlab="Genomic Position (Nt)", cex.lab=1.5, cex.main=1.5)
		axis(1, at=seq(1,31000,2500), labels=seq(0,31000,2500),cex=1.5)
		axis(2, las=1, at=seq(0,1,.1),cex.axis=1.2)
		Genome <- diff(c(1, 266, 13468, 21556, 21563, 25385, 25393, 26221, 26245, 26473, 26523, 27192, 27202, 27388, 27394, 27756, 27888, 27894, 28260, 28274, 29534, 29558, 29675, 29903)) # This vector holds each interval (all ORFs). Endings were adjusted with +1 for calculations
		cols_bar = c("01 5-UTR"="cornflowerblue", "02 ORF1a"="coral1","03 ORF1b"="turquoise3", "Sep1"=NA, "04 S"="chartreuse2", "Sep2"=NA, "05 ORF3a"="purple", "Sep3"=NA, "06 E"="firebrick", "Sep4"=NA, "07 M"="gold2", "Sep5"=NA, "08 ORF6"="hotpink", "Sep6"=NA, "09 ORF7a"="forestgreen", "10 ORF7b"="darkblue", "Sep7"=NA, "11 ORF8"="darkslategray", "Sep8"=NA, "12 N"="darkorange2", "Sep9"=NA, "14 ORF10"="brown2", "15 3-UTR"="mediumpurple")
		sub_cols <- cols_bar[!is.na(cols_bar)]
		par(oma = c(6, 0, 43, 6)) # Adjust for plotting the genes
		barplot(las=1,cbind(Genome), horiz = TRUE, beside = FALSE, col=cols_bar, border=NA, xaxt='n',add=T)
		par(fig = c(0, 1, 0, 1), oma = c(1.5, 0, 1, 1.5), mar = c(8.1, 0, 7.1, 0), new = TRUE)
		plot(0,0,type = "n", bty = "n", xaxt = "n", yaxt = "n",xlab="", ylab="")
		legend("topright",legend=c("Sub Fix","Sub VAMIP","Del Fix","Del VAMIP","Ins Fix","Ins VAMIP"), lty=ltys, col=cols, title="Lines", lwd=3, bg="white")
		legend("right",legend=names(sub_cols), pch=15, col=sub_cols, pt.cex=1.5, title="Genomic positions", bg="white")
		dev.off()
	}
}
save.image("Genomes-VAMIPsVsGenomes_chkpt2.Rdata") # Update the image
load("Genomes-VAMIPsVsGenomes_chkpt2.Rdata") # Reload to flush memory

# UPDATE 2022-10-06: We now want to evaluate the whole set of mutations
load("Genomes-VAMIPsVsGenomes_chkpt1.Rdata")
meta <- read.table("07_metadata/2023-05-01_Metadata_all_in_VAMIP_set_that_matched-derep_sortByTime_xbbFix.tsv", header=T, sep ='\t',stringsAsFactors = FALSE, check.names=F, fill=T, quote="") # Load metadata with acc num for xref. IDs should be unique
dict <- meta[,"New_Lineage"]; names(dict) <- meta[,"Linaje Pangolin"] # We'll use the full name for this analysis, so we need a dictionary to change them correctly
nt2AA <- read.table("07_metadata/MutNtxAA_complete.tsv", sep="\t",header=T, skip=0, comment.char='',quote="",fill=F,row.names=1, check.names=FALSE) # Load nt to AA name table
nt2AA[nt2AA[,"MutAA"]=="NC","MutAA"] <- paste0(nt2AA[nt2AA[,"MutAA"]=="NC","Gene"],":NC") # add a category to flag UTRs
mut_VOC <- read.table("10_Analysis_per_genome/genome_VOC_completeness/2023-05-03_Main_mutations_per_VOC-filtered-rel.tsv", sep="\t",header=T, skip=0, comment.char='',quote="",fill=F,row.names=1, check.names=FALSE) # Load input table
# mut_VOC <- mut_VOC[,c("222", "519", "Alpha", "Gamma", "Delta", "Omicron_BA.1.x", "Omicron_BA.2.x", "Omicron_BA.4.x", "Omicron_BA.5.x", "XB", "XBB", "XAS", "MultLin", "Others")] # This fix was required for older an older version
# write.table(mut_VOC, "10_Analysis_per_genome/genome_VOC_completeness/2023-05-03_Main_mutations_per_VOC-filtered-rel.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)
mut_Lin <- read.table("10_Analysis_per_genome/genome_variant_completeness/2023-05-03_Main_mutations_per_variant-filtered-rel.tsv", sep="\t",header=T, skip=0, comment.char='',quote="",fill=F,row.names=1, check.names=FALSE) # Load input table
colnames(mut_Lin) <- dict[colnames(mut_Lin)] # Replace abbreviated names with full names
mut_Lin <- mut_Lin[,sort(colnames(mut_Lin))] # Sort variants by name
mut_Lin <- mut_Lin[,-grep("A\\.2\\.5$",colnames(mut_Lin))] # Remove lineage A.2.5
mut_Lin <- mut_Lin[,-grep("B\\.1$",colnames(mut_Lin))] # Remove lineage B.1
mut_Lin <- mut_Lin[,-grep("B\\.1\\.1$",colnames(mut_Lin))] # and B.1.1
rownames(mut_VOC) <- sub("\\|.*","",rownames(mut_VOC))
rownames(mut_VOC) <- paste(rownames(mut_VOC),nt2AA[rownames(mut_VOC),"MutAA"], sep='|') # and append this to the mutation names
rownames(mut_Lin) <- sub("\\|.*","",rownames(mut_Lin))
rownames(mut_Lin) <- paste(rownames(mut_Lin),nt2AA[rownames(mut_Lin),"MutAA"], sep='|')

main_mut_map_VOC <- apply(mut_VOC,2,function(x){out=x[x>=0.6];out[order(names(out))]}) # Extract a map of mutations per variant (at least 60% of all genomes from each variant should have them). This holds how common that mutation is in the variant. The last part just orders by position. This is set at 60 because they are broader groups
main_mut_map_Lin <- apply(mut_Lin,2,function(x){out=x[x>=0.7];out[order(names(out))]}) # These are set at 70 as they are less variable groups
df[is.na(df)] <- 0 # We need to remove NAs before any subsetting
rownames(df) <- sub("\\|.*","",rownames(df)) # Remove tags for delta and omicron as we are expanding the VOCs we are observing (this is a trailing fix for older versions)
rownames(df) <- paste(rownames(df),nt2AA[rownames(df),"MutAA"], sep='|') # and append this to the mutation names
dfVOC <- df[sort(unique(unlist(lapply(main_mut_map_VOC,names)))),] # For processing, we will create two subsets of the main mutation table, one for VOCs (163 unique items for any VOCs at 70%)
# dim(dfVOC)
# [1]   235 17926
dfLin <- df[sort(unique(unlist(lapply(main_mut_map_Lin,names)))),] # and one for all variants with >=20 genomes in the set (unique mutations seen in at least 70% in any variant)
# dim(dfLin)
# [1]   587 29661
rm(df)
# save.image("Genomes-VAMIPsVsGenomes_chkpt3.Rdata")
# load("Genomes-VAMIPsVsGenomes_chkpt3.Rdata")

strip <- function(string){ # Takes a string object and strips rare characters
	string <- iconv(string,from="UTF-8",to="ASCII//TRANSLIT") # Remove accents and rare chars
# 	string <- tolower(string) # change to lower
	string <- gsub("(?<=[\\s])\\s*|^\\s+|\\s+$", "", string, perl=TRUE) # Remove multiple spaces
	string <- gsub("\\+","",string) # remove + characters
	string <- gsub("\\|", "_-_",string)
	string <- gsub(" ", "_",string)
	string <- gsub("\\(", "_",string)
	string <- gsub("\\)", "_",string)
	string <- gsub("\\/", "_",string)
	string <- gsub("__", "_",string)
	return(string)
}
FreqMut_perGenome <- function(i){ # For each genome
	dir.create("10_Analysis_per_genome/Main_mut_Freq-VOC", showWarnings=FALSE)
	dir.create("10_Analysis_per_genome/Main_mut_Freq-Lineage", showWarnings=FALSE)
	print(i) # Start processing the ith genome
	name <- colnames(dfVOC)[i] # This holds the name of the genome
	allFreqsVOC <- dfVOC[i] # Exctact the frequencies of all mutations detected in that genome (it is important to preserve it as a dataframe to keep rownames
	allFreqsLin <- dfLin[i]
	freqVOC <- lapply(lapply(main_mut_map_VOC,names), function(x){round(allFreqsVOC[x,1],2)}) # Now, use the VOC/lineage map to extract the corresponding values
	freqLin <- lapply(lapply(main_mut_map_Lin,names), function(x){round(allFreqsLin[x,1],2)})
	out_tableVOC <- matrix("",ncol=length(main_mut_map_VOC), nrow=nrow(dfVOC)); colnames(out_tableVOC) <- names(main_mut_map_VOC); rownames(out_tableVOC) <- rownames(dfVOC) # and prepare an output tables containing the values in a square matrices
	out_tableLin <- matrix("",ncol=length(main_mut_map_Lin), nrow=nrow(dfLin)); colnames(out_tableLin) <- names(main_mut_map_Lin); rownames(out_tableLin) <- rownames(dfLin)
	for(var in names(main_mut_map_VOC)){out_tableVOC[names(main_mut_map_VOC[[var]]),var] <- freqVOC[[var]]} # fill the matrices with the corresponding items
# 	out_tableVOC <- out_tableVOC[,c("222","519","Alpha","Gamma","XB","Delta","Omicron_BA.1.x", "Omicron_BA.2.x", "Omicron_BA.4.x", "Omicron_BA.5.x", "XAS", "XBB.1", "XBB.1.2", "XBB.1.4", "XBB.1.5", "XBB.2")] # Order by VOC appearance and remove others
	out_tableVOC <- out_tableVOC[,c("222","519","Alpha","Gamma","Delta","Omicron_BA.1.x", "Omicron_BA.2.x", "Omicron_BA.4.x", "Omicron_BA.5.x","XB", "XBB", "XAS")] # Order by VOC appearance and remove others
	for(var in names(main_mut_map_Lin)){out_tableLin[names(main_mut_map_Lin[[var]]),var] <- freqLin[[var]]}
	out_tableLin <- out_tableLin[,sort(colnames(out_tableLin))]
	exp_wFreqVOC <- apply(out_tableVOC, 2, function(x){vect=as.numeric(x); vect=vect[!is.na(vect)]}) # Get the expected mutations per VOC/Variant with the corresponding frequencies
	exp_wFreqLin <- apply(out_tableLin, 2, function(x){vect=as.numeric(x); vect=vect[!is.na(vect)]})
	expVOC <- unlist(lapply(exp_wFreqVOC,length)) # Get the expected items
	expLin <- unlist(lapply(exp_wFreqLin,length))
	obsVOC <- unlist(lapply(exp_wFreqVOC,function(x){sum(x>0)})) # And the observed items > 0
	obsLin <- unlist(lapply(exp_wFreqLin,function(x){sum(x>0)}))
	completenessVOC <- round(obsVOC/expVOC,2) # Calculate how complete each variant is
	completenessLin <- round(obsLin/expLin,2)
	obsVOC_VAMIPs <- unlist(lapply(exp_wFreqVOC,function(x){sum((x>0)*(x<0.5))})) # add how many ot the observed items are from VAMIPs and SNPs and their separate completeness
	obsVOC_SNPs <- unlist(lapply(exp_wFreqVOC,function(x){sum(x>=0.5)}))
	completenessVOC_VAMIPs <- round(obsVOC_VAMIPs/expVOC,2)
	completenessVOC_SNPs <- round(obsVOC_SNPs/expVOC,2)
	obsLin_VAMIPs <- unlist(lapply(exp_wFreqLin,function(x){sum((x>0)*(x<0.5))})) # add how many ot the observed items are from VAMIPs and SNPs and their separate completeness
	obsLin_SNPs <- unlist(lapply(exp_wFreqLin,function(x){sum(x>=0.5)}))
	completenessLin_VAMIPs <- round(obsLin_VAMIPs/expLin,2)
	completenessLin_SNPs <- round(obsLin_SNPs/expLin,2)

write.table(rbind(out_tableVOC,expVOC,obsVOC, obsVOC_SNPs, obsVOC_VAMIPs, completenessVOC, completenessVOC_SNPs, completenessVOC_VAMIPs), paste0("10_Analysis_per_genome/Main_mut_Freq-VOC/", strip(name), ".tsv"), sep="\t", quote=FALSE, row.names=T, col.names=NA)
write.table(rbind(out_tableLin,expLin,obsLin, obsLin_SNPs, obsLin_VAMIPs, completenessLin, completenessLin_SNPs, completenessLin_VAMIPs), paste0("10_Analysis_per_genome/Main_mut_Freq-Lineage/", strip(name), ".tsv"), sep="\t", quote=FALSE, row.names=T, col.names=NA)
}

save.image("Genomes-VAMIPsVsGenomes_chkpt3.Rdata")
load("Genomes-VAMIPsVsGenomes_chkpt3.Rdata")
# i <- 4527 # test
for(i in 1:ncol(dfVOC)){
	print(i)
	FreqMut_perGenome(i)
}


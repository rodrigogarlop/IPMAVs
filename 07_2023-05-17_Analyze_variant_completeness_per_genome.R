# Started: 2023-05-17 by Rodrigo García for Profr. Carlos Arias's GAL group at IBt, UNAM.
# R version 4.3.0 (2023-04-21) -- "Already Tomorrow"
# This script was is intended to detect putative coinfections using a propietary completeness table built for each genome in the IPMAV set. This version takes into account when each two variants co-circulated in Mexico and the expected muations for each of them in the Mexican IPMAV set.
# The main input table is a map of frequencies with only the key mutations per variant in the set. The name of the table includes the GISAID ID separated by _-_ characters (these must be present for it to work)
# Secondary tables include:
# meta (Metadata table): Table containing the variants, dates, sex and other metadata, compiled for the whole IPMAV set
# tab_occur (Occurrence table): Square matrix containing how many days each pair of variants co-occurred. This was created by analyzing the whole mexican set in GISAID and filtering outliers with the script /home/rod/Documents/01_Projects/SARS/VAMIP/11_Epi_Genome/Variant_Analysis_VAMIPs_cooccurrence.R. The resulting tables were stored at /home/rod/Documents/01_Projects/SARS/VAMIP/07_metadata/02_Imported_from_Vigilancia_2023-03-16_recomb/ It was created with the script /home/rod/Documents/01_Projects/SARS/Vigilancia/2023-03-16_recomb/Variant_Analysis_extra_cooccurrence_info_for_VAMIPs.R, not just for the vamip set.
# tab_freq (Frequency table): This table was created to establish the expected frequency observed in other variants. It was created with the script /home/rod/Documents/01_Projects/SARS/VAMIP/Extract_putative_coinfections.R

# Run as follows:
# Rscript 2023-015-17_Analyze_VOC_completeness_per_genome.R <input_genome_table> <prefix_output> <expected_variant_coocurrence_table> <expected_frequency_table_for_other_variants>
# Paths may be included in arguments
# Test run:
# input="/home/rod/Documents/01_Projects/SARS/VAMIP/10_Analysis_per_genome/Main_mut_Freq-Lineage/2022-02-22_-_EPI_ISL_11110841_-_B.1.1.529_20B_Omicron_BA.1.x_-_047_-_F_-_Yucatan.tsv"
# input="/home/rod/Documents/01_Projects/SARS/VAMIP/10_Analysis_per_genome/Main_mut_Freq-Lineage/2023-01-27_-_EPI_ISL_16968378_-_B.1.1.529.5.3.1.1.1.1.1.1.2_20B_Omicron_BA.5.x_-_043_-_F_-_Nuevo_Leon.tsv"
# input="/home/rod/Documents/01_Projects/SARS/VAMIP/10_Analysis_per_genome/Main_mut_Freq-Lineage/2021-11-29_-_EPI_ISL_7970117_-_B.1.617.2.26_20A_Delta_-_025_-_F_-_Mexico_City.tsv"
input="/home/rod/Documents/01_Projects/SARS/VAMIP/10_Analysis_per_genome/Main_mut_Freq-Lineage/2021-06-07_-_EPI_ISL_2801669_-_B.1.1.28.1_20J_501Y.V3_Gamma_-_065_-_F_-_Nuevo_Leon.tsv"
prefix="/home/rod/Documents/01_Projects/SARS/VAMIP/10_Analysis_per_genome/Main_mut_Freq-Lineage/Freq_analysis/test"
meta="/home/rod/Documents/01_Projects/SARS/VAMIP/07_metadata/2023-05-01_Metadata_all_in_VAMIP_set_that_matched-derep_sortByTime_xbbFix.tsv"
tab_occur="/home/rod/Documents/01_Projects/SARS/VAMIP/07_metadata/02_Imported_from_Vigilancia_2023-03-16_recomb/cooccurring_variants-days-extended.tsv"
tab_freq="/home/rod/Documents/01_Projects/SARS/VAMIP/10_Analysis_per_genome/genome_variant_completeness/2023-05-08_Vars-mutIntersection_0.7rel_fullLineage.tsv"


###################### LOAD INPUT, FUNCTIONS AND LIBRARIES ######################
library("vioplot")
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 5) { # at least, 5 arguments should be included: <input_genome_table> <prefix_output> <expected_variant_coocurrence_table> <expected_frequency_table_for_other_variants>
  stop("A minimum of 5 arguments are mandatory, please run as follows: cat table.tsv|Rscript 2023-015-17_Analyze_VOC_completeness_per_genome.R <prefix_output> <metadata_table> <expected_variant_coocurrence_table> <expected_frequency_table_for_other_variants>", call.=FALSE)
}
input <- as.character(args[1]) # Get a string handle to read the input table (per genome with completeness)
prefix <- as.character(args[2]) # Get a string handle to create output names
meta <- as.character(args[3]) # Get the location of the metadata table
tab_occur <- as.character(args[4])  # Get the location of the main coocurrence matrix, showing which variants coexisted
tab_freq <- as.character(args[5]) # Get the location of the main frequency matrix, showing how common each mutation is expected to be in the set.

# print("Loaded parameters:")
print(paste("Input:",input))
# print(paste("Prefix:",prefix))
# print(paste("Table 1, metadata:",meta))
# print(paste("Table 2, expected completeness:",tab_occur))
# print(paste("Table 3, expected frequency:",tab_freq))

df <- read.table(input, sep="\t",header=T, skip=0, comment.char='',quote="",fill=F, row.names=1, check.names=F) # Load the input
meta <- read.table(meta, header=TRUE, sep ='\t',stringsAsFactors = FALSE, check.names=F, fill=T, quote="", skip=0) # Load metadata with id folio for xref. IDs should be unique
rownames(meta) <- meta[,"Accession ID"] # Use the epi Acc_num as id
tab_occur <- read.table(tab_occur, header=TRUE, sep ='\t',stringsAsFactors = FALSE, check.names=F, fill=TRUE, quote="", skip=0, row.names=1) # now, load the cooccurrence table
tab_freq <- read.table(tab_freq, header=TRUE, sep ='\t',stringsAsFactors = FALSE, check.names=F, fill=TRUE, quote="", skip=0, row.names=1) # and the frequency table

#################### MAIN ######################
# First, extract some important metadata
genome <- strsplit(input, "_-_")[[1]][2] # start by extracting the ID from the input name (the table should include it)
voc <- meta[genome, "VOCs_three"] # Then extract which VOC it is in order to compare it
pango <- meta[genome, "Linaje Pangolin"]
variant <- meta[genome, "New_Lineage"]
date <- meta[genome, "Fecha de recoleccion"]
if(sum(unlist(df["completenessLin",]))==0){ # there are a few variants that share no mutations with known VOCs, these are ignored
  stop("variants with that share no key mutations with known other variants will be ignored in this analysis, ABORTED", call.=FALSE)
}
# Then, compare to expected frequencies.
if(length(grep(paste0("^",variant,"$"), colnames(tab_freq)))==0){ # variants with very few observations are unreliable as we don't have enough data to determine the key mutations correctly and are thus a reason to abort execution
  stop("variants with very few observations are unreliable as we don't have enough data to determine the key mutations correctly, ABORTED", call.=FALSE)
}
expfreq <- tab_freq[,variant]; names(expfreq) <- rownames(tab_freq) # Extract the expected values for all other variants
pearson <- round(cor(expfreq[colnames(df)],unlist(df["completenessLin",])),4) # and store the pearson's correlation value with the completenessLin
pearson2 <- round(cor(expfreq[colnames(df)],unlist(df["completenessLin_SNPs",])),4) # and store the pearson's correlation value with the completenessLin_SNPs
# write.table(t(c(date,genome,voc,pango,variant, pearson, pearson2)), paste0(prefix,"_variant_result.tsv"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
if(pearson>=0.997){ # abort if no unexpected distribution is observed
	stop("Like this one, most genomes have no evidence of putative coinfection. IGNORING", call.=FALSE)
}
# If it was indeed unexpected, then extract the rest of the important metadata
state <- meta[genome, "Estado"]
age <- meta[genome, "Edad"]
sex <- meta[genome, "Genero"]
status <- meta[genome, "Estatus del paciente"]
newName <- paste(date,pango,voc,genome, state, age, sex, status, sep="|")
# Now, filter items that did not co-occurr in the same period
coocurrence <- tab_occur[,variant]; names(coocurrence) <- rownames(tab_occur)
df <- df[,coocurrence[colnames(df)]>0] # drop all that were not co-circulating
# Next, compare to expected frequencies
expVsObs <- as.data.frame(cbind("Expected"=expfreq[colnames(df)],"Observed"=unlist(df["completenessLin",]))) # For this, create a small matrix fo variants vs frequencies (observed and expected)
expVsObs <- cbind(expVsObs,"Diff"=apply(expVsObs,1,function(x){abs(x[2]-x[1])})) # and append the difference as a last column
# We want to detect items that are widely deviating from what was expected. For the cutoff, I'm using the mean+1sd
# meanPlus1SD <- mean(expVsObs[,"Diff"])+sd(expVsObs[,"Diff"])
meanPlus.5SD <- mean(expVsObs[,"Diff"])+0.5*sd(expVsObs[,"Diff"]) # This alternative version considers halft an s.d., which captures around 39.2% cases instead of 68.2% in a normal distro
assigned <- expVsObs[rownames(expVsObs)==variant,] # temporarily store the expected_vs_observed for the assigned variant
rownames(assigned) <- paste0("-",rownames(assigned),"-")
expVsObs <- expVsObs[rownames(expVsObs)!=variant,] # Remove self hit (the actual variant that was assigned)
# Finally, identify items largely deviating from the expected values. For min completeness, I'll be requesting at least 0.5 and for the deviation I'm using the calculation above
putative_coinfections <- expVsObs[apply(expVsObs, 1, function(x){as.logical((x[2]>=0.5)*(x[3]>meanPlus.5SD))}),]
if(nrow(putative_coinfections)==0){
	stop("No alternative variant has an unexpected completeness for this genome. IGNORING", call.=FALSE)
}
out_tab <- df[,c(variant,rownames(putative_coinfections))] # The first item will be the official assignation
colnames(out_tab)[1] <- paste0("-",colnames(out_tab)[1],"-")
out_tab[is.na(out_tab)] <- -1; out_tab <- out_tab[rowSums(out_tab)!=(-1*ncol(out_tab)),];out_tab[out_tab==-1] <- NA # Remove mutations that are not observed in either the official variant or the secondary items
temp <- t(rbind(assigned,putative_coinfections)); temp <- temp[-2,]; rownames(temp) <- paste0("completeness",rownames(temp))
temp <- round(temp,4)
eot <-grep("expLin",rownames(out_tab)) # Detect where the end of the table is
sub_full <-sub <- out_tab[1:(eot-1),]# subset mutations
sub <- sub[is.na(sub[,1]),] # from the rest, extract which are not seen in the assigned variant
extra <- apply(sub,2,function(x){set=x[!is.na(x)];c("extraObserved"=sum(set>0),"extraExpected"=length(set))}) # counth the expected and observed items within the remaining differences
sub_VAMIP <- sub; sub_VAMIP[sub_VAMIP>=0.5] <- 0  # and create a copy to calculate only vamips (those extras where they are seen as vamips only)
extra <- rbind(extra, "extraObservedVAMIP"=apply(sub_VAMIP,2,function(x){set=x[!is.na(x)];sum(set>0)})) # and append the extra observed.
extra <- rbind(extra,"extraObservedSNP"=(extra["extraObserved",]-extra["extraObservedVAMIP",]))
out_tab <- rbind(out_tab,temp,extra)
# Now test for rare cases. The first filter is for those not having at least 3 extra items from the putative coinfecting variants
passing <- extra["extraObserved",]>4; passing[1] <- TRUE # Flag those passing and prevent the first item from being removed (those that have at least this many extra mutations
if(sum(passing)<2){ # 2 is the minimum for an extra variant
	stop("No alternative variant has an enough mutations from other variants to consider coinfections. IGNORING", call.=FALSE)
}
# The second filter is for those having 50% or less of the expected extra mutations for those variants
ratio1_obs <- extra["extraObserved",]/extra["extraExpected",];ratio1_obs[1] <-1; ratio1_obs <- ratio1_obs>0.5
if(sum(ratio1_obs)<2){ # 2 is the minimum for an extra variant
	stop("No alternative variant has an enough mutations from other variants to consider coinfections. IGNORING", call.=FALSE)
}
ratio2_obsVAMIPs <- extra["extraObservedVAMIP",]/extra["extraExpected",];ratio2_obsVAMIPs[1] <-1; ratio2_obsVAMIPs <- ratio2_obsVAMIPs>0.25
if(max(ratio2_obsVAMIPs)<=0.25){
	stop("No alternative variant has an enough mutations from other variants to consider coinfections. IGNORING", call.=FALSE)
}
logic_vect <- as.logical(passing*ratio1_obs*ratio2_obsVAMIPs) # Update the passing vector
if(sum(logic_vect==TRUE)==1){
	stop("No alternative variant recognized as a putative coinfection. IGNORING", call.=FALSE)
}
out_tab <- out_tab[,logic_vect] # and subset to keep only those passing
# coinfect <- paste(names(out_tab)[2:ncol(out_tab)],sep="|") # Store which were putative coinfections
coinfect <- paste0(names(out_tab)[2:ncol(out_tab)]," Mut:", out_tab["obsLin",2:ncol(out_tab)], "/", out_tab["expLin",2:ncol(out_tab)], " Uni:", "(", out_tab["extraObservedVAMIP",2:ncol(out_tab)], ")", out_tab["extraObserved",2:ncol(out_tab)], "/", out_tab["extraExpected",2:ncol(out_tab)])
out_tab <- out_tab[rowSums(!is.na(out_tab))>0,] # Remove mutations that are not observed in either the official variant or the secondary items
sub_full <- sub_full[,logic_vect] # filter the frequencies as well
sub_full <- sub_full[rowSums(!is.na(sub_full))>0,] # and filter those not seen in the table as well
reset <- !is.na(sub_full[,1]) # create a boolean vector to flag which mutations are seen in the target variant
temp <- sub_full[,1] # and temporarily save the frequencies for the assigned variant
sub_full <- apply(sub_full, 2, function(x){x[reset]<-NA;x}) # Remove those seen in the main variant
sub_full[,1] <- temp # Revert the first column (that of the actual assigned variant) because it should be preserved
firstlast <- apply(sub_full, 2, function(x){x <- x[!is.na(x)];c("PosFirst"=as.numeric(sub(":.*","",names(x)[1])),"PosLast"=as.numeric(sub(":.*","",names(x)[length(x)])))}) # get the first and last mutation per variant
firstlast <- rbind(firstlast,"Diff"=(firstlast[2,]-firstlast[1,])) # append the difference (I did not use function deff() in order to use a custom name)
firstlast <- rbind(firstlast, "Partial"=(firstlast["Diff",]<20000)*99999) # Prints 99999 if the span of key mutations per variant is shorter than 20K in lenght. The number was arbitrary and was selected to obtain a number that wouldn't appear otherwise. We'd expect for early variants to be flagged here as well as recombinations
out_tab <- rbind(out_tab,firstlast) # Append the extra columns to the output table
plotData <- apply(sub_full,2,function(x){x[!is.na(x)]},simplify = FALSE) # Create a list object with no NAs for each item

# ### OUTPUT RESULTS ###
temp1 <- out_tab["extraObserved",]; temp1[1] <- out_tab["obsLin",1]
temp2 <- out_tab["extraExpected",]; temp2[1] <- out_tab["expLin",1]
pdf(paste0(prefix,"_Lin_result.pdf"),width=length(plotData)+3.5)
	par(mar=c(7.1,5.1,4.1,2.1))
	vioplot(las=1, plotData, boxwex=as.numeric(out_tab["completenessVOC",]), main=newName, cex.main=1, ylab="Key Mutations' Frequency          ", col="white", ylim=c(0,1.2), frame=FALSE, xaxt='n', yaxt='n', notch=TRUE, border="purple",boxlwd = 3, cex=2, lwd=2,lineCol="purple", rectCol="mediumpurple1",colMed="violet")
	stripchart(plotData, vertical = TRUE, method = "jitter",pch = 2, add = TRUE, col = "gray30", cex=0.8, jitter = 0.1)
	mtext("Box width depicts variant completeness",cex=1.2)
# 	par(mar=c(3.1,2.1,2.1,2.1))
	axis(1, las=2, at=1:length(plotData), labels=sub("Omicron_","",names(plotData)),lwd.ticks = FALSE, lwd=NA, cex.axis=0.5)
	axis(2, las=1, at=c(seq(0,1,0.1),1.05,1.15), labels=c(seq(0,1,0.1),"ExMut","Comp"), cex.axis=1)
	text(1:length(plotData), 1.1, paste0(out_tab["completenessLin",]*100,"%\n",temp1,"/",temp2),cex=1.5)
dev.off()
write.table(t(c(date,genome,voc,pango,variant, pearson, pearson2, coinfect)), paste0(prefix,"_Lin_result.tsv"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE) # Output a table with useful parameters for deciding on coinfections (including the correlations and which are the putative coinfections)
out_tab[is.na(out_tab)] <- ""
write.table(out_tab, paste0(prefix,"_Lin_full_MutTable.tsv"), sep="\t", quote=FALSE, row.names=TRUE, col.names=NA) # Output the updated table





















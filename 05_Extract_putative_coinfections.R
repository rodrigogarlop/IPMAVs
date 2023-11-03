# Started 2023-05-08 by Rodrigo García-López for Carlos Arias's Virology Group at IBt, UNAM, Cuernavaca, Mexico as part of the CoViGen-Mex SARS-CoV-2 survillance in Mexico.
# R version 4.2.3 (2023-03-15) -- "Shortstop Beagle"
# Under GNU GPLv3 license
# Disclamer: This script was written for use with specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.
# This script was created to take the large collection of SARS-CoV-2 mutations including IPMAVs and the expected mutations per variant to identify which genomes are putative coinfections where more than one variant is sequenced.
# The first input table is a table of expected mutations. This is a contingency table containing all mutations found (including VAMIPs) in the set vs all the variants. The script 2023-05-03_VAMIPS_per_genome_comparison.R was used to create this table. There is a relative and absolute version, we'll use the former. There are also two of these, one calculated for VOCs, one for separate variants, we'll use both but with different cutoffs as the former has a higher variability within the groups.
# An example of the table is as follows:
#              Delta Omicron_BA.5.x Omicron_BA.1.x   519 Omicron_BA.2.x Others
# 00004:A->G  0.0004              0              0 3e-04              0 0.0013
# 00004:A->T  0.0015              0              0 3e-04              0 0.0000

# ### LOAD LIBRARIES/FUNCTIONS ###
calculate_intersection <- function(matrix, cutoff=0.6){ # This function uses a collection of all mutations vs specific variants and a cutoff and returns three matrices: 1.- a raw list of mutations passing the cutoff for frequency, 2.- An intersection triangular matrix of items shared among each pair of variants. 3.- Same as 2 but relative (by row).
	expMut <- (matrix>=cutoff)*1 # Convert to boolean matrix,
	expMut <- expMut[(rowSums(expMut)>0),] # and clean to remove empty items (not a central part of the actual variants in the set)
	sumExpMut <- colSums(expMut) # get an additional vector having the total items per variant
	simi_Mut <- apply(expMut, 2, function(x){colSums(x*expMut)}) # Calculate the intersection for each item vs all variants
	simi_Mut_rel <- simi_Mut/sumExpMut # Addionally calculate the frequency based on the expected items per variant
	out <- list(raw=expMut,abs=simi_Mut, rel=simi_Mut_rel) # Export a list of three matrices
	return(out)
}

# ### MAIN ###
# First, load the desired table of expected mutations (see documentation above)
df <- read.table("/home/rod/Documents/01_Projects/SARS/VAMIP/10_Analysis_per_genome/genome_VOC_completeness/2023-05-03_Main_mutations_per_VOC-filtered-rel.tsv", sep="\t",header=T, skip=0, comment.char='',quote="",fill=F,row.names=1, check.names=FALSE) # Load input table
VOCs <- calculate_intersection(df,0.6)
write.table(VOCs$raw, "/home/rod/Documents/01_Projects/SARS/VAMIP/10_Analysis_per_genome/genome_VOC_completeness/2023-05-08_VOCs-representative_mutations_0.6Freq.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)
write.table(VOCs$abs, "/home/rod/Documents/01_Projects/SARS/VAMIP/10_Analysis_per_genome/genome_VOC_completeness/2023-05-08_VOCs-mutIntersection_0.6Abs.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)
write.table(VOCs$rel, "/home/rod/Documents/01_Projects/SARS/VAMIP/10_Analysis_per_genome/genome_VOC_completeness/2023-05-08_VOCs-mutIntersection_0.6rel.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)
# Again, for the full lineages
df <- read.table("/home/rod/Documents/01_Projects/SARS/VAMIP/10_Analysis_per_genome/genome_variant_completeness/2023-05-03_Main_mutations_per_variant-filtered-rel.tsv", sep="\t",header=T, skip=0, comment.char='',quote="",fill=F,row.names=1, check.names=FALSE) # Load input table
Vars <- calculate_intersection(df,0.7)
write.table(Vars$raw, "/home/rod/Documents/01_Projects/SARS/VAMIP/10_Analysis_per_genome/genome_variant_completeness/2023-05-08_Vars-representative_mutations_0.7Freq.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)
write.table(Vars$abs, "/home/rod/Documents/01_Projects/SARS/VAMIP/10_Analysis_per_genome/genome_variant_completeness/2023-05-08_Vars-mutIntersection_0.7Abs.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)
write.table(Vars$rel, "/home/rod/Documents/01_Projects/SARS/VAMIP/10_Analysis_per_genome/genome_variant_completeness/2023-05-08_Vars-mutIntersection_0.7rel.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)

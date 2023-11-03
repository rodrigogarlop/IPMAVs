# Updated 2023-07-28 Added a 5 obs filter
# Updated 2023-05-02 from version 2022-08-25_Analyze_vamip_contingency_Bar_n_boxplots.R at /home/rod/Documents/01_Projects/SARS/VAMIP/Resp_2023-04-11/
# Started: 2022-08-25
# by Rodrigo García-López for Carlos Arias's Virology Group at IBt, UNAM, Cuernavaca, Mexico as part of the CoViGen-Mex SARS-CoV-2 survillance in Mexico.
# R version 4.2.3 (2023-03-15) -- "Shortstop Beagle"
# Under GNU GPLv3 license
# Disclamer: This script was written for use with specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.
# This script was written to study each mutation (allele) separatelly, based on the allelic map of the whole CoViGen-Mex dataset
# The script is a modified version of the analysis scripts ML_sequences.R and 2022-06-26_Analyze_vamip_contingency.R, and is intended to analyze VAMIPs in longitudinal data from a contingency table containing genomes (columns) x particular mutations (rows). The former include some metadata in their names, which can be parsed previously. Old update notes older than this script were preseved for documentation purposes

# IMPORTANT NOTE 1: We referred here to fixated items as "snps", please note these are the ones with freq 0.5 or higher
# IMPORTANT NOTE 2: Please note most results were not used in the final version of the published article, but these can be used to study mutations longitudinally.

# ### LOAD LIBRARIES AND FUNCTIONS ###
library("pheatmap") # This is used for heatmaps (better than R base's
library("stringr") # This is required to extract the metadata easily
plot_single_mutation <- function(vect, map, item_name="Item", comparison_name="Comparison", min_nonzero=2){ # Input items are a frequency vector with NAs (vect), a list containing positions per group (e.g. months, age, etc.; map), names for the actual vector (mutation are expected; item_name) and the category that is explored (comparison_name), and minimum expected subcategories (groups) that should have the item (min_nonzero).
	list_totals <- unlist(lapply(map,function(x){temp <- vect[x]; temp <- temp[!is.na(temp)]; length(temp)}))
	if(sum(list_totals > 0) >= min_nonzero){
		if(max(list_totals)<10){return()} # only consider mutations where the max genomes per category (e.g. any given month) equal at least this
		if(sum(list_totals)<30){return()} # only consider mutations that have at least 30 total items
		list_out <- lapply(map,function(x){temp <- vect[x]; temp <- temp[!is.na(temp)]; if(length(temp)==0){temp <- 0}; temp}) # for each item in the map, extract non-NA values from the vector into a list output
		dir.create(comparison_name,showWarnings = FALSE)
		ref_totals <- unlist(lapply(map, length)) # get the totals per subcategory (e.g. month)
		item_name_new <- sub("\\|", "_", item_name) # Fix characters for file names
		item_name_new <- gsub(":", "_", item_name_new)
		item_name_new <- gsub(" ", "", item_name_new)
		item_name_new <- gsub(",", "n", item_name_new)
		item_name_new <- gsub(">", "", item_name_new)
		item_name_new <- gsub("-", "t", item_name_new)
		pdf(paste0(comparison_name, "/", comparison_name,"_",item_name_new,".pdf"), width=14)
		par(oma = c(3, 1, 1, 10)) # This is just a creative fix to plot the legend outside
		# 		par(mar=c(5,4,2,5))
		bp <- barplot(list_totals/ref_totals, las=2, main=paste(item_name, "by", comparison_name), ylab="Frequency", ylim=c(-0.1,1.1), col="coral1", border=NA, yaxt='n')
		boxplot(list_out, add=T, las=2, at=bp, border="navyblue", col=NA, outline = T, yaxt='n')
		axis(2,at=seq(0,1,0.1), las=2)
		axis(4,at=seq(0,1,0.1), labels=seq(0,1,0.1)*max(list_totals), las=2, col="chartreuse4",col.axis="chartreuse4")
		points(bp,list_totals/max(list_totals), col="chartreuse3", pch=17, cex=1.5)
		par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0.5), mar = c(0, 0, 0, 0), new = TRUE)
		plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
		legend("right", legend = c("Mut Freq","Allele Freq","Total with Mut"), xpd = TRUE, horiz = FALSE, inset = c(0,0), bty = "n", col = c("coral1", "navyblue", "chartreuse3"), lty = c(NA,1,NA), cex = 0.8, pch=c(15,NA,17))
		dev.off()
		return(paste0("Created: ",comparison_name,"_",item_name_new,".pdf"))
	}
}
plot_single_mutation_long <- function(vect, map, item_name="Item", comparison_name="Comparison", min_nonzero=2){ # Input items are a frequency vector with NAs (vect), a list containing positions per group (e.g. months, age, etc.; map), names for the actual vector (mutation are expected; item_name) and the category that is explored (comparison_name), and minimum expected subcategories (groups) that should have the item (min_nonzero).
	list_totals <- unlist(lapply(map,function(x){temp <- vect[x]; temp <- temp[!is.na(temp)]; length(temp)}))
	if(sum(list_totals > 0) >= min_nonzero){
		if(max(list_totals)<10){return()} # only consider mutations where the max genomes per category (e.g. any given month) equal at least this
		if(sum(list_totals)<30){return()} # only consider mutations that have at least 30 total items
		list_out <- lapply(map,function(x){temp <- vect[x]; temp <- temp[!is.na(temp)]; if(length(temp)==0){temp <- 0}; temp}) # for each item in the map, extract non-NA values from the vector into a list output
		dir.create(comparison_name,showWarnings = FALSE)
		ref_totals <- unlist(lapply(map, length)) # get the totals per subcategory (e.g. month)
		item_name_new <- sub("\\|", "_", item_name) # Fix characters for file names
		item_name_new <- gsub(":", "_", item_name_new)
		item_name_new <- gsub(" ", "", item_name_new)
		item_name_new <- gsub(",", "n", item_name_new)
		item_name_new <- gsub(">", "", item_name_new)
		item_name_new <- gsub("-", "t", item_name_new)
		pdf(paste0(comparison_name, "/", comparison_name,"_",item_name_new,".pdf"), width=28)
		par(oma = c(2, 1, 1, 10)) # This is just a creative fix to plot the legend outside
		# 		par(mar=c(5,4,2,5))
		bp <- barplot(list_totals/ref_totals, las=2, main=paste(item_name, "by", comparison_name), ylab="Frequency", ylim=c(-0.1,1.1), col="coral1", border=NA, yaxt='n')
		boxplot(list_out, add=T, las=2, at=bp, border="navyblue", col=NA, outline = T, yaxt='n')
		axis(2,at=seq(0,1,0.1), las=2)
		axis(4,at=seq(0,1,0.1), labels=seq(0,1,0.1)*max(list_totals), las=2, col="chartreuse4",col.axis="chartreuse4")
		points(bp,list_totals/max(list_totals), col="chartreuse3", pch=17, cex=1.5)
		par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0.5), mar = c(0, 0, 0, 0), new = TRUE)
		plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
		legend("right", legend = c("Mut Freq","Allele Freq","Total with Mut"), xpd = TRUE, horiz = FALSE, inset = c(0,0), bty = "n", col = c("coral1", "navyblue", "chartreuse3"), lty = c(NA,1,NA), cex = 0.8, pch=c(15,NA,17))
		dev.off()
		return(paste0("Created: ",comparison_name,"_",item_name_new,".pdf"))
	}
}
	plot_longitudinal_distro_per_category(unlist(df[mut,]), longitudinal_map_weeks, all_mainvar, sub, "Variants_per_week",4)
# UPDATE 2022-08-25: A new version of the plot_single_mutation to classify items was added, this will get an additional genome classification for subclassifying per category (e.g. variant composition per month per mutation)
plot_longitudinal_distro_per_category <- function(vect, map, genome_classif, item_name="Item", comparison_name="Comparison", min_nonzero=2){ # Input items are a frequency vector with NAs (vect), a list containing positions per group (e.g. months, age, etc.; map), names for the actual vector (mutation are expected; item_name) and the category that is explored (comparison_name), and minimum expected subcategories (groups) that should have the item (min_nonzero). genome_classif is now included to hold a classification of the genomes and should correspond exactly to the same positions found in map
	names(vect) <- NULL # This is only required for testing, really
	sub_map <- lapply(map,function(x){temp <- vect[x];mask <- !is.na(temp); out <- x[mask]; if(length(out)==0){OUT <- NA}; out}) # Identify those that have items with nonNA values
	list_totals <- unlist(lapply(sub_map,function(x){out <- length(x);if(out==1){if(is.na(x)){out <- 0}};out}))
	if(sum(list_totals > 0) >= min_nonzero){
		if(max(list_totals)<10){return()} # only consider mutations where the max genomes per category (e.g. any given month) equal at least this
		if(sum(list_totals)<30){return()} # only consider mutations that have at least 30 total items
		list_out <- lapply(sub_map,function(x){ifelse(!is.na(x),vect[x],0)}) # Extract the frequencies where present
		names(genome_classif) <- NULL # Since the genome_classif has the same order, we require no labels, just the cat contents
		classif_genome_vars <- names(table(genome_classif)) # get a list of all available genome classifications (e.g. variants)
		classif_map_cats <- names(list_totals) # and one for the actual map (categories used for separating groups; e.g. months)
		map_vs_cats <- matrix(0, nrow = length(classif_genome_vars), ncol = length(classif_map_cats)) # output a zero matrix with total classifications for genomes and groups
		colnames(map_vs_cats) <- classif_map_cats
		rownames(map_vs_cats) <- classif_genome_vars
		genome_cat_out <- lapply(sub_map,function(x){if(sum(is.na(x))>0){NA}else{table(genome_classif[as.numeric(x)])}})
		# genome_cat_out <- lapply(sub_map,function(x){table(genome_classif[as.numeric(x)])}) # get the corresponding classifications based on genome_classif
		for(i in classif_map_cats){
			items <- genome_cat_out[[i]]
			if(sum(is.na(items))>0){next}
			temp <- names(items)
			for(j in temp){
				# print(genome_cat_out[[i]][j])
				map_vs_cats[j,i] <- genome_cat_out[[i]][j]
			}
		}
		dir.create(comparison_name,showWarnings = FALSE)
		dir.create(paste0(comparison_name,"/tab"),showWarnings = FALSE)
		ref_totals <- unlist(lapply(map, length)) # get the totals per subcategory (e.g. month)
		item_name_new <- sub("\\|", "_", item_name) # Fix characters for file names
		item_name_new <- gsub(":", "_", item_name_new)
		item_name_new <- gsub(" ", "", item_name_new)
		item_name_new <- gsub(",", "n", item_name_new)
		item_name_new <- gsub(">", "", item_name_new)
		item_name_new <- gsub("-", "t", item_name_new)
		pdf(paste0(comparison_name, "/", comparison_name,"_",item_name_new,".pdf"), width=14)
		par(oma = c(3, 1, 1, 10)) # This is just a creative fix to plot the legend outside
		# 		par(mar=c(5,4,2,5))
		col <- c('chartreuse3', 'cornflowerblue', 'darkgoldenrod1', 'peachpuff3','mediumorchid2', 'turquoise3', 'wheat4','slategray2',"coral1","aquamarine2","blue2","violetred2","palegreen3","purple3","magenta1","limegreen","darkorange2","darkgray")
		# bp <- barplot(map_vs_cats/max(list_totals), las=2, main=paste(item_name, "by", comparison_name), ylab="Frequency", col=col, border=NA, yaxt='n',cha ylim=c(-0.1,1.1)) # Use this to plot totals
		bp <- barplot(	t(apply(map_vs_cats,1,function(x){x/ref_totals})), las=2, main=paste(item_name, "by", comparison_name), ylab="Frequency", col=col, border=NA, yaxt='n', ylim=c(-0.1,1.1)) # Use this to scale to actual total prevalence
		tabout <- cbind(t(map_vs_cats), t(as.data.frame(lapply(list_out,function(x){quantile(x, seq(0,1,.25))}))))
		write.table(tabout, paste0(comparison_name, "/tab/", comparison_name,"_",item_name_new,".tsv"), sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)
		boxplot(list_out, add=T, las=2, at=bp, border="navyblue", col=NA, outline = T, xaxt='n', yaxt='n')
		axis(2,at=seq(0,1,0.1), las=2, col="navyblue",col.axis="navyblue")
		axis(4,at=seq(0,1,0.1), labels=seq(0,1,0.1)*max(list_totals), las=2, col="black",col.axis="black")
		points(bp,list_totals/max(list_totals), col="black", pch=17, cex=0.5)
		par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0.5), mar = c(0, 0, 0, 0), new = TRUE)
		plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
		legend("right", legend = classif_genome_vars, xpd = TRUE, horiz = FALSE, inset = c(0,0), bty = "n", col = col[1:length(classif_genome_vars)], cex = 0.8, pch=15)
		# legend("right", legend = c("Mut Freq","Allele Freq","Total with Mut"), xpd = TRUE, horiz = FALSE, inset = c(0,0), bty = "n", col = c("coral1", "navyblue", "chartreuse3"), lty = c(NA,1,NA), cex = 0.8, pch=c(15,NA,17))
		dev.off()
		return(paste0("Created: ",comparison_name,"_",item_name_new,".pdf"))
	}
}
save_pheatmap_pdf <- function(x, filename, width, height) { # This is a printing function for pheatmaps with custom filename, width and height. X is the input heatmap made with the pheatmap package.
	pdf(filename,width,height)
	grid::grid.newpage()
	grid::grid.draw(x$gtable)
	dev.off()
}
strip <- function(string){ # Takes a string object and strips rare characters
	string <- iconv(string,from="UTF-8",to="ASCII//TRANSLIT") # Remove accents and rare chars
# 	string <- tolower(string) # change to lower
	string <- gsub("(?<=[\\s])\\s*|^\\s+|\\s+$", "", string, perl=TRUE) # Remove multiple spaces
	string <- gsub("\\+","",string) # remove + characters
	string <- gsub("\\("," ",string)
	string <- gsub("\\)"," ",string)
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
group_freq_mean <- function(gmap,inmat) { # Get the group map from the corresponding function and get the prevalence per row item. This is calculated based on the total non-NAs.
	nam <- names(gmap) # Save names for indices
	out_raw <- matrix(NA, nrow=nrow(inmat), ncol=length(nam)); rownames(out_raw) <- rownames(inmat); colnames(out_raw) <- nam # Create a void matrix container
	for(i in nam){ # now, for each table
		# 		print(i)
		single_grp <- inmat[,gmap[[i]]] # First, subset the input table using the group map for a single group
		temp <- apply(single_grp, 1, function(x){mean(x[which(!is.na(x))])})  # and calculate the total items
		temp[is.nan(temp)] <- NA
		out_raw[,i] <- temp
		# 		out_rel[,i] <- out_raw[,i]*100/ncol(single_grp)
	}
	return(out_raw)
}
group_freq_centile <- function(gmap,inmat,cent=50) { # Get the group map from the corresponding function and get the prevalence per row item. This is calculated based on the total non-NAs.
	cent <- cent*0.01
	nam <- names(gmap) # Save names for indices
	out_raw <- matrix(NA, nrow=nrow(inmat), ncol=length(nam)); rownames(out_raw) <- rownames(inmat); colnames(out_raw) <- nam # Create a void matrix container
	for(i in nam){ # now, for each table
		# 		print(i)
		single_grp <- inmat[,gmap[[i]]] # First, subset the input table using the group map for a single group
		out_raw[,i] <- apply(single_grp, 1, function(x){quantile(x[which(!is.na(x))],cent)})  # and calculate the quantile (centile 0-100)
		# 		temp[is.nan(temp)] <- NA
		# 		out_raw[,i] <- temp
		# 		out_rel[,i] <- out_raw[,i]*100/ncol(single_grp)
	}
	return(out_raw)
}


# ### LOAD INPUTS ###
df <- read.table("Alt_variant_Freq.tsv", sep="\t",header=T, skip=0, comment.char='',quote="",fill=F,row.names=1, check.names=FALSE) # Load input table
save.image("Analyze_vamip_contingency_Bar_n_boxplots_chkpt1.Rdata")
# load("Analyze_vamip_contingency_Bar_n_boxplots_chkpt1.Rdata")
# Each column was named as follows: "2021-03-04|L004_202101083013|B.1.1.222(20B_Others)|046|M|Puebla" and they have been sorted before
dim(df)
# [1]  25053 29661
sum(!is.na(df))
# [1] 1866675

# PREPARE METADATA
# Process the name to extract basic info and the Folio ID
# names are processed into a table: Date|Folio ID|Pango(clade)|age|gender|State
# Example: "2021-03-04|L004_202101083013|B.1.1.222(20B_Others)|046|M|Puebla"
# minimetatab <- t(as.data.frame(str_split(colnames(df), "\\|"))) # Extract a mini metadata table from the names (including
minimetatab <- t(as.data.frame(sapply(colnames(df),function(x){strsplit(x,"\\|")[[1]]}))) #UPDATE 2023-05-02: After updating to R.4.2.3, the stringr, the str_split function stopped working so I changed it back to R base
rownames(minimetatab) <- NULL # row names are removed
meta <- read.table("07_metadata/2023-05-01_Metadata_all_in_VAMIP_set_that_matched.tsv", header=T, sep ='\t',stringsAsFactors = FALSE, check.names=F, fill=T, quote="", row.names=1) # Load metadata with id folio for xref. IDs should be unique
temp <- table(meta[,"Accession ID"]) # In case these are not unique, use the following solution to keep only the first unique item
temp <- sapply(names(which(temp>1)),function(x){grep(x,meta[,"Accession ID"])[-1]}) # identify all repetead items
meta <- meta[-temp,]
rownames(meta) <- meta[,"Accession ID"]
dim(meta)
# [1] 29661    60
# Now, explore basic dataset composition per genome
# Extract dates from the column names (each genome name includes the sample collection date)
# The rest of the metadata needs to be extracted from the external metadata table using the ID Folio for xref
meta <- meta[unique(minimetatab[,2]),]
dim(meta)
# [1] 29661    60
write.table(meta, "07_metadata/2023-05-01_Metadata_all_in_VAMIP_set_that_matched-derep_sortByTime.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
# write.table(meta, "07_metadata/Metadata_matched_only.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
# This table may be a larger than the mutations table (columns) but only matching items will be used

# ### MAIN ###
# Create two copies to compare only vamips (<0.05) and only major variants
df_vamip <- df; df_vamip[df_vamip>0.5] <- NA # this will hold vamips (>0.5 are now NAs)
df_novamip <- df; df_novamip[df_novamip<=0.5] <- NA # this will hold vamips (>0.5 are now NAs)

# Some metadata are already found in the genome's name and will be processed directly, others must be extracted from the external meta table
# DATE
# Some metadata are already found in the genome's name and will be processed directly, others must be extracted from the external meta table
genomes <- find_items(minimetatab[,2], meta, "File basename") # Extract the genome names (this should include their corresponding batches
# DATE
# all_Days <- as.Date(minimetatab[,1]) # First, the day # These were originally extracted from the name in the table but there is no need to calculate it again since it is already in the metadata table
# all_Months <- format(all_Days,"%Y-%m") # Then, the month
all_Days <- as.Date(find_items(minimetatab[,2], meta, "Fecha de recoleccion"))
all_Months <- find_items(minimetatab[,2], meta, "Month")
all_Weeks <- find_items(minimetatab[,2], meta, "Week")
# head(cbind("Full"=colnames(df),"Genome"=genomes, "Date"=as.character(all_Days), "Month"=all_Months, "Week"=all_Weeks))

# CLADES & VARIANTS
all_lineages <- find_items(minimetatab[,2],meta, "New_Lineage") # Start with the pango lineages
all_mainvar <- find_items(minimetatab[,2],meta, "VOCs_two") # Next, the main variants (VOCs and 519)
all_clade <- find_items(minimetatab[,2],meta, "Clado Nexstrain") # Next, the main clades
# AGE
all_age <- find_items(minimetatab[,2],meta, "Edad") # Start with the raw age
all_age_vac <- find_items(minimetatab[,2],meta, "Age_vac") # Now with vaccination group
# GENDER
all_sex <- find_items(minimetatab[,2],meta, "Genero")
# PATIENT STATUS
all_status <- find_items(minimetatab[,2],meta, "Estatus del paciente")
all_status[all_status=="Asymptomatic"] <- "Amb"
all_status[all_status=="Symptomatic"] <- "Hosp"
all_status[all_status=="Fatal"] <- "Dec"
all_status[all_status=="Released"] <- "Amb"
# LOCATION
all_region7 <- find_items(minimetatab[,2],meta, "Region_7")
# CT # These are no longer working on the 2023-03-16-based table
# all_ctrdrp <- round(as.numeric(find_items(minimetatab[,2],meta, "RdRP")))
# all_cte <- round(as.numeric(find_items(minimetatab[,2],meta, "E")))
# Additional items (no longer needed as they only have info on fastas)
# xtra_mut_Nt <- find_items(minimetatab[,2],meta, "Mutaciones Nucleotidos")
# xtra_mut_AA <- find_items(minimetatab[,2],meta, "Mutaciones Aminoacido")

out <- cbind("Full"=colnames(df),"Genome"=genomes, "Date"=as.character(all_Days), "Month"=all_Months, "Week"=all_Weeks, "Lineage"=all_lineages, "Variant"=all_mainvar, "Clade"=all_clade, "Age"=all_age, "Age_vac"=all_age_vac, "Sex"=all_sex, "Status"=all_status, "Region"=all_region7)
write.table(out,"07_metadata/2023-05-02_metadata_VAMIP_set-Mar2020-Jan2023.tsv",sep='\t', row.names=FALSE, col.names=TRUE)

# BASIC PLOTS
pdf("08_Analyses/Total_genomes-month.pdf", width=14)
par(oma=c(4,1,1,0))
barplot(las=2,table(all_Months),border=NA,col=NA,yaxt='n',ylim=c(0,2000), main="All Genomes - Monthly Total",ylab="Total Genomes")
axis(2,las=1, at=seq(0,2500,250))
abline(h=seq(0,2500,250), col="gray80", lty=2)
barplot(las=2,table(all_Months),border=NA,col="cornflowerblue", yaxt='n', xaxt='n', ylim=c(0,2000), add=TRUE)
dev.off()
pdf("08_Analyses/Total_genomes-week.pdf",width=28)
par(oma=c(4,1,1,0))
barplot(las=2,table(all_Weeks),border=NA,col=NA,yaxt='n', main="All Genomes - Weekly Total",ylab="Total Genomes")
axis(2,las=1, at=seq(0,2000,100))
abline(h=seq(0,2000,100), col="gray80", lty=2)
barplot(las=2,table(all_Weeks),border=NA,col="cornflowerblue", xaxt='n', yaxt='n', add=TRUE)
dev.off()
pdf("08_Analyses/Total_genomes-variants.pdf")
par(oma=c(4,1,1,0))
barplot(las=2,table(all_mainvar),border=NA,col=NA,yaxt='n',ylim=c(0,9000), main="All Genomes - Main VOCs Total", ylab="Total Genomes")
axis(2,las=1, at=seq(0,9000,1000))
abline(h=seq(0,9000,1000), col="gray80", lty=2)
barplot(las=2,table(all_mainvar),border=NA, col="coral1",xaxt='n', yaxt='n', ylim=c(0,9000), add=TRUE)
dev.off()
pdf("08_Analyses/Total_genomes-clades.pdf",width=14)
par(oma=c(2,1,1,0))
barplot(las=2,table(all_clade),border=NA,col="coral1",yaxt='n', main="All Genomes - Main VOCs Total", ylab="Total Genomes")
axis(2,las=1, at=seq(0,20000,2500))
dev.off()
temp <- table(all_age)
temp <- temp[-length(temp)]
# temp <- temp[as.character(0:101)]
pdf("08_Analyses/Total_genomes-age.pdf",width=14)
par(oma=c(2,1,1,0))
barplot(las=2,temp,border=NA,col=NA,yaxt='n', main="All Genomes - Age Total", ylab="Total Genomes", xlab="Ages")
axis(2,las=1, at=seq(0,700,100))
abline(h=seq(0,700,100), col="gray80", lty=2)
barplot(las=2,temp,border=NA,col="forestgreen",xaxt='n',yaxt='n', add=TRUE)
dev.off()
pdf("08_Analyses/Total_genomes-age_vac.pdf")
par(oma=c(2,1,1,0))
barplot(las=2,table(all_age_vac),border=NA,col=NA,yaxt='n', main="All Genomes - Age Vac Total", ylab="Total Genomes", xlab="Age Vac")
axis(2,las=1, at=seq(0,10000,1000))
abline(h=seq(0,10000,1000), col="gray80", lty=2)
barplot(las=2,table(all_age_vac),border=NA,col="forestgreen", xaxt='n', yaxt='n', add=TRUE)
dev.off()
pdf("08_Analyses/Total_genomes-sex.pdf",width=4)
par(oma=c(1,1,1,0))
barplot(las=2,table(all_sex),border=NA,col="purple",yaxt='n', main="All Genomes - Sex", ylab="Total Genomes", xlab="Sex")
axis(2,las=1, at=seq(0,18000,1000))
dev.off()
pdf("08_Analyses/Total_genomes-region7.pdf",width=7)
par(oma=c(1,1,1,0))
barplot(las=2,table(all_region7),border=NA,col="firebrick",yaxt='n', main="All Genomes - Region", ylab="Total Genomes", xlab="Region")
axis(2,las=1, at=seq(0,8000,500))
dev.off()
pdf("08_Analyses/Total_genomes-status.pdf",width=4)
par(oma=c(1,1,1,0))
barplot(las=2,table(all_status),border=NA,col="gold2",yaxt='n', main="All Genomes - Patient Status", ylab="Total Genomes", xlab="Status")
axis(2,las=1, at=seq(0,20000,2000))
dev.off()
save.image("Analyze_vamip_contingency_Bar_n_boxplots_chkpt2.Rdata")
load("Analyze_vamip_contingency_Bar_n_boxplots_chkpt2.Rdata")

# COMPARE PREVALENCE
# First, we want to compare a raw view of all totals, with and without the vamips:
All=rowSums(group_prevalence(group_map(all_sex),df)$raw)
only_VAMIPs=rowSums(group_prevalence(group_map(all_sex),df_vamip)$raw)
only_Consensus=rowSums(group_prevalence(group_map(all_sex),df_novamip)$raw)
compare_all <- cbind(All, only_VAMIPs, only_Consensus)
compare_all_total <- cbind(compare_all, compare_all/ncol(df)*100)
colnames(compare_all_total) <- c(colnames(compare_all_total)[1:3],paste(colnames(compare_all_total)[1:3],"%"))
write.table(compare_all_total, "08_Analyses/VAMIP_vs_Consensus.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)

pdf("08_Analyses/Histogram_mutations.pdf")
hist(log10(compare_all_total[,1]),main="Histograma: Genomas con cada mutación (todas)",las=1, col="turquoise3",border="white", xaxt='n', xlab='Cuántos genomas tienen esa mutación (10^log10)', ylab='Mutaciones puntuales', breaks=50,xlim=range(log10(compare_all_total[,1])))
temp <- c(5,2,10,20,50,100,200,500,1000,2000,5000,10000,max(compare_all_total[,1]))
axis(1,at=log10(temp),labels=(temp),las=2)
hist(log10(compare_all_total[,2]),main="Histograma: Genomas con cada mutación (sólo VAMIP)",las=1, col="cornflowerblue",border="white", xaxt='n', xlab='Cuántos genomas tienen esa mutación (10^log10)', ylab='Mutaciones puntuales', breaks=50,xlim=range(log10(compare_all_total[,1])))
temp <- c(5,2,10,20,50,100,200,500,1000,2000,5000,10000,max(compare_all_total[,1]))
axis(1,at=log10(temp),labels=(temp),las=2)
hist(log10(compare_all_total[,3]),main="Histograma: Genomas con cada mutación (sólo en Consensos)",las=1, col="coral1",border="white", xaxt='n', xlab='Cuántos genomas tienen esa mutación (10^log10)', ylab='Mutaciones puntuales', breaks=50,xlim=range(log10(compare_all_total[,1])))
temp <- c(5,2,10,20,50,100,200,500,1000,2000,5000,10000,max(compare_all_total[,1]))
axis(1,at=log10(temp),labels=(temp),las=2)
dev.off()
data <- rev(c("All Mut"=nrow(df),"All VAMIP"=sum(compare_all_total[,2]>0),"Consensus and VAMIP"=sum((compare_all_total[,3]>0)*(compare_all_total[,2]>0)),"All in Consensus"=sum(compare_all_total[,3]>0),"Only VAMIP"=sum(compare_all_total[,2]==compare_all_total[,1]),"Only in Consensus"=sum(compare_all_total[,3]==compare_all_total[,1])))
pdf("08_Analyses/VAMIP_consensus_compare.pdf")
par(oma=c(1,7,1,0))
barplot(las=2,data,horiz=T, col=c("cornflowerblue","coral1","coral1","cornflowerblue","coral1","chartreuse3"), border=NA)
dev.off()
# Create a screening of all mutations by month prevalence
g_month <- group_prevalence(group_map(sub("\\(.*","",all_Months)),df)
g_month_vamip <- group_prevalence(group_map(sub("\\(.*","",all_Months)),df_vamip)
g_month_novamip <- group_prevalence(group_map(sub("\\(.*","",all_Months)),df_novamip)
sub <- head(g_month$rel[names(sort(rowSums(g_month$rel),decreasing=TRUE)),],2000) # get the 2000 most common mutations only
sub[sub==0] <- NA
sub_vamip <- head(g_month_vamip$rel[names(sort(rowSums(g_month_vamip$rel),decreasing=TRUE)),],2000)
sub_vamip[sub_vamip==0] <- NA
sub_novamip <- head(g_month_novamip$rel[names(sort(rowSums(g_month_novamip$rel),decreasing=TRUE)),],2000)
# with these, we can get an overhead view of all the most common mutations
sub_novamip[sub_novamip==0] <- NA
out_heatmap <- pheatmap(sub, cluster_cols=FALSE, cluster_rows=FALSE)
save_pheatmap_pdf(out_heatmap, "08_Analyses/Month_heatmap_2000TopMut.pdf", 9, 230)
out_heatmap <- pheatmap(sub_vamip, cluster_cols=FALSE, cluster_rows=FALSE)
save_pheatmap_pdf(out_heatmap, "08_Analyses/Month_vamip_heatmap_2000TopMut.pdf", 9, 230)
out_heatmap <- pheatmap(sub_novamip, cluster_cols=FALSE, cluster_rows=FALSE)
save_pheatmap_pdf(out_heatmap, "08_Analyses/Month_novamip_heatmap_2000TopMut.pdf", 9, 230)
# Now, use the 2000 most abundant items in sub_vamip to extract them from the novamip set
sub_novamip_same_as_sub_vamip <- g_month_novamip$rel[rownames(sub_vamip),]
sub_novamip_same_as_sub_vamip[sub_novamip_same_as_sub_vamip==0] <- NA
out_heatmap <- pheatmap(sub_novamip_same_as_sub_vamip, cluster_cols=FALSE, cluster_rows=FALSE)
save_pheatmap_pdf(out_heatmap, "08_Analyses/Month_novamip_same_set_as_vamip_heatmap_2000TopMut.pdf", 9, 230)
# Compare those with vamips and novamips:
rownames(sub_novamip_same_as_sub_vamip) <- paste(rownames(sub_novamip_same_as_sub_vamip),"X novamip") # First, rename both sets
rownames(sub_novamip_same_as_sub_vamip) <- paste(sprintf("%04d", 1:2000),rownames(sub_novamip_same_as_sub_vamip))
rownames(sub_vamip) <- paste(rownames(sub_vamip),"VAMIP")
rownames(sub_vamip) <- paste(sprintf("%04d", 1:2000),rownames(sub_vamip))
compare <- rbind(sub_vamip,sub_novamip_same_as_sub_vamip)
compare <- compare[order(rownames(compare)),]
out_heatmap <- pheatmap(compare, cluster_cols=FALSE, cluster_rows=FALSE)
save_pheatmap_pdf(out_heatmap, "08_Analyses/Month_heatmap_compare_2000TopMut.pdf", 9, 460)

save.image("Analyze_vamip_contingency_Bar_n_boxplots_chkpt3.Rdata")
load("Analyze_vamip_contingency_Bar_n_boxplots_chkpt3.Rdata")

# Update 2023-07-28: The output of vamip comparison was changed to show only those with at least 5 observations
vectx <- apply(df,1, function(x){sum(x[!is.na(x)]>0)}) # Get total mutations with non-zero values
dfb <- df[vectx>4,] # Of these, subset those with at least 5 items
dim(dfb)
# [1] 14973 29661
has_vamips <- apply(dfb,2, function(x){sum(x[!is.na(x)]<0.5)}) # Get how many vamips were found per genome
sum(has_vamips<=5)*100/length(has_vamips)
# [1] 61.25889
sum(has_vamips<=5)
# [1] 18170
has_snps <- apply(dfb,2, function(x){sum(x[!is.na(x)]>=0.5)}) # Likewise, get total snps
sum(has_snps<=5)*100/length(has_snps)
# [1] 0.04045717
sum(has_snps<=5)
# [1] 12 # Most snps have over 5 items

pdf("08_Analyses/IMPORTANT_Total_IPMAVs_SNPs_per_Genome.pdf",width=10)
hist(las=1,has_vamips,breaks=200, main="IPMAVs and SNPs distribution per genome", xlab="Total Mutations", ylab="Total Genomes", col="coral1", border="white",xlim=c(0,100), ylim=c(0,9000),xaxt='n', yaxt='n')
hist(las=1,has_snps,breaks=50, col=rgb(0,1,1,0.5), border="white",add=TRUE)
legend("topright", legend=c("IPMAVs","SNPs"), pch=15, col=c("coral1", rgb(0,1,1)))
axis(las=1,1, seq(0,100, 10))
axis(las=1,2, seq(0,9000, 1000))
# abline(v=c(27,44,65,75,95)) # These are the modes
text(11, 350, labels="11")
text(27, 1050, labels="27")
text(43, 2280, labels="43")
text(65, 1100, labels="65")
text(75, 1800, labels="75")
text(95, 400, labels="95")
dev.off()

get_composition <- function(intA, intB){ # from intA mutations to intB-1 mutations, or else [intA, intB)
	mutations <- names(has_snps[as.logical((has_snps>=intA)*(has_snps<intB))]) # Get the genos where snps observed are between 40 and 50
	mutations <- sapply(mutations,function(x){strsplit(x,"\\|")[[1]][2]}) # extract only the items names
	out <- round(table(find_items(mutations, meta, "VOCs_three"))*100/length(mutations),2) # Get the % of Deltas in ranges 40-50
	return(out)
}
get_composition(10,20) # Check different total mutation ranges
# B.1.1.222 B.1.1.519    Others
#     32.39      0.75     66.87
get_composition(20,30)
# B.1.1.222 B.1.1.519     Delta     Gamma    Others
#      8.70     76.37      0.48      0.07     14.38
get_composition(30,40)
#          Alpha      B.1.1.222      B.1.1.519          Delta          Gamma
#           4.62           0.56          29.41          48.30           5.89
#        MultLin Omicron_BA.1.x         Others             XB
#           0.26           0.07          10.47           0.41
get_composition(40,50)
#          Alpha      B.1.1.519          Delta          Gamma        MultLin
#           4.55           0.12          79.91           7.19           0.14
# Omicron_BA.1.x         Others             XB
#           5.33           1.69           1.07
get_composition(50,60)
#          Alpha          Delta          Gamma        MultLin Omicron_BA.1.x
#           0.67          42.03           0.21           0.47          53.55
# Omicron_BA.2.x Omicron_BA.4.x Omicron_BA.5.x         Others             XB
#           1.25           0.05           0.52           1.14           0.10
get_composition(60,70)
#          Delta        MultLin Omicron_BA.1.x Omicron_BA.2.x Omicron_BA.4.x
#           0.14           0.10          78.79           8.94           0.48
# Omicron_BA.5.x            XAH            XAP            XAS          XBB.1
#          11.14           0.10           0.03           0.24           0.03
get_composition(70,80)
#        MultLin Omicron_BA.1.x Omicron_BA.2.x Omicron_BA.4.x Omicron_BA.5.x
#           0.05           1.21          25.43           6.37          66.15
#            XAF            XAG            XAH            XAM            XAS
#           0.02           0.02           0.02           0.02           0.65
#          XBB.1          XBB.3
#           0.06           0.02
get_composition(80,90)
# Omicron_BA.1.x Omicron_BA.2.x Omicron_BA.4.x Omicron_BA.5.x          XBB.1
#           0.04          15.21           9.96          73.97           0.58
#        XBB.1.2        XBB.1.5
#           0.07           0.18
get_composition(90,100)
# Omicron_BA.2.x Omicron_BA.4.x Omicron_BA.5.x            XBB          XBB.1
#          15.62           0.28           2.70           0.71          45.45
#        XBB.1.2        XBB.1.4        XBB.1.5        XBB.1.6        XBB.1.9
#           6.53           1.99          22.59           0.57           0.99
#          XBB.2          XBB.3
#           2.27           0.28

test <- has_snps/has_vamips # Test the ratio between both sets
length(test)
# [1] 29661
test <- test[!is.infinite(test)] # Remove those with 0 vamips (yields infinite)
length(test)
# [1] 28623


vectxcols <- apply(df,2, function(x){sum(x[!is.na(x)]>0)}) # and again, get mutations with non-zero values
sum(vectxcols==0) # Check total columns with
# [1] 0 # This is as expected, none have only the lowest frequency mutations in the set.

compare_all_total <- compare_all_total[compare_all_total[,1]>4,] # subset to keep only those >5
dim(compare_all_total)
# [1] 14973     6
write.table(compare_all_total, "08_Analyses/VAMIP_vs_Consensus.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)


# Identify type of mutations
nt2AA <- read.table("07_metadata/MutNtxAA_complete.tsv", sep="\t",header=T, skip=0, comment.char='',quote="",fill=F,row.names=1, check.names=FALSE) # Load nt to AA name table
nt2AA[nt2AA[,"MutAA"]=="NC","MutAA"] <- paste0(nt2AA[nt2AA[,"MutAA"]=="NC","Gene"],":NC") # add a category to flag UTRs
vamipSnps <- cbind(compare_all_total, nt2AA[rownames(compare_all_total),]) # Append AA info
write.table(vamipSnps, "08_Analyses/VAMIP_vs_Consensus_withAA.tsv", sep="\t", quote=FALSE, row.names=T, col.names=NA)

save.image("Analyze_vamip_contingency_Bar_n_boxplots_chkpt3.1.Rdata")
round(table(vamipSnps[,"Gene"])*100/nrow(vamipSnps),2)
# 3pUTR 5pUTR     E     M     N    NC ORF10 ORF1a ORF1b ORF3a  ORF6 ORF7a ORF7b
#  2.00  2.09  0.97  2.17  5.16  0.79  0.78 42.10 18.57  3.39  0.57  1.71  0.59
#  ORF8     S
#  1.59 17.52

table(vamipSnps[,"MutType"])

#  Del_FS Del_NFS  Ins_FS Ins_NFS      NC  NonSyn     Syn
#     383      92     110      14     730    8275    5369
# 2690
Coding_SNPs <- 8275 + 5369
round(table(vamipSnps[,"MutType"])*100/Coding_SNPs,2)[6:7]
# NonSyn    Syn
#  60.65  39.35

sort(vamipSnps[vamipSnps[,"only_VAMIPs"]==0,"only_Consensus"]) # Get the number of items in those that had no vamips
sort(vamipSnps[vamipSnps[,"only_Consensus"]==0,"only_VAMIPs"]) # and viceversa

cor.test(table(vamipSnps[vamipSnps[,"only_Consensus"]>0,"Gene"]), table(vamipSnps[vamipSnps[,"only_VAMIPs"]>0,"Gene"]), method="spearman") # This will test the gene distribution for vamip vs consensus (expected to be high)
# Do a wilcoxon test but first normalize to make them comparable
wilcox.test(prop.table(table(vamipSnps[vamipSnps[,"only_Consensus"]>0,"Gene"])),prop.table(table(vamipSnps[vamipSnps[,"only_VAMIPs"]>0,"Gene"])),paired=TRUE)

vamipSnps <- cbind(vamipSnps,"Location"=sub(":.*", "", rownames(vamipSnps)))

# Now, split snps and indels
fill_missing_items <- function(mat,vector){ # Gets a table and a vector with all rows that should be present, outputs an expanded row collection with missing dates for complete calendar in that range. rownames should have date format as %Y-%m-%d
	xtable <- as.data.frame(matrix(0,nrow=length(vector),ncol=ncol(mat)), stringsAsFactors = FALSE) # Create empty vessel for output
	rownames(xtable) <- vector # use the vector as rownames
	colnames(xtable) <- colnames(mat) # and inherit the names
	invisible(sapply(rownames(mat),function(x) {xtable[x,] <<- mat[x,]}))	# append the original values in the corresponding places (write to higher env variable
	return(xtable)
}
all_pos <- sprintf("%05d", 1:29903) # Add a vector for all available positions
# Start with indels
vamipSnps_dels <- vamipSnps[grep(">\\-", rownames(vamipSnps)),]
vamipSnps_dels <- rowsum(vamipSnps_dels[,2:3], group=vamipSnps_dels[,"Location"])
vamipSnps_dels <- fill_missing_items(vamipSnps_dels,all_pos) # create a full table containing 0s for missing items, then rename the columns
colnames(vamipSnps_dels) <- paste0(colnames(vamipSnps_dels),"_del")
# Then insertions
vamipSnps_ins <- vamipSnps[grep(">\\+", rownames(vamipSnps)),]
vamipSnps_ins <- rowsum(vamipSnps_ins[,2:3], group=vamipSnps_ins[,"Location"])
vamipSnps_ins <- fill_missing_items(vamipSnps_ins,all_pos)
colnames(vamipSnps_ins) <- paste0(colnames(vamipSnps_ins),"_ins")
# Now get the SNPs
temp <- vamipSnps[grep(">\\-", rownames(vamipSnps),invert=TRUE),]
vamipSnps_snps <- temp[grep(">\\+", rownames(temp),invert=TRUE),]
vamipSnps_snps <- rowsum(vamipSnps_snps[,2:3], group=vamipSnps_snps[,"Location"])
vamipSnps_snps <- fill_missing_items(vamipSnps_snps,all_pos)
colnames(vamipSnps_snps) <- paste0(colnames(vamipSnps_snps),"_SNP")
# Merge them
vamipSnps_all <- cbind(vamipSnps_snps,vamipSnps_dels, vamipSnps_ins)
# I'd prefer to have Consensus first
vamipSnps_all <- vamipSnps_all[,c("only_Consensus_SNP", "only_Consensus_del", "only_Consensus_ins", "only_VAMIPs_del", "only_VAMIPs_ins", "only_VAMIPs_SNP")]
# vamipSnps_all <- log10(vamipSnps_all)
# vamipSnps_all[is.infinite(as.matrix(vamipSnps_all))] <- 0
pdf("07_metadata/Mutations_Genome_VAMIPs_vs_Consensus.pdf",width=18, height=10)
par(oma = c(1.5, 0, 1, 6)) # This is just a creative fix to plot the legend outside
par(mar=c(5.1, 7.1, 4.1, 2.1))
ltys <- c(1,2,3,2,3,1)
cols <- c("turquoise","pink", "chartreuse3", "violetred", "forestgreen", "blue")
matplot(vamipSnps_all[,1:6],type='l', col=cols,lty=ltys,lwd=2, xlim=c(1,29903),ylim=c(0,30000), cex.lab=1.5, cex.main=1.5, main="Genomic Distribution of Mutations (IPMAV and Consensus)", xaxt='n', yaxt='n',frame.plot=FALSE, xlab="Genomic Position (Nt)", ylab="")
axis(1, at=seq(1,31000,2500), labels=seq(0,31000,2500),cex.axis=1)
axis(2, las=1, at=seq(0,30000,2500),cex.axis=1.2)
title(ylab="Samples with mutation in this position", line=5, cex.lab=2)
Genome <- diff(c(1, 266, 13468, 21556, 21563, 25385, 25393, 26221, 26245, 26473, 26523, 27192, 27202, 27388, 27394, 27756, 27888, 27894, 28260, 28274, 29534, 29558, 29675, 29903)) # This vector holds each interval (all ORFs). Endings were adjusted with +1 for calculations
cols_bar = c("01 5-UTR"="cornflowerblue", "02 ORF1a"="coral1","03 ORF1b"="turquoise3", "Sep1"=NA, "04 S"="chartreuse2", "Sep2"=NA, "05 ORF3a"="purple", "Sep3"=NA, "06 E"="firebrick", "Sep4"=NA, "07 M"="gold2", "Sep5"=NA, "08 ORF6"="hotpink", "Sep6"=NA, "09 ORF7a"="forestgreen", "10 ORF7b"="darkblue", "Sep7"=NA, "11 ORF8"="darkslategray", "Sep8"=NA, "12 N"="darkorange2", "Sep9"=NA, "14 ORF10"="brown2", "15 3-UTR"="mediumpurple")
sub_cols <- cols_bar[!is.na(cols_bar)]
par(oma = c(6, 0, 43, 6)) # Adjust for plotting the genes
barplot(las=1,cbind(Genome), horiz = TRUE, beside = FALSE, col=cols_bar, border=NA, xaxt='n', width=30000, add=T)
par(fig = c(0, 1, 0, 1), oma = c(1.5, 0, 1, 1.5), mar = c(8.1, 0, 7.1, 0), new = TRUE)
plot(0,0,type = "n", bty = "n", xaxt = "n", yaxt = "n",xlab="", ylab="")
legend("topright",legend=c("Fixated SNPs","IPMAV SNPs","Fixated Del","IPMAV Del","Fixated Ins","IPMAV Ins"), lty=c(1,1,2,2,3,3), col=c("turquoise","blue", "pink", "violetred", "chartreuse3", "forestgreen"), title="Lines", lwd=3, bg="white")
legend("right",legend=names(sub_cols), pch=15, col=sub_cols, pt.cex=1.5, title="Genomic positions", bg="white")
dev.off()

# Get prevalence of mutations having any IPMAVs
sum(vamipSnps[,4]>0)
# [1] 14973
sum(vamipSnps[,2]>0)/nrow(vamipSnps)
# [1] 0.9587257


# Update 2023-07-29: We will now filter those in NC regions, those having less than 5 mutations and those with no vamips
vectx <- apply(df,1, function(x){sum(x[!is.na(x)]>0)}) # Get total mutations with non-zero values
dfb <- df[vectx>4,] # Of these, subset those with at least 5 items
dim(dfb)
# [1] 14973 29661
has_snps <- apply(dfb,1, function(x){sum(x[!is.na(x)]>=0.5)}) # Remove mutations with no fixed snps
dfb <- dfb[has_snps>0,]
dim(dfb)
# [1] 11096 29661
# Identify type of mutations
nt2AA <- read.table("07_metadata/MutNtxAA_complete.tsv", sep="\t",header=T, skip=0, comment.char='',quote="",fill=F,row.names=1, check.names=FALSE) # Load nt to AA name table
nt2AA[nt2AA[,"MutAA"]=="NC","MutAA"] <- paste0(nt2AA[nt2AA[,"MutAA"]=="NC","Gene"],":NC") # add a category to flag UTRs
remove <- grep("NC$",nt2AA[rownames(df),"MutAA"]) # store a vector of items to remove
dfb <- dfb[-remove,]
dim(dfb[-remove,])
# [1] 10712 29661

df_vamip <- dfb; df_vamip[df_vamip>0.5] <- NA # this will hold vamips (>0.5 are now NAs)
df_novamip <- dfb; df_novamip[df_novamip<=0.5] <- NA # this will hold vamips (>0.5 are now NAs)


# Create a screening of all mutations by month prevalence
g_month <- group_prevalence(group_map(all_Months),dfb)
g_month_vamip <- group_prevalence(group_map(all_Months),df_vamip)
g_month_novamip <- group_prevalence(group_map(all_Months),df_novamip)
sub <- head(g_month$rel[names(sort(rowSums(g_month$rel),decreasing=TRUE)),],2000) # get the 2000 most common mutations only
sub[sub==0] <- NA
sub_vamip <- head(g_month_vamip$rel[names(sort(rowSums(g_month_vamip$rel),decreasing=TRUE)),],2000)
sub_vamip[sub_vamip==0] <- NA
sub_novamip <- head(g_month_novamip$rel[names(sort(rowSums(g_month_novamip$rel),decreasing=TRUE)),],2000)
# with these, we can get an overhead view of all the most common mutations
sub_novamip[sub_novamip==0] <- NA
out_heatmap <- pheatmap(sub, cluster_cols=FALSE, cluster_rows=FALSE)
save_pheatmap_pdf(out_heatmap, "08_Analyses/Month_heatmap_2000TopMut.pdf", 9, 230)
out_heatmap <- pheatmap(sub_vamip, cluster_cols=FALSE, cluster_rows=FALSE)
save_pheatmap_pdf(out_heatmap, "08_Analyses/Month_vamip_heatmap_2000TopMut.pdf", 9, 230)
out_heatmap <- pheatmap(sub_novamip, cluster_cols=FALSE, cluster_rows=FALSE)
save_pheatmap_pdf(out_heatmap, "08_Analyses/Month_novamip_heatmap_2000TopMut.pdf", 9, 230)
# Now, use the 2000 most abundant items in sub_vamip to extract them from the novamip set
sub_novamip_same_as_sub_vamip <- g_month_novamip$rel[rownames(sub_vamip),]
sub_novamip_same_as_sub_vamip[sub_novamip_same_as_sub_vamip==0] <- NA
out_heatmap <- pheatmap(sub_novamip_same_as_sub_vamip, cluster_cols=FALSE, cluster_rows=FALSE)
save_pheatmap_pdf(out_heatmap, "08_Analyses/Month_novamip_same_set_as_vamip_heatmap_2000TopMut.pdf", 9, 230)
# Compare those with vamips and novamips:
rownames(sub_novamip_same_as_sub_vamip) <- paste(rownames(sub_novamip_same_as_sub_vamip),"X novamip") # First, rename both sets
rownames(sub_novamip_same_as_sub_vamip) <- paste(sprintf("%04d", 1:2000),rownames(sub_novamip_same_as_sub_vamip))
rownames(sub_vamip) <- paste(rownames(sub_vamip),"VAMIP")
rownames(sub_vamip) <- paste(sprintf("%04d", 1:2000),rownames(sub_vamip))
compare <- rbind(sub_vamip,sub_novamip_same_as_sub_vamip)
compare <- compare[order(rownames(compare)),]
out_heatmap <- pheatmap(compare, cluster_cols=FALSE, cluster_rows=FALSE)
save_pheatmap_pdf(out_heatmap, "08_Analyses/Month_heatmap_compare_2000TopMut.pdf", 9, 460)

save.image("Analyze_vamip_contingency_Bar_n_boxplots_chkpt3.Rdata")
load("Analyze_vamip_contingency_Bar_n_boxplots_chkpt3.Rdata")

# Extract mean frequencies per month
freq_mean_month <- group_freq_mean(group_map(all_Months),df) # Create a colection of means (frequency) per  month
freq_mean_month <- freq_mean_month[apply(freq_mean_month, 1, function(x){sum(!is.na(x))})>=5,] # keep only those appearing in 5 months (these may not be contiguous)
out_heatmap <- pheatmap(freq_mean_month, cluster_cols=FALSE, cluster_rows=FALSE) #
save_pheatmap_pdf(out_heatmap, "08_Analyses/freq_mean_month.pdf", 9, 1200)
# Repeat but with row clustering
sub <- freq_mean_month
sub[is.na(sub)]=-10
# sub <- sub[,colSums(sub)>0] # Filters columns not having any
dim(sub)
# [1] 4827   15
clust_row <- hclust(dist(sub), method = "mcquitty") # These versions create clusters (WPGMA)
clust_col <- hclust(dist(t(sub)), method = "mcquitty")
out_heatmap <- pheatmap(freq_mean_month, cluster_cols=FALSE, cluster_rows=clust_row) # Repeat but use row clustering
save_pheatmap_pdf(out_heatmap, "08_Analyses/freq_mean_month_crows.pdf", 9, 1200)

freq_q50_month <- group_freq_centile(group_map(all_Months), df, 50) # Create a colection of centiles (of the frequency) per  month (no considering NAs). In this case, we use quantile 50.
freq_q50_month <- freq_q50_month[apply(freq_q50_month, 1, function(x){sum(!is.na(x))})>=5,] # keep only those appearing in 5 months (these may not be contiguous)
out_heatmap <- pheatmap(freq_q50_month, cluster_cols=FALSE, cluster_rows=FALSE) #
save_pheatmap_pdf(out_heatmap, "08_Analyses/freq_q50_month.pdf", 9, 1200)

# Now, for variants (general)
g_var <- group_prevalence(group_map(all_mainvar),df)
g_var_vamip <- group_prevalence(group_map(all_mainvar),df_vamip)
g_var_novamip <- group_prevalence(group_map(all_mainvar),df_novamip)
sub <- head(g_var$rel[names(sort(rowSums(g_var$rel),decreasing=TRUE)),],2000) # get the 2000 most common mutations only
sub[sub==0] <- NA
sub_vamip <- head(g_var_vamip$rel[names(sort(rowSums(g_var_vamip$rel),decreasing=TRUE)),],2000)
sub_vamip[sub_vamip==0] <- NA
sub_novamip <- head(g_var_novamip$rel[names(sort(rowSums(g_var_novamip$rel),decreasing=TRUE)),],2000)
# with these, we can get an overhead view of all the most common mutations
sub_novamip[sub_novamip==0] <- NA
out_heatmap <- pheatmap(sub, cluster_cols=FALSE, cluster_rows=FALSE)
save_pheatmap_pdf(out_heatmap, "08_Analyses/Variant_heatmap_2000TopMut.pdf", 9, 230)
out_heatmap <- pheatmap(sub_vamip, cluster_cols=FALSE, cluster_rows=FALSE)
save_pheatmap_pdf(out_heatmap, "08_Analyses/Variant_vamip_heatmap_2000TopMut.pdf", 9, 230)
out_heatmap <- pheatmap(sub_novamip, cluster_cols=FALSE, cluster_rows=FALSE)
save_pheatmap_pdf(out_heatmap, "08_Analyses/Variant_novamip_heatmap_2000TopMut.pdf", 9, 230)


# Compare those with vamips and novamips:
rownames(sub_novamip_same_as_sub_vamip) <- paste(rownames(sub_novamip_same_as_sub_vamip),"X novamip") # First, rename both sets
rownames(sub_novamip_same_as_sub_vamip) <- paste(sprintf("%04d", 1:2000),rownames(sub_novamip_same_as_sub_vamip))
rownames(sub_vamip) <- paste(rownames(sub_vamip),"VAMIP")
rownames(sub_vamip) <- paste(sprintf("%04d", 1:2000),rownames(sub_vamip))
compare <- rbind(sub_vamip,sub_novamip_same_as_sub_vamip)
compare <- compare[order(rownames(compare)),]
out_heatmap <- pheatmap(compare, cluster_cols=FALSE, cluster_rows=FALSE)
save_pheatmap_pdf(out_heatmap, "08_Analyses/Variant_heatmap_compare_2000TopMut.pdf", 9, 460)
save.image("Analyze_vamip_contingency_Bar_n_boxplots_chkpt4.Rdata")

# ### LONGITUDINAL MAPS ###
load("Analyze_vamip_contingency_Bar_n_boxplots_chkpt2.Rdata")
meta[meta[,"Estatus del paciente"]=="Asymptomatic","Estatus del paciente"] <- "Amb"
meta[meta[,"Estatus del paciente"]=="Fatal","Estatus del paciente"] <- "Dec"
meta[meta[,"Estatus del paciente"]=="Released","Estatus del paciente"] <- "Hosp"
meta[meta[,"Estatus del paciente"]=="Symptomatic","Estatus del paciente"] <- "Hosp"
write.table(meta, "07_metadata/2023-07-24_Metadata_for_longitudinal_analyses.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
tot_per_mut <- apply(df,1,function(x){sum(!is.na(x))}) # Calculate the total number of observations per mutations (not considering NAs)
# Item 10029:C->T|Omicr:ORF1a:T3255I has 13443
# grep("10029", rownames(df)) # and it is found in position 2404
# We now want to create a map for tracing each mutations throughout months
longitudinal_map_months <- group_map(all_Months)
test_mutation <- grep("23604:C->A", rownames(df)) # mutation 23604:C->A|Omicr:S:P681H has 7628 observations in the table, not only from omicron but throughout months
test_mut_freqs <- df[test_mutation,]
lapply(longitudinal_map_months,function(x){temp <- test_mut_freqs[x]; temp <- temp[!is.na(temp)]; if(length(temp)==0){temp <- 0}; quantile(temp,seq(0,1,0.1))})
test <- lapply(longitudinal_map_months,function(x){temp <- test_mut_freqs[x]; temp <- temp[!is.na(temp)]; if(length(temp)==0){temp <- 0}; quantile(temp,seq(0,1,0.1))})
test <- lapply(longitudinal_map_months,function(x){temp <- test_mut_freqs[x]; temp <- temp[!is.na(temp)]; if(length(temp)==0){temp <- 0}; temp})
test <- unlist(lapply(longitudinal_map_months,function(x){temp <- test_mut_freqs[x]; temp <- temp[!is.na(temp)]; length(temp)}))
test

save.image("Analyze_vamip_contingency_Bar_n_boxplots_chkpt5.Rdata")
load("Analyze_vamip_contingency_Bar_n_boxplots_chkpt5.Rdata")

for(mut in rownames(df)){
	print(mut)
	plot_single_mutation(df[mut,], longitudinal_map_months, mut, "Month_14_fullname",14)
}

for(mut in rownames(df)){
	print(mut)
	plot_single_mutation_long(df[mut,], group_map(all_Weeks), mut, "Weeks_4",4)
}

list_ages <- group_map(all_age)
new_names <- sprintf("%03d", (as.numeric(names(list_ages))))
new_names[length(new_names)] <- "u"
# new_names <- paste0("a",new_names)
names(list_ages) <- new_names
list_ages <- list_ages[order(names(list_ages))]
for(mut in rownames(df)){
	print(mut)
	plot_single_mutation(df[mut,], list_ages, mut, "Age",2)
}

for(mut in rownames(df)){
	print(mut)
	plot_single_mutation(df[mut,], group_map(all_age_vac), mut, "Age_vac",1)
}

list_var <- group_map(all_mainvar)
names(list_var)
# [1] "519"     "Alpha"   "Delta"   "Gamma"   "Omicron" "Others"  "Recomb"
# list_var <- list_var[c(1,2,4,3,5,7,6)]
for(mut in rownames(df)){
	print(mut)
	plot_single_mutation(df[mut,], list_var, mut, "Variants",1)
}

for(mut in rownames(df)){
	print(mut)
	plot_single_mutation(df[mut,], group_map(all_sex), mut, "Sex",1)
}

for(mut in rownames(df)){
	print(mut)
	plot_single_mutation(df[mut,], group_map(all_region7), mut, "Region_7",1)
}

for(mut in rownames(df)){
	print(mut)
	plot_single_mutation(df[mut,], group_map(all_status), mut, "Status",1)
}

deltas <- group_map(all_mainvar)$Delta
deltas <- deltas[deltas > 100]
deltas <- df[,deltas]
temp <- deltas
temp[is.na(temp)] <- 0
deltas <- deltas[rowSums(temp)>0,]

omicrons <- group_map(all_mainvar)$Omicron
omicrons <- df[,omicrons]
temp <- omicrons
temp[is.na(temp)] <- 0
omicrons <- omicrons[rowSums(temp)>0,]
rm("temp")
rm("df")
save.image("checkpoint2.Rdata")

# NEW ATTEMPT: Describe the actual variant composition per mutation (min 30 observations, at least 10 obs in any category)
mut_nt2AA <- read.table("07_metadata/MutNtxAA_complete.tsv", sep="\t",header=T, skip=0, comment.char='',quote="",fill=F,row.names=1, check.names=FALSE) # Load input table
mut_nt2AA[,"New_name"] <- paste(mut_nt2AA[,"MutAA"],mut_nt2AA[,"MutType"], sep='_')

# IMPORTANT: This should have the same order of the columns in df for the following lines to work correctly
# minimetatab <- t(as.data.frame(str_split(colnames(df), "\\|"))) # Extract a mini metadata table from the names (including
minimetatab <- t(as.data.frame(sapply(colnames(df),function(x){strsplit(x,"\\|")[[1]]}))) #UPDATE 2023-05-02: After updating to R.4.2.3, the stringr, the str_split function stopped working so I changed it back to R base
rownames(minimetatab) <- NULL # row names are removed
all_mainvar <- find_items(minimetatab[,2],meta, "VOCs_two") # If it does not exist, create a map of the main variants (VOCs and 519).

save.image("Analyze_vamip_contingency_Bar_n_boxplots_chkpt6.Rdata")
load("Analyze_vamip_contingency_Bar_n_boxplots_chkpt6.Rdata")

for(mut in rownames(df)){
# 	print(mut)
	sub <- sub("\\|.*","",mut)
	# print(unlist(df[mut,]))
	sub <- paste(sub, mut_nt2AA[sub,"New_name"],sep="_")
	print(sub)
	plot_longitudinal_distro_per_category(unlist(df[mut,]), longitudinal_map_months, all_mainvar, sub, "Variants_per_month",1)
}

# Next, test with VOCs
all_mainvar <- find_items(minimetatab[,2],meta, "VOCs_two") # Next, the main variants (VOCs and 519)
list_var <- group_map(all_mainvar)
for(mut in rownames(df)){
# 	print(mut)
	sub <- sub("\\|.*","",mut)
	# print(unlist(df[mut,]))
	sub <- paste(sub, mut_nt2AA[sub,"New_name"],sep="_")
	print(sub)
	plot_longitudinal_distro_per_category(unlist(df[mut,]), list_var, all_mainvar, sub, "Variants_per_mutation",2)
}

all_age[all_age=="u"] <- NA
list_ages <- group_map(all_age)
new_names <- sprintf("%03d", (as.numeric(names(list_ages))))
new_names[length(new_names)] <- "u"
# new_names <- paste0("a",new_names)
names(list_ages) <- new_names
list_ages <- list_ages[order(names(list_ages))]
for(mut in rownames(df)){
	sub <- sub("\\|.*","",mut)
	# print(unlist(df[mut,]))
	sub <- paste(sub, mut_nt2AA[sub,"New_name"],sep="_")
	print(sub)
	plot_longitudinal_distro_per_category(unlist(df[mut,]), list_ages, all_mainvar, sub, "Variants_per_age",4)
}

# Now, test for all mutants but with weeks
longitudinal_map_weeks <- group_map(all_Weeks)
for(mut in rownames(df)){
# 	print(mut)
	sub <- sub("\\|.*","",mut)
	# print(unlist(df[mut,]))
	sub <- paste(sub, mut_nt2AA[sub,"New_name"],sep="_")
	print(sub)
	plot_longitudinal_distro_per_category(unlist(df[mut,]), longitudinal_map_weeks, all_mainvar, sub, "Variants_per_week",4)
}

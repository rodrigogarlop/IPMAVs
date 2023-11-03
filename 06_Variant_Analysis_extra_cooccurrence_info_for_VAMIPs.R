# Updated 2023-08-05 Added an exception for XBB recombinants
# Started 2023-05-09 by Rodrigo García-López for Carlos Arias's Virology Group at IBt, UNAM, Cuernavaca, Mexico as part of the CoViGen-Mex SARS-CoV-2 survillance in Mexico. Based on the full protocol to analyze the cummulative CoViGen-Mex data
# R version 4.2.3 (2023-03-15) -- "Shortstop Beagle"
# Under GNU GPLv3 license
# Disclamer: This script was written for use with specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.
# This script is used to explore the CoViGen-Mex report table and get the total items that co-occur in time

### Institute COMPARISON ###
# ### Single variable comparisons (split by variant)
# ### DATE
# Several dates are missing and the two variants have two different histories. Thus, we should be careful with the empty categories. In this case, we'll first get the actual extant dates in a 2-way date vs lineage table
### FUNCTIONS ###
remove_1item_cols <- function(mat){ # Check if any column has <2 items and remove it if it does
	mat <- mat[,apply(mat,2,function(x) length(unique(x)))>1]
	return(mat)
}
twoway_table <- function(mat){ # Cross two variables from a 2 column matrix (var1, var2), return a table of var2 x var1 with n columns depending on total var2 items
	xtable <- xtabs(rep(1,nrow(mat))~mat[,1]+mat[,2], data=mat)
	return(xtable)
}
fill_missing_items <- function(mat,vector){ # Gets a table and a vector with all rows that should be present, outputs an expanded row collection with missing dates for complete calendar in that range. rownames should have date format as %Y-%m-%d
	xtable <- as.data.frame(matrix(0,nrow=length(vector),ncol=ncol(mat)), stringsAsFactors = FALSE) # Create empty vessel for output
	rownames(xtable) <- vector # use the vector as rownames
	colnames(xtable) <- colnames(mat) # and inherit the names
	invisible(sapply(rownames(mat),function(x) {xtable[x,] <<- mat[x,]}))	# append the original values in the corresponding places (write to higher env variable
	return(xtable)
}
rolling_N_avg <- function(mat,int){ # Input should be a table with continuous data at rows and the desired interval for the mean. If an even number is provided, the average will be placed one position to the right of as there is no single middle number
	before <- after <- trunc(int/2) # initialize with same range before and after each position
	if(int%%2==0){after <- after-1} # If even, shift the upper half of the range by 1 position (the mean will be calculated for the next value next to the middle as there is no exact number in it)
	roll <- sapply((before+1):(nrow(mat)-after),function(y) {apply(mat,2, function(x) mean(as.numeric(x[(y-before):(y+after)])))})
	colnames(roll) <- row.names(mat)[(before+1):(nrow(mat)-after)]
	roll <- t(roll)
	return(roll)
}
# ### Define functions for crossing variables
date_vs_X <- function(mat,string,int) { # Input should have at least a "Date" column with format as %Y-%m-%d, a target vector of dates that should be present (passed as a string vector) and an integer for the days in the rolling average
	dateX <- twoway_table(mat[,c("Date",string)])
	dates <- seq(as.Date(rownames(dateX)[1]),as.Date(rownames(dateX)[nrow(dateX)]),by="day") # create the complete date range
	dateX <- fill_missing_items(dateX,dates) # use the predicted missing days to get the whole date spectrum (adding 0s when required)
	dateX <- rolling_N_avg(dateX,int) # smoothen with a N-day average
	return(dateX)
}
# Define functions for specific and general-purpose plotting
define_plot_scheme <- function(){
	lty <- c("solid","dashed","dotted","dotdash","longdash","twodash")
	col <- c('chartreuse3', 'cornflowerblue', 'darkgoldenrod1', 'peachpuff3','mediumorchid2', 'turquoise3', 'wheat4','slategray2',"black","coral1","aquamarine2","blue2","violetred2","palegreen3","purple3","magenta1","limegreen","darkorange2","darkgray")
	lwd <- c(1.3,2)
	pars <- expand.grid(col = col, lty = lty, lwd = lwd, stringsAsFactors = FALSE) # This will create all
	return(pars)
}
# ### MAIN ###
df <- read.table("01_data/All_Variants_noMuts.tsv",header=T, sep='\t', skip=0, comment.char="",fill=T, check.names=FALSE, stringsAsFactors = FALSE) # columns regarding mutations were removed as some EOF error arose
# names(df)[which(names(df)=="Fecha de recoleccion")] <- "Date" # In case only some labels need updating
dim(df)
# [1] 86468    56
# The original column names are in Spanish
df <- df[,c("Accession ID", "Nombre", "Fecha de recoleccion", "Fecha de envio", "Estado", "Municipio", "Genero", "Edad","E","RdRP","Estatus del paciente","Clado Nexstrain", "New_Lineage","Linaje Pangolin", "Delta", "Delta_sub", "Omicron","Omicron_sub", "Omicron_two", "Omicron_three", "Omicron_four", "Omicron_BWart", "Recomb", "Recomb_two", "VOCs", "VOCs_two", "VOCs_three", "VOCs_four", "Region_7", "Region_5", "Region_4", "Region_3", "1999_Inegi", "CEEY", "Regiones_operativas", "Om2_A", "Om2_B", "Om2_C", "Month", "Year", "Week", "Age_range_by10", "Age_range_by5","Age_range_manual", "Age_vac")]
dim(df)
# [1] 86468    45
df[,"VOCs_three"] <- sub("XBB\\..*","XBB",df[,"VOCs_three"]) # Update 2023-08-05: This was fixed as to include all XBB together (no longer separated by subvariant)
# write.table(cbind("Old_name"=colnames(df),"New_name"=colnames(df)), "01_data/dummyCols_forAnalysis.tsv", sep="\t", col.names=TRUE, quote = FALSE, row.names=FALSE)
rename <- read.table("01_data/rename_fields_v4.tsv",header=T, sep='\t', skip=0, comment.char='',fill=F, check.names=FALSE, stringsAsFactors = FALSE, row.names=1) # columns are renamed
names(df) <- rename[names(df),]
df[is.na(df)]="u" # Fill missing items
# Remove known bad registries
# df <- df[!(df[,"Delta_or_not"]=="Delta" & (as.Date(df[,"Date"]) < as.Date("2021-04-21", "%Y-%m-%d"))),] # Delta samples before May 2020
df <- remove_1item_cols(df) # if any column doesn't have at least 2 items, remove it (first "official" date was Apr 21th, 2021)
df[,grep("Age$", colnames(df))] <- fix_age(df[,grep("Age$", colnames(df))]) # Transform u in ages to NAs (by coercion while casting to numeric)
df[df=="u"]=NA # change "u"s to NAs to prevent warnings # UNCOMMENT if required
df[df=="Unk"]=NA

# UPDATE 2023-05-09: We calculate a date x all lineages (full names)
dict <- sapply(names(table(df[,"Lineage"])),function(x){df[grep(paste0(x,"$"),df[,"Lineage"])[1],"Pango"]}) # Get the correspondance obs per lineage per date
dateVsLin <- date_vs_X(df,"Lineage",1)
write.table(dateVsLin,"Raw_tables/01_date_VS_AllVariants_RAW.tsv", sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)
mat <- dateVsLin>=1

# dummy <- dateVsLin[,15:16] # Test for existing two different cases
dummy <- dateVsLin[,c("B.1.1.519","B.1.1.322","XAP","B.1.617.2.20")] # Some other difficult cases
vector <- dummy[,4]
first_and_last_date_no_outliers <- function(vector,days_cutoff=60){ # input date vector and days cutoff after which to consider an outlier)
	dates_obs <- as.Date(sort(names(vector[vector>0]))) # Extract unique dates where variant was detected
# 	if(length(dates_obs) < 2){return(c(dates_obs,dates_obs))} # If only one item is present, start=end and exit
	first_seen <- dates_obs[1] # Store the first date
	last_seen <- dates_obs[length(dates_obs)] # Store the last date
	dif <- c(0,diff(dates_obs)) # Create a vector for differences
	temp <- data.frame(cbind("dates"=as.character(dates_obs), dif))
	target <- which(dif>=days_cutoff) # Detect item number of those that might be outliers (two months or more apart)
	place <- target*100/length(dif) # get the position of these putative outliers (1%=first to 100%=last)
	down <- target[place<40] # Cases in the middle (40-60%) will be ignored, first process those in the first items
	if(length(down)>=1){ # If there are any outliers at the beggining
		first_seen <- dates_obs[down[length(down)]] # then update the new start to the last item in the list (those that came first are now ignored)
	}
	up <- target[place>=60] # Next, process those in the last items (larger than 60%)
	if(length(up)>=1){ # If there are any outliers at the beggining
		last_seen <- dates_obs[up[1]-1] # then update the new start to the item just before the first in this set (those that came afterwards are now ignored)
	}
	out <- c(as.character(first_seen), as.character(last_seen), as.character((last_seen-first_seen)+1), as.character(length(dates_obs)))
	dates_obs <- dates_obs[grep(first_seen, dates_obs):grep(last_seen, dates_obs)] # The next 3 lines are for testing the results
	dif <- c(0,diff(dates_obs))
# 	return(dif)
	return(out)
}
# test <- apply(dateVsLin, 2, first_and_last_date_no_outliers, 15) # These are for testing results
# test[!sapply(test, is.null)] # Filter null elements
# sort(names(table(df[df[,"Lineage"]=="B.1.1.222","Date"]))) # Test dates for a single variant
period <- as.data.frame(t(apply(dateVsLin, 2, first_and_last_date_no_outliers))) # Create a date map for the appearance of each variant but removing outliers
colnames(period) <- c("first_seen", "last_seen", "days_passed", "actual_days_observed") # Rename the columns
period <- cbind("Variant"=dict[rownames(period)],"Full_name"=rownames(period),period) # and append the PANGO names for variants
# first_seen <- apply(mat,2,function(x){rownames(mat)[which(x==TRUE)[1]]})
# last_seen <- apply(mat,2,function(x){rownames(mat)[rev(which(x==TRUE))[1]]})
# period <- t(rbind(Variant=dict[names(first_seen)], Full_name=names(first_seen),first_seen,last_seen))
write.table(period,"variants_period.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

out <- apply(mat, 2, function(x){colSums(x*mat)}) # Create a new matrix with how many days a variant appears with each other
order <- period[,3]; names(order) <- rownames(period)
out <- out[,names(sort(as.Date(order)))]
out <- out[names(sort(as.Date(order))),]
write.table(out,"cooccurring_variants-days-observed.tsv", sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)


period[,3] <- as.Date(period[,3])
period[,4] <- as.Date(period[,4])
# Next, calculate which ones coexisted based on their periods
intersect_Days <- function(period1, period2){ # Each period consists of a 2-items vector, including a starting and finishing date in the same string and separated by a |. These are passed as lists
	days <- seq.Date(as.Date(as.Date(period1[1][[1]])),as.Date(period1[2][[1]]),1) # Next, expand all dates in the first period
	days <- c(days,seq.Date(as.Date(period2[1][[1]]),as.Date(period2[2][[1]]),1)) # and append those in the second
	matching_days <- names(which(table(days)>1)) # Finally, keep only the ones seen twice
	return(length(matching_days))
}

matching_days_exp <- matrix("",nrow(period),nrow(period)) # Create a void matrix with all variant combinations
rownames(matching_days_exp) <- rownames(period)
colnames(matching_days_exp) <- rownames(period)
for(i in 1:nrow(period)){ # To fill it, calculate each pairwise combo and how many days are shared
	for(j in 1:nrow(period)){
		print(c(i,j))
		matching_days_exp[i,j] <- intersect_Days(period[i,3:4], period[j,3:4]) # use the previously defined function to calculate the overlapping days
	}
}
write.table(matching_days_exp,"cooccurring_variants-days-predicted.tsv", sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)

# Finally, create an extended period and recalculate the matrix (first seen will be backtracked 7 days earlier and last seen will be forwarded until 7 days later)
extended <- period # Work with a copy
extended[,"first_seen"] <- as.character(as.Date(extended[,"first_seen"])-7) # Add 7 days before the initial day
extended[,"last_seen"] <- as.character(as.Date(extended[,"last_seen"])+7) # Add 7 days after the last day
matching_days_ext <- matrix("",nrow(extended),nrow(extended)) # Create a void matrix with all variant combinations
rownames(matching_days_ext) <- rownames(extended)
colnames(matching_days_ext) <- rownames(extended)
for(i in 1:nrow(extended)){ # To fill it, calculate each pairwise combo and how many days are shared
	for(j in 1:nrow(extended)){
		print(c(i,j))
		matching_days_ext[i,j] <- intersect_Days(extended[i,3:4], extended[j,3:4]) # use the previously defined function to calculate the overlapping days
	}
}
write.table(matching_days_ext,"cooccurring_variants-days-extended.tsv", sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)

save.image("chkpt2.Rdata")

# We'll repeat this with the VOCs only
dateVsLin <- date_vs_X(df,"VOCs_three",1)
write.table(dateVsLin,"Raw_tables/01_date_VS_VOCs_two_RAW.tsv", sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)
mat <- dateVsLin>=1

period <- as.data.frame(t(apply(dateVsLin, 2, first_and_last_date_no_outliers))) # Create a date map for the appearance of each variant but removing outliers
colnames(period) <- c("first_seen", "last_seen", "days_passed", "actual_days_observed") # Rename the columns
period <- cbind("Variant"=dict[rownames(period)],"Full_name"=rownames(period),period) # and append the PANGO names for variants
write.table(period,"VOCs_period.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

out <- apply(mat, 2, function(x){colSums(x*mat)}) # Create a new matrix with how many days a variant appears with each other
order <- period[,3]; names(order) <- rownames(period)
out <- out[,names(sort(as.Date(order)))]
out <- out[names(sort(as.Date(order))),]
write.table(out,"cooccurring_VOCs-days-observed.tsv", sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)

# Next, calculate which ones coexisted based on their periods
matching_days_exp <- matrix("",nrow(period),nrow(period)) # Create a void matrix with all variant combinations
rownames(matching_days_exp) <- rownames(period)
colnames(matching_days_exp) <- rownames(period)
for(i in 1:nrow(period)){ # To fill it, calculate each pairwise combo and how many days are shared
	for(j in 1:nrow(period)){
		print(c(i,j))
		matching_days_exp[i,j] <- intersect_Days(period[i,3:4], period[j,3:4]) # use the previously defined function to calculate the overlapping days
	}
}
write.table(matching_days_exp,"cooccurring_VOCs-days-predicted.tsv", sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)

# Finally, create an extended period and recalculate the matrix (first seen will be backtracked 7 days earlier and last seen will be forwarded until 7 days later)
extended <- period # Work with a copy
extended[,"first_seen"] <- as.character(as.Date(extended[,"first_seen"])-7) # Add 7 days before the initial day
extended[,"last_seen"] <- as.character(as.Date(extended[,"last_seen"])+7) # Add 7 days after the last day
matching_days_ext <- matrix("",nrow(extended),nrow(extended)) # Create a void matrix with all variant combinations
rownames(matching_days_ext) <- rownames(extended)
colnames(matching_days_ext) <- rownames(extended)
for(i in 1:nrow(extended)){ # To fill it, calculate each pairwise combo and how many days are shared
	for(j in 1:nrow(extended)){
		print(c(i,j))
		matching_days_ext[i,j] <- intersect_Days(extended[i,3:4], extended[j,3:4]) # use the previously defined function to calculate the overlapping days
	}
}
write.table(matching_days_ext,"cooccurring_VOCs-days-extended.tsv", sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)
matching_days_extb <- matching_days_ext[,c("B.1.1.222", "B.1.1.519", "Alpha", "Gamma", "Delta", "Omicron_BA.1.x", "Omicron_BA.2.x", "Omicron_BA.4.x", "Omicron_BA.5.x", "XB", "XBB", "XAS")]
matching_days_extb <- matching_days_extb[c("B.1.1.222", "B.1.1.519", "Alpha", "Gamma", "Delta", "Omicron_BA.1.x", "Omicron_BA.2.x", "Omicron_BA.4.x", "Omicron_BA.5.x", "XB", "XBB", "XAS"),]
write.table(matching_days_extb,"cooccurring_VOCs-days-extended_sort.tsv", sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)

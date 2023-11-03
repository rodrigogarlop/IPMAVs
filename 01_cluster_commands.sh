# Updated on 2023-05-01
# Updated on 2023-04-11
# Updated on 2023-03-14
# Started on 2022-11-25
# by Rodrigo García-López for Carlos Arias's Virology Group at IBt, UNAM, Cuernavaca, Mexico as part of the CoViGen-Mex SARS-CoV-2 survillance in Mexico.
# Under GNU GPLv3 license
# This script was created to process iVar tables produced alongside consensus sequences for SARS-CoV-2 viruses. They first require estimation of Ns in each sample followed by reencoding of the genomes and filters
# Disclamer: This script was written for use with specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.
# Files have been deposited in an hdd at /share/DiskCoViGen/VAMIP in our local cluster, and were separated by batch (lote)
# data folder: /scratch/btaboada/Data/CoV2019/

# Gzip all fastas and tables
cd /share/DiskCoViGen/VAMIP/05_variants
find . -name "*.fasta" -exec gzip "{}" \;
find . -name "*_Q_Dp.ivar.tsv" -exec gzip "{}" \;

# Summarize contents
# We now need a list of all files having an ivar table
cd /share/DiskCoViGen/VAMIP/05_variants
# This will exclude all items whose names have not been processed (those flagged with string "test")
find ~+ -name "*.ivar.tsv.gz"|grep -v test >2023-03-14_ivar_tables.list # The ~+ expands the current path to full path
mkdir /share/DiskCoViGen/VAMIP/resp_2023-03-14_ivar_tables
cat 2023-03-14_ivar_tables.list| grep CIAD >2023-03-14_CIAD_only.tsv # This is latter used for a filter
# We now need the nucleotide contents
cd /share/DiskCoViGen/VAMIP/05_variants
find ~+ -name "*.ivar2.fasta.gz"|grep -v test >2023-03-14_fastas.list # List all target consensus (ignore those flagged tests)
# Using the list of files, we get which ones were the original (raw, just after joinning lanes) files
find ~+ -name "*fastq.gz"|grep -v "01_Preprocess\|02_derep\|03_assembly" >2023-03-14_raw_fastqs.list
# We can now extract out of the original files, how many sequences they had
grep -v "/L0" 2023-03-14_raw_fastqs.list|grep -v 00_Extra >2023-03-14_Names_wo_batchIDs.txt # Get those without batch identifier
cut -d "/" -f 6-7 2023-03-14_Names_wo_batchIDs.txt|sort|uniq|grep Lote # Now print the folders so that samples can be renamed

for i in $(cat 2023-03-14_raw_fastqs.list);do echo ls -lh $i;done|bash >2023-03-14_rawfastq_size.tsv
# This was carried fo R1 and R2 because there are very rare cases where they don't match in total seqs
grep _R1.fastq.gz 2023-03-14_rawfastq_size.tsv >2023-03-14_sample_sizeR1.tsv
# For the raw files, we'll create a catalog with batch, Institute and size. Id is the folio IMSS ID (e.g. 202201000977) which is appendend as the las column (8)
paste <(cut -f 2 2023-03-14_sample_sizeR1.tsv) <(sed -e 's/\t.*academico /\t/' -e 's/\//\t/g' -e 's/ /\t/' -e 's/\t\(202.*\)\(_S.*\)/\t\1\2\t\1/' 2023-03-14_sample_sizeR1.tsv) >2023-03-14_sample_ok.tsv

cp 2023-03-14_sample_ok.tsv 2023-03-14_sample_ok_target.list.tsv

# Summarize genomic nucleotide contents:
#NOTE: Avoid divisions by 0 as these produce empty rows with only the filename
for i in $(cat 2023-03-14_fastas.list);do echo 'printf "'$i'\t" >>2023-03-14_nucleotide_contents.tsv;perl -ne '\''{next if $_=~ /^>/; chomp($_); $length=length($_); $Ns=tr/[N]//;$As=tr/[A]//; $Ts=tr/[T]//; $Cs=tr/[C]//; $Gs=tr/[G]//; $ATs=tr/[AT]//; $CGs=tr/[CG]//; $nonNs=tr/[ATCG]//; $Np=$Ns*100/$length; $nonNp=$nonNs*100/$length; $Ap=$As*100/$length; $Tp=$Ts*100/$length; $Cp=$Cs*100/$length; $Gp=$Gs*100/$length; $ATp=$ATs*100/$nonNs; $CGp=$CGs*100/$nonNs; $out=sprintf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f", $length, $nonNs, $Ns, $As, $Ts, $Cs, $Gs, $ATs, $CGs, $nonNp, $Np, $Ap, $Tp, $Cp, $Gp, $ATp, $CGp); print $out}'\'' <(zcat '$i') >>2023-03-14_nucleotide_contents.tsv; printf "\n" >>2023-03-14_nucleotide_contents.tsv';done >summarize_contents.sh
nohup bash summarize_contents.sh &

# Extract only those with non-empty items
awk -F'\t' '$2!=""' 2023-03-14_nucleotide_contents.tsv >2023-03-14_nucleotide_contents_nonEmtpy.tsv
# Now, append the IMSS folio as an extra column
paste 2023-03-14_nucleotide_contents_nonEmtpy.tsv <(sed -e "s/\//\t/g" 2023-03-14_nucleotide_contents_nonEmtpy.tsv -e 's/_S.*//' -e 's/_Q_Dp.*//'|cut -f 9) >2023-03-14_nucleotide_contents_wID.tsv
paste 2023-03-14_nucleotide_contents_nonEmtpy.tsv <(sed -e "s/\//\t/g" 2023-03-14_nucleotide_contents_nonEmtpy.tsv|cut -f 9|sed -e 's/_Q_Dp.*//' -e 's/_S[0-9]*//') >2023-03-14_nucleotide_contents_wID.tsv
# Additionally, we create a list of those having >25%, >15% and >10% Ns that are targeted for removal later on
awk -F'\t' '$12>=25' 2023-03-14_nucleotide_contents_wID.tsv|cut -f 1|sed -e 's/.*03_assembly\///' -e 's/\(_Q_Dp\).*//' >2023-03-14_for_removal_gt_25perc_Ns.tsv
awk -F'\t' '$12>=15' 2023-03-14_nucleotide_contents_wID.tsv|cut -f 1|sed -e 's/.*03_assembly\///' -e 's/\(_Q_Dp\).*//' >2023-03-14_for_removal_gt_15perc_Ns.tsv
awk -F'\t' '$12>=10' 2023-03-14_nucleotide_contents_wID.tsv|cut -f 1|sed -e 's/.*03_assembly\///' -e 's/\(_Q_Dp\).*//' >2023-03-14_for_removal_gt_10perc_Ns.tsv
awk -F'\t' '$12>=5' 2023-03-14_nucleotide_contents_wID.tsv|cut -f 1|sed -e 's/.*03_assembly\///' -e 's/\(_Q_Dp\).*//' >2023-03-14_for_removal_gt_05perc_Ns.tsv
for i in $(seq 0 100);do echo "awk -F'\t' '\$12>=$i' 2023-03-14_nucleotide_contents_wID.tsv|cut -f 1|sed -e 's/.*03_assembly\///' -e 's/\(_Q_Dp\).*//'|wc -l";done|bash >temp
printf "MaxNsCutoff\tBadGenomes\n" >2023-03-14_for_removal_gt_perc_Ns_ALL_percentages.tsv
paste <(seq 0 100) temp >>2023-03-14_for_removal_gt_perc_Ns_ALL_percentages.tsv

# Now, from the list of target samples, create a sublist of files to ignore (accept all from main IMSS batches [L0X], CIAD, INC, P1000 and RETRO"
cut -f 10 2023-03-14_sample_ok_target.list.tsv|grep -v "^L0\|^CIAD\|^P1000\|^RETRO\|^INC" >2023-03-14_IGNORED_samples.list
sed -e 's/_R1.fastq.*/_Q_Dp/' -e 's/^/\//' 2023-03-14_IGNORED_samples.list|grep -Ff - 2023-03-14_nucleotide_contents_wID.tsv >2023-03-14_IGNORED_samples_nucleotide_contents.tsv
# With this, we create a list of all available genomes that are part of the surveillance batches and other target samples (samples only, no controls) # IMPORTANT NOTE: This list contains the actual working samples for the downstream analyses
sed -e 's/_R1.fastq.*/_Q_Dp/' -e 's/^/\//' 2023-03-14_IGNORED_samples.list|grep -Fvf - 2023-03-14_nucleotide_contents_wID.tsv >2023-03-14_Main_mapped-Vigilancia.list
# wc -l 2023-03-14_Main_mapped-Vigilancia.list
# 32824 2023-03-14_Main_mapped-Vigilancia.list
cut -f 1 2023-03-14_Main_mapped-Vigilancia.list|sed 's/\.ivar2.*/.ivar.tsv.gz/' >2023-03-14_filtered_ivar_tables.list

# With both the Main_mapped_-Vigilancia.list and the 2023-03-14_sample_ok_target.list.tsv tables, we can now append size information to the nucleotide contents table using 2023-03-14_nucleotide_contents_wID.tsv as xref table
printf "File\tgenome_size\tNonNs\tNs\tAs\tTs\tCs\tGs\tATs\tCGs\tnonN%%\tN%%\tA%%\tT%%\tC%%\tG%%\tAT%%\tCG%%\tGenome\tPerm\tFile_size\n" >2023-03-14_nucleotide_contents_wSize.tsv
for i in $(cat 2023-03-14_Main_mapped-Vigilancia.list);do paste <(grep -m 1 "/$i\_S" 2023-03-14_nucleotide_contents_wID.tsv) <(grep -m 1 "/$i\_S" 2023-03-14_sample_ok_target.list.tsv|cut -f 2,3|cut -d " " -f 4);done >>2023-03-14_nucleotide_contents_wSize.tsv
# To decide upon this matter, we graph some an Rscript to graph some features


# ### Summarize ivar tables ###
# IMPORTANT NOTE: we must first make sure all ivar tables are in gzip format for the following lines to work correctly
find . -name "*ivar.tsv" -exec gzip "{}" \;
cd /share/DiskCoViGen/VAMIP/05_variants
mkdir -p /scratch/rodrigog/01_Projects/VAMIP/06_VAMIP_tables/01_per_sample

# again, but this time, create a single file per ivar table:
outdir="/scratch/rodrigog/01_Projects/VAMIP/06_VAMIP_tables/01_per_sample";for i in $(cat 2023-03-14_filtered_ivar_tables.list);do base=$(basename $i);echo 'printf "'$(echo $i|sed -e "s/.*05_variants\///" -e "s/03_assembly\///" -e "s/\//\\\\t/g" -e "s/_Q_Dp.ivar.tsv.gz//" -e "s/_Q_Dp_genome.ivar.tsv.gz//")'\t" >'$outdir'/'${base%_Q_Dp.ivar.tsv.gz}'_vamip.tsv; awk '\''NR>=2 {new_var=">"$2"M"$3"N"$4"F"$11"X"$12"R"$5"L"$9" "; print new_var}'\'' <(zcat '$i')|tr . d|tr - m|tr + p|tr -d "\n" >>'$outdir/${base%_Q_Dp.ivar.tsv.gz}'_vamip.tsv;printf "\n" >>'$outdir/${base%_Q_Dp.ivar.tsv.gz}'_vamip.tsv';done >create_single_sample_vamip_tables.sh
# qsub script qsub_command_per_line.sh was created to run all lines in parallel (any bash script for one command per line)
mkdir -p $outdir/00_qsubout
qsub -R y -l h_rt=23:59:59 -l h_vmem=100M -N vamip_tables -o $outdir/00_qsubout -t 1-$(cat create_single_sample_vamip_tables.sh|wc -l) qsub_command_per_line.sh /share/DiskCoViGen/VAMIP/05_variants/create_single_sample_vamip_tables.sh
mv 01_per_sample/00_qsubout/ $outdir/.. # I decided to move the logs afterwards
# ls $outdir|grep "tsv$"|wc -l
# 32815
cd $outdir
rename "_Q_Dp_genome.ivar.tsv.gz_vamip" "vamip" * # Some items carry out an unwanted suffix
cd /share/DiskCoViGen/VAMIP/05_variants
ls $outdir|grep "tsv$" >2022-03-28_All_current_VAMIP_coded_tables.tsv
# Repeated files were also excluded

# We now filter items with too many Ns
cd /scratch/rodrigog/01_Projects/VAMIP/06_VAMIP_tables
mkdir 01_per_sample-ignored 01_per_sample_gt_25perc_Ns 01_per_sample_gt_15perc_Ns 01_per_sample_gt_10perc_Ns
# Now, filter those that were selected
curdir="/share/DiskCoViGen/VAMIP/05_variants"
for i in $(cat $curdir/2023-03-14_bad_vamip.list);do echo 'mv 01_per_sample/'$i' 01_per_sample-ignored';done|bash
# And those with more than 25%, 15% and 10% Ns # I kept the errors on purpose to track which are no longer present
for i in $(cat $curdir/2023-03-14_for_removal_gt_25perc_Ns.tsv);do echo 'mv 01_per_sample/'$i'_vamip.tsv 01_per_sample_gt_25perc_Ns';done|bash 2>gt25_no_longer_present.tsv # Use the first list (>25%Ns to filter items with the lowest quality)
for i in $(cat $curdir/2023-03-14_for_removal_gt_25perc_Ns.tsv|grep -wvFf - $curdir/2023-03-14_for_removal_gt_15perc_Ns.tsv);do echo 'mv 01_per_sample/'$i'_vamip.tsv 01_per_sample_gt_15perc_Ns';done|bash 2>gt15_no_longer_present.tsv # use the remainder for the second tier (gt >15%Ns)
for i in $(cat $curdir/2023-03-14_for_removal_gt_15perc_Ns.tsv|grep -wvFf - $curdir/2023-03-14_for_removal_gt_10perc_Ns.tsv);do echo 'mv 01_per_sample/'$i'_vamip.tsv 01_per_sample_gt_10perc_Ns';done|bash 2>gt10_no_longer_present.tsv # Finally, do the same for those with >10%Ns
# IMPORTANT NOTE: Only genomes bearing <10% Ns were kept for dowstream analyses

# Now survey each position
cd /scratch/rodrigog/01_Projects/VAMIP/06_VAMIP_tables/
mkdir -p /scratch/rodrigog/01_Projects/VAMIP/06_VAMIP_tables/02_mut_by_pos
for i in {1..30000};do echo 'grep -m 1 -o ">'$i'M\w*" 01_per_sample/*|sed -e "s/.*01_per_sample\///" -e "s/_R2_/_/" -e "s/_vamip.tsv:>/\t/" >02_mut_by_pos/pos'$i.tsv'';done >mutations_by_position.sh

# For use with a SGE queue system:
qsub -R y -l h_rt=23:59:59 -l h_vmem=100M -pe thread 1 -tc 60 -N SNP_positions -o /scratch/rodrigog/01_Projects/VAMIP/06_VAMIP_tables/00_qsubout -t 1-$(cat mutations_by_position.sh|wc -l) ~/bin/qsub_command_per_line.sh /scratch/rodrigog/01_Projects/VAMIP/06_VAMIP_tables/mutations_by_position.sh /scratch/rodrigog/01_Projects/VAMIP/06_VAMIP_tables/

# Alt for local run:
cd /home/rod/Documents/01_Projects/SARS/VAMIP/06_VAMIP_tables/
split -a 2 -d -l 5000 --additional-suffix .sh mutations_by_position.sh run_parallel
for i in run_parallel*;do echo nohup bash $i \&\> nohup_$i.out \&;done
# UPDATE 2022-04-29: The largest genome size in the current set is 29,909, so we delete non-extant items
for i in {29910..30000};do echo rm 02_mut_by_pos/pos$i.tsv;done|bash
# the refence genome, Wuhan-Hu-1 (https://www.ncbi.nlm.nih.gov/nuccore/1798174254) is 29,903 nt long but the last item having any mutation is pos29901.tsv, so the total genome length will be considered as the maximum genome size for percentage calculations.
wc -l 02_mut_by_pos/*|tr -s " "|sort -h|sed -e 's/^ //' -e 's/ /\t/' -e 's/02_mut_by_pos\/pos//' -e 's/.tsv//' >total_mutations_by_position.tsv

# Remove those with no mutations
find 02_mut_by_pos/ -type f -empty -delete
ls 02_mut_by_pos/|wc -l
# 29774
# As of 2022-03-14, there were 28,684 (95.92%) positions that report any mutation
# # UPDATE 2023-04-11: These are now 29,774 (96.34%)
awk -F '\t' '$1==0' total_mutations_by_position.tsv |wc -l
# Only 135 (3.66%) positions had no evidence of mutations
# These empty files are thus removed # This was changed to remove them entirely
# Now remove those with only 1 observation (only one genome has them
mkdir 02_mut_by_pos/singletons
for i in $(wc -l 02_mut_by_pos/*.tsv|tr -s " "|sort -h|sed -e 's/^ //' -e 's/ /\t/'|awk '$1==1'|cut -f 2);do echo mv $i 02_mut_by_pos/singletons/;done|bash
ls 02_mut_by_pos/singletons/|wc -l
# 2,341 (7.83%) items were singletons # UPDATE 2023-04-11: 414 (1.38%)
ls 02_mut_by_pos/*.tsv|wc -l
# There are 26,343 (88.09%) positions with more than 1 observation # UPDATE 2023-04-11: 29360
awk -F '\t' '$1<=2' total_mutations_by_position.tsv |wc -l
# from the remaining items, 2,881 (9.63%) had at least 2 observations (IMPORTANT NOTE: are filtered below). # UPDATE 2023-04-11: 1212
awk -F '\t' '$1>=5' total_mutations_by_position.tsv |wc -l
# and 18,083 (60.47%) had at least 5 observations # UPDATE 2023-04-11: 26650
awk -F '\t' '$1<5' total_mutations_by_position.tsv |wc -l
# whereas 11,827 (39.55%) had less than 5 (and will not be considered for downstream analyses) # UPDATE 2023-04-11: 3260
mkdir 02_mut_by_pos/items_2-4 # I further filtered it with to those not in at least 5 observations
for i in $(wc -l 02_mut_by_pos/*.tsv|tr -s " "|sort -h|sed -e 's/^ //' -e 's/ /\t/'|awk '$1<5'|cut -f 2);do echo mv $i 02_mut_by_pos/items_2-4;done|bash
ls 02_mut_by_pos/*.tsv|wc -l
# UPDATE 2023-05-01: 26654

# Now, create a global table and decode signs
cat 02_mut_by_pos/*.tsv >temp; paste <(cut -f 1 temp) <(cut -f 2 temp|sed -e 's/d/\./' -e 's/p/+/' -e 's/m/-/' -e 's/[MNFXRL]/\t/g') >All_mut_NoSingl.tsv;rm temp
wc -l All_mut_NoSingl.tsv
# 1,212,141 mutations will be considered for the current set. # UPDATE 2023-04-11: 2627183

# This will be used to parse mutations with R (commands found locally @ /home/rod/Documents/01_Projects/SARS/VAMIP/explore_vamips.R)). Two extra tables are required:
# 1.- xref table for nt mutations to aa mutations @ UNIQUE_def_mutations_DeltaOmicron.tsv
# 2.- xref table for sample names to lineage and other tags @ 01_data/IMPORT_All_found_xref_with_names_fixed.tsv

# Get AA into ivar tables with iVar:
in_dir="/share/DiskCoViGen/VAMIP/05_variants/05_ivar_tables/bam"
out_dir="/share/DiskCoViGen/VAMIP/05_variants/05_ivar_tables/ivar"
mkdir -p $in_dir $out_dir
cd /share/DiskCoViGen/VAMIP/05_variants/05_ivar_tables/bam
for i in $(find /share/DiskCoViGen/VAMIP/05_variants/ -name "*.bam");do echo 'ln -s '$i' .';done|bash
ls >bams.list
cd ..
echo 'qsub -R y -l h_rt=23:59:59 -l h_vmem=1G -pe thread 2 -tc 30 -N ivar_AA -o /scratch/rodrigog/qsub/ -t 1-'$(cat bam/bams.list|wc -l)' /share/DiskCoViGen/VAMIP/Programs/ivar_variant_table_AA.sh '$in_dir' bams.list /share/DiskCoViGen/VAMIP/DB/MN908947.fasta /share/DiskCoViGen/VAMIP/DB/Wuhan-Hu-1_MN908947_fix.gff3 '$out_dir''|bash

echo 'qsub -R y -l h_rt=23:59:59 -l h_vmem=1G -pe thread 2 -tc 30 -N ivar_AA -o /scratch/rodrigog/qsub/ -t 1-'$(cat bam/bams_INC.list|wc -l)' /share/DiskCoViGen/VAMIP/Programs/ivar_variant_table_AA.sh '$in_dir' bams_INC.list /share/DiskCoViGen/VAMIP/DB/MN908947.fasta /share/DiskCoViGen/VAMIP/DB/Wuhan-Hu-1_MN908947_fix.gff3 '$out_dir''|bash

cd /run/media/rod/Bicho/Tables_ivar/ivar/
tar -cvf ivar.tar ivar/ # This was then copied to my local host
gzip *.tsv; cd ..
for i in $(ls ivar/*.tsv.gz);do base=$(basename $i);echo 'zcat '$i'|sed "s/^/'$base'\t/"';done|bash > single_ivar_table.tsv

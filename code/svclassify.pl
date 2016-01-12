#!/usr/bin/perl

####
#### Usage: perl svclassify.pl [flanking_size] [coverage_cutoff] [sv_size_cutoff] [time_cutoff] [median_deletion_cutoff] [median_insertion_cutoff] [Is_Moleculo_data] [Is_PacBio_data] -b <BAM> -f <Genome.fa> -sv <BED> -hom <BED> -het <BED> -sine <BED> -simple <BED> -o <TAB>
#### Created as on 05-07-2015 by Hemang Parikh, Ph.D.
#### hemang.parikh@nist.gov and parikhhemangm@gmail.com
#### Requires samtools and bedtools to be in a path
####

#To load following perl modules
use strict;
use warnings;
use List::Util qw(min max sum);
use Getopt::Long;
use Statistics::Descriptive;

#Input parameters
my (@INPUT_BAM, $INPUT_FASTA, $INPUT_REGION, $INPUT_HOMOZYGOUS, $INPUT_HETEROZYGOUS, $INPUT_REPEAT_SINE, $INPUT_REPEAT_SIMPLE, $OUTPUT_FILE, $bedtools, $samtools);
my ($head, $flanking, $cutoff, $sv_cutoff, $time_cutoff, $median_del_cutoff, $median_ins_cutoff, $is_moleculo, $is_pacbio, $verbose);

#Declare variables
GetOptions(
	'b=s{1,}' => \@INPUT_BAM,
	'f|genome_reference:s' => \$INPUT_FASTA,
	'sv|genome_region:s' => \$INPUT_REGION,
	'hom|bed_homozygous:s' => \$INPUT_HOMOZYGOUS,
	'het|bed_heterozygous:s' => \$INPUT_HETEROZYGOUS,
	'sine|repeatmasker_SINE:s' => \$INPUT_REPEAT_SINE,
	'simple|repeatmasker_Simple:s' => \$INPUT_REPEAT_SIMPLE,
	'hd|header:s' => \$head,
	's|flanking_size:s' => \$flanking,
	'c|coverage_cutoff:s' => \$cutoff,
	'sv_c|sv_size_cutoff:s' => \$sv_cutoff,
	't|sv_time_cutoff:s' => \$time_cutoff,
	'del|median_deletion_cutoff:s' => \$median_del_cutoff,
	'ins|median_insertion_cutoff:s' => \$median_ins_cutoff,
	'mol|Is_Moleculo_data:s' => \$is_moleculo,
	'pac|Is_PacBio_data:s' => \$is_pacbio,
	'o|output:s' => \$OUTPUT_FILE,
	'v' => \$verbose,
	"h|help|?" => \&usage);

#To check input files
my $INPUT_BAM = $INPUT_BAM[0];
my $INPUT_INDEX = $INPUT_BAM . ".bai";

if(defined($INPUT_BAM)){$INPUT_BAM = $INPUT_BAM} else {print usage();die "\n\nWhere is the BAM file?; We can't find BAM file.\n\n"}
if(defined($INPUT_INDEX)){$INPUT_INDEX = $INPUT_INDEX} else {print usage();die "\n\nWhere is the BAM Index file?; We can't find BAM index file.\n\n"}
if(defined($INPUT_REGION)){$INPUT_REGION = $INPUT_REGION} else {print usage();die "\n\nWhere is the genomic region file?; We can't find genomic region file.\n\n"}
if(defined($OUTPUT_FILE)){$OUTPUT_FILE = $OUTPUT_FILE} else {print usage();die "\n\nPlease specify the name of output file.\n\n"}

#To make following input files as optional
if(defined($INPUT_FASTA)){$INPUT_FASTA = $INPUT_FASTA} else {print "\n\nThe reference fasta file is not provided hence the percentage of GC content will not be calculated.\n\n"}
if(defined($INPUT_HOMOZYGOUS)){$INPUT_HOMOZYGOUS = $INPUT_HOMOZYGOUS} else {print "\n\nThe SNP homozygous file is not provided hence the number of homozygous SNP genotype calls will not be calculated.\n\n"}
if(defined($INPUT_HETEROZYGOUS)){$INPUT_HETEROZYGOUS = $INPUT_HETEROZYGOUS} else {print "\n\nThe SNP heterozygous file is not provided hence the number of heterozygous SNP genotype calls will not be calculated.\n\n"}
if(defined($INPUT_REPEAT_SINE)){$INPUT_REPEAT_SINE = $INPUT_REPEAT_SINE} else {print "\n\nThe UCSC's RepeatMasker file with short interspersed nuclear elements (SINE), long interspersed nuclear elements (LINE) and long terminal repeat elements (LTR) is not provided hence the percentage of SINE, LINE and LTR will not be calculated.\n\n"}
if(defined($INPUT_REPEAT_SIMPLE)){$INPUT_REPEAT_SIMPLE = $INPUT_REPEAT_SIMPLE} else {print "\n\nThe UCSC's RepeatMasker file with simple repeats, low complexity repeats and satellite repeats is not provided hence the percentage of simple, low complexity and satellite repeats will not be calculated.\n\n"}

#To chek for samtools and bedtools to be in the path
$bedtools=`which intersectBed` ;
if(!defined($bedtools)){die "\nError:\n\tno bedtools. Please install bedtools and add to the system path.\n";}
$samtools=`which samtools`;
if($samtools !~ /(samtools)/i){die "\nError:\n\tno samtools. Please install samtools and add to the system path.\n";}

#To define usage of the svclassify
print "Usage = perl svclassify.pl -b $INPUT_BAM -sv $INPUT_REGION ";
if(defined($INPUT_FASTA)) {print "-f $INPUT_FASTA "}; 
if(defined($INPUT_HOMOZYGOUS)) {print "-hom $INPUT_HOMOZYGOUS "}; 
if(defined($INPUT_HETEROZYGOUS)) {print "-het $INPUT_HETEROZYGOUS "};
if(defined($INPUT_REPEAT_SINE)) {print "-sine $INPUT_REPEAT_SINE "};
if(defined($INPUT_REPEAT_SIMPLE)) {print "-simple $INPUT_REPEAT_SIMPLE "};
print "-o $OUTPUT_FILE \n\n";

sub usage {
	print "\nusage: perl svclassify.pl -b <BAM> -f <Genome.fa> -sv <BED> -hom <BED> -het <BED> -sine <BED> -simple <BED> -o <TAB> \n";
	print "\t-b \t\t\tThe input BAM file. This bam file must be coordinate-sorted and have an index file in the same directory\n";
	print "\t-sv, --genome_region\t\t\tTab separated genomic regions in BED format\n";
	print "\t-f, --genome_reference \t\t\tThe input genome fasta file. This file is optional\n";
	print "\t-hom, --bed_homozygous\t\t\tHomozygous SNP genotype calls in BED format. Homozygous SNP genotype calls must be generated from the same individual as the aligned BAM file. This file is optional\n";
	print "\t-het, --bed_heterozygous\t\tHeterozygous SNP genotype calls in BED format. Heterozygous SNP genotype calls must be generated from the same individual as the aligned BAM file. This file is optional\n";
	print "\t-sine, --repeatmasker_SINE\t\tRepeatMasker file with short interspersed nuclear elements (SINE), long interspersed nuclear elements (LINE) and long terminal repeat elements (LTR) in BED format. This file is generated using the UCSC Genome Browser’s RepeatMasker Track. This file is optional\n";
	print "\t-simple, --repeatmasker_Simple\t\tRepeatMasker file with simple repeats, low complexity repeats and satellite repeats in BED format. This file is generated using the UCSC Genome Browser’s RepeatMasker Track. This file is optional\n";
	print "\t-hd, --header\t\t\tDoes the genomic region file have a header line? [Input values are either 0 or 1, Default = 1]\n";
	print "\t-s, --flanking_size\t\t\tFlanking size of genomic region [Default = 350]\n";
	print "\t-c, --coverage_cutoff\t\t\tCutoff of the coverage for categorizing SVs based on average coverage. This value should be half of the mean coverage of the BAM file [Default = 25]\n";	
	print "\t-sv_c, --sv_size_cutoff\t\t\tCutoff of the structural variant size for which annotations will not be calculated for any SV with more than this cutoff size [Default = 1000000]\n";
	print "\t-t, --sv_time_cutoff\t\t\tCutoff of the time for which structural variants annotations will not be calculated for any SV taking more than this cutoff time (in seconds) for calcuations [Default = 900]\n";
	print "\t-del, --median_deletion_cutoff\t\tCutoff of the median deletion based on random genome. This paramter is useful for PacBio data because PacBio reads have high deletion error rates, Del (the mean of deleted bases of the reads) can be normalized for PacBio by subtracting the mean Del per read length of random regions [Default = 0]\n";
	print "\t-ins, --median_insertion_cutoff\t\tCutoff of the median insertion based on random genome. This paramter is useful for PacBio data because PacBio reads have high insertion error rates, Ins (the mean of inserted bases of the reads) can be normalized for PacBio by subtracting the mean Ins per read length of random regions [Default = 0]\n";
	print "\t-mol, --Is_Moleculo_data\t\tIs data generated using Moleculo sequencing technology? [Input values are either 0 or 1, Default = 0]\n";
	print "\t-pac, --Is_PacBio_data\t\t\tIs data generated using Pacific Biosciences sequencing technology? [Input values are either 0 or 1, Default = 0]\n";
   	print "\t-o, --output\t\t\t\tOutput file name. Output file format is a TAB format\n";
	exit 1;
}

#To define flanking regions, coverage cutoff and other parameters' default values
if(defined($head)){$head = int($head)} else {$head = 1}
if(defined($flanking)){$flanking = int($flanking)} else {$flanking = 350}
if(defined($cutoff)){$cutoff = $cutoff} else {$cutoff = 25}
if(defined($sv_cutoff)){$sv_cutoff = $sv_cutoff} else {$sv_cutoff = 1000000}
if(defined($time_cutoff)){$time_cutoff = $time_cutoff} else {$time_cutoff = 900}
if(defined($median_del_cutoff)){$median_del_cutoff = $median_del_cutoff} else {$median_del_cutoff = 0}
if(defined($median_ins_cutoff)){$median_ins_cutoff = $median_ins_cutoff} else {$median_ins_cutoff = 0}
if(defined($is_moleculo)){$is_moleculo = $is_moleculo} else {$is_moleculo = 0}
if(defined($is_pacbio)){$is_pacbio = $is_pacbio} else {$is_pacbio = 0}

##################################################################
# To create temporary variable name for all the temporary files
##################################################################
srand (time ^ $$ ^ unpack "%L*", `ps axww | gzip -f`);
my $random_name = join "", map { ("a".."z")[rand 26] } 1..20;
my $command = "";

#############################################################
# To create BED and BAM files based on genomic regions
#############################################################
my $temp_bed = $random_name . "_temp_region.bed";
my $temp_large_SV = $OUTPUT_FILE . "_Not_Calculated_SV";
open (FLANKING, "<", $INPUT_REGION) or die "Can't open the genomic region file. Please check if the file exists\n\n";
open (OUT_REGION, ">", $temp_bed) or die $!;
open (OUT_INFO, ">", $OUTPUT_FILE) or die $!;
open (OUT_LARGE_SV, ">", $temp_large_SV) or die $!;

#To define numbers of columns for BED file
my $BED_Size = 0;
my @header_list = ();

print "Making BED file of the genomic and flanking regions\n";

#To remove the header line
if ($head == 1){
	my $header_line = <FLANKING>;
	chomp($header_line);
	@header_list = split("\t", $header_line);
}

while (my $line_flanking = <FLANKING>) {
	chomp($line_flanking);
	my @region_flanking = ();
	@region_flanking = split("\t", $line_flanking);

	#To get numbers of columns in the BED file
	$BED_Size = @region_flanking;
	
	if((0+@region_flanking)<3){
		die "\nError:\n\tGenomic region BED file is not in proper format. It should be tab separated values in BED format\n";
	}
	
	#To check if the start position is less than the defined flanking size
	if($region_flanking[1] <= $flanking){
		$flanking = 0;	
	}

	#To print out large SVs
	if(($region_flanking[2]-$region_flanking[1]) > $sv_cutoff){
		print OUT_LARGE_SV $line_flanking . "\n";
	}

	else{
		if(int($region_flanking[1]-$flanking) >= 0){
			if((0+@region_flanking) == 3){
				print OUT_REGION $region_flanking[0] . "\t" . ($region_flanking[1]-int($flanking)) . "\t" . ($region_flanking[2]+int($flanking)) . "\n";
			}
			else{
				for my $i (0 .. ((0+@region_flanking)-1)){
					if($i==0){
						print OUT_REGION $region_flanking[0] . "\t";
					}
					elsif($i==1){
						print OUT_REGION ($region_flanking[1]-int($flanking)) . "\t";
					}
					elsif($i==2){
						print OUT_REGION ($region_flanking[2]+int($flanking)) . "\t";
					}
					elsif($i==((0+@region_flanking)-1)){
						print OUT_REGION $region_flanking[$i] . "\n";
					}
					else{
						print OUT_REGION $region_flanking[$i] . "\t";					
					}				
				}				
			}
		}
	}
}


#To print header for output file
if ($head == 1){
	print OUT_INFO $header_list[0] . "\t". $header_list[1] . "\t". $header_list[2] . "\t". "SV_size\tSV_Cat\t";
}
else{
	print OUT_INFO "Chr\tStart\tEnd\tSV_size\tSV_Cat\t";
}

#To print header for additional input parameters
if($BED_Size > 3){
	for my $j (3 .. (($BED_Size)-1)){
		if ($head == 1){
			print OUT_INFO $header_list[$j] . "\t";
		}
		else{
			print OUT_INFO "Info_" . int($j-2) . "\t";	
		}	
	}
}

#To print entire header if data are not generated using PacBio or Moleculo technologies
if(($is_moleculo ne "1") && ($is_pacbio ne "1")){
	print OUT_INFO "L_Cov\tL_Cov_sd\tL_Cov_pro\tL_Insert\tL_Insert_sd\tL_Insert_10_percentile\tL_Insert_90_percentile\tL_Dis_unmap\tL_Dis_map\tL_Dis_all\tL_Dis_unmap_ratio\tL_Dis_map_ratio\tL_Mapping_q\tL_Mapping_q_sd\tL_Mapping_pro\tL_Mapping_10_percentile\tL_Mapping_90_percentile\tL_Soft\tL_Soft_sd\tL_Soft_pro\tL_Soft_10_percentile\tL_Soft_90_percentile\tL_Del\tL_Del_sd\tL_Del_10_percentile\tL_Del_90_percentile\tL_Ins\tL_Ins_sd\tL_Ins_10_percentile\tL_Ins_90_percentile\tL_Diff\tL_Diff_sd\tL_Diff_10_percentile\tL_Diff_90_percentile\tLM_Cov\tLM_Cov_sd\tLM_Cov_pro\tLM_Insert\tLM_Insert_sd\tLM_Insert_10_percentile\tLM_Insert_90_percentile\tLM_Dis_unmap\tLM_Dis_map\tLM_Dis_all\tLM_Dis_unmap_ratio\tLM_Dis_map_ratio\tLM_Mapping_q\tLM_Mapping_q_sd\tLM_Mapping_pro\tLM_Mapping_10_percentile\tLM_Mapping_90_percentile\tLM_Soft\tLM_Soft_sd\tLM_Soft_pro\tLM_Soft_10_percentile\tLM_Soft_90_percentile\tLM_Del\tLM_Del_sd\tLM_Del_10_percentile\tLM_Del_90_percentile\tLM_Ins\tLM_Ins_sd\tLM_Ins_10_percentile\tLM_Ins_90_percentile\tLM_Diff\tLM_Diff_sd\tLM_Diff_10_percentile\tLM_Diff_90_percentile\tM_Cov\tM_Cov_sd\tM_Cov_pro\tM_Insert\tM_Insert_sd\tM_Insert_10_percentile\tM_Insert_90_percentile\tM_Dis_unmap\tM_Dis_map\tM_Dis_all\tM_Dis_unmap_ratio\tM_Dis_map_ratio\tM_Mapping_q\tM_Mapping_q_sd\tM_Mapping_pro\tM_Mapping_10_percentile\tM_Mapping_90_percentile\tM_Soft\tM_Soft_sd\tM_Soft_pro\tM_Soft_10_percentile\tM_Soft_90_percentile\tM_Del\tM_Del_sd\tM_Del_10_percentile\tM_Del_90_percentile\tM_Ins\tM_Ins_sd\tM_Ins_10_percentile\tM_Ins_90_percentile\tM_Diff\tM_Diff_sd\tM_Diff_10_percentile\tM_Diff_90_percentile\tM_Cov_Cat\t";

	if(defined($INPUT_HOMOZYGOUS)){
		print OUT_INFO "M_Homvar\tM_Homvar_SV\t";
	}

	if(defined($INPUT_HETEROZYGOUS)){
		print OUT_INFO "M_Hetvar\tM_Hetvar_SV\t";
	}

	if(defined($INPUT_FASTA)){
		print OUT_INFO "M_GCcontent\t";
	}

	if(defined($INPUT_REPEAT_SINE)){
		print OUT_INFO "M_Sine_Line_Ltr_SV\t";
	}

	if(defined($INPUT_REPEAT_SIMPLE)){
		print OUT_INFO "M_Simple_Low_Satellite_SV\t";
	}

	print OUT_INFO "RM_Cov\tRM_Cov_sd\tRM_Cov_pro\tRM_Insert\tRM_Insert_sd\tRM_Insert_10_percentile\tRM_Insert_90_percentile\tRM_Dis_unmap\tRM_Dis_map\tRM_Dis_all\tRM_Dis_unmap_ratio\tRM_Dis_map_ratio\tRM_Mapping_q\tRM_Mapping_q_sd\tRM_Mapping_pro\tRM_Mapping_10_percentile\tRM_Mapping_90_percentile\tRM_Soft\tRM_Soft_sd\tRM_Soft_pro\tRM_Soft_10_percentile\tRM_Soft_90_percentile\tRM_Del\tRM_Del_sd\tRM_Del_10_percentile\tRM_Del_90_percentile\tRM_Ins\tRM_Ins_sd\tRM_Ins_10_percentile\tRM_Ins_90_percentile\tRM_Diff\tRM_Diff_sd\tRM_Diff_10_percentile\tRM_Diff_90_percentile\tR_Cov\tR_Cov_sd\tR_Cov_pro\tR_Insert\tR_Insert_sd\tR_Insert_10_percentile\tR_Insert_90_percentile\tR_Dis_unmap\tR_Dis_map\tR_Dis_all\tR_Dis_unmap_ratio\tR_Dis_map_ratio\tR_Mapping_q\tR_Mapping_q_sd\tR_Mapping_pro\tR_Mapping_10_percentile\tR_Mapping_90_percentile\tR_Soft\tR_Soft_sd\tR_Soft_pro\tR_Soft_10_percentile\tR_Soft_90_percentile\tR_Del\tR_Del_sd\tR_Del_10_percentile\tR_Del_90_percentile\tR_Ins\tR_Ins_sd\tR_Ins_10_percentile\tR_Ins_90_percentile\tR_Diff\tR_Diff_sd\tR_Diff_10_percentile\tR_Diff_90_percentile\n";

}

if($is_moleculo eq "1"){
	print OUT_INFO "L_Cov\tL_Cov_sd\tL_Cov_pro\tL_Mapping_q\tL_Mapping_q_sd\tL_Mapping_pro\tL_Mapping_10_percentile\tL_Mapping_90_percentile\tL_Soft\tL_Soft_sd\tL_Soft_pro\tL_Soft_10_percentile\tL_Soft_90_percentile\tL_Del\tL_Del_sd\tL_Del_10_percentile\tL_Del_90_percentile\tL_Ins\tL_Ins_sd\tL_Ins_10_percentile\tL_Ins_90_percentile\tL_Diff\tL_Diff_sd\tL_Diff_10_percentile\tL_Diff_90_percentile\tLM_Cov\tLM_Cov_sd\tLM_Cov_pro\tLM_Mapping_q\tLM_Mapping_q_sd\tLM_Mapping_pro\tLM_Mapping_10_percentile\tLM_Mapping_90_percentile\tLM_Soft\tLM_Soft_sd\tLM_Soft_pro\tLM_Soft_10_percentile\tLM_Soft_90_percentile\tLM_Del\tLM_Del_sd\tLM_Del_10_percentile\tLM_Del_90_percentile\tLM_Ins\tLM_Ins_sd\tLM_Ins_10_percentile\tLM_Ins_90_percentile\tLM_Diff\tLM_Diff_sd\tLM_Diff_10_percentile\tLM_Diff_90_percentile\tM_Cov\tM_Cov_sd\tM_Cov_pro\tM_Mapping_q\tM_Mapping_q_sd\tM_Mapping_pro\tM_Mapping_10_percentile\tM_Mapping_90_percentile\tM_Soft\tM_Soft_sd\tM_Soft_pro\tM_Soft_10_percentile\tM_Soft_90_percentile\tM_Del\tM_Del_sd\tM_Del_10_percentile\tM_Del_90_percentile\tM_Ins\tM_Ins_sd\tM_Ins_10_percentile\tM_Ins_90_percentile\tM_Diff\tM_Diff_sd\tM_Diff_10_percentile\tM_Diff_90_percentile\tM_Cov_Cat\t";

	if(defined($INPUT_HOMOZYGOUS)){
		print OUT_INFO "M_Homvar\tM_Homvar_SV\t";
	}

	if(defined($INPUT_HETEROZYGOUS)){
		print OUT_INFO "M_Hetvar\tM_Hetvar_SV\t";
	}

	if(defined($INPUT_FASTA)){
		print OUT_INFO "M_GCcontent\t";
	}

	if(defined($INPUT_REPEAT_SINE)){
		print OUT_INFO "M_Sine_Line_Ltr_SV\t";
	}

	if(defined($INPUT_REPEAT_SIMPLE)){
		print OUT_INFO "M_Simple_Low_Satellite_SV\t";
	}

	print OUT_INFO "RM_Cov\tRM_Cov_sd\tRM_Cov_pro\tRM_Mapping_q\tRM_Mapping_q_sd\tRM_Mapping_pro\tRM_Mapping_10_percentile\tRM_Mapping_90_percentile\tRM_Soft\tRM_Soft_sd\tRM_Soft_pro\tRM_Soft_10_percentile\tRM_Soft_90_percentile\tRM_Del\tRM_Del_sd\tRM_Del_10_percentile\tRM_Del_90_percentile\tRM_Ins\tRM_Ins_sd\tRM_Ins_10_percentile\tRM_Ins_90_percentile\tRM_Diff\tRM_Diff_sd\tRM_Diff_10_percentile\tRM_Diff_90_percentile\tR_Cov\tR_Cov_sd\tR_Cov_pro\tR_Mapping_q\tR_Mapping_q_sd\tR_Mapping_pro\tR_Mapping_10_percentile\tR_Mapping_90_percentile\tR_Soft\tR_Soft_sd\tR_Soft_pro\tR_Soft_10_percentile\tR_Soft_90_percentile\tR_Del\tR_Del_sd\tR_Del_10_percentile\tR_Del_90_percentile\tR_Ins\tR_Ins_sd\tR_Ins_10_percentile\tR_Ins_90_percentile\tR_Diff\tR_Diff_sd\tR_Diff_10_percentile\tR_Diff_90_percentile\n";	
}

if($is_pacbio eq "1"){
	print OUT_INFO "L_Cov\tL_Cov_sd\tL_Cov_pro\tL_Del\tL_Del_sd\tL_Del_10_percentile\tL_Del_90_percentile\tL_Ins\tL_Ins_sd\tL_Ins_10_percentile\tL_Ins_90_percentile\tL_Diff\tL_Diff_sd\tL_Diff_10_percentile\tL_Diff_90_percentile\tLM_Cov\tLM_Cov_sd\tLM_Cov_pro\tLM_Del\tLM_Del_sd\tLM_Del_10_percentile\tLM_Del_90_percentile\tLM_Ins\tLM_Ins_sd\tLM_Ins_10_percentile\tLM_Ins_90_percentile\tLM_Diff\tLM_Diff_sd\tLM_Diff_10_percentile\tLM_Diff_90_percentile\tM_Cov\tM_Cov_sd\tM_Cov_pro\tM_Del\tM_Del_sd\tM_Del_10_percentile\tM_Del_90_percentile\tM_Ins\tM_Ins_sd\tM_Ins_10_percentile\tM_Ins_90_percentile\tM_Diff\tM_Diff_sd\tM_Diff_10_percentile\tM_Diff_90_percentile\tM_Cov_Cat\t";

	if(defined($INPUT_HOMOZYGOUS)){
		print OUT_INFO "M_Homvar\tM_Homvar_SV\t";
	}

	if(defined($INPUT_HETEROZYGOUS)){
		print OUT_INFO "M_Hetvar\tM_Hetvar_SV\t";
	}

	if(defined($INPUT_FASTA)){
		print OUT_INFO "M_GCcontent\t";
	}

	if(defined($INPUT_REPEAT_SINE)){
		print OUT_INFO "M_Sine_Line_Ltr_SV\t";
	}

	if(defined($INPUT_REPEAT_SIMPLE)){
		print OUT_INFO "M_Simple_Low_Satellite_SV\t";
	}

	print OUT_INFO "RM_Cov\tRM_Cov_sd\tRM_Cov_pro\tRM_Del\tRM_Del_sd\tRM_Del_10_percentile\tRM_Del_90_percentile\tRM_Ins\tRM_Ins_sd\tRM_Ins_10_percentile\tRM_Ins_90_percentile\tRM_Diff\tRM_Diff_sd\tRM_Diff_10_percentile\tRM_Diff_90_percentile\tR_Cov\tR_Cov_sd\tR_Cov_pro\tR_Del\tR_Del_sd\tR_Del_10_percentile\tR_Del_90_percentile\tR_Ins\tR_Ins_sd\tR_Ins_10_percentile\tR_Ins_90_percentile\tR_Diff\tR_Diff_sd\tR_Diff_10_percentile\tR_Diff_90_percentile\n";
}


close FLANKING or warn $! ? "Error closing genomic regions file: $!": "Exit status $? from genomic regions file";
close OUT_REGION or warn $! ? "Error closing temporary  genomic regions file: $!": "Exit status $? from temporary genomic regions file";

#To make a BAM file of the genomic regions
my $temp_unsorted_bam = $random_name . "_temp_unsorted_bam.bam";
$command = "samtools view -b -L " . $temp_bed . " " . $INPUT_BAM . " >> " . $temp_unsorted_bam;
print "Making a BAM file of the genomic regions\n";
if($verbose){print "$command\n";}
system("$command");

#To make a sorted BAM file
my $temp_bam = $random_name . "_temp_bam";
$command = "samtools sort " . $temp_unsorted_bam . " " . $temp_bam;
print "Sorting a BAM file of the genomic regions\n";
if($verbose){print "$command\n";}
system("$command");

#To make an index file from BAM file
$command = "samtools index " . $temp_bam . ".bam";
print "Indexing a BAM file of the genomic regions\n";
if($verbose){print "$command\n";}
system("$command");

print "Calculating parameters for each structural variant\n";


###########################################################################################
# Create subroutines for
# Calculating coverage distribution for the genomic regions
# Calculating insert size distribution
# Calculating discordant read pairs
# Calculating mapping quality
# Calculating soft clipped reads information
# Calculating numbers of heterozygous and homozygous SNP genotype calls
# Calculating % of repeatmaskers' SINE, LINE, LTR, simple repeats and satellite repeats
###########################################################################################
#To store output of subroutines
my $sub_percentile_output = "";
my $sub_average_output = "";
my $sub_info_output = "";
my $sub_sv_output = "";

#To check if the SV is running longer than sv_time_cutoff
my $in_eval = 0;

sub sub_percentile{

	my ($INPUT_PERCENTILE, $percentile_value) = @_;	
	open(IN_PERCENTILE, "<", $INPUT_PERCENTILE) or warn "";
	my $count_percentile_line = 0;
	my @window = ();
	
	while(<IN_PERCENTILE>){
		$count_percentile_line++;
	}

	#Calculate the size of the sliding window
	my $remember_count = 1 + (100 - $percentile_value) * $count_percentile_line/100;
	
	if($count_percentile_line!=0){	
		seek IN_PERCENTILE, 0, 0;

		@window = sort { $a <=> $b }
		map scalar <IN_PERCENTILE>,			
		1 .. $remember_count;
		chomp @window;

		while (<IN_PERCENTILE>) {
			chomp;
			next if $_ < $window[0];
			shift @window;
			my $k = 0;
			$k++ while $k <= $#window and $window[$k] <= $_;
			splice @window, $k, 0, $_;
		}
		$sub_percentile_output = $window[0] . "\t";
	}
	else{
		$sub_percentile_output = "0" . "\t";
	}
	close IN_PERCENTILE or warn $! ? "Error closing input percentile file: $!": "Exit status $? from input percentile file";

	return $sub_percentile_output;
}


sub sub_average{
	my ($INPUT_AVERAGE, $min_value, $length_details) = @_;	
	my $depth_count = 0;
	my $proportion_value = 0;
	my $line_count = 0;
	my @average_data = ();
	my $analysis_file;

	open(IN_AVERAGE_FILE, "<", $INPUT_AVERAGE) or warn "";
	my $count_line = 0;
	while(<IN_AVERAGE_FILE>){
		$count_line++;
	}
	
	close IN_AVERAGE_FILE or warn $! ? "Error closing input file: $!": "Exit status $? from input file";
	
	if($count_line!=0){

		open (INPUT_AVERAGE_FILE, "<", $INPUT_AVERAGE) or warn "";

		if (defined $length_details) {

			if(($length_details>$count_line) && ($count_line<=15000)){ 
	      		$line_count = $length_details - $count_line;
      		}
    	}

	
		while (my $input_average_line=<INPUT_AVERAGE_FILE>) {
  		
			chomp($input_average_line);
			$line_count += 1;
			
			push(@average_data, $input_average_line);
			
			
				if (defined $min_value) {
					if ($input_average_line<=$min_value){
						$depth_count += 1;		
					}
				}
		}

		my $stat_average_data=Statistics::Descriptive::Full->new();
		$stat_average_data->add_data(@average_data);
		$sub_average_output = $stat_average_data->mean() . "\t" . $stat_average_data->standard_deviation() . "\t";
						
		if (defined $min_value) {
			if (defined $length_details) {
				if(($length_details>$count_line) && ($count_line<=15000)){
					$proportion_value = (($length_details - $count_line) + ($depth_count))/$line_count;
				}
				else{
					$proportion_value = $depth_count/$line_count;				
				}
			}
			else{	
				$proportion_value = $depth_count/$line_count;
			}
			$sub_average_output = $stat_average_data->mean() . "\t" . $stat_average_data->standard_deviation() . "\t" . $proportion_value . "\t";
		}

		else{
			$sub_average_output = $stat_average_data->mean() . "\t" . $stat_average_data->standard_deviation() . "\t";
		}
		
		close INPUT_AVERAGE_FILE or warn $! ? "Error closing input file: $!": "Exit status $? from input file";
	}
	
	else{
		
		if (defined $min_value) {
			$sub_average_output = "0" . "\t" . "0" . "\t" . "0" . "\t";
		}
		else {
			$sub_average_output = "0" . "\t" . "0" . "\t";
		}
	}

	return $sub_average_output;
}

######################################################################
#To create subroutine to generate all the parameters from a bam file
######################################################################
sub sub_info{
	
	#To get input parameters from each SV positions
	my($temp_line_info, $temp_sorted_bam) = @_;
	
	###############################################################
	#To obtain depth of coverage for each SV
	###############################################################
	my $temp_line_depth = $temp_line_info;
	my @temp_depth = split(":", $temp_line_depth);
	my $chr = $temp_depth[0];
	my ($start, $end) = split("\-",$temp_depth[1]);
	my $sv_length = int($end)-int($start);

	my $temp_depth_file = $random_name . "_temp_depth";
	my $temp_depth_small_file = $temp_depth_file . "_small";
			
	open (OUT_DEPTH_FILE, ">", $temp_depth_file) or die $!;
	
	$command="samtools depth -r " . $temp_line_info . " " . $temp_sorted_bam . " >>" . $temp_depth_file;

	system("$command");
	
		
	close OUT_DEPTH_FILE or warn $! ? "Error closing the depth file: $!": "Exit status $? from the depth file";
	
	my $ln_temp_depth_small_file = `awk \'END { print NR }\'  $temp_depth_file `;

	#To check if the depth file is too large
	if (chomp($ln_temp_depth_small_file) > 15000)
	{
		my $rand_depth_number = 10000/chomp($ln_temp_depth_small_file);
		`awk -F\'\\t\' \' BEGIN {srand()} !/^\$/ {if (rand() <= $rand_depth_number) print \$3}\' $temp_depth_file > $temp_depth_small_file `;
	}
	else
	{
		`awk -F\'\\t\' \'{print \$3}\' $temp_depth_file > $temp_depth_small_file `;
	}

	my $depth_average = &sub_average($temp_depth_small_file, 5, $sv_length);
	
	$sub_info_output = $sub_info_output . $depth_average;

	
	###############################################################
	#To obtain insert size for each SV
	###############################################################
	my $temp_insert_file = $random_name . "_temp_insert";
	my $temp_insert_small_file = $temp_insert_file . "_small";
	my $temp_insert_mini_file = $temp_insert_file . "_mini";
		
	open (OUT_INSERT_FILE, ">", $temp_insert_file) or warn "";
	
	if (($sv_length>100000) && ($is_moleculo ne "1") && ($is_pacbio ne "1"))
	{
		$command="samtools view -f2 -s 0.1 " . $temp_sorted_bam . " " . $temp_line_info . " >>" . $temp_insert_file;
		system("$command");
	}
	else
	{
		$command="samtools view -f2 " . $temp_sorted_bam . " " . $temp_line_info . " >>" . $temp_insert_file;
		system("$command");
	}
	close OUT_INSERT_FILE or warn $! ? "Error closing the insert file: $!": "Exit status $? from the insert file";
	
	my $ln_temp_insert_small_file = `awk \'END { print NR }\'  $temp_insert_file `;
	if (chomp($ln_temp_insert_small_file)>15000)
	{
		my $rand_insert_number = 10000/chomp($ln_temp_insert_small_file);
		`awk -F\'\\t\' \' BEGIN {srand()} !/^\$/ {if (rand() <= $rand_insert_number) print \$9}\' $temp_insert_file > $temp_insert_mini_file `;
		`awk -F\'\\t\' \'{if (\$0<0){\$0=-\$0}else{\$0=\$0}print \$0}\' $temp_insert_mini_file > $temp_insert_small_file `;
	}
	else
	{
		`awk -F\'\\t\' \'{if (\$9<0){\$9=-\$9}else{\$9=\$9}print \$9}\' $temp_insert_file > $temp_insert_small_file `;
	}

	#To skip insert size calculations if data are generated using PacBio or Moleculo technologies
	if(($is_moleculo ne "1") && ($is_pacbio ne "1")){
		my $insert_average = &sub_average($temp_insert_small_file);
		my $insert_10_percentile = &sub_percentile($temp_insert_small_file, 10);
		my $insert_90_percentile = &sub_percentile($temp_insert_small_file, 90);

		$sub_info_output = $sub_info_output . $insert_average . $insert_10_percentile . $insert_90_percentile;
	}
		
	###############################################################

	###############################################################
	#To obtain discordant read counts for each SV
	###############################################################
	
	my $temp_discordant_unmap = $random_name . "_temp_discordant_unmap";
		
	open (OUT_DISCORDANT_UNMAP, ">", $temp_discordant_unmap) or warn "";

	if (($sv_length>100000) && ($is_moleculo ne "1") && ($is_pacbio ne "1"))
	{
		$command = "samtools view -f9 -F 1792 -s 0.1 " . $temp_sorted_bam . " " . $temp_line_info . " >>" . $temp_discordant_unmap;
		system("$command");
	}
	else
	{
		$command = "samtools view -f9 -F 1792 " . $temp_sorted_bam . " " . $temp_line_info . " >>" . $temp_discordant_unmap;
		system("$command");
	}
	
	
	close OUT_DISCORDANT_UNMAP or warn $! ? "Error closing the discordant unmap file: $!": "Exit status $? from the discordant unmap file";

	open(IN_DISCORDANT_UNMAP, "<", $temp_discordant_unmap) or warn "";
	my $discordant_unmap = 0;
	while(<IN_DISCORDANT_UNMAP>){
		$discordant_unmap++;
	}
	close IN_DISCORDANT_UNMAP or warn $! ? "Error closing the discordant unmap file: $!": "Exit status $? from the discordant unmap file";
	
	my $temp_discordant_map = $random_name . "_temp_discordant_map";
		
	open (OUT_DISCORDANT_MAP, ">", $temp_discordant_map) or warn "";

	if (($sv_length>100000) && ($is_moleculo ne "1") && ($is_pacbio ne "1"))
	{
		$command = "samtools view -f1 -F 1802 -s 0.1 " . $temp_sorted_bam . " " . $temp_line_info . " >>" . $temp_discordant_map;
		system("$command");
	}
	else
	{
		$command = "samtools view -f1 -F 1802 " . $temp_sorted_bam . " " . $temp_line_info . " >>" . $temp_discordant_map;
		system("$command");
	}
	
	close OUT_DISCORDANT_MAP or warn $! ? "Error closing the discordant map file: $!": "Exit status $? from the discordant map file";

	open(IN_DISCORDANT_MAP, "<", $temp_discordant_map) or warn "";
	my $discordant_map = 0;
	while(<IN_DISCORDANT_MAP>){
		$discordant_map++;
	}
	close IN_DISCORDANT_MAP or warn $! ? "Error closing the discordant map file: $!": "Exit status $? from the discordant map file";
	
	my $temp_discordant_all = $random_name . "_temp_discordant_all";
		
	open (OUT_DISCORDANT_ALL, ">", $temp_discordant_all) or warn "";

	if (($sv_length>100000) && ($is_moleculo ne "1") && ($is_pacbio ne "1"))
	{
		$command = "samtools view -f2 -s 0.1 " . $temp_sorted_bam . " " . $temp_line_info . " >>" . $temp_discordant_all;
		system("$command");
	}
	else
	{
		$command = "samtools view -f2 " . $temp_sorted_bam . " " . $temp_line_info . " >>" . $temp_discordant_all;
		system("$command");
	}
	
	close OUT_DISCORDANT_ALL or warn $! ? "Error closing the discordant all file: $!": "Exit status $? from the discordant all file";

	open(IN_DISCORDANT_ALL, "<", $temp_discordant_all) or warn "";
	my $discordant_all = 0;
	while(<IN_DISCORDANT_ALL>){
		$discordant_all++;
	}
	close IN_DISCORDANT_ALL or warn $! ? "Error closing the discordant all file: $!": "Exit status $? from the discordant all file";
	
	
	if(int($discordant_unmap)!=0){	
		chomp($discordant_unmap);
	}
	else{
		$discordant_unmap = 0;
	}
	if(int($discordant_map)!=0){	
		chomp($discordant_map);
	}
	else{
		$discordant_map = 0;
	}
	
	if(int($discordant_all)!=0){	
		chomp($discordant_all);
	}
	else{
		$discordant_all = 0;
	}
	
	my $discordant_unmap_ratio = 0;
	my $discordant_map_ratio = 0;
	
	#To get ratio to total discordant read pairs
	if(int($discordant_all)!=0)
	{
		$discordant_unmap_ratio=int($discordant_unmap)/	int($discordant_all);
		$discordant_map_ratio=int($discordant_map)/int($discordant_all);
				
		#To skip discordant details if data are generated using PacBio or Moleculo technologies
		if(($is_moleculo ne "1") && ($is_pacbio ne "1")){
			$sub_info_output = $sub_info_output . $discordant_unmap . "\t" . $discordant_map . "\t" . $discordant_all . "\t" . $discordant_unmap_ratio . "\t" . $discordant_map_ratio . "\t";
		}
	}	
	else
	{
		#To skip discordant details if data are generated using PacBio or Moleculo technologies
		if(($is_moleculo ne "1") && ($is_pacbio ne "1")){
			$sub_info_output = $sub_info_output . $discordant_unmap . "\t" . $discordant_map . "\t0\t" . $discordant_unmap_ratio . "\t" . $discordant_map_ratio . "\t";
		}
	}
		
	my $temp_quality_cigar_file = $random_name . "_temp_quality_cigar";
	my $temp_quality_small_file = $temp_quality_cigar_file . "_small";
	my $temp_cigar_small_file = $temp_quality_cigar_file . "cigar_small";
			
	open (OUT_QUALITY_CIGAR_FILE, ">", $temp_quality_cigar_file) or warn "";
	
	if (($sv_length>100000) && ($is_moleculo ne "1") && ($is_pacbio ne "1"))
	{
		$command = "samtools view -s 0.1 " . $temp_sorted_bam . " " . $temp_line_info . " >>" . $temp_quality_cigar_file;
		system("$command");
	}
	else
	{
		$command = "samtools view " . $temp_sorted_bam . " " . $temp_line_info . " >>" . $temp_quality_cigar_file;
		system("$command");
	}
	
	close OUT_QUALITY_CIGAR_FILE or warn $! ? "Error closing the mapping quality and CIGAR files: $!": "Exit status $? from the mapping quality and CIGAR files";

	my $ln_temp_quality_small_file = `awk \'END { print NR }\'  $temp_quality_cigar_file `;
	if (chomp($ln_temp_quality_small_file)>15000)
	{
		my $rand_quality_number=10000/chomp($ln_temp_quality_small_file);
		`awk -F\'\\t\' \' BEGIN {srand()} !/^\$/ {if (rand() <= $rand_quality_number) print \$5}\' $temp_quality_cigar_file > $temp_quality_small_file `;
		
	}
	else
	{
		`awk -F\'\\t\' \'{print \$5}\' $temp_quality_cigar_file > $temp_quality_small_file `;
	}
	

	
	#To skip mapping quality calculations if data are generated using PacBio
	if($is_pacbio ne "1"){
		my $mapping_average = &sub_average($temp_quality_small_file, 0);
		my $mapping_10_percentile = &sub_percentile($temp_quality_small_file, 10);
		my $mapping_90_percentile = &sub_percentile($temp_quality_small_file, 90);

		$sub_info_output = $sub_info_output . $mapping_average . $mapping_10_percentile . $mapping_90_percentile;
		
		
	}
	
	my $ln_temp_cigar_small_file = `awk \'END { print NR }\'  $temp_quality_cigar_file `;
	if (chomp($ln_temp_cigar_small_file)>15000)
	{
		my $rand_cigar_number=10000/chomp($ln_temp_cigar_small_file);
		`awk -F\'\\t\' \' BEGIN {srand()} !/^\$/ {if (rand() <= $rand_cigar_number) print \$6}\' $temp_quality_cigar_file > $temp_cigar_small_file `;
	}
	else
	{
		`awk -F\'\\t\' \'{print \$6}\' $temp_quality_cigar_file > $temp_cigar_small_file `;
	}
		
	
	
	###############################################################
	#To obtain cigar values for each SV
	###############################################################
 	my @soft_cigar=();
	my @cigar_del=();
	my @cigar_ins=();
	my @cigar_diff=();
	
	my @cigar=();
	
	my $temp_cigar_file=$random_name . "_temp_cigar";
	my $temp_soft_cigar_file=$random_name . "_temp_soft_cigar";
	my $temp_del_file=$random_name . "_temp_del";
	my $temp_ins_file=$random_name . "_temp_ins";
	my $temp_diff_file=$random_name . "_temp_diff";
	
	open (OUT_CIGAR_FILE, ">", $temp_cigar_file) or warn "";
	open (OUT_SOFT_CIGAR_FILE, ">", $temp_soft_cigar_file) or warn "";
	open (OUT_DEL_FILE, ">", $temp_del_file) or warn "";
	open (OUT_INS_FILE, ">", $temp_ins_file) or warn "";
	open (OUT_DIFF_FILE, ">", $temp_diff_file) or warn "";

	
	open(IN_CIGAR_FILE, "<", $temp_cigar_small_file) or warn "";
	my $cigar_line = 0;
	while(<IN_CIGAR_FILE>){
		$cigar_line++;
	}
	close IN_CIGAR_FILE or warn $! ? "Error closing the mapping quality and CIGAR files: $!": "Exit status $? from the mapping quality and CIGAR files";
	
	if(int($cigar_line)==0){
		#To skip soft clipped calculations if data are generated using PacBio
		if($is_pacbio ne "1"){
			$sub_info_output = $sub_info_output . "0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t";
		}
		else{
			$sub_info_output = $sub_info_output . "0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t";		
		}
		
	}
	else{
		
		open(INPUT_CIGAR_FILE, "<", $temp_cigar_small_file) or die $!;
		while(my $row_cigar=<INPUT_CIGAR_FILE>){
			chomp $row_cigar;

			print OUT_CIGAR_FILE $row_cigar . "\n";	

			my @temp_cigar=();
			my @temp_cigar_del=();
			my @temp_cigar_ins=();
			my @temp_cigar_diff=();
			my @read_length=();
			my $sum_cigar_del=0;
			my $sum_cigar_ins=0;
			my $sum_cigar_diff=0;
			
			#To get length of each read
			while($row_cigar=~m/(\d+)([MI=X])/g){
				push(@read_length, $1);
			}
			
			while($row_cigar=~m/(\d+)S/g){
				push(@temp_cigar, $1);
			}
			if (@temp_cigar){
				my $max_cigar = max(@temp_cigar);
				printf OUT_SOFT_CIGAR_FILE "%f\n", $max_cigar;
			}

			while($row_cigar=~m/(\d+)D/g){
				push(@temp_cigar_del, $1);				
			}
			while($row_cigar=~m/(\d+)I/g){
				push(@temp_cigar_ins, $1);				
			}
			
			if (int(0+@read_length)==0){
				();			
			} 			
			else{			
				if (@temp_cigar_del){
					$sum_cigar_del = (sum(@temp_cigar_del))-($median_del_cutoff * sum(@read_length));
					printf OUT_DEL_FILE "%f\n", $sum_cigar_del;	
				}
				if (@temp_cigar_ins){
					$sum_cigar_ins = (sum(@temp_cigar_ins))-($median_ins_cutoff * sum(@read_length));
					printf OUT_INS_FILE "%f\n", $sum_cigar_ins;	
				}
						
				if (@temp_cigar_ins){
					if(@temp_cigar_del){
						$sum_cigar_diff = ($sum_cigar_ins - $sum_cigar_del);
						printf OUT_DIFF_FILE "%f\n", $sum_cigar_diff;		
					}
					else{
						printf OUT_DIFF_FILE "%f\n", $sum_cigar_ins;
					}
				}
				else{
					if(@temp_cigar_del){
						printf OUT_DIFF_FILE "%f\n", -($sum_cigar_del);
					}
					else{
						printf OUT_DIFF_FILE "%f\n", $sum_cigar_diff;
					}			
				}		
		}			
	}
	


		close OUT_CIGAR_FILE or warn $! ? "Error closing the cigar file: $!": "Exit status $? from the cigar file";
		close OUT_SOFT_CIGAR_FILE or warn $! ? "Error closing the soft cigar file: $!": "Exit status $? from the soft cigar file";
		close OUT_DEL_FILE or warn $! ? "Error closing the deletion file: $!": "Exit status $? from the deletion file";
		close OUT_INS_FILE or warn $! ? "Error closing the insertion file: $!": "Exit status $? from the insertion file";
		close OUT_DIFF_FILE or warn $! ? "Error closing the difference file: $!": "Exit status $? from the difference file";
		
		open(INPUT_CIGAR_FILE, "<", $temp_cigar_file) or die $!;
		my $count_cigar_line = 0;
		while(<INPUT_CIGAR_FILE>){
			$count_cigar_line++;
		}
		close INPUT_CIGAR_FILE or warn $! ? "Error closing the cigar file: $!": "Exit status $? from the cigar file";
	
		open(INPUT_SOFT_CIGAR_FILE, "<", $temp_soft_cigar_file) or die $!;
		my $count_soft_cigar_line = 0;
		while(<INPUT_SOFT_CIGAR_FILE>){
			$count_soft_cigar_line++;
		}
		close INPUT_SOFT_CIGAR_FILE or warn $! ? "Error closing the soft cigar file: $!": "Exit status $? from the soft cigar file";

		open(INPUT_DEL_FILE, "<", $temp_del_file) or die $!;
		my $count_del_line = 0;
		while(<INPUT_DEL_FILE>){
			$count_del_line++;
		}
		close INPUT_DEL_FILE or warn $! ? "Error closing the deletion file: $!": "Exit status $? from the deletion file";

		open(INPUT_INS_FILE, "<", $temp_ins_file) or die $!;
		my $count_ins_line = 0;
		while(<INPUT_INS_FILE>){
			$count_ins_line++;
		}
		close INPUT_INS_FILE or warn $! ? "Error closing the insertion file: $!": "Exit status $? from the insertion file";

		open(INPUT_DIFF_FILE, "<", $temp_diff_file) or die $!;
		my $count_diff_line = 0;
		while(<INPUT_DIFF_FILE>){
			$count_diff_line++;
		}
		close INPUT_DIFF_FILE or warn $! ? "Error closing the difference file: $!": "Exit status $? from the difference file";


		open (OPEN_SOFT_CIGAR_FILE, ">>", $temp_soft_cigar_file) or die $!;
		open (OPEN_DEL_FILE, ">>", $temp_del_file) or die $!;
		open (OPEN_INS_FILE, ">>", $temp_ins_file) or die $!;
		open (OPEN_DIFF_FILE, ">>", $temp_diff_file) or die $!;

	
		if(int($count_soft_cigar_line)<int($count_cigar_line)){
			my @zero_soft_cigar = (0) x (int($count_cigar_line)-int($count_soft_cigar_line));
			foreach my $each_zero_soft_cigar(@zero_soft_cigar){
				printf OPEN_SOFT_CIGAR_FILE "%f\n", (int($each_zero_soft_cigar));
			}
		}

		if(int($count_del_line)<int($count_cigar_line)){
			my @zero_cigar_del = (0) x (int($count_cigar_line)-int($count_del_line));
			foreach my $each_zero_cigar_del(@zero_cigar_del){
				printf OPEN_DEL_FILE "%f\n", (int($each_zero_cigar_del));
			}
		}	

		if(int($count_ins_line)<int($count_cigar_line)){
			my @zero_cigar_ins = (0) x (int($count_cigar_line)-int($count_ins_line));
			foreach my $each_zero_cigar_ins(@zero_cigar_ins){
				printf OPEN_INS_FILE "%f\n", (int($each_zero_cigar_ins));
			}	
		}

		if(int($count_diff_line)<int($count_cigar_line)){
			my @zero_cigar_diff = (0) x (int($count_cigar_line)-int($count_diff_line));
			foreach my $each_zero_cigar_diff(@zero_cigar_diff){
				printf OPEN_DIFF_FILE "%f\n", (int($each_zero_cigar_diff));
			}	
		}

		close OPEN_SOFT_CIGAR_FILE or warn $! ? "Error closing the soft cigar file: $!": "Exit status $? from the soft cigar file";
		close OPEN_DEL_FILE or warn $! ? "Error closing the deletion file: $!": "Exit status $? from the deletion file";
		close OPEN_INS_FILE or warn $! ? "Error closing the insertion file: $!": "Exit status $? from the insertion file";
		close OPEN_DIFF_FILE or warn $! ? "Error closing the difference file: $!": "Exit status $? from the difference file";

		#To skip soft clipped calculations if data are generated using PacBio
		if($is_pacbio ne "1"){

			my @cigar_values = split("\t", &sub_average($temp_soft_cigar_file, 5));
			my $cigar_10_percentile = &sub_percentile($temp_soft_cigar_file, 10);
			my $cigar_90_percentile = &sub_percentile($temp_soft_cigar_file, 90);

			$sub_info_output = $sub_info_output . $cigar_values[0] . "\t" . $cigar_values[1] . "\t" . (1-$cigar_values[2]) . "\t" . $cigar_10_percentile . $cigar_90_percentile;
		
		}		

		my $del_average = &sub_average($temp_del_file);
		my $del_10_percentile = &sub_percentile($temp_del_file, 10);
		my $del_90_percentile = &sub_percentile($temp_del_file, 90);

		$sub_info_output = $sub_info_output . $del_average . $del_10_percentile . $del_90_percentile;

		my $ins_average = &sub_average($temp_ins_file);
		my $ins_10_percentile = &sub_percentile($temp_ins_file, 10);
		my $ins_90_percentile = &sub_percentile($temp_ins_file, 90);

		$sub_info_output = $sub_info_output . $ins_average . $ins_10_percentile . $ins_90_percentile;
			
		my $diff_average = &sub_average($temp_diff_file);
		my $diff_10_percentile = &sub_percentile($temp_diff_file, 10);
		my $diff_90_percentile = &sub_percentile($temp_diff_file, 90);

		$sub_info_output = $sub_info_output . $diff_average . $diff_10_percentile . $diff_90_percentile;
		
	}

	return $sub_info_output;
	
}
###############################################################


#To get some characteristics for SV 
sub SV_info{

	#To get input parameters from each SV positions
	my($temp_SV_line_info, $temp_SV_sorted_bam) = @_;
	
	#To create temp file for each SV
	my $temp_sv_bed=$random_name . "_temp_SV.bed";

	###############################################################
	#To obtain homozygous and heterozygous SNP counts as well as GCcontent for each SV
	###############################################################
	my @temp_SV_line = split(":", $temp_SV_line_info);
	my $chr_SV = $temp_SV_line[0];
	my ($start_SV, $end_SV) = split("\-",$temp_SV_line[1]);
	my $SV_total_length = int($end_SV)-int($start_SV);

	open (OUT_SV_BED, ">", $temp_sv_bed) or warn "";
	print OUT_SV_BED $chr_SV . "\t" . $start_SV . "\t" . $end_SV . "\n";	
	close OUT_SV_BED or warn $! ? "Error closing the temp SV bed file: $!": "Exit status $? from the SV bed file";
	
	my $temp_SV_depth_file=$random_name . "_temp_SV_depth";
		
	open (OUT_SV_DEPTH_FILE, ">", $temp_SV_depth_file) or warn "";
	
	$command="samtools depth -r " . $temp_SV_line_info . " " . $temp_SV_sorted_bam . " >>" . $temp_SV_depth_file;
	system("$command");
	
	close OUT_SV_DEPTH_FILE or warn $! ? "Error closing the SV depth file: $!": "Exit status $? from the SV depth file";
	
	open(IN_SV_DEPTH_FILE, "<", $temp_SV_depth_file) or warn "";
	my $SV_average_depth_line = 0;
	my $SV_average_depth = 0;
	while(<IN_SV_DEPTH_FILE>){
		$SV_average_depth_line++;
	}
	close IN_SV_DEPTH_FILE or warn $! ? "Error closing the SV depth file: $!": "Exit status $? from the SV depth file";

	
	if(int($SV_average_depth_line)!=0){
		$SV_average_depth=`awk -F\'\\t\' \'{count_depth+=\$3} END {print count_depth/(int($SV_total_length)+1)}\' $temp_SV_depth_file `;
	}
	else{
		$SV_average_depth=0;	
	}
	
	#To print categories of sizes of SVs		
	if($SV_average_depth<=$cutoff)
	{
		$sub_sv_output = $sub_sv_output . "0" . "\t";		
	}
	elsif(($SV_average_depth)>($cutoff) && ($SV_average_depth)<=(2*($cutoff)))
	{
		$sub_sv_output = $sub_sv_output . "1" . "\t";		
	}
	else
	{
		$sub_sv_output = $sub_sv_output . "2" . "\t";
	}
	
	if(defined($INPUT_HOMOZYGOUS) && $in_eval==0){
		my $temp_homozygous_file=$random_name . "_temp_homozygous";
			
		open (OUT_HOMOZYGOUS_FILE, ">", $temp_homozygous_file) or warn "";

		$command="bedtools intersect -a " . $INPUT_HOMOZYGOUS . " -b " . $temp_sv_bed . " >>" . $temp_homozygous_file;
		system("$command");
		close OUT_HOMOZYGOUS_FILE or warn $! ? "Error closing homozygous file: $!": "Exit status $? from homozygous file";
		
		open(IN_HOMOZYGOUS_FILE, "<", $temp_homozygous_file) or warn "";
		my $count_homozygous_line = 0;
		
		while(<IN_HOMOZYGOUS_FILE>){
			$count_homozygous_line++;
		}
		close IN_HOMOZYGOUS_FILE or warn $! ? "Error closing homozygous file: $!": "Exit status $? from homozygous file";
				
		if(int($count_homozygous_line)!=0){	
			chomp($count_homozygous_line);
		}
		else{
			$count_homozygous_line=0;
		}

		if(int($count_homozygous_line)==0){
			$sub_sv_output = $sub_sv_output . "0" . "\t" . "0" . "\t";	
		}
		else{
			my $count_homozygous_line_sv=int($count_homozygous_line)/int($SV_total_length);
			$sub_sv_output = $sub_sv_output . $count_homozygous_line . "\t" .  $count_homozygous_line_sv . "\t";
		}
	}

	if(defined($INPUT_HETEROZYGOUS) && $in_eval==0){
	
		my $temp_heterozygous_file=$random_name . "_temp_heterozygous";
			
		open (OUT_HETEROZYGOUS_FILE, ">", $temp_heterozygous_file) or warn "";

		$command="bedtools intersect -a " . $INPUT_HETEROZYGOUS . " -b " . $temp_sv_bed . " >>" . $temp_heterozygous_file;
		system("$command");
		close OUT_HETEROZYGOUS_FILE or warn $! ? "Error closing heterozygous file: $!": "Exit status $? from heterozygous file";

		
		open(IN_HETEROZYGOUS_FILE, "<", $temp_heterozygous_file) or warn "";
		my $count_heterozygous_line = 0;
		
		while(<IN_HETEROZYGOUS_FILE>){
			$count_heterozygous_line++;
		}
		close IN_HETEROZYGOUS_FILE or warn $! ? "Error closing heterozygous file: $!": "Exit status $? from heterozygous file";
		
		if(int($count_heterozygous_line)!=0){	
			chomp($count_heterozygous_line);
		}
		else{
			$count_heterozygous_line=0;
		}

		if(int($count_heterozygous_line)==0){
			$sub_sv_output = $sub_sv_output . "0" . "\t" . "0" . "\t";	
		}
		else{
			my $count_heterozygous_line_sv=int($count_heterozygous_line)/int($SV_total_length);		
			$sub_sv_output = $sub_sv_output . $count_heterozygous_line . "\t" . $count_heterozygous_line_sv . "\t";
		}
	}
	

	if(defined($INPUT_FASTA) && $in_eval==0){

		my $temp_gccontent_file=$random_name . "_temp_gccontent";
			
		open (OUT_GCCONTENT_FILE, ">", $temp_gccontent_file) or warn "";

		$command="bedtools nuc -fi " . $INPUT_FASTA . " -bed " . $temp_sv_bed . " >>" . $temp_gccontent_file;
		system("$command");
		close OUT_GCCONTENT_FILE or warn $! ? "Error closing GCContent file: $!": "Exit status $? from GCContent file";


		open(IN_GCCONTENT_FILE, "<", $temp_gccontent_file) or warn "";
		my $GC_content_line = 0;
		
		while(<IN_GCCONTENT_FILE>){
			$GC_content_line++;
		}
		close IN_GCCONTENT_FILE or warn $! ? "Error closing GCContent file: $!": "Exit status $? from GCContent file";
		
		my $GC_content=0;
		
		if(int($GC_content_line)!=0){
			$GC_content=`awk -F\'\\t\' \' NR!=1 {print \$5}\' $temp_gccontent_file `;
			chomp($GC_content);
			if (length $GC_content == 0){
				$GC_content="Genome region is outside of reference genome";
			}
		}
		else{
			$GC_content=0;
		}

		$sub_sv_output = $sub_sv_output . $GC_content . "\t";
	}
	
	
		
	###############################################################
	#To obtain overlap with RepeatMasker for each SV
	###############################################################
	if(defined($INPUT_REPEAT_SINE) && $in_eval==0){

		#To create temp file for each SV
		my $temp_repeat_sine=$random_name . "_temp_repeat_sine.bed";
		my $temp_repeat_sine_sort=$random_name . "_temp_repeat_sine_sort.bed";
		my $temp_repeat_sine_merge=$random_name . "_temp_repeat_sine_merge.bed";
			
		open (OUT_REPEAT_SINE, ">", $temp_repeat_sine) or warn "";
		
		$command="bedtools intersect -a " . $temp_sv_bed . " -b " . $INPUT_REPEAT_SINE . " >>" . $temp_repeat_sine;
		system("$command");
		
		close OUT_REPEAT_SINE or warn $! ? "Error closing the RepeatMasker SINE file: $!": "Exit status $? from the RepeatMasker SINE file";
		
		open(IN_REPEAT_SINE, "<", $temp_repeat_sine) or warn "";
		my $count_SINE_line = 0;
		
		while(<IN_REPEAT_SINE>){
			$count_SINE_line++;
		}
		close IN_REPEAT_SINE or warn $! ? "Error closing RepeatMasker SINE file: $!": "Exit status $? from RepeatMasker SINE file";
			
			
		my $count_SINE=0;
		
			
		if(int($count_SINE_line)!=0 && $in_eval==0){
			open (OUT_REPEAT_SINE_SORT, ">", $temp_repeat_sine_sort) or warn "";
				
			$command="sort -k1,1 -k2,2g -o " . $temp_repeat_sine_sort . " " . $temp_repeat_sine;
			system("$command");
			
			close OUT_REPEAT_SINE_SORT or warn $! ? "Error closing RepeatMasker SINE sort file: $!": "Exit status $? from RepeatMasker SINE file";
			open (OUT_REPEAT_SINE_MERGE, ">", $temp_repeat_sine_merge) or warn "";
			
			$command="bedtools merge -i " . $temp_repeat_sine_sort . ">>" . $temp_repeat_sine_merge;
			system("$command");
			close OUT_REPEAT_SINE_MERGE or warn $! ? "Error closing RepeatMasker SINE merge file: $!": "Exit status $? from RepeatMasker Simple file";

			$count_SINE=`awk -F \'\\t\' \'BEGIN{SUM=0}{SUM+=int(\$3)-int(\$2)}END{print SUM}\' $temp_repeat_sine_merge `;
			chomp($count_SINE);
		}

		$sub_sv_output = $sub_sv_output . (int($count_SINE)/int($SV_total_length)) . "\t";
	}

	if(defined($INPUT_REPEAT_SIMPLE) && $in_eval==0){
		my $temp_repeat_simple=$random_name . "_temp_repeat_simple.bed";
		my $temp_repeat_simple_sort=$random_name . "_temp_repeat_simple_sort.bed";
		my $temp_repeat_simple_merge=$random_name . "_temp_repeat_simple_merge.bed";

		open (OUT_REPEAT_SIMPLE, ">", $temp_repeat_simple) or warn "";
				
		$command="bedtools intersect -a " . $temp_sv_bed . " -b " . $INPUT_REPEAT_SIMPLE . " >>" . $temp_repeat_simple;
		system("$command");

		close OUT_REPEAT_SIMPLE or warn $! ? "Error closing the RepeatMasker SIMPLE file: $!": "Exit status $? from the RepeatMasker SIMPLE file";

		open(IN_REPEAT_SIMPLE, "<", $temp_repeat_simple) or warn "";
		my $count_SIMPLE_line = 0;
			
		while(<IN_REPEAT_SIMPLE>){
			$count_SIMPLE_line++;
		}
		close IN_REPEAT_SIMPLE or warn $! ? "Error closing RepeatMasker SIMPLE file: $!": "Exit status $? from RepeatMasker SIMPLE file";

		my $count_SIMPLE=0;

		if(int($count_SIMPLE_line)!=0 && $in_eval==0){
			open (OUT_REPEAT_SIMPLE_SORT, ">", $temp_repeat_simple_sort) or warn "";
					
			$command="sort -k1,1 -k2,2g -o " . $temp_repeat_simple_sort . " " . $temp_repeat_simple;
			system("$command");

			close OUT_REPEAT_SIMPLE_SORT or warn $! ? "Error closing RepeatMasker SIMPLE sort file: $!": "Exit status $? from RepeatMasker SINE file";
			open (OUT_REPEAT_SIMPLE_MERGE, ">", $temp_repeat_simple_merge) or warn "";

			$command="bedtools merge -i " . $temp_repeat_simple_sort . ">>" . $temp_repeat_simple_merge;
			system("$command");
			close OUT_REPEAT_SIMPLE_MERGE or warn $! ? "Error closing RepeatMasker SIMPLE merge file: $!": "Exit status $? from RepeatMasker Simple file";

			$count_SIMPLE=`awk -F \'\\t\' \'BEGIN{SUM=0}{SUM+=int(\$3)-int(\$2)}END{print SUM}\' $temp_repeat_simple_merge `;
			chomp($count_SIMPLE);
		}

		$sub_sv_output = $sub_sv_output . (int($count_SIMPLE)/int($SV_total_length)) . "\t";
	}

	return $sub_sv_output;
	
}		
###############################################################	####

###################################################################
# To calculate various parameters for each structural variant based on bam file
###################################################################
open (REGION, "<", $temp_bed) or die $!;

#To store outputs
my $M_sub_info_out="";
my $M_sub_sv_out="";
my $R_sub_info_out="";

#To read each line of structural variant and obtain different parameters
while (my $line_region = <REGION>) {

	my $chr_info_out="";

	#To define boundaries for each structural variant
	my $left_middle_temp_line = "";
	my $right_middle_temp_line = ""; 
	chomp($line_region);	
	my @region_array = ();
	@region_array = split("\t", $line_region);


	###################################################################
	# To calculate various parameters for three column input file
	###################################################################
	if((0+@region_array) == 3){
		$chr_info_out = $chr_info_out . $region_array[0] . "\t" . ($region_array[1]+$flanking) . "\t" . ($region_array[2]-$flanking) . "\t" . (int($region_array[2]-$flanking)-int($region_array[1]+$flanking)) . "\t";

		#To print categories of sizes of SVs		
		if(int($region_array[2]-$flanking)-int($region_array[1]+$flanking) <= 100)
		{
			$chr_info_out = $chr_info_out . "0" . "\t";		
		}
		elsif(int($region_array[2]-$flanking)-int($region_array[1]+$flanking) > 100 && int($region_array[2]-$flanking)-int($region_array[1]+$flanking) <= 1000)
		{
			$chr_info_out = $chr_info_out . "1" . "\t";		
		}
		elsif(int($region_array[2]-$flanking)-int($region_array[1]+$flanking) > 1000 && int($region_array[2]-$flanking)-int($region_array[1]+$flanking) <= 10000)
		{
			$chr_info_out = $chr_info_out . "2" . "\t";		
		}
		else
		{
			$chr_info_out = $chr_info_out . "3" . "\t";
		}
		
		#For left flanking region of structural variant 		
		my $left_temp_line = $region_array[0] . ":" . int($region_array[1]) . "-" . int($region_array[1]+$flanking-1);
		
		#For left flanking region and length of structural variant
		if((int($region_array[1]+$flanking-1)+int($flanking)) >= int($region_array[2]-$flanking)){
			$left_middle_temp_line = $region_array[0] . ":" . int($region_array[1]+$flanking-1) . "-" . int($region_array[2]-$flanking);
		}
		else{
			$left_middle_temp_line = $region_array[0] . ":" . int($region_array[1]+$flanking-1) . "-" . int(int($region_array[1]+$flanking-1)+int($flanking));
		}
		
		#For structural variant itself		
		my $temp_line = $region_array[0] . ":" . int($region_array[1]+$flanking-1) . "-" . int($region_array[2]-$flanking);
		
		#For right flanking region and length of structural variant
		if((int($region_array[2]-$flanking)-int($flanking)) <= int($region_array[1]+$flanking)){
			$right_middle_temp_line = $region_array[0] . ":" . int($region_array[1]+$flanking-1) . "-" . int($region_array[2]-$flanking);
		}
		else{
			$right_middle_temp_line = $region_array[0] . ":" . int(int($region_array[2]-$flanking)-int($flanking)) . "-" . int(int($region_array[2]-$flanking));
		}
		
		#For right flanking region of structural variant
		my $right_temp_line = $region_array[0] . ":" . int($region_array[2]-$flanking) . "-" . int($region_array[2]);
		my $bam =  $temp_bam .".bam";

   		eval {
   			local $SIG{ALRM} = sub { print OUT_LARGE_SV $line_region . "\n"; $in_eval = 1; goto end_eval; };
					
			#To schedule alarm in $time_cutoff seconds
   			alarm $time_cutoff;

       		#To run sub_info subroutine
   			&sub_info($left_temp_line, $bam);
   			&sub_info($left_middle_temp_line, $bam);
   			$M_sub_info_out = &sub_info($temp_line, $bam);
   			$sub_info_output = "";
   			$M_sub_sv_out = &SV_info($temp_line, $bam);
			$sub_sv_output = "";
   			&sub_info($right_middle_temp_line, $bam);
   			$R_sub_info_out = &sub_info($right_temp_line, $bam);
   			$R_sub_info_out =~ s/\s+$//;
			$sub_info_output = "";       				
				
			#To cancel the alarm
   			alarm 0;                    
   						
		};
   		
              
   		if ( $@ && $@ !~ /alarm clock restart/ ) {
   			$in_eval = 1;
   			goto end_eval;
		}

		else{
			print OUT_INFO $chr_info_out . $M_sub_info_out . $M_sub_sv_out . $R_sub_info_out . "\n";
		}
		end_eval:
			next;
			
}
	
	###################################################################
	# To calculate various parameters for more than three column input file
	###################################################################
	else{
		for my $j (0 .. ((0+@region_array)-1)){

			if($j==0){
				$chr_info_out = $chr_info_out . $region_array[0] . "\t";
			}
			elsif($j==1){
				$chr_info_out = $chr_info_out . ($region_array[1]+$flanking) . "\t";
			}
			elsif($j==2){
				$chr_info_out = $chr_info_out . ($region_array[2]-$flanking) . "\t" . (int($region_array[2]-$flanking)-int($region_array[1]+$flanking)) . "\t";
				
				#To print categories of sizes of SVs		
				if(int($region_array[2]-$flanking)-int($region_array[1]+$flanking)<=100)
				{
					$chr_info_out = $chr_info_out . "0" . "\t";		
				}
				elsif(int($region_array[2]-$flanking)-int($region_array[1]+$flanking)>100 && int($region_array[2]-$flanking)-int($region_array[1]+$flanking)<=1000)
				{
					$chr_info_out = $chr_info_out . "1" . "\t";		
				}
				elsif(int($region_array[2]-$flanking)-int($region_array[1]+$flanking)>1000 && int($region_array[2]-$flanking)-int($region_array[1]+$flanking)<=10000)
				{
					$chr_info_out = $chr_info_out . "2" . "\t";		
				}
				else
				{
					$chr_info_out = $chr_info_out . "3" . "\t";
				}			
			}
			elsif($j==((0+@region_array)-1)){

				#To print all the additional columns of input file
				$chr_info_out = $chr_info_out . $region_array[$j] . "\t";
							
				#For left flanking region of structural variant 		
				my $left_temp_line = $region_array[0] . ":" . int($region_array[1]) . "-" . int($region_array[1]+$flanking-1);
		
				#For left flanking region and length of structural variant				
				if((int($region_array[1]+$flanking-1)+int($flanking))>=int($region_array[2]-$flanking)){
					$left_middle_temp_line = $region_array[0] . ":" . int($region_array[1]+$flanking-1) . "-" . int($region_array[2]-$flanking);
				}
				else{
					$left_middle_temp_line = $region_array[0] . ":" . int($region_array[1]+$flanking-1) . "-" . int(int($region_array[1]+$flanking-1)+int($flanking));
				}
				
				#For structural variant itself	
				my $temp_line = $region_array[0] . ":" . int($region_array[1]+$flanking-1) . "-" . int($region_array[2]-$flanking);
		
				#For right flanking region and length of structural variant
				if((int($region_array[2]-$flanking)-int($flanking))<=int($region_array[1]+$flanking)){
					$right_middle_temp_line = $region_array[0] . ":" . int($region_array[1]+$flanking-1) . "-" . int($region_array[2]-$flanking);
				}
				else{
					$right_middle_temp_line = $region_array[0] . ":" . int(int($region_array[2]-$flanking)-int($flanking)) . "-" . int(int($region_array[2]-$flanking));
				}
				
				#For right flanking region of structural variant
				my $right_temp_line = $region_array[0] . ":" . int($region_array[2]-$flanking) . "-" . int($region_array[2]);
				my $bam =  $temp_bam .".bam";				
	
				eval { 
									
					local $SIG{ALRM} = sub { print OUT_LARGE_SV $line_region . "\n"; $in_eval = 1; goto end_eval; };
					
					#To schedule alarm in $time_cutoff seconds
   					alarm $time_cutoff;
			 					
   					#To run sub_info subroutine
   					&sub_info($left_temp_line, $bam);
   					&sub_info($left_middle_temp_line, $bam);
   					$M_sub_info_out = &sub_info($temp_line, $bam);
   					$sub_info_output = "";
   					$M_sub_sv_out = &SV_info($temp_line, $bam);
					$sub_sv_output = "";
   					&sub_info($right_middle_temp_line, $bam);
   					$R_sub_info_out = &sub_info($right_temp_line, $bam);
   					$R_sub_info_out =~ s/\s+$//;
					$sub_info_output = "";       				
				
					#To cancel the alarm
   					alarm 0;                    
   					
				};
   				
   				if ($@ && $@ !~ /alarm clock restart/ ) {
   					$in_eval = 1;
   					goto end_eval;
				}

				else{
					print OUT_INFO $chr_info_out . $M_sub_info_out . $M_sub_sv_out . $R_sub_info_out . "\n";
				}
				
				end_eval:
					next;							
		}
		else{
			#To print all the additional columns of input file
			$chr_info_out = $chr_info_out . $region_array[$j] . "\t";			
		}
		}
	}
}
close REGION or warn $! ? "Error closing the temp genomic region file: $!": "Exit status $? from the genomic region file";
close OUT_INFO or warn $! ? "Error closing the output file: $!": "Exit status $? from the output file";

my $large_SV_filesize = -s OUT_LARGE_SV;
close OUT_LARGE_SV or warn $! ? "Error closing the output file with large SV genomic regions: $!": "Exit status $? from the large SV file";

#To remove the file if the size is zero
if ($large_SV_filesize==0){
	unlink glob $temp_large_SV;
} 

#Deleting temp files
unlink glob $random_name . "*";


#!/usr/bin/perl

# use additional perl packages
use strict						;
use warnings					;
use File::Path					; # for deleting non-empty directories
use File::Copy					; # to use copy and move
use File::Find::Rule			;
use File::Path qw(rmtree)		;
use Data::Dumper				;
use Statistics::R				;
use Cwd							; #use getcwd
use POSIX qw(strftime)			;
use bignum						; # avoid floating point issues in while and for loops



#------------------- DEFAULT PARAMETER SETUP ---------------------------------------------------------

######################################################################################################
# Parameter Setting								                                                     #
######################################################################################################
my	%value_of_scr_par	=	(

	'seq_algorithm'		=>	{
	
							'nap'						=>	1		,	# SeaLion parameter: 1-> analysis of split support NAP  (using standard correction factor)
	}	,
	
	'seq_filterRisk'	=>	{
		
							'risk1'						=>	1		,	# SeaLion parameter: 1-> quartet filtering if quartet is in 'risk1' defined Nc/Na (uncorrected) ratio-range 
							'risk2'						=>	0		,	# SeaLion parameter: 1-> quartet filtering if quartet is in 'risk2' defined Nc/Na (uncorrected) ratio-range 
							'tlrisk'					=>	0.70	,	# SeaLion parameter: lower threshold of Nc/Na risk range
							'turisk'					=>	1		,	# SeaLion parameter: upper threshold of Nc/Na risk range
							'rscale'					=>	0.01	,	# SeaLion parameter: scale steps for optimum Nc(Na filtering
							'ralpha'					=>	3		,	# SeaLion parameter: alpha value of f(P) function for risk filter optimization
	}	,

	'seq_filterDist'	=>	{
		
							'qdist1'					=>	1		,	# SeaLion parameter: 1-> quartet filtering if split support distance between best and 2nd best rooted quartet tree < threshold qdist1st2nd
							'tldist'					=>	0.0		,	# SeaLion parameter: lower threshold of support distance range
							'tudist'					=>	0.1		,	# SeaLion parameter: upper threshold of support range
							'dscale'					=>	0.01	,	# SeaLion parameter: scale steps for optimum distance filtering 
							'dalpha'					=>	3		,	# SeaLion parameter: alpha value of f(P) function for dist filter optimization
	}	,
	
	'pdf_prtResult'		=>	{
							
							'main'						=>	1		,	# SeaLion parameter: main result files (e.g. ternär plots about tree success, lineplots about quaret filtering success...)
	}	,
	
	'sys_parameter'		=>	{
							
							'sysstart'					=>	0							,	# SeaLion parameter: 1 -> command for process start, 0 -> open config menu
							'msa_nopath'				=>	''							,	# SeaLion parameter: name of the clanfile corresponding alignment infile without path dir (only defined with script internal msa definition)
							'msafile'					=>	''							,	# SeaLion parameter: name of the clanfile corresponding alignment infile
							'clfile'					=>	''							,	# SeaLion parameter: name of the alignment infile corresponding clan-infile (.txt) with possible path prefix
							'clf_nopath'				=>	''							,	# SeaLion parameter: name of the alignment infile corresponding clan-infile (.txt) without possible path prefix (only defined with script internal msa definition)
							'outgroup'					=>	''							,	# SeaLion parameter: specified code of defined outgroup-clan in clan-infile
							'l_prt'						=>	0							,	# SeaLion parameter: if defined '1', generate latex tables in pdf and svg. Blocked if '0'
							'r_prt'						=>	0							,	# SeaLion parameter: if defined '1', generate R-result plots in pdf and svg. Blocked if '0'
							'screen_prt'				=>	0							,	# SeaLion parameter: if defined '1', print detailed calculation steps on the terminal. Blocked if '0'
							't_qtrPrt'					=>	50							,	# SeaLion parameter: threshold of 4clan tree corresponding newick-string characters until the newick string print in grafic plots is replaced by 'qtree'-coding to reduce print-space in grafics
							'write_spd'					=>	0							,	# SeaLion parameter: if defined '1', printing of split-pattern distribution files under -restart option
							'i_main'					=>	''							,	# SeaLion parameter: constant name part of penguin-analysis main resultfolder (unnecessary after complete script programming)
							'o_main'					=>	'SeaLion_results'			,	# SeaLion parameter: name of the main output-resultfolder
							'pen_bin'					=>	'sealion1.pl'				,	# SeaLion parameter: actual penguin script name
							'ice_bin'					=>	'icebreaker.o'				,	# SeaLion parameter: name of the terminal-executable icebreaker tree-fishing script
							'mean'						=>	1							,	# SeaLion parameter: if defined '1', calculation and print of mean support results. Blocked if '0'
							'median'					=>	1							,	# SeaLion parameter: if defined '1', calculation and print of median support results. Blocked if '0'
							'request'					=>	1							,	# SeaLion parameter: if defined '1', no user request if maximum limit of quartet analyses is lower than total number of available quartets of a 4clan combination
							'restart'					=>	1							,	# SeaLion parameter: if defined '1', using already existing P4 files from defined filepath 'i_main'
							'p4data'					=>	0							,	# SeaLion parameter: if defined '1', keep p4 results in an extra folder, otherwise delete p4 output after data processing
							'rycoding'					=>	0							,	# SeaLion parameter: if defined '1', recode sequence characters to RY, do not otherwise
							'spdfname'					=>	"*_split_calc_q*.txt"		,	# SeaLion parameter: Filename_prt for split-pattern distribution outfiles. Pattern is needed for path finding and input reading
	}	,
	
	'p4_parameter'		=>	{
							
							'alpha'						=>	1.0			,	# SeaLion parameter: P4 start value alpha
							'pinv'						=>	0.3			,	# SeaLion parameter: P4 start value proportion invariable sites (I)
							'model_aa'					=>	'lg'		,	# SeaLion parameter: ML model for aa expected quartet split-pattern distribution, available: d78, jtt, wag, mtrev24, lg, blosum62
							'model_nu'					=>	'GTR'		,	# SeaLion parameter: ML model for nu expected quartet split-pattern distribution, available: JC, F81, K2P, HKY, GTR
	}	,

	'qlimit'			=>	{
							
							'qlength'					=>	10000		,	# SeaLion parameter: Minimum allowed sequence length of each quaretet
							'qnumber'					=>	20000		,	# SeaLion parameter: Maximum number of single quartet analyses per 4clan combination
	}	,
	
	'author_info'		=>	{
		
							'name'						=>	'Patrick Kueck'				,	
							'institute'					=>	'LIB Bonn, GER'				,	
							'email'						=>	'pkueck[at]leibniz-lib.de'	,	
							'homepage'					=>	'www.bonn.leibniz-lib.de'	,
							'version'					=>	'1.0'						,
							'language'					=>	'Perl'						,
							'lastUpdate'				=>	'March, 2024'				,
							'p4'						=>	'http://p4.nhm.ac.uk/'		
	}	,
);
######################################################################################################
######################################################################################################



#------------------- R PLOT PRINT DEFINITIONS --------------------------------------------------------

######################################################################################################
# Print control of single R plots: if hash value == 1 -> print, otherwise -> print de-activated      #
######################################################################################################
my	%print_Rplot_quartets	=	(
	
	# MAIN
	'ggtree_rooted_clade_trees'					=>	1 , # R1*	treeplot:		MAIN	UserVersion	
	'pline_taxsupport_clan4'					=>	1 , # T2*	Lineplot:		MAIN	UserVersion	R lineplot of single, mean, and median taxon related unfiltered quartet support
	'ptern_taxsupport_clan4'					=>	1 , # T1*	Ternaryplot:	MAIN	UserVersion	R ternary plot of single, mean, and median taxon related unfiltered quartet support
	'pbarplot_treesuccess_mclan'				=>	1 , # MQ3*	Barplot:		MAIN	UserVersion	adapted	R ternary plot of single, mean, and median unfiltered quartet support
	'ptern_seqsupport_clan4'					=>	1 , # MQ2*	Ternaryplot:	MAIN	UserVersion	adapted	R ternary plot of single, mean, and median unfiltered quartet support
	'ptern_seqfilter_clan4_fil_unfil'			=>	1 , # MQ4*	Ternaryplot:	MAIN	UserVersion	adapted	R ternary plot of single, mean, and median support in respect of fillter approaches, highlighting potentially rNa/Nc risky quartets
	'ptern_seqfilter_clan4'						=>	1 , # MQ2*	Ternaryplot:	MAIN	UserVersion adapted R ternary plot of single, mean, and median support in respect of fillter approaches, highlighting potentially rNa/Nc risky quartets
	'pline_clan4_MeanMedian'					=>	1 , # MQ1*	Lineplot: 		MAIN	UserVersion	adapted	R lineplot of Penguin overall mean and median support success of the three rooted quartet-trees using original CF site info (NAP)
	'ppoint_rNcNaOrgRisk_FN'					=>	1 , # Q3*	Pointplot:		MAIN	DevVersion	adapted	ratio Nc/Na (unreduced), related to 4clan corresponding file-number for each rooted tree
	'ppoint_QseqLength_best_tree'				=>	1 , # Q7*	Pointplot:		MAIN	DevVersion	R pointplot, showing best quartet-tree support for each 4clan combination, sorted by the quartet related number of split-pattern (seq length)
	'pline_taxreject_clan4'						=>	1 , # Q5*	Lineplot:		MAIN	UserVersion	adapted	R lineplot of N rejected quartets in 4clan combinations using 'risk' quartet-filtering
	'pline_taxreject_total'						=>	1 , # Q5*	Lineplot:		MAIN	UserVersion	adapted	R lineplot of N rejected quartets overall 4clan combinations using 'risk' quartet-filtering
	'pline_taxremain_clan4'						=>	1 , # Q5*	Lineplot:		MAIN	UserVersion	adapted	R lineplot of N remaining quartets in 4clan combinations using 'risk' quartet-filtering
	'pline_taxremain_total'						=>	1 , # Q5*	Lineplot:		MAIN	UserVersion	adapted	R lineplot of N remaining quartets overall 4clan combinations using 'risk' quartet-filtering
	'pline_qfilter_optimzation'					=>	1 , # MQ6*	Lineplot:		MAIN	UserVersion	adapted	R lineplot of N remaining quartets overall 4clan combinations using 'risk' quartet-filtering
	'pline_qnumber_optimzation'					=>	1 , # MQ6*	Lineplot:		MAIN	UserVersion	adapted	R lineplot of N remaining quartets overall 4clan combinations using 'risk' quartet-filtering
	'pbarplot_filtered_Nquartets'				=>	1 , # Q6*	Barplot			MAIN	UserVersion	adapted	R barplot of rejected and total number of 4clan rejected quartets
	'ppoint_rNcNaOrg_CFobsCFexp_unfilt'			=>	1 , # Q3*	Pointplot:		MAIN	UserVersion	adapted	Proportion of Nc/Na split signal in relation to CF exp/obs differences for each of the three rooted quartet trees given different lengths BL2
	'ppoint_rNcNa_CFobsCFexp_filt'				=>	1 , # Q3*	Pointplot:		MAIN	UserVersion	adapted	Proportion of Nc/Na split signal in relation to CF exp/obs differences for each of the three rooted quartet trees given different lengths BL2
);
######################################################################################################
######################################################################################################



#------------------- LATEX PRINT DEFINITIONS ---------------------------------------------------------

######################################################################################################
# Print control of single Latex tables: if hash value == 1 -> print, otherwise -> print deactivated  #
######################################################################################################
my	%print_latex	=	(
		
	'table_4clancodes_sealion'					=>	1 , # LP3*	Table:	MAIN	List of analysed 4clan combinations and internal coding for 4clan result prints
	'table_clancodes_sealion'					=>	1 , # LP2*	Table:	MAIN	List of identified clans and internal coding for 4clan result prints
	'table_parameter_sealion'					=>	1 , # LP1*	Table:	MAIN	List of specified SeaLion analysis parameter
	'table_treesupport_4clan'					=>	1 , # LQ1*	Table:	MAIN	List of filtered and unfiltered 4clan-tree support of each of the three possible 4clan trees
	'table_besttreesupport_4clan'				=>	1 , # LQ2*	Table:	MAIN	List of filtered and unfiltered best 4clan-tree support
	'table_4seq_associated_taxa_filtered'		=>	1 , # LQ3*	Table:	MAIN	List of taxon participations in filtered quartets (in total and per 4clan combination)
	'table_4seq_associated_taxa_unfiltered'		=>	1 ,	# LQ4*	Table:	MAIN	List of taxon participations in unfiltered quartets (in total and per 4clan combination)
	'table_4seq_unassociated_taxa_unfiltered'	=>	1 ,	# LQ5*	Table:	MAIN	List of taxa who have not been analysed in any unfiltered quartet
	'table_optimized_4seq_filter_threshold'		=>	1 , # LQ6*	Table:	MAIN	List of 4seq-filter optimzed threshold values
	'table_number_analysed_quartets'			=>	1 , # LQ7*	Table:	MAIN	List of unfiltered and filtered analysed quartets
	'table_besttreesupport_mclan'				=>	1 , # LRC1*	Table:	MAIN	List of multiclan best-tree suppport
	'table_singletreesupport_mclan'				=>	1 , # LRC2*	Table:	MAIN	List of single multiclan tree suppport
	'table_4clan_quartet_trees'					=>	1 ,	# LQT1*	Table:	MAIN	List of taxon related triangle Q1,Q2,Q3 coded 4-clan quartet relationships
);

#########################################
# key1	: Latex table unique printcode (e.g. 'latexC')
# key2	: *.tex filename; 
# value	: @list of assigned table lines
my	%hoh_linedata_of_latexfilename_of_printcode ;	
#########################################

#########################################
# key1	: textfile unique printcode (e.g. 'textA')
# value	: filename for print
my	%textfilename_of_printcode = (	
	
	'txtA' 		=> 'LP1_SeaLion_parameter_info.txt'   ,
	'txtLRC1'	=> 'LRC1_SeaLion_rooted-clades_best-tree_support.txt',
	'txtLRC2'	=> 'LRC2_SeaLion_rooted-clades_single-tree_support.txt',
	'txtLQ1'	=> 'LQ1_SeaLion_clade-quartet_single-tree_support.txt',
	'txtLQ2'	=> 'LQ2_SeaLion_clade-quartet_best-tree_support.txt',
	'txtLQ3'	=> 'LQ3_SeaLion_clade-quartet_species_participations_filtered.txt',
	'txtLQ4'	=> 'LQ4_SeaLion_clade-quartet_species_participations_unfiltered.txt',
	'txtLQ5'	=> 'LQ5_SeaLion_species_unassignable.txt',
	'txtLQ6'	=> 'LQ6_SeaLion_clade-quartet_optimized_RISK_filter_threshold.txt',
	'txtLQ7'	=> 'LQ7_SeaLion_clade-quartet_analysed_number_of_species-quartets.txt',
	'txtLQT1'	=> 'LQT1_SeaLion_clade-quartet_tree_codes.txt',
) ;
######################################################################################################
######################################################################################################



#------------------- ARGV handling -------------------------------------------------------------------

######################################################################################################
# READ COMMAND LINE                                                                                  #
######################################################################################################

#########################################
# All option values are stored in %value_of_scr_par
&argv_handling(

	\%value_of_scr_par								, # key1: sealion user available parameter options; value: optional parameter setting									IN: defined; OUT: changed
) ;
######################################################################################################

START:
#print "\n\texit\n\n;";exit;
#------------------- Build Result Folder -------------------------------------------------------------

######################################################################################################
# BUILD MAIN-FOLDER and SUBFOLDER FOR OUTPUT STORING                                                 #
######################################################################################################

#########################################
# key1	: subfolder code ('pdf', 'tsv'... )
# value	: key corresponding resultfolder path
my	%subfolder_q					= () ;
#########################################

#########################################
# key1	: subsubfolder code ('pdf_single', 'pdf_sumed' )
# value	: key corresponding resultfolder path
my	%subsubfolder_q					= () ;
#########################################

#########################################
# build resultfolder
&build_resultfolder(
	
	\%value_of_scr_par	,	# key1: method_code; key2: method_option; value: float or integer		IN: defined; OUT: unchanged
	\%subfolder_q		,	# key1: subfolder code ('pdf', 'tsv'... ); value: resultfolder path		IN: undefined; OUT: defined
	\%subsubfolder_q	,	# key1: subfolder code ('pdf_single', 'pdf_sumed' ); value: path		IN: undefined; OUT: defined
) ;
######################################################################################################



#------------------- Assign Calculation & Filter Method Parameter ------------------------------------

######################################################################################################
# Assign Specified Calculation & Filter Parameter                                                    #
######################################################################################################

#########################################
# key1	: type ("ftaxon, ctaxon, tsuppo, original, expanded")
# value	: list of type corresponding parameter specifications 
my	%algorithm_def_methods			= () ;
################

################
# assign parameter
&assign_selected_methods(
	
	\%value_of_scr_par		,	# key1: method_code; key2: method_option; value: float or integer		IN: defined; OUT: unchanged
	\%algorithm_def_methods	,	# key1: methode type; key2: method code; value: integer					IN: undefined; OUT: defined
);
#########################################

#########################################
# key1	: type of supportfilter ('risk', 'dist')
# value : list of ftype corresponding threshold scales
my	%hol_scales_of_filtertype ;
################
	
################	
&assign_filter_scales(
	
	\%value_of_scr_par 			, # key1: method_code; key2: method_option; value: float or integer										IN: defined; OUT: unchanged
	\%hol_scales_of_filtertype	, # key1: type of supportfilter ('risk', 'dist'); value: list of ftype corresponding threshold scales	IN: undefined; OUT: defined
) ;
#########################################

######################################################################################################



#------------------- Screen/File print of parameter --------------------------------------------------

######################################################################################################
# PRINT PARAMETER-SETUP to SCREEN amd FILE                                                           #
######################################################################################################

#########################################
# Print specified parameters on terminal screen
&print_script_settings(

	\%value_of_scr_par								, # key1: sealion user available parameter options; value: optional parameter setting									IN: defined; OUT: unchanged
	\%subfolder_q									, # key1: subfolder code ('pdf_single', 'pdf_sumed' ); value: path														IN: defined; OUT: unchanged
	\%hoh_linedata_of_latexfilename_of_printcode	, # key1: Latex table unique printcode (e.g. 'latexC'); key2: *.tex filename; value: @list of assigned table lines		IN: defined; OUT: changed	
	\%textfilename_of_printcode						, # key1: textfile unique printcode (e.g. 'textA'); value: filename for print											IN: defined; OUT: unchanged
) ;
######################################################################################################



#------------------- CLANFILE handling ---------------------------------------------------------------

######################################################################################################
# READ AND CHECK INFILE DATA: CLANFILE                                                               #
######################################################################################################

#########################################
# key1	:	sequence-name; 
# value	:	clan definition
my	%clan_of_seq					= () ;
#########################################

#########################################
# key1	:	clanname; 
# value	:	list of clan assigned taxa
my	%hol_taxa_of_clan				= () ;	
#########################################

#########################################
# outgroup clan
my	$outgroup_clan	= $value_of_scr_par{sys_parameter}{outgroup} ;	
#########################################

#########################################
# clanfile without path
my $clanfile	= $value_of_scr_par{sys_parameter}{clf_nopath} ;
#########################################

#########################################
# Clan Infile Check (incl. Data Read IN)
&check_clan_infile(

	\%value_of_scr_par						,	# key1: sealion user available parameter options; value: optional parameter setting			IN: defined; OUT: unchanged
	\$outgroup_clan							,	# user specified outgroup clan																IN: defined; OUT: unchanged
	\%hol_taxa_of_clan						,	# key1: defined clan name; value: list of clan assigned taxa								IN: undefined; OUT: defined
	\%clan_of_seq							,	# key1: sequence-name; value: clan definition												IN: undefined, OUT: defined
) ;
######################################################################################################
#exit;


#------------------- ALIGNMENT + P4 analysis ---------------------------------------------------------

######################################################################################################
# CHECK MSA, BUILD 4CLAN+4TAXON COMBINATIONS, AND PROCESS P4 ANALYSIS                                #
######################################################################################################
unless ( $value_of_scr_par{sys_parameter}{restart} == 1 ){
	
	
	#------------------- ALIGNMENT handling --------------------------------------------------------------
	
	######################################################################################################
	# READ AND CHECK INFILE DATA: MSA                                                                    #
	######################################################################################################
	
	#########################################
	# key1	: type of sequence condition 
	# value : corresponding value
	# $data_of_infile_property{type}		= sequence tpye ('nuc' || 'aa')
	# $data_of_infile_property{length}		= original sequence length (e.g.: 10000)
	# $data_of_infile_property{remain_pos}	= Number of remaining positions after exclusion of defined sites
	my %data_of_infile_property 		= () ;			
	#########################################
	
	#########################################
	# key1	: taxonname;
	# value	: sequence
	my %unreduced_sequence_of_taxa		= () ;
	#########################################
	
	#########################################
	# key1	: taxonname;
	# key2	: seqposition number
	# value	: site character
	my %hol_sequence_states_of_taxon	= () ;
	#########################################
	
	#########################################
	# MSA Infile Check (incl. Data Read IN)
	&check_msa_infile (
	
		\%value_of_scr_par						,	# key1: sealion user available parameter options; value: optional parameter setting			IN: defined; OUT: unchanged
		\%data_of_infile_property				,	# key1: property; value : property value													IN: undefined; OUT: defined
		\%unreduced_sequence_of_taxa			,	# key1: taxonname; value sequence															IN: undefined; OUT: defined
		\%hol_sequence_states_of_taxon			,	# key1: taxon name; key2: state position number; value: sequence state at that position		IN: undefined; OUT: defined
		\%clan_of_seq							,	# key1: sequence-name; value: clan definition												IN: defined; OUT: changed
		\%hol_taxa_of_clan						,	# key1: defined clan name; value: list of clan assigned taxa								IN: defined; OUT: changed
	) ;
	######################################################################################################
	
	
	#------------------- Build 4clan combinations --------------------------------------------------------
	
	######################################################################################################
	# Build Clan-Quartets                                                                                #
	######################################################################################################
	
	#########################################
	# list of 4clan combinations
	my @clan_quartets	= () ;
	#########################################
	
	&generate_clan_quartets(
	
		\%hol_taxa_of_clan						,	# key1: clan name; value: list of clan assigned taxa										IN: defined; OUT: unchanged
		\@clan_quartets							,	# list of generated 4-clan combinations														IN: undefined; OUT: defined
		\%value_of_scr_par						,	# key1: sealion user available parameter options; value: optional parameter setting			IN: defined; OUT: unchanged
	) ;
	######################################################################################################
	
	
	#------------------- Build 4Clan Quartets & P4 Analysis ----------------------------------------------
	
	######################################################################################################
	# FOR EACH 4CLAN BUILD & P4-ANALYSIS of SINGLE QUARTET COMBINTIONS                                   #
	######################################################################################################
	&quartet_building_and_p4_main(
		
		\%value_of_scr_par				, # key1: method_code; key2: method_option; value: float or integer		IN: defined; OUT: unchanged
		\%subfolder_q					, # key: output-type (e.g. 'pdf'); value: pathdir						IN: defined; OUT: unchanged
		\%hol_taxa_of_clan				, # key1: clan; value: list of assigned taxa							IN: defined; OUT: unchanged
		\@clan_quartets					, # list of generated 4-clan combinations								IN: defined; OUT: unchanged
		\%unreduced_sequence_of_taxa	, # key1: taxonname; value sequence										IN: defined; OUT: unchanged
		\%hol_sequence_states_of_taxon	, # key1: taxon name; key2: state position number; value: state 		IN: defined; OUT: undefined
		\%data_of_infile_property		, # key1: property; value : property value								IN: defined; OUT: unchanged
		\%clan_of_seq					, # key1: sequence-name; value: clan definition							IN: defined; OUT: unchanged
	) ;
	
	%unreduced_sequence_of_taxa		= () ;
	%hol_sequence_states_of_taxon	= () ;
	%data_of_infile_property		= () ;
	
	############################################
	# set spd file-print command to zero to avoid de-novo printing during spd reading and analysis process
	$value_of_scr_par{sys_parameter}{write_spd}	= 0 ;
	######################################################################################################
}
######################################################################################################

#exit;

#------------------- Sample SPD-DATA from defined DIR ------------------------------------------------

######################################################################################################
# CHECK defined DIR for SPD input files                                                              #
######################################################################################################

###################################################
# constant string code of penguin main result folder
my	$main_resultfolder_code		=	$value_of_scr_par{sys_parameter}{i_main} ;
###################################################

###################################################
# store each penguin calculation-file '*_detailed_split_calc_*' file with path as list in @found
my	$search_file = $value_of_scr_par{sys_parameter}{spdfname} ;

print	"\n\t---------------------------------------------------\n",
		"\tReading SPD-files ", $search_file, " of DIR ", $main_resultfolder_code, "/*\n"	;

my	@files = &sampling_subdirectories( 
	
	\$main_resultfolder_code 	,	# main directory for filepath search 		IN:defined; OUT: unchanged
	\$search_file				,	# file to search for in search directory	IN:defined; OUT: unchanged
); 
#for ( sort @files ){ print "\nfound ", $_ } exit;
######################################################################################################



#------------------- Sample 4clan combinations from defined DIR --------------------------------------

######################################################################################################
## 4clan extraction & assignment due to quartet calculation filepath(s)  '*_detailed_split_calc_*'   #
######################################################################################################
#print	"\t4clan extraction, re-coding and assignment from SPD filepaths\n";

#########################################
# key1	:	4clan combination code; 
# value	:	list of quartet caluclation filepaths 
my	%hol_fp_qcal_of_4clan_code	;	
#########################################

#########################################
# key1	:	clan name; 
# value	:	assigned code
my	%code_of_clan	;	
#########################################

#########################################
# key1	:	4clan assigned clan names; 
# value	:	assigned code
my	%code_of_4clan	;	
#########################################

#########################################
# extraction of clan names from filepath and recoding
print	"\tAssigning SPD related 4Clan combinations from subfolders of DIR ", $main_resultfolder_code, "/*\n"	;

&clan4_extraction_from_filepaths(
	
	\@files								, # list of quartet individual split-claculation filepaths 										IN: defined; OUT: unchanged
	\$main_resultfolder_code			, # constant string of the main penguin resultfolder name										IN: defined; OUT: unchanged
	\%hol_fp_qcal_of_4clan_code			, # list of 4clan corresponding qcalc filepaths													IN: undefined; OUT: defined
	\%value_of_scr_par					, # key1: sealion user available parameter options; value: optional parameter setting			IN: defined; OUT: unchanged
	\%code_of_clan						, # key1: clanname; value assigned code															IN: undefined; OUT: defined
	\%code_of_4clan						, # key1: 4clan assigned clan names; value assigned code										IN: undefined; OUT: defined
	\%subfolder_q						, # key1: output-type (e.g. 'pdf'); value: pathdir												IN: defined; OUT: unchanged
);
#########################################

#########################################
# latex LP2 print of identified clan name and recode
&latex_data_sample_LP2( 
					
	\%hoh_linedata_of_latexfilename_of_printcode	, # key1: Latex table unique printcode (e.g. 'latexC'); key2: *.tex filename; value: @list of assigned table lines		IN: defined; OUT: changed	
	\%code_of_clan									, # key1: clanname; value assigned code																					IN: defined; OUT: unchanged
	\%textfilename_of_printcode						, # key1: textfile unique printcode (e.g. 'textA'); value: filename for print											IN: defined; OUT: unchanged
	\%subfolder_q									, # key1: output-type (e.g. 'pdf'); value: pathdir																		IN: defined; OUT: unchanged
) ;
#########################################

#########################################
# latex LP3 print of identified 4clan combinations and recode
&latex_data_sample_LP3( 
					
	\%hoh_linedata_of_latexfilename_of_printcode	, # key1: Latex table unique printcode (e.g. 'latexC'); key2: *.tex filename; value: @list of assigned table lines		IN: defined; OUT: changed	
	\%code_of_4clan									, # key1: 4clan assigned clan names; value assigned code																IN: defined; OUT: unchanged
	\%textfilename_of_printcode						, # key1: textfile unique printcode (e.g. 'textA'); value: filename for print											IN: defined; OUT: unchanged
	\%subfolder_q									, # key1: output-type (e.g. 'pdf'); value: pathdir																		IN: defined; OUT: unchanged
) ;
######################################################################################################



#------------------- 4CLAN SPECIFIC ANALYSIS FROM SINGLE-QUARTET SPLIT-PATTERN DATA ------------------

######################################################################################################
##  Analyse each 4-clan combination due to quartet-calculation files '*_detailed_split_calc_*'       #
######################################################################################################

#########################################
# defined @pdf_files where all quartet split-pattern generated pdf files are stored for final pdf summarization.
my	@pdf_files				;

# key: result type; value: list of associated pdf's in order of construction date
my	%hol_pdf_files_of_type ;
#########################################

#########################################
# storage of printed metafiles, used for R plotting
# key1	: R grafic unique printcode (e.g. 'metafileT') 
# key2	: *.tsv filename 
# value	: *tsv filename data addings (e.g. clanfilename, type of support...)
my	%hoh_tsvdata_of_tsvfilename_of_printcode ;
#########################################

#########################################
# storage of final pdf's for pdf summarization
# key: name of the main penguin resultfolder ($presult_folder)
# value: list of associated pdf's for resultfolder corresponding result summarization
my	%hol_pdf_final_sum_of_presult_folder ;
#########################################https://www.dropbox.com/scl/fi/9t9nnay57w5k38ad2syew/aligroove_1_08_GUI_version_windows.zip?rlkey=jeugg5e2t5dnqvulc4c9xillv&dl=0

#########################################
# multiple hash list, listing for each of the three rooted quartet-trees and each 4clan combination
# single quartet inferred support (NAP or NAPV) in analysed order
# key1	: 4clan combination code
# key2	: type of support ('nap' or 'napv')
# key3	: rooted-quartet tree 
# value	: list of single quartet-tree support (normed2)
my	%hol_sup_of_clan4_of_calc_of_qtree ;
#########################################

#########################################
# multiple hash, storing the final support
# (nap or napv) related to each rooted quartet-tree, seed,
# and 4clan combination, used for ternary plot underlying tsv table file
# key1	: 4clan combination code
# key2	: type of support ('nap' or 'napv')
# key3	: rooted-quartet tree 
# key4	: quartet-number
# value	: quartet-tree support value (normed2)	
my	%hoh_sup_of_clan4_of_calc_of_qtree_of_fc ;
#########################################

#########################################
# key1	: 4clan combination code 
# key2	: type of support ('nap' or 'napv')
# key3	: quartet-number 
# key4	: type of quartet-tree support distance ('first2nd' or 'first3rd')
my	%hoh_sup_dist_of_clan4_type_of_treepair ;
#########################################

#########################################
# key1	: 4clan combination code 
# key2	: type of supportfilter ('nap_risk1', 'napr_risk1', or 'napv_risk1' ...)
# key3	: quartet-number
# value	: quartet-filter corresponding score value (1 -> reject quartet)
my	%hoh_seen_rejected_of_4clan_of_filtertype_of_fn ;
#########################################
	
#########################################
# file counter overall 4clan combinations
my	$fc = 0 ; 
#########################################

#########################################
# key1	: type of supportfilter ('nap_risk1', 'napr_risk1', or 'napv_risk1' ...)
# key2	: mean|median
# value : list of 'tree\tmean|medianScore'
my	%hol_sup_of_type_of_mfilter ;
#########################################

#########################################
# key1	: quartet-number ($fc)
# value : quartet assigned taxa (alphabetically ordered, joined by ':'
my	%taxa4_of_fc ;
#########################################

#########################################
# key1	: 4clan combination code 
# value : N analysed quartets
my	%seen_clan4 ;
#########################################

#########################################
# key1	: 4clan combination code 
# key2	: taxon
# value : N analysed quartets
my	%hoh_seen_quartets_of_taxon_of_clan4 ;
#########################################

#########################################
# key1	: 4seq filtermethod (e.g. nap_risk1)
# key2	: 4clan combination code 
# value : N rejected quartets*3 (counted three times, each time for each quartet tree)
my	%hoh_N_qrej_of_clan4_of_ft ;
#########################################

#########################################
# key1	: 4seq filtermethod (e.g. nap_risk1)
# key2	: taxon
# key3	: 4clan combination code + 'total N' over all 4clans
# value : list of rejected quartet-numbers
my	%hol_qn_of_qfiltered_otu_of_4clan_of_ft ;
#########################################

#########################################
# key1	: 4seq filtermethod (e.g. nap_risk1)
# key2	: taxon
# key3	: 4clan combination code + 'total N' over all 4clans
# value : N remaining quartets
my	%hoh_qremaining_of_otu_of_4clan_of_ft ;
#########################################

#########################################
# key1	: quartet-number ($fc)
# value : quartet related sequence length (corrected)
my	%seqL_of_fc ;
#########################################

#########################################
# key1	: quartet-number ($fc)
# key2	: calculation method (nap, napv, napr)
# value : quartet related best supported tree
my	%bestTree_of_fc_of_cm ;
#########################################

#########################################
# key1	: calculation method + filter method (nap_risk1, napv_dist2, napr_dist3) 
# key2	: quartet-number ($fc)
# value : seen +1
my	%seen_remain_of_cmfm_of_fc ;
#########################################

#########################################
# key1	: 4clan combination code 
# value	: list of assigned quartet-numbers ($fc)
my	%hol_fc_of_clan4 ;
#########################################

#########################################
# key1	: 4clan tree
# value	: 4clan combination code 
my	%clan4_of_tree ;
#########################################

#########################################
# key1	: filter method incl OTU clipping approaches
# value	: seen +1
my	%seen_filter_method ;
#########################################

#########################################
# key1	: quartet-number ($fc)
# key2	: calculation method + filter method (nap_risk1, napv_dist2, napr_dist3) 
# key3	: quartet-tree
# value	: rNc/Na
my	%rNcNa_of_cmfm_of_fc_of_tree ;
#########################################

#########################################
for my $clan4 ( sort keys %hol_fp_qcal_of_4clan_code ){
	
	
	###################################################
	# store all 4clan corresponding qcalc filepaths as reference-list in @$aref_fp
	my	$aref_fp = \@{$hol_fp_qcal_of_4clan_code{$clan4}} ;
	my	$Nq_4clan_total = @$aref_fp ;
	###################################################
	
	###################################################
	# store actual 4clan quartet numbers ($fc) as hashkey
	my	%fc_of_4clan = () ;
	###################################################
	
	###################################################
	# key1	: type of supportfilter ('nap_risk1', 'napr_risk1', or 'napv_risk1' ...)
	# key2	: threshold scale
	# key3	: quartet-number
	# value	: quartet-filter corresponding score value (1 -> reject quartet)
	my	%hoh_qrej_of_filtertype_of_scale_of_fn ;
	###################################################
	
	
	
	#------------------- Analysis of single-quartet split pattern frequencies ----------------------------
	
	######################################################################################################
	# Analyse SPD distribution of each single quartet-calculation file                                   #
	######################################################################################################
	for my $p_qcalc_file ( sort @$aref_fp ){
			
		
		###################################################
		# Set Filecounter +1
		$fc++ ;
		$fc_of_4clan{$fc}++  ;
		push @{$hol_fc_of_clan4{$clan4}}, $fc ;
		
		# count total number of clan4 analysed quartets
		$seen_clan4{$clan4}++ ;
		###################################################
		
		
		
		######################################################################################################
		## ReadIN '*_detailed_split_calc_*'                                                                 ##
		######################################################################################################
		
		###################################################
		# read IN qcalc-file '*_detailed_split_calc_*' and store each line in @l_qcalcfile
		my @l_qcalcfile = &readIN_file( \$p_qcalc_file ) ;
		#for my $l (@l_qcalcfile){print "\n", $l} exit;
		######################################################################################################
		
		
		
		######################################################################################################
		## extract taxon combination from penguin '*_detailed_split_calc_*'                                 ##
		######################################################################################################
		
		###################################################
		# Extract the quartet corresponding 4taxon combination from the 1st line of '*_detailed_split_calc_*'
		# Store all four taxa in @qtaxa (same order as in the fileline)
		print 	"\n\t----------------------------------\n"	,
				"\tanalysing ", $clan4, " quartet ", $fc, "\t<", $l_qcalcfile[0] , ">\n" ; 
		
		#########################################
		# list of quartet related taxon-combination (taxa are ordered as they appear in the first line of the infile)
		my	@qtaxa	= () ;
		#########################################
		
		&extract_4taxon_combination_from_fileline(
			
			\$l_qcalcfile[0]						, # fileline with 4-taxon info (must be the first line in '*_detailed_split_calc_*'						IN: defined; OUT: unchanged
			\@qtaxa									, # list of quartet depending taxa in the same order as in the infile									IN: undefined; OUT: defined
			\%taxa4_of_fc							, # key1: quartet-number ($fc); value : quartet assigned taxa (alphabetically ordered, joined by ':'	IN: undefined; OUT: defined
			\$fc									, # actual quartet-number																				IN: defined; OUT: unchanged
			\%hoh_seen_quartets_of_taxon_of_clan4	, # key1: 4clan combination code; key2: taxon; value: N analysed quartets								IN: undefined; OUT: defined
			\$clan4									, # Clan4 of given quartet combination																	IN: defined; OUT: unchanged
			\%clan_of_seq							, # hash with sequence name (key) assigned clan (value)													IN: defined; OUT: unchanged
			\%value_of_scr_par						, # key1: sealion user available parameter options; value: optional parameter setting					IN: defined; OUT: unchanged
			\$p_qcalc_file							, # actual spd infile																					IN: defined; OUT: unchanged
		);
		######################################################################################################
		
		
		
		######################################################################################################
		## generate split-pattern code                                                                      ##
		######################################################################################################
			
		###################################################
		# assigning split pattern due to penguin defined split pattern codes for terminal print of pattern assignment
		# key: pattern code
		# value: correponding taxon split
		# e.g. key: A  value: T1,T2||T3,T4 
		if ( $value_of_scr_par{sys_parameter}{screen_prt} == 1 ){
			
			&assign_split_pattern_to_code(
				
				\@qtaxa								, # list of quartet depending taxa in the same order as in the infile	IN: defined; OUT: unchanged
				\%clan_of_seq						, # hash with sequence name (key) assigned clan (value)					IN: defined; OUT: unchanged
			) ;
		}
		######################################################################################################
		
		
		
		######################################################################################################
		## assign split-pattern frequencies in '*_detailed_split_calc_*' to pattern code                    ##
		######################################################################################################
		
		###################################################
		# pattern distribution (n) observed and expected is always located in infile line 7 to 21, tabstop delimited:
		my	(
				%obs_n_of_pattern 					, # key: pattern code (e.g. 'A'); value number of observed pattern
				%exp_n_of_tree1_pattern 			, # key: pattern code (e.g. 'A'); value number of expected pattern given tree 1
				%exp_n_of_tree2_pattern 			, # key: pattern code (e.g. 'A'); value number of expected pattern given tree 2
				%exp_n_of_tree3_pattern 			, # key: pattern code (e.g. 'A'); value number of expected pattern given tree 3
				%obs_p_of_pattern 					, # key: pattern code (e.g. 'A'); value percentage of observed pattern
				%exp_p_of_tree1_pattern 			, # key: pattern code (e.g. 'A'); value percentage of expected pattern given tree 1
				%exp_p_of_tree2_pattern 			, # key: pattern code (e.g. 'A'); value percentage of expected pattern given tree 2
				%exp_p_of_tree3_pattern 			, # key: pattern code (e.g. 'A'); value percentage of expected pattern given tree 3
				$Npat_total_obs			  			, # Total number of different split pattern
		) ;
		
		&read_write_spd_files(
			
			\@l_qcalcfile							, # list of infile-lines of '*_detailed_split_calc_*'												IN: defined; OUT: unchanged
			\%obs_n_of_pattern 						, # key: pattern code (e.g. 'A'); value number of observed pattern frequency						IN: undefined; OUT: defined
			\%exp_n_of_tree1_pattern 				, # key: pattern code (e.g. 'A'); value number of expected pattern frequency given tree 1			IN: undefined; OUT: defined
			\%exp_n_of_tree2_pattern 				, # key: pattern code (e.g. 'A'); value number of expected pattern frequency given tree 2			IN: undefined; OUT: defined
			\%exp_n_of_tree3_pattern 				, # key: pattern code (e.g. 'A'); value number of expected pattern frequency given tree 3			IN: undefined; OUT: defined
			\%obs_p_of_pattern 						, # key: pattern code (e.g. 'A'); value percentage of observed pattern frequency					IN: undefined; OUT: defined
			\%exp_p_of_tree1_pattern 				, # key: pattern code (e.g. 'A'); value percentage of expected pattern frequency given tree 1		IN: undefined; OUT: defined
			\%exp_p_of_tree2_pattern 				, # key: pattern code (e.g. 'A'); value percentage of expected pattern frequency given tree 2		IN: undefined; OUT: defined
			\%exp_p_of_tree3_pattern 				, # key: pattern code (e.g. 'A'); value percentage of expected pattern frequency given tree 3		IN: undefined; OUT: defined
			\$Npat_total_obs						, # Total number of different split pattern															IN: undefined; OUT: defined
			\%value_of_scr_par						, # key1: sealion user available parameter options; value: optional parameter setting				IN: defined; OUT: unchanged
			\%subfolder_q							, # key1: output-type (e.g. 'pdf'); value: pathdir													IN: defined; OUT: unchanged
			\@qtaxa									, # list of quartet depending taxa in the same order as in the infile								IN: defined; OUT: unchanged
			\%clan_of_seq							, # key1: sequence name; value: assigned clan 														IN: defined; OUT: unchanged
			\$fc									, # actual quartet-number																			IN: defined; OUT: unchanged
			\%seqL_of_fc							, # key1: quartet-number ($fc); value: quartet related sequence length (corrected)					IN: undefined; OUT: defined
			\$p_qcalc_file							, # actual spd infile																				IN: defined; OUT: unchanged
		) ;
		######################################################################################################
		
		
		
		######################################################################################################
		## Generating outgroup rooted & non-rooted quartet-trees and list sampling in pattern signal order  ##
		######################################################################################################
		
		###################################################
		# produce all three outgroup polarized 4clan trees based on specified outgroup clan $outgroup_clan
		# NOTICE! The order of the rooted quartet-trees is in the same as order as the quartet-information!
		# DON'T CHANGE ORDER OF @newick_rooted
		
		#########################################
		# list of polarized, outgroup-clan rooted quartet trees with names in alphabetical order
		# List order of trees corresponds with tree position linked split-pattern evaluation
		my	@newick_rooted = () ;
		#########################################
		
		#########################################
		# list of polarized, unrooted quartet trees with names in alphabetical order, used for box plots of ML expected scores...
		# ...'pbox_CFexpTr1_BL2_besttrees'
		# ...'pbox_CFexpTr2_BL2_besttrees', and
		# ...'pbox_CFexpTr3_BL2_besttrees'
		my	@newick_unrooted = () ;
		#########################################
		
		#########################################
		# generate the three outgroup polarized, unrooted 4-clan trees
		&build_quartet_trees(
			
			\@qtaxa									, # list of quartet depending taxa in the same order as in the infile	IN: defined; OUT: unchanged
			\%clan_of_seq							, # hash with sequence name (key) assigned clan (value)					IN: defined; OUT: unchanged
			\$outgroup_clan							, # clan name of the outgroup 											IN: defined; OUT: unchanged
			\@newick_unrooted						, # list of unrooted quartet trees										IN: undefined; OUT: defined
			\@newick_rooted							, # list of rooted quartet trees										IN: undefined; OUT: defined
			\%value_of_scr_par						, # key1: sealion user available parameter options; value: setting		IN: defined; OUT: unchanged
			\%clan4_of_tree							, # key1: 4clan tree; value: 4clan code									IN: undefined; OUT: defined
			\$clan4									, # Clan4 of given quartet combination									IN: defined; OUT: unchanged
		) ;
		#########################################
		
		#########################################
		# check if length of 4clan newick strings is acceptable in grafic prints
		&check_newick_chars(
			
			\%value_of_scr_par						, # key1: sealion user available parameter options; value: optional parameter setting	IN: defined, OUT -> changed
			\$newick_unrooted[0]					, # list of unrooted quartet trees														IN: defined, OUT -> unchanged
			\$clan4									, # Clan4 of given quartet combination													IN: defined, OUT -> unchanged
		) ;
		######################################################################################################

		
		
		######################################################################################################
		## Split-pattern evaluation and Support calculation of actual quartet combination                   ##
		######################################################################################################
		
		#########################################
		# 1) Split pattern evaluation of actual 4-clan depending 4taxon combination
		# Splits are evaluated due @newick_rooted tree order
		# Further evaluated data are stored in corresponding in-variables
		&evaluation_quartet_pattern_main(
			
			\@newick_rooted										, # list of rooted quartet trees																																				IN: defined, OUT: unchanged
			\@qtaxa												, # list of quartet depending taxa in the same order as in the infile																											IN: defined, OUT: unchanged
			\%obs_n_of_pattern									, # observed pattern frequencies (key: pattern code; value: frequency number																									IN: defined, OUT: unchanged
			\%exp_n_of_tree1_pattern							, # pattern frequencies expected given tree1 (key: pattern code; value: frequency number																						IN: defined, OUT: unchanged
			\%exp_n_of_tree2_pattern							, # pattern frequencies expected given tree2 (key: pattern code; value: frequency number																						IN: defined, OUT: unchanged
			\%exp_n_of_tree3_pattern							, # pattern frequencies expected given tree3 (key: pattern code; value: frequency number																						IN: defined, OUT: unchanged
			\%obs_p_of_pattern									, # observed pattern percentages (key: pattern code; value: frequency number																									IN: defined, OUT: unchanged
			\%exp_p_of_tree1_pattern							, # pattern percentages expected given tree1 (key: pattern code; value: frequency number																						IN: defined, OUT: unchanged
			\%exp_p_of_tree2_pattern							, # pattern percentages expected given tree2 (key: pattern code; value: frequency number																						IN: defined, OUT: unchanged
			\%exp_p_of_tree3_pattern							, # pattern percentages expected given tree3 (key: pattern code; value: frequency number																						IN: defined, OUT: unchanged
			\$Npat_total_obs									, # total number of (gap excluded) site positions																																IN: defined, OUT: unchanged
			\%clan_of_seq										, # hash with sequence name (key) assigned clan (value)																															IN: defined, OUT: unchanged
			\$outgroup_clan										, # clan name of the outgroup																																					IN: defined, OUT: unchanged
			\$clan4												, # Clan4 of given quartet combination																																			IN: defined, OUT: unchanged
			\$fc												, # actual sealion quartet number																																				IN: defined, OUT: unchanged
			\%hoh_sup_of_clan4_of_calc_of_qtree_of_fc			, # multiple hash, storing the final support (nap or napv) related to each rooted quartet-tree																					IN: defined, OUT: changed
			\%hol_sup_of_clan4_of_calc_of_qtree					, # multiple hash, listing for each of the three rooted quartet-trees of given 4clan combination all single quartet inferred support (NAP or NAPV)								IN: defined, OUT: changed
			\%hoh_tsvdata_of_tsvfilename_of_printcode			, # key1: R grafic unique printcode (e.g. 'metafileT'); key2: *.tsv filename; value: *tsv filename data addings (e.g. clanfilename, type of support...)							IN: defined, OUT: changed
			\%subfolder_q										, # hashpath to output-type specific subdirs of single quartet-pattern analyses; key: output-type (e.g. 'pdf'); value: pathdir													IN: defined, OUT: unchanged
			\%value_of_scr_par									, # key1: sealion user available parameter options; value: optional parameter setting																 							IN: defined, OUT: unchanged
			\%hoh_sup_dist_of_clan4_type_of_treepair			, # key1: clan4 code; key2: type of support ('nap', 'napr', or 'napv'); key3: tree-pair code ('first2nd' or 'first3rd'); value: distance Nap									IN: undefined, OUT: defined 
			\%hoh_qrej_of_filtertype_of_scale_of_fn				, # key1: type of supportfilter ('nap_risk1', 'napr_risk1', or 'napv_risk1'); key2: scale value; key3: quartet number; value: filter corresponding score-value					IN: undefined, OUT: defined
			\%bestTree_of_fc_of_cm								, # key1: quartet-number ($fc), key2: calculation method (nap, napv, napr); value: quartet related best supported tree															IN: undefined, OUT: defined 
			\%rNcNa_of_cmfm_of_fc_of_tree						, # key1: quartet-number ($fc), key2: calculation+filter method (nap_risk1, napv_dist2, napr_dist3); key3: tree; value: rNc/na													IN: undefined, OUT: defined 
			\%algorithm_def_methods								, # list of defined support methods (nap, napv ...)																																IN: defined; OUT: unchanged
			\%hol_scales_of_filtertype							, # key1: type of supportfilter ('risk', 'qdist'); value: list of ftype corresponding threshold scales																			IN: defined; OUT: unchanged
			\%seen_filter_method								, # key1: filter method incl OTU clipping approaches; value: seen +1																											IN: undefined; OUT: defined
		) ;
		######################################################################################################
	}
	######################################################################################################
	
	print 	"\n\t----------------------------------\n"	,
			"\tanalysing single quartet support for 4clan trees of ", $clan4 ,
			"\n\n\t\t-all quartets\n" 				; 
	
	
	if ( ( $value_of_scr_par{seq_filterRisk}{risk1}  == 1 ) ||
		 ( $value_of_scr_par{seq_filterRisk}{risk2}  == 1 ) ||
		 ( $value_of_scr_par{seq_filterDist}{qdist1} == 1 ) ){
		 	 
		&filter_threshold_optimization(
			
			\%subfolder_q											, # hashpath to output-type specific subdirs of single quartet-pattern analyses; key: output-type (e.g. 'pdf'); value: pathdir													IN: defined; OUT: unchanged
			\%hoh_tsvdata_of_tsvfilename_of_printcode				, # key1: R grafic unique printcode (e.g. 'metafileT'); key2: *.tsv filename; value: *tsv filename data addings (e.g. clanfilename, type of support...)							IN: defined; OUT: changed
			\%hoh_linedata_of_latexfilename_of_printcode			, # key1: Latex table unique printcode (e.g. 'latexC'); key2: *.tex filename; value: @list of assigned table lines																IN: defined; OUT: changed
			\%textfilename_of_printcode								, # key1: textfile unique printcode (e.g. 'textA'); value: filename for print																									IN: defined; OUT: unchanged
			\%value_of_scr_par										, # key1: sealion user available parameter options; value: optional parameter setting																 							IN: defined, OUT: unchanged
			\$clan4													, # Clan4 of given quartet combination																																			IN: defined, OUT: unchanged
			\%hol_fc_of_clan4										, # key1: 4clan code; value: list of 4clan analysed quartet numbers																												IN: defined, OUT: unchanged
			\%hoh_qrej_of_filtertype_of_scale_of_fn					, # key1: type of supportfilter ('nap_risk1', 'napr_risk1', or 'napv_risk1'); key2: scale value; key3: quartet number; value: filter corresponding score-value					IN: defined, OUT: undefined
			\%hoh_seen_rejected_of_4clan_of_filtertype_of_fn		, # key1: 4clan code; key2: type of supportfilter ('nap_risk1', 'napr_risk1', or 'napv_risk1'); key3: actual sealion quartet number; value: filter corresponding score-value	IN: undefined; OUT: defined 
			\%hol_scales_of_filtertype								, # key1: type of supportfilter ('risk', 'qdist'); value: list of ftype corresponding threshold scales																			IN: defined; OUT: unchanged
			\%seen_filter_method									, # key1: filter method approaches; value: seen +1																																IN: defined; OUT: unchanged
		) ;
		
		&print_tsv_Q3b(
			
			\%subfolder_q											, # hashpath to output-type specific subdirs of single quartet-pattern analyses; key: output-type (e.g. 'pdf'); value: pathdir													IN: defined; OUT: unchanged
			\%hoh_tsvdata_of_tsvfilename_of_printcode				, # key1: R grafic unique printcode (e.g. 'metafileT'); key2: *.tsv filename; value: *tsv filename data addings (e.g. clanfilename, type of support...)							IN: defined; OUT: changed
			\$clan4													, # Clan4 of given quartet combination																																			IN: defined; OUT: unchanged
			\%fc_of_4clan											, # key: quartet numbers of current 4clan combination																															IN: defined; OUT: unchanged
			\%hoh_seen_rejected_of_4clan_of_filtertype_of_fn		, # key1: 4clan code; key2: type of supportfilter ('nap_risk1', 'napr_risk1', or 'napv_risk1'); key3: actual sealion quartet number; value: filter corresponding score-value	IN: defined; OUT: unchanged
			\%seen_filter_method									, # key1: filter method approaches; value: seen +1																																IN: defined; OUT: unchanged
		) ;
	}
	
	#------------------- Quartet Filtering ---------------------------------------------------------------
	
	######################################################################################################
	# Filter 'risk' or 'dist' evaluated support for each type (napv, nap, napr)      					 #
	######################################################################################################
	
	#####################################################################################
	# Filter 'risk' evaluated support for each type (napv, nap, napr)                   # 
	#####################################################################################
	my	(
			%hoh_seen_filter_of_type 	# key1: filter ('dist' or 'risk'); value: list of combis method+filter ('nap_dist1', 'nap_risk2'...)
	) ;
	
	if ( ( $value_of_scr_par{seq_filterRisk}{risk1} == 1 ) ||
		 ( $value_of_scr_par{seq_filterRisk}{risk2} == 1 ) ){
			 
		print	"\t\t-filtered quartets ('risk' option)\n" ;
	
		&filter_quartets_risk(
			
			\$clan4												, # Clan4 of given quartet combination																																			IN: defined; OUT: unchanged
			\%hol_sup_of_clan4_of_calc_of_qtree					, # multiple hash, listing for each of the three rooted quartet-trees of given 4clan combination all single quartet inferred support (NAP or NAPV)								IN: defined; OUT: changed
			\%hoh_sup_of_clan4_of_calc_of_qtree_of_fc			, # multiple hash, storing the final support (nap or napv) related to each rooted quartet-tree																					IN: defined; OUT: changed
			\%hoh_seen_rejected_of_4clan_of_filtertype_of_fn	, # key1: 4clan code; key2: type of supportfilter ('nap_risk1', 'napr_risk1', or 'napv_risk1'); key3: actual sealion quartet number; value: filter corresponding score-value	IN: defined; OUT: unchanged 
			\%value_of_scr_par									, # key1: sealion user available parameter options; value: optional parameter setting																 							IN: defined; OUT: unchanged
			\%algorithm_def_methods								, # list of defined support methods (nap, napv ...)																																IN: defined; OUT: changed
			\%hoh_seen_filter_of_type							, # key1: filter ('dist' or 'risk'); key2: combi-code method+filter ('nap_dist1', 'nap_dist2'...); value: number of filter passing quartets										IN: defined; OUT: changed
			\%hol_qn_of_qfiltered_otu_of_4clan_of_ft			, # key1: 4seq filtermethod (e.g. nap_risk1); key2: taxon; key3: 4clan|total; value: list of rejected qn's																		IN: undefined; OUT: defined
			\%hoh_qremaining_of_otu_of_4clan_of_ft				, # key1: 4seq filtermethod (e.g. nap_risk1); key2: taxon; key3: 4clan|total; value: N remaining quartets																		IN: undefined; OUT: defined
			\%taxa4_of_fc										, # key1: quartet-number ($fc); value : quartet assigned taxa (alphabetically ordered, joined by ':')																			IN: defined; OUT: unchanged
			\%hoh_N_qrej_of_clan4_of_ft							, # key1: 4seq filtermethod (e.g. nap_risk1); key2: 4clan combination code; value:  N rejected quartets*3 (counted three times, each time for each quartet tree)				IN: undefined; OUT: defined
			\%seen_remain_of_cmfm_of_fc							, # key1: calculation method + filter method (nap_risk1, napv_dist2, napr_dist3); key2: quartet-number ($fc); value: seen +1													IN: undefined; OUT: defined
		) ;
	}
	#####################################################################################
	
	
	#####################################################################################
	# Filter 'dist' evaluated support for each type (napv, nap, napr)                   # 
	#####################################################################################
	if ( ( $value_of_scr_par{seq_filterDist}{qdist1} == 1 ) ){
		 	 
		print	"\t\t-filtered quartets ('dist' option)\n" ;
				
		&filter_quartets_dist(
			
			\$clan4												, # Clan4 of given quartet combination																																			IN: defined; OUT: unchanged
			\%hol_sup_of_clan4_of_calc_of_qtree					, # multiple hash, listing for each of the three rooted quartet-trees of given 4clan combination all single quartet inferred support (NAP or NAPV)								IN: defined; OUT: changed
			\%hoh_sup_of_clan4_of_calc_of_qtree_of_fc			, # multiple hash, storing the final support (nap or napv) related to each rooted quartet-tree																					IN: defined; OUT: changed
			\%hoh_seen_rejected_of_4clan_of_filtertype_of_fn	, # key1: 4clan code; key2: type of supportfilter ('nap_risk1', 'napr_risk1', or 'napv_risk1'); key3: actual sealion quartet number; value: filter corresponding score-value	IN: defined; OUT: unchanged 
			\%value_of_scr_par									, # key1: sealion user available parameter options; value: optional parameter setting																 							IN: defined; OUT: unchanged
			\%algorithm_def_methods								, # list of defined support methods (nap, napv ...)																																IN: defined; OUT: changed
			\%hoh_seen_filter_of_type							, # key1: filter ('dist' or 'risk'); key2: combi-code method+filter ('nap_dist1', 'nap_dist2'...); value: number of filter passing quartets										IN: defined; OUT: changed
			\%hol_qn_of_qfiltered_otu_of_4clan_of_ft			, # key1: 4seq filtermethod (e.g. nap_risk1); key2: taxon; key3: 4clan|total; value: list of rejected qn's																		IN: undefined; OUT: defined
			\%hoh_qremaining_of_otu_of_4clan_of_ft				, # key1: 4seq filtermethod (e.g. nap_risk1); key2: taxon; key3: 4clan|total; value: N remaining quartets																		IN: undefined; OUT: defined
			\%taxa4_of_fc										, # key1: quartet-number ($fc); value : quartet assigned taxa (alphabetically ordered, joined by ':')																			IN: defined; OUT: unchanged
			\%hoh_N_qrej_of_clan4_of_ft							, # key1: 4seq filtermethod (e.g. nap_risk1); key2: 4clan combination code; value:  N rejected quartets*3 (counted three times, each time for each quartet tree)				IN: undefined; OUT: defined
			\%seen_remain_of_cmfm_of_fc							, # key1: calculation method + filter method (nap_risk1, napv_dist2, napr_dist3); key2: quartet-number ($fc); value: seen +1													IN: undefined; OUT: defined
		) ;
	}
	#####################################################################################
	
	#####################################################################################
	# Filter both, 'risk' + 'dist', evaluated support for each type (napv, nap, napr)   # 
	#####################################################################################
	if ( ( $hoh_seen_filter_of_type{risk} ) &&
		 ( $hoh_seen_filter_of_type{dist} ) ){
			 
		print	"\t\t-filtered quartets ('dist' + 'risk' option)\n" ;
				
		&filter_quartets_riskdist(
			
			\$clan4												, # Clan4 of given quartet combination																																			IN: defined; OUT: unchanged
			\%hol_sup_of_clan4_of_calc_of_qtree					, # multiple hash, listing for each of the three rooted quartet-trees of given 4clan combination all single quartet inferred support (NAP or NAPV)								IN: defined; OUT: changed
			\%hoh_sup_of_clan4_of_calc_of_qtree_of_fc			, # multiple hash, storing the final support (nap or napv) related to each rooted quartet-tree																					IN: defined; OUT: changed
			\%hoh_seen_rejected_of_4clan_of_filtertype_of_fn	, # key1: 4clan code; key2: type of supportfilter ('nap_risk1', 'napr_risk1', or 'napv_risk1'); key3: actual sealion quartet number; value: filter corresponding score-value	IN: defined; OUT: unchanged 
			\%value_of_scr_par									, # key1: sealion user available parameter options; value: optional parameter setting																 							IN: defined; OUT: unchanged
			\%algorithm_def_methods								, # list of defined support methods (nap, napv ...)																																IN: defined; OUT: changed
			\%hoh_seen_filter_of_type							, # key1: filter ('dist' or 'risk'); key2: combi-code method+filter ('nap_dist1', 'nap_dist2'...); value: number of filter passing quartets										IN: defined; OUT: changed
			\%hol_qn_of_qfiltered_otu_of_4clan_of_ft			, # key1: 4seq filtermethod (e.g. nap_risk1); key2: taxon; key3: 4clan|total; value: list of rejected qn's																		IN: defined; OUT: changed
			\%hoh_qremaining_of_otu_of_4clan_of_ft				, # key1: 4seq filtermethod (e.g. nap_risk1); key2: taxon; key3: 4clan|total; value: N remaining quartets																		IN: defined; OUT: changed
			\%taxa4_of_fc										, # key1: quartet-number ($fc); value : quartet assigned taxa (alphabetically ordered, joined by ':')																			IN: defined; OUT: unchanged
			\%hoh_N_qrej_of_clan4_of_ft							, # key1: 4seq filtermethod (e.g. nap_risk1); key2: 4clan combination code; value:  N rejected quartets*3 (counted three times, each time for each quartet tree)				IN: defined; OUT: changed
			\%seen_remain_of_cmfm_of_fc							, # key1: calculation method + filter method (nap_risk1, napv_dist2, napr_dist3); key2: quartet-number ($fc); value: seen +1													IN: undefined; OUT: defined
			\%seen_filter_method								, # key1: filter method; value: seen +1																																			IN: defined; OUT: changed
		) ;
	}
	#####################################################################################
	
	######################################################################################################
	
	
	
	#------------------- Calculate 4clan tree support and table print ------------------------------------
	
	######################################################################################################
	## Calculate & print overall mean & median 4-clan tree support for each type (napv, nap, napr)      ##
	######################################################################################################
	
	#####################################################################################
	# calculate, print, and store mean and median support for each 4clan ccombination   # 
	#####################################################################################
	&calc_4clan_tree_support(	
		
		\$clan4												, # Clan4 of given quartet combination																													IN: defined; OUT: unchanged
		\$Nq_4clan_total									, # total number of single 4clan quartet analyses (unfiltered)																							IN: defined; OUT: unchanged
		\%hol_sup_of_clan4_of_calc_of_qtree					, # multiple hash, listing for each of the three rooted quartet-trees of given 4clan combination all single quartet inferred support (NAP or NAPV)		IN: defined; OUT: unchanged
		\%hoh_sup_of_clan4_of_calc_of_qtree_of_fc			, # multiple hash, storing the final support (nap or napv) related to each rooted quartet-tree															IN: defined; OUT: changed
		\%hol_sup_of_type_of_mfilter						, # key1: type of supportfilter; key2: mean|median; value: list of 'tree\tmean|medianScore'																IN: undefined; OUT: defined
		\%hoh_linedata_of_latexfilename_of_printcode		, # key1: Latex table unique printcode (e.g. 'latexC'); key2: *.tex filename; value: @list of assigned table lines										IN: defined; OUT: changed
		\%value_of_scr_par									, # key1: sealion user available parameter options; value: optional parameter setting																 	IN: defined; OUT: unchanged
		\%subfolder_q										, # hashpath to output-type specific subdirs of single quartet-pattern analyses; key: output-type (e.g. 'pdf'); value: pathdir							IN: defined; OUT: unchanged
		\%textfilename_of_printcode							, # key1: textfile unique printcode (e.g. 'textA'); value: filename for print																			IN: defined; OUT: unchanged
	);	
	#####################################################################################
	
	
	#####################################################################################
	# print single quartet and average 4-clan tree support (napv, nap, napr)            #
	#####################################################################################
	&print_tsv_MQ2(
		
		\%subfolder_q										, # hashpath to output-type specific subdirs of single quartet-pattern analyses; key: output-type (e.g. 'pdf'); value: pathdir													IN: defined; OUT: unchanged
		\%hoh_tsvdata_of_tsvfilename_of_printcode			, # key1: R grafic unique printcode (e.g. 'metafileT'); key2: *.tsv filename; value: *tsv filename data addings (e.g. clanfilename, type of support...)							IN: defined; OUT: changed
		\$clan4												, # Clan4 of given quartet combination																																			IN: defined; OUT: unchanged
		\%fc_of_4clan										, # key: quartet numbers of current 4clan combination																															IN: defined; OUT: unchanged
		\%hoh_sup_of_clan4_of_calc_of_qtree_of_fc			, # multiple hash, storing the final support (nap or napv) related to each rooted quartet-tree																					IN: defined; OUT: unchanged
		\%value_of_scr_par									, # key1: sealion user available parameter options; value: optional parameter setting																 							IN: defined; OUT: unchanged
		\%hoh_seen_rejected_of_4clan_of_filtertype_of_fn	, # key1: 4clan code; key2: type of supportfilter ('nap_risk1', 'napr_risk1', or 'napv_risk1'); key3: actual sealion quartet number; value: filter corresponding score-value	IN: defined; OUT: unchanged 
		\%algorithm_def_methods								, # list of defined support methods (nap, napv ...)																																IN: defined; OUT: unchanged 
		\%hoh_seen_filter_of_type							, # key1: filter ('dist' or 'risk'); key2: combi-code method+filter ('nap_dist1', 'nap_dist2'...); value: number of filter passing quartets										IN: defined; OUT: unchanged
	);
	#####################################################################################
	
	
	#####################################################################################
	# print single quartet and average 4-clan tree support (napv, nap, napr)            #
	#####################################################################################
	&print_tsv_MQ4(
		
		\%subfolder_q										, # hashpath to output-type specific subdirs of single quartet-pattern analyses; key: output-type (e.g. 'pdf'); value: pathdir													IN: defined; OUT: unchanged
		\%hoh_tsvdata_of_tsvfilename_of_printcode			, # key1: R grafic unique printcode (e.g. 'metafileT'); key2: *.tsv filename; value: *tsv filename data addings (e.g. clanfilename, type of support...)							IN: defined; OUT: changed
		\$clan4												, # Clan4 of given quartet combination																																			IN: defined; OUT: unchanged
		\%fc_of_4clan										, # key: quartet numbers of current 4clan combination																															IN: defined; OUT: unchanged
		\%hoh_sup_of_clan4_of_calc_of_qtree_of_fc			, # multiple hash, storing the final support (nap or napv) related to each rooted quartet-tree																					IN: defined; OUT: unchanged
		\%value_of_scr_par									, # key1: sealion user available parameter options; value: optional parameter setting																 							IN: defined; OUT: unchanged
		\%hoh_seen_rejected_of_4clan_of_filtertype_of_fn	, # key1: 4clan code; key2: type of supportfilter ('nap_risk1', 'napr_risk1', or 'napv_risk1'); key3: actual sealion quartet number; value: filter corresponding score-value	IN: defined; OUT: unchanged 
		\%algorithm_def_methods								, # list of defined support methods (nap, napv ...)																																IN: defined; OUT: unchanged 
		\%hoh_seen_filter_of_type							, # key1: filter ('dist' or 'risk'); key2: combi-code method+filter ('nap_dist1', 'nap_dist2'...); value: number of filter passing quartets										IN: defined; OUT: unchanged
	);
	#####################################################################################
	
	######################################################################################################
}
######################################################################################################

#exit;


#------------------- Global table prints of SPD calculations -----------------------------------------

######################################################################################################
# Global table prints of 4SEQ split signal calculations                                              #
######################################################################################################

#####################################################################################
# sampling latex data lines LQ3 (N quartet particpations of each taxon, filtered)   #
#####################################################################################
&latex_data_sample_LQ3(
						
		\%hoh_linedata_of_latexfilename_of_printcode		, # key1: Latex table unique printcode (e.g. 'latexC'); key2: *.tex filename; value: @list of assigned table lines																	IN: defined; OUT: changed
		\%hol_qn_of_qfiltered_otu_of_4clan_of_ft			, # key1: 4seq filtermethod (e.g. nap_risk1); key2: taxon; key3: 4clan|total; value: list of rejected qn's																			IN: defined; OUT: unchanged
		\%hoh_qremaining_of_otu_of_4clan_of_ft				, # key1: 4seq filtermethod (e.g. nap_risk1); key2: taxon; key3: 4clan|total; value: N remaining quartets																			IN: defined; OUT: unchanged
		\%clan_of_seq										, # key1: sequence-name; value: clan definition																																		IN: defined; OUT: unchanged
		\%seen_clan4										, # key1: clan4; value: N analyseed quartets																																		IN: defined; OUT: unchanged
		\%subfolder_q										, # hashpath to output-type specific subdirs of single quartet-pattern analyses; key: output-type (e.g. 'pdf'); value: pathdir														IN: defined; OUT: unchanged
		\%textfilename_of_printcode							, # key1: textfile unique printcode (e.g. 'textA'); value: filename for print																										IN: defined; OUT: unchanged
		\%value_of_scr_par									, # key1: sealion user available parameter options; value: optional parameter setting																 								IN: defined; OUT: unchanged
		\%seen_filter_method								, # key1: filter method incl OTU clipping approaches; value: seen +1																												IN: defined; OUT: unchanged
) ;
#####################################################################################


#####################################################################################
# sampling latex data lines LQ4 (N quartet particpations of each taxon, unfiltered) # 
#####################################################################################
&latex_data_sample_LQ4(
						
		\%hoh_linedata_of_latexfilename_of_printcode		, # key1: Latex table unique printcode (e.g. 'latexC'); key2: *.tex filename; value: @list of assigned table lines																	IN: defined; OUT: changed
		\%hoh_seen_quartets_of_taxon_of_clan4				, # key1: 4clan combination code; key2: taxon; value: N analysed quartets																											IN: defined; OUT: unchanged
		\%clan_of_seq										, # key1: sequence-name; value: clan definition																																		IN: defined; OUT: unchanged
		\%algorithm_def_methods								, # key1: methode type; key2: method code; value: integer																															IN: defined; OUT: unchanged
		\%subfolder_q										, # hashpath to output-type specific subdirs of single quartet-pattern analyses; key: output-type (e.g. 'pdf'); value: pathdir														IN: defined; OUT: unchanged
		\%textfilename_of_printcode							, # key1: textfile unique printcode (e.g. 'textA'); value: filename for print																										IN: defined; OUT: unchanged
) ;
#####################################################################################


#####################################################################################
# sampling latex data lines LQ5 (un-analysed taxa, unfiltered)                      # 
#####################################################################################
&latex_data_sample_LQ5(
						
		\%hoh_linedata_of_latexfilename_of_printcode		, # key1: Latex table unique printcode (e.g. 'latexC'); key2: *.tex filename; value: @list of assigned table lines																	IN: defined; OUT: changed
		\%hoh_seen_quartets_of_taxon_of_clan4				, # key1: 4clan combination code; key2: taxon; value: N analysed quartets																											IN: defined; OUT: unchanged
		\%clan_of_seq										, # key1: sequence-name; value: clan definition																																		IN: defined; OUT: unchanged
		\%algorithm_def_methods								, # key1: methode type; key2: method code; value: integer																															IN: defined; OUT: unchanged
		\%subfolder_q										, # hashpath to output-type specific subdirs of single quartet-pattern analyses; key: output-type (e.g. 'pdf'); value: pathdir														IN: defined; OUT: unchanged
		\%textfilename_of_printcode							, # key1: textfile unique printcode (e.g. 'textA'); value: filename for print																										IN: defined; OUT: unchanged
) ;
#####################################################################################


#####################################################################################
# sampling tsv data lines Q5 (N rejected quartets for each taxon)                   # 
#####################################################################################
&print_tsv_Q5(
						
		\%subfolder_q										, # hashpath to output-type specific subdirs of single quartet-pattern analyses; key: output-type (e.g. 'pdf'); value: pathdir														IN: defined; OUT: unchanged
		\%hoh_tsvdata_of_tsvfilename_of_printcode			, # key1: R grafic unique printcode (e.g. 'metafileT'); key2: *.tsv filename; value: *tsv filename data addings (e.g. clanfilename, type of support...)								IN: defined; OUT: changed
		\%hol_qn_of_qfiltered_otu_of_4clan_of_ft			, # key1: 4seq filtermethod (e.g. nap_risk1); key2: taxon; key3: 4clan|total; value: list of rejected qn's																			IN: defined; OUT: unchanged
		\%hoh_qremaining_of_otu_of_4clan_of_ft				, # key1: 4seq filtermethod (e.g. nap_risk1); key2: taxon; key3: 4clan|total; value: N remaining quartets																			IN: defined; OUT: unchanged
		\%clan_of_seq										, # key1: sequence-name; value: clan definition																																		IN: defined; OUT: unchanged
		\%seen_clan4										, # key1: clan4; value: N analyseed quartets																																		IN: defined; OUT: unchanged
		\%seen_filter_method								, # key1: filter method incl OTU clipping approaches; value: seen +1																												IN: defined; OUT: unchanged
) ;
#####################################################################################


#####################################################################################
# sampling tsv data lines Q6 (N rejected and remaining quartets)                    # 
#####################################################################################
&print_tsv_Q6(
						
		\%subfolder_q										, # hashpath to output-type specific subdirs of single quartet-pattern analyses; key: output-type (e.g. 'pdf'); value: pathdir														IN: defined; OUT: unchanged
		\%hoh_tsvdata_of_tsvfilename_of_printcode			, # key1: R grafic unique printcode (e.g. 'metafileT'); key2: *.tsv filename; value: *tsv filename data addings (e.g. clanfilename, type of support...)								IN: defined; OUT: changed
		\%hoh_N_qrej_of_clan4_of_ft							, # key1: 4seq filtermethod (e.g. nap_risk1); key2: 4clan combination code; value:  N rejected quartets																				IN: defined; OUT: unchanged
		\%seen_clan4										, # key1: clan4; value: N analyseed quartets																																		IN: defined; OUT: unchanged
		\%seen_filter_method								, # key1: filter method incl OTU clipping approaches; value: seen +1																												IN: defined; OUT: unchanged
) ;
#####################################################################################


#####################################################################################
# sampling tex data lines LQ7 (N rejected and remaining quartets)                   # 
#####################################################################################
&latex_data_sample_LQ7(
						
		\%hoh_linedata_of_latexfilename_of_printcode		, # key1: Latex table unique printcode (e.g. 'latexC'); key2: *.tex filename; value: @list of assigned table lines																	IN: defined; OUT: changed
		\%textfilename_of_printcode							, # key1: textfile unique printcode (e.g. 'textA'); value: filename for print																										IN: defined; OUT: unchanged
		\%subfolder_q										, # hashpath to output-type specific subdirs of single quartet-pattern analyses; key: output-type (e.g. 'pdf'); value: pathdir														IN: defined; OUT: unchanged
		\%hoh_N_qrej_of_clan4_of_ft							, # key1: 4seq filtermethod (e.g. nap_risk1); key2: 4clan combination code; value:  N rejected quartets																				IN: defined; OUT: unchanged
		\%seen_clan4										, # key1: clan4; value: N analyseed quartets																																		IN: defined; OUT: unchanged
		\%seen_filter_method								, # key1: filter method incl OTU clipping approaches; value: seen +1																												IN: defined; OUT: unchanged
) ;
######################################################################################################

#exit;

#------------------- OTU ANALYSIS ------------------------------------------------------------------
	
######################################################################################################
## single OTU analysis for each analysed quartet calculation and filter-method (single and combined) #
######################################################################################################
my	@seq_filters = sort keys %hol_sup_of_type_of_mfilter ;

for my $ftype ( @seq_filters ){
	
	###################################################
	# step by step processing of quartet-calc method and assigned filter-methods
	# $ftype = 'calc_filter', e.g. 'nap_risk1'; 
	&evaluation_otu_support_main(
		
		\%subfolder_q										, # hashpath to output-type specific subdirs of single quartet-pattern analyses; key: output-type (e.g. 'pdf'); value: pathdir													IN: defined, OUT: unchanged
		\%hoh_tsvdata_of_tsvfilename_of_printcode			, # key1: R grafic unique printcode (e.g. 'metafileT'); key2: *.tsv filename; value: *tsv filename data addings (e.g. clanfilename, type of support...)							IN: defined, OUT: changed
		\%hoh_sup_of_clan4_of_calc_of_qtree_of_fc			, # multiple hash, storing the final support (nap or napv) related to each rooted quartet-tree																					IN: defined; OUT: unchanged
		\%hoh_seen_rejected_of_4clan_of_filtertype_of_fn	, # key1: 4clan code; key2: type of supportfilter ('nap_risk1', 'napr_risk1', or 'napv_risk1'); key3: actual sealion quartet number; value: filter corresponding score-value	IN: defined; OUT: unchanged 
		\%taxa4_of_fc										, # key1: quartet-number ($fc); value : quartet assigned taxa (alphabetically ordered, joined by ':'																			IN: defined; OUT: unchanged 
		\$ftype												, # type of calculation (nap, napv...) and or filter (rsik_1, risk1_dist1...) method																							IN: defined; OUT: unchanged 
		\%clan_of_seq										, # key1: taxon; value: assigned clan																																			IN: defined; OUT: unchanged 
		\%hol_sup_of_type_of_mfilter						, # key1: type of supportfilter ('nap_risk1', 'napr_risk1', or 'napv_risk1' ...); key2: mean|median; value : list of 'tree\tmean|medianScore'									IN: defined; OUT: changed
		\%algorithm_def_methods								, # key1: 'ftaxon'; value: list of defined taxon filter																															IN: defined; OUT: unchanged
		\%value_of_scr_par									, # key1: sealion user available parameter options; value: optional parameter setting																 							IN: defined; OUT: unchanged
		\%hoh_linedata_of_latexfilename_of_printcode		, # key1: Latex table unique printcode (e.g. 'latexC'); key2: *.tex filename; value: @list of assigned table lines																IN: defined; OUT: changed
		\%seen_remain_of_cmfm_of_fc							, # key1: calculation method + filter method (nap_risk1, napv_dist2, napr_dist3); key2: quartet-number ($fc); value: seen +1													IN: defined; OUT: changed
		\%seen_filter_method								, # key1: filter method incl OTU clipping approaches; value: seen +1																											IN: defined; OUT: changed
		\%textfilename_of_printcode							, # key1: textfile unique printcode (e.g. 'textA'); value: filename for print																									IN: defined; OUT: unchanged
	) ;
	###################################################
}
######################################################################################################



#####################################################################################
# sampling tsv data lines Q7 (N sites of best quartet-tree support)                 # 
#####################################################################################
&print_tsv_Q7(
	
	\%subfolder_q											, # hashpath to output-type specific subdirs of single quartet-pattern analyses; key: output-type (e.g. 'pdf'); value: pathdir														IN: defined; OUT: unchanged
	\%hoh_tsvdata_of_tsvfilename_of_printcode				, # key1: R grafic unique printcode (e.g. 'metafileT'); key2: *.tsv filename; value: *tsv filename data addings (e.g. clanfilename, type of support...)								IN: defined; OUT: changed
	\%seqL_of_fc											, # key1: quartet-number ($fc); value: quartet related sequence length (corrected)																									IN: defined; OUT: unchanged
	\%bestTree_of_fc_of_cm									, # key1: quartet-number ($fc), key2: calculation method (nap, napv, napr); value: quartet related best supported tree																IN: defined; OUT: unchanged
	\%seen_remain_of_cmfm_of_fc								, # key1: calculation method + filter method (nap_risk1, napv_dist2, napr_dist3); key2: quartet-number ($fc); value: seen +1														IN: defined; OUT: unchanged
	\%taxa4_of_fc											, # key1: quartet-number ($fc); value : quartet assigned taxa (alphabetically ordered, joined by ':'																				IN: defined; OUT: unchanged
	\%hol_fc_of_clan4										, # key1: 4clan combination code; value: list of assigned quartet-numbers ($fc)																										IN: defined; OUT: unchanged
	\%algorithm_def_methods									, # key1: methode type; key2: method code; value: integer																															IN: defined; OUT: unchanged
	\%seen_filter_method									, # key1: filter method incl OTU clipping approaches; value: seen +1																												IN: defined; OUT: unchanged
) ;
#####################################################################################

#####################################################################################
# sampling tsv data lines MQ1 (4clan support of each calculation and filter method) # 
#####################################################################################
&print_tsv_MQ1(
	
	\%subfolder_q											, # hashpath to output-type specific subdirs of single quartet-pattern analyses; key: output-type (e.g. 'pdf'); value: pathdir														IN: defined; OUT: unchanged
	\%hoh_tsvdata_of_tsvfilename_of_printcode				, # key1: R grafic unique printcode (e.g. 'metafileT'); key2: *.tsv filename; value: *tsv filename data addings (e.g. clanfilename, type of support...)								IN: defined; OUT: changed
	\%hol_fc_of_clan4										, # key1: 4clan combination code; value: list of assigned quartet-numbers ($fc)																										IN: defined; OUT: unchanged
	\%hol_sup_of_type_of_mfilter							, # key1: type of supportfilter ('nap_risk1', 'napr_risk1', or 'napv_risk1' ...); key2: mean|median; value: list of 'tree\tmean|medianScore'										IN: defined; OUT: unchanged
	\%value_of_scr_par										, # key1: sealion user available parameter options; value: optional parameter setting																								IN: defined; OUT: unchanged
	\%clan4_of_tree											, # key1: 4clan tree; value: 4clan code																																				IN: defined; OUT: unchanged
	\%algorithm_def_methods									, # key1: methode type; key2: method code; value: integer																															IN: defined; OUT: unchanged
) ;
#####################################################################################



#------------------- ICEBREAKER MATRIX PREPARATION ---------------------------------------------------
	
######################################################################################################
# preparation of support matrix of clan relationships (Icebreaker input)                             #
######################################################################################################
	
#####################################################################################
# step by step processing of quartet-calc method and assigned filter-methods
for my $ftype ( sort keys %hol_sup_of_type_of_mfilter ){
	
	print 	"\n\t----------------------------------\n"	,
			"\tIcebreaker multi-clan analysis of ", $ftype, "\n" ;
	
	for my $stype ( sort keys %{$hol_sup_of_type_of_mfilter{$ftype}} ){
		
		####################################
		# only conducted if all trees of each 4clan combination are present
		# in some cases, filtering might lead to a complete rejection of 4clan support values
		my $N_4clan_combinations = keys %seen_clan4 ;
		my $N_trees_total		 = @{$hol_sup_of_type_of_mfilter{$ftype}{$stype}} ;
		
		if ( $N_trees_total == $N_4clan_combinations * 3 ){
						
			if ( $value_of_scr_par{sys_parameter}{$stype} == 1 ){
			
				######################################################################
				# Building multiclan support matrix (icebreaker-input)               #
				######################################################################
				
				####################################
				my %hoh_score_of_rooted_clan_of_sister_pair ;
				# key 1: rooted clan (y-axis e.g. RO)
				# key 2: sister pair of rooted clan (x-axis e.g. AB)
				# value: phyquart identified support score
				
				my %seen_sister_pair ;
				# key : sister_pair combination (e.g. 'AB')
				# value: counter of appearance in multiple clan analyses
				
				
				my	$sf_type = $stype."_".$ftype ;
				####################################
				
				####################################
				# build icebreaker infile-mcmatrix
				&build_multiclan_support_matrix(
					
					\%hol_sup_of_type_of_mfilter				, # key1: type of supportfilter ('nap_risk1', 'napr_risk1', or 'napv_risk1' ...); key2: mean|median; value: list of 'tree\tmean|medianScore'	IN: defined; OUT: unchanged
					\%hoh_score_of_rooted_clan_of_sister_pair	, # key1: rooted clan (y-axis e.g. RO), key 2: sister pair of rooted clan (x-axis e.g. AB), value: phyquart identified support score			IN: undefined/defined; OUT: defined/changed
					\%seen_sister_pair							, # key1: sister_pair combination (e.g. 'AB', value: counter of appearance in multiple clan analyses											IN: undefined/defined; OUT: defined/changed
					\%value_of_scr_par							, # key1: sealion user available parameter options; value: optional parameter setting															IN: defined; OUT: unchanged
					\$stype										, # mean|median																																	IN: defined; OUT: unchanged
					\$ftype										, # supportfilter ('nap_risk1', 'napr_risk1', or 'napv_risk1' ...)																				IN: defined; OUT: unchanged
				) ;
				####################################
				
				####################################
				# write and execute icebreaker infile-mcmatrix
				my	$outfile_prefix ;
				
				&print_and_processing_multiclan_support_matrix(
					
					\%hoh_score_of_rooted_clan_of_sister_pair	, # key 1: rooted clan (y-axis e.g. RO), key 2: sister pair of rooted clan (x-axis e.g. AB), value: phyquart identified support score		IN: defined; OUT: unchanged
					\%seen_sister_pair							, # key : sister_pair combination (e.g. 'AB', value: counter of appearance in multiple clan analyses										IN: defined; OUT: unchanged
					\%subfolder_q								, # hashpath to output-type specific subdirs of single quartet-pattern analyses; key: output-type (e.g. 'pdf'); value: pathdir				IN: defined; OUT: unchanged
					\$sf_type									, # code of quartet-calculation and filter-method	(e.g. 'nap_risk1')																		IN: defined; OUT: unchanged
					\$clanfile									, # name of given clanfile																													IN: defined; OUT: unchanged
					\$outgroup_clan								, # outgroup clan-code																														IN: defined; OUT: unchanged
					\%value_of_scr_par							, # key1: sealion user available parameter options; value: optional parameter setting														IN: defined; OUT: unchanged
					\$outfile_prefix							, # filename-prefix of tree-fishing resultfiles																								IN: undefined; OUT: defined
				);
				######################################################################
				
				
				######################################################################
				# Icebreaker outfile processing                                      #
				######################################################################
				&treefishing_result_processing(
					
					\%value_of_scr_par								, # key1: sealion user available parameter options; value: optional parameter setting																				IN: defined; OUT: unchanged
					\%subfolder_q									, # hashpath to output-type specific subdirs of single quartet-pattern analyses; key: output-type (e.g. 'pdf'); value: pathdir										IN: defined; OUT: unchanged
					\%hoh_tsvdata_of_tsvfilename_of_printcode		, # key1: R grafic unique printcode (e.g. 'metafileT'); key2: *.tsv filename; value: *tsv filename data addings (e.g. clanfilename, type of support...)				IN: defined; OUT: changed
					\$sf_type										, # code of quartet-calculation and filter-method	(e.g. 'nap_risk1')																								IN: defined; OUT: unchanged
					\$clanfile										, # name of given clanfile																																			IN: defined; OUT: unchanged
					\$outfile_prefix								, # filename-prefix of tree-fishing resultfiles																														IN: defined; OUT: unchanged
					\%hoh_linedata_of_latexfilename_of_printcode	, # key1: Latex table unique printcode (e.g. 'latexC'); key2: *.tex filename; value: @list of assigned table lines													IN: defined; OUT: changed
					\%textfilename_of_printcode						, # key1: textfile unique printcode (e.g. 'textA'); value: filename for print																						IN: defined; OUT: unchanged
					\%algorithm_def_methods							, # key1: methode type; key2: method code; value: integer																											IN: defined; OUT: unchanged
				) ;
				######################################################################
			}
			####################################
		}
	}
}
######################################################################################################



#------------------- Latex ---------------------------------------------------------------------------
	
######################################################################################################
# Latex table prints                                                                                 #
######################################################################################################
if ( $value_of_scr_par{sys_parameter}{l_prt} == 1 ){
	
	&latex_quartet_analyses(
		
		\%print_latex									, # Print control of single latex tables																						IN: defined; OUT: unchanged
		\%hoh_linedata_of_latexfilename_of_printcode	, # key1: Latex table unique printcode (e.g. 'latexC'); key2: *.tex filename; value: @list of assigned table lines				IN: defined; OUT: unchanged
		\%hol_pdf_files_of_type							, # key: result type; value: list of associated pdf's in order of construction date												IN: defined; OUT: changed
		\@pdf_files										, # list variable of single quartet inference denerated pdf's																	IN: defined; OUT: changed
		\%subfolder_q									, # key1: file type ('svg', 'pdf', 'tex', 'tsv') ; value: printdir																IN: defined; OUT: unchanged
		\%subsubfolder_q								, # key1: pdf related subfolder codes ('pdf_single', 'pdf_sumed'...)															IN: defined; OUT: unchanged	
		\%value_of_scr_par								, # key1: sealion user available parameter options; value: optional parameter setting											IN: defined; OUT: unchanged
		\%algorithm_def_methods							, # list of defined support methods (nap, napv ...)																				IN: defined; OUT: unchanged 
	)
}
######################################################################################################



#------------------- R PLOTTING ----------------------------------------------------------------------
	
######################################################################################################
# R plot processings                                                                                 #
######################################################################################################
if ( $value_of_scr_par{sys_parameter}{r_prt} == 1 ){
	
	&plotsR_quartet_analyses(
		
		\%print_Rplot_quartets						, # Print control of single R plots due to single split-pattern quartet analyses																			IN: defined; OUT: unchanged
		\%hoh_tsvdata_of_tsvfilename_of_printcode	, # key1: R grafic unique printcode (e.g. 'metafileT'); key2: *.tsv filename; value: *tsv filename data addings (e.g. clanfilename, type of support...)  	IN: defined; OUT: unchanged
		\%hol_pdf_files_of_type						, # key: result type; value: list of associated pdf's in order of construction date																			IN: defined; OUT: changed
		\@pdf_files									, # list variable of single quartet inference denerated pdf's																								IN: defined; OUT: changed
		\%subfolder_q								, # key1: file type ('svg', 'pdf', 'tex', 'tsv') ; value: printdir																							IN: defined; OUT: unchanged
		\%subsubfolder_q							, # key1: pdf related subfolder codes ('pdf_single', 'pdf_sumed'...)																						IN: defined; OUT: unchanged	
		\%value_of_scr_par							, # key1: sealion user available parameter options; value: optional parameter setting																 		IN: defined; OUT: unchanged
		\%algorithm_def_methods						, # list of defined support methods (nap, napv ...)																											IN: defined; OUT: unchanged 
	)
}
######################################################################################################
#exit;


#------------------- FINAL PDF SUMMARIZATION ---------------------------------------------------------

######################################################################################################
# pdf summarization with pdfunite                                                                    #
######################################################################################################
if ( ( $value_of_scr_par{sys_parameter}{r_prt} == 1 )	||
	 ( $value_of_scr_par{sys_parameter}{l_prt} == 1 )	){
	
	print	"\n\n\t---------------------------------------------------\n",
			"\tSummarising PDF's...\n" ;
	
	for my $method ( sort keys %hol_pdf_files_of_type ){
		
		my	@prts			=	split "_", $method ;
		my	$cm				=	$prts[1] ;
		
		my	$qfilter_type	=	&recode_filter_printout( \$method ) ;
		my	$new_pdf 		=	$subsubfolder_q{pdf_sumed}."/SeaLion_results_sumed_".$prts[0]."_".$qfilter_type."_".$prts[1].".pdf" ;
		
		my	$command_line	=	join " ", @{$hol_pdf_files_of_type{$method}}	;
			$command_line	=	"pdfunite ".$command_line." $new_pdf"			;	# extend commandline with converter scriptname ('pdfunite'), followed by all (blank separated) infiles, followed by outut filename
				
		print	"\n\t\t...", $new_pdf ;
		
		system($command_line)
	}
}
######################################################################################################



#------------------- END SCREEN ----------------------------------------------------------------------

######################################################################################################
# close script by printing end-screen                                                                #
######################################################################################################
&print_endscreen( \%value_of_scr_par ) ;

exit ;
######################################################################################################



#------------------- SOUBROUTINES --------------------------------------------------------------------

######################################################################################################
## SUBROUTINES                                                                                       #
######################################################################################################
sub argv_handling{

	my	$href_val_of_scr_par	= $_[0] ;	# key1: sealion user available parameter options; value: optional parameter setting										IN: defined; OUT: changed
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# &assign_menu_parameter
	# &check_user_parameter
	# &config_menu
	######################################################################
	
	######################################################################
	# Read, check, assign, and print user specified command line options #
	######################################################################
	
	##############################################
	# assign parameter to %$href_val_of_scr_par
	&assign_menu_parameter(
		
			\%$href_val_of_scr_par	,	# key1: parameter option; value: parameter setting			IN: defined; OUT: changed
	) ;
	##############################################

	##############################################
	# Check command line defined parameter
	my ( $commandline ) = join " ", @ARGV ; @ARGV = () ;
	#print "\ncmd line ", $commandline, "\n\n"; exit;

	if ( $commandline ){

		&check_user_parameter(
			
			\%$href_val_of_scr_par	,	# key1: sealion user available parameter options; value: optional parameter setting			IN: defined; OUT: changed
			\$commandline			,	# command line of one or multiple parameters												IN: defined; OUT: changed
		) ;
	}
	##############################################
	
	##############################################
	# assign parameter to %$href_val_of_scr_par
	&assign_menu_parameter(
		
			\%$href_val_of_scr_par	,	# key1: parameter option; value: parameter setting			IN: defined; OUT: changed
	) ;
	##############################################
	
	######################################################################
	
	
	######################################################################
	# Unless -s start command is defined. open config menu               #
	######################################################################
	
	##############################################
	# unless defined start command -s... 
	# ...open config-menu & parameter-request, 
	# ...new parameter check & assignment
	until( $value_of_scr_par{sys_parameter}{sysstart} ){
		
		my $config_line_new ;
		
		&config_menu(
	
			\%$href_val_of_scr_par	,	# key1: parameter option; value: parameter setting			IN: defined; OUT: unchanged
			\$config_line_new		,	# via config menu specified parameter line					IN: undefined; OUT: defined
		);
		
		&check_user_parameter(
			
			\%$href_val_of_scr_par	,	# key1: options; value: optional parameter setting			IN: defined; OUT: changed
			\$config_line_new		,	# command line of one or multiple parameters				IN: defined; OUT: changed
		) ;
		
		&assign_menu_parameter(
		
			\%$href_val_of_scr_par	,	# key1: parameter option; value: parameter setting			IN: defined; OUT: changed
		) ;
	}
	##############################################
	
	##############################################
	# unless restart is activated when msa and clanfile is not command line defined, extra check for script internally defined inputs
	unless ( $href_val_of_scr_par->{sys_parameter}{restart} == 1 ){ 
		
		unless ( ( $href_val_of_scr_par->{sys_parameter}{msafile}  ) && 
				 ( $href_val_of_scr_par->{sys_parameter}{clfile}   ) ){ 
			
			print "\nPARAMETER-ERROR: undefined alignmentfile and/or clanfile!\n\n" ;
			
			&help( \'h', \%$href_val_of_scr_par ) 
		}
	}
	##############################################
	
	##############################################
	# if restart is activated, extra check for script internally defined input spd-file path
	if ( $href_val_of_scr_par->{sys_parameter}{restart} == 1 ){ 
		
		unless ( $href_val_of_scr_par->{sys_parameter}{i_main} ){ 
			
			print "\nPARAMETER-ERROR: no inputpath defined for -restart option!\n\n" ;
			
			&help( \'imain', \%$href_val_of_scr_par ) 
		}
	}
	##############################################
	
	##############################################
	# if system start '-s' defined goto 'START' (goto command avoids possible back jumps to the help menu)
	if ( $value_of_scr_par{sys_parameter}{sysstart} ){ goto START }
	else{ die "\nBUG-ERROR: Cannot assign 'sysstart' in subroutine &argv_handling!\nPlease, report BUG to software developer!\n\n" }
	##############################################
	
	######################################################################
}

sub assign_filter_scales{
	
	my	$href_value_of_scr_par 			= $_[0] ; # key1: method_code; key2: method_option; value: float or integer										IN: defined; OUT: unchanged
	my	$href_hol_scales_of_filtertype	= $_[1] ; # key1: type of supportfilter ('risk', 'qdist'); value: list of ftype corresponding threshold scales	IN: undefined; OUT: defined
	
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# &calc_scales
	######################################################################
	
	######################################################################
	## generation of possible sclae values between defined limit-range  ##
	######################################################################
	
	#########################################
	# if one RISK-filter is defined at least, calculate and assign list of risk scales using &calc_scales
	if ( ( $value_of_scr_par{seq_filterRisk}{risk1} == 1 ) ||
		 ( $value_of_scr_par{seq_filterRisk}{risk2} == 1 ) ||
		 ( $value_of_scr_par{seq_filterRisk}{risk3} == 1 ) ||
		 ( $value_of_scr_par{seq_filterRisk}{risk4} == 1 ) ){
		
		@{$href_hol_scales_of_filtertype->{risk}} = &calc_scales(
			
			\$href_value_of_scr_par->{seq_filterRisk}{turisk},	# upper risk-scale threshold 	IN: defined; OUT: unchanged
			\$href_value_of_scr_par->{seq_filterRisk}{tlrisk},	# lower risk-scale threshold 	IN: defined; OUT: unchanged
			\$href_value_of_scr_par->{seq_filterRisk}{rscale},	# scale steps				 	IN: defined; OUT: unchanged
			\'risk'												# filtertype (risk or qdist)	IN: defined; OUT: unchanged
		) ;
	}
	#########################################
	
	#########################################
	# if one DIST-filter is defined at least, calculate and assign list of dist scales using &calc_scales
	if ( $value_of_scr_par{seq_filterDist}{qdist1} == 1 ){
		
		@{$href_hol_scales_of_filtertype->{qdist}} = &calc_scales(
			
			\$href_value_of_scr_par->{seq_filterDist}{tudist},	# upper risk-scale threshold 	IN: defined; OUT: unchanged
			\$href_value_of_scr_par->{seq_filterDist}{tldist},	# lower risk-scale threshold 	IN: defined; OUT: unchanged
			\$href_value_of_scr_par->{seq_filterDist}{dscale},	# scale steps				 	IN: defined; OUT: unchanged
			\'qdist'											# filtertype (risk or qdist)	IN: defined; OUT: unchanged
		) ;
	}
	#########################################
	
	######################################################################
	#exit;
}

sub assign_menu_parameter{
	
	my	$href_value_of_scr_par		= $_[0] ; # key1: available options; value: parameter setting		IN: defined; OUT: unchanged
	
	######################################################################
	# called from subroutine(s)...
	# &argv_handling
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	############################################################
	# alignmenfile
	$href_value_of_scr_par->{menu_parameter}{msa_infile} = 'none' ;
	if 		( $href_value_of_scr_par->{sys_parameter}{msa_nopath} 	){ $href_value_of_scr_par->{menu_parameter}{msa_infile} = $href_value_of_scr_par->{sys_parameter}{msa_nopath} }
	elsif 	( $href_value_of_scr_par->{sys_parameter}{restart} 		){ $href_value_of_scr_par->{menu_parameter}{msa_infile} = 'restart' }
	############################################################
	
	
	############################################################
	# outgroup
	$href_value_of_scr_par->{menu_parameter}{outgroup} = 'none' ;
	if 	( $href_value_of_scr_par->{sys_parameter}{outgroup} ){ $href_value_of_scr_par->{menu_parameter}{outgroup} = $href_value_of_scr_par->{sys_parameter}{outgroup} }
	############################################################
	
	
	############################################################
	# Quartet-Filter 'RISK' setup
	$href_value_of_scr_par->{menu_parameter}{risk_default} = 'disabled';
	my	@risk_setting ;
	
	$href_value_of_scr_par->{menu_parameter}{risk1} = 'no' ;
	if 	( $href_value_of_scr_par->{seq_filterRisk}{risk1} ){ $href_value_of_scr_par->{menu_parameter}{risk_default2} = 'activated'; $href_value_of_scr_par->{menu_parameter}{risk1} = 'yes'; push @risk_setting, 'risk1' }
	
	$href_value_of_scr_par->{menu_parameter}{risk2} = 'no' ;
	if 	( $href_value_of_scr_par->{seq_filterRisk}{risk2} ){ $href_value_of_scr_par->{menu_parameter}{risk_default2} = 'activated'; $href_value_of_scr_par->{menu_parameter}{risk2} = 'yes'; push @risk_setting, 'risk2' }
	
	if ( @risk_setting ){ $href_value_of_scr_par->{menu_parameter}{risk_default} = join ",", @risk_setting }
	############################################################
	
	
	############################################################
	# Quartet-Filter 'DIST' setup
	$href_value_of_scr_par->{menu_parameter}{qdist_default}	= 'disabled';
	$href_value_of_scr_par->{menu_parameter}{dist1}			= 'no' ;
	
	if 	( $href_value_of_scr_par->{seq_filterDist}{qdist1} ){ 
		
		$href_value_of_scr_par->{menu_parameter}{dist1}			= 'yes'; 
		$href_value_of_scr_par->{menu_parameter}{qdist_default} = 'activated' 
	}
	############################################################
	
	
	############################################################
	# OTU-Filter selected 4clan classification(s)
	$href_value_of_scr_par->{menu_parameter}{tdist_default}		 = 'disabled';
	$href_value_of_scr_par->{menu_parameter}{otu_strict_setting} = 'disabled';
	my @otu_str_setting 	;
	
	$href_value_of_scr_par->{menu_parameter}{strict_out} = 'no' ;
	if 	( $href_value_of_scr_par->{tax_filterDist}{strict_out} ){ 
		
		$href_value_of_scr_par->{menu_parameter}{strict_out} = 'yes'; 
		push @otu_str_setting, 'outgroup'	; 
		$href_value_of_scr_par->{menu_parameter}{tdist_default} = 'activated' 
	}
	
	$href_value_of_scr_par->{menu_parameter}{strict_in} = 'no' ;
	if 	( $href_value_of_scr_par->{tax_filterDist}{strict_in}  ){ 
		
		$href_value_of_scr_par->{menu_parameter}{strict_in}  = 'yes'; 
		push @otu_str_setting, 'ingroup'	; 
		$href_value_of_scr_par->{menu_parameter}{tdist_default} = 'activated' 
	}
	
	$href_value_of_scr_par->{menu_parameter}{strict_all} = 'no' ;
	if 	( $href_value_of_scr_par->{tax_filterDist}{strict_all} ){ 
		
		$href_value_of_scr_par->{menu_parameter}{strict_all} = 'yes'; 
		push @otu_str_setting, 'all' ;
		$href_value_of_scr_par->{menu_parameter}{tdist_default} = 'activated' 
	}
	
	if ( @otu_str_setting ){ $href_value_of_scr_par->{menu_parameter}{otu_strict_setting} = join ",", @otu_str_setting }
	############################################################
	
	
	############################################################
	# Restart option
	$href_value_of_scr_par->{menu_parameter}{restart_setting} = 'disabled';
	if ( $href_value_of_scr_par->{sys_parameter}{restart} == 1 ){ $href_value_of_scr_par->{menu_parameter}{restart_setting} = 'activated' }
	############################################################
	
	
	############################################################
	# SPD infile path
	$href_value_of_scr_par->{menu_parameter}{spd_folder} = 'none' ;
	if 	( $href_value_of_scr_par->{sys_parameter}{i_main} ){ $href_value_of_scr_par->{menu_parameter}{spd_folder} = $href_value_of_scr_par->{sys_parameter}{i_main} }
	############################################################
	
	
	############################################################
	# Additional Fileprint(s)
	$href_value_of_scr_par->{menu_parameter}{prt_setting} = 'none';
	my @print_setting ;
	
	$href_value_of_scr_par->{menu_parameter}{prt_latex} = 'no' ;
	if 	( $href_value_of_scr_par->{sys_parameter}{l_prt} ){ 
		
		$href_value_of_scr_par->{menu_parameter}{prt_latex} = 'yes'; 
		
		push @print_setting, 'LatexTable' 
	}
	
	$href_value_of_scr_par->{menu_parameter}{prt_Rplot} = 'no' ;
	if 	( $href_value_of_scr_par->{sys_parameter}{r_prt} ){ 
		
		$href_value_of_scr_par->{menu_parameter}{prt_Rplot} = 'yes'; 
		
		push @print_setting, 'Rplot' 
	}
	
	if ( @print_setting ){ $href_value_of_scr_par->{menu_parameter}{prt_setting} = join ",", @print_setting }
	############################################################
	
	
	############################################################
	# RY coding
	$href_value_of_scr_par->{menu_parameter}{ry_setting} = 'disabled';
	if ( $href_value_of_scr_par->{sys_parameter}{rycoding} == 1 ){ $href_value_of_scr_par->{menu_parameter}{ry_setting} = 'activated' }
	############################################################
	
	
	############################################################
	# Script query
	$href_value_of_scr_par->{menu_parameter}{query_setting} = 'disabled';
	if ( $href_value_of_scr_par->{sys_parameter}{request} == 1 ){ $href_value_of_scr_par->{menu_parameter}{query_setting} = 'activated' }
	############################################################
	
	
	############################################################
	# Assign defined support averages
	$href_value_of_scr_par->{menu_parameter}{average} = ();
	if		( ( $href_value_of_scr_par->{sys_parameter}{median} == 1 ) &&
			  ( $href_value_of_scr_par->{sys_parameter}{mean}   == 1 ) ){ $href_value_of_scr_par->{menu_parameter}{average} = 'mean+median' }
	elsif 	( ( $href_value_of_scr_par->{sys_parameter}{median} == 1 ) &&
			  ( $href_value_of_scr_par->{sys_parameter}{mean}   == 0 ) ){ $href_value_of_scr_par->{menu_parameter}{average} = 'median' 		}	
	else	 															{ $href_value_of_scr_par->{menu_parameter}{average} = 'mean' 		}
	############################################################
	
	
	############################################################
	# Single parameter options
	$href_value_of_scr_par->{menu_parameter}{model_aa} 		= $href_value_of_scr_par->{p4_parameter}{model_aa} 		;
	$href_value_of_scr_par->{menu_parameter}{model_nu} 		= $href_value_of_scr_par->{p4_parameter}{model_nu} 		;
	$href_value_of_scr_par->{menu_parameter}{start_alpha}  	= $href_value_of_scr_par->{p4_parameter}{alpha} 		;
	$href_value_of_scr_par->{menu_parameter}{start_pinv}   	= $href_value_of_scr_par->{p4_parameter}{pinv} 			;
	$href_value_of_scr_par->{menu_parameter}{output_folder} = $href_value_of_scr_par->{sys_parameter}{o_main} 		;
	$href_value_of_scr_par->{menu_parameter}{clfile}		= $href_value_of_scr_par->{sys_parameter}{clf_nopath}	;
	$href_value_of_scr_par->{menu_parameter}{maxquartets}	= $href_value_of_scr_par->{qlimit}{qnumber}				;
	$href_value_of_scr_par->{menu_parameter}{minlengths}	= $href_value_of_scr_par->{qlimit}{qlength}				;
	$href_value_of_scr_par->{menu_parameter}{upperlimit}	= $href_value_of_scr_par->{seq_filterRisk}{turisk}		;
	$href_value_of_scr_par->{menu_parameter}{lowerlimit}	= $href_value_of_scr_par->{seq_filterRisk}{tlrisk}		;
	$href_value_of_scr_par->{menu_parameter}{rscale}		= $href_value_of_scr_par->{seq_filterRisk}{rscale}		;
	$href_value_of_scr_par->{menu_parameter}{ralpha}		= $href_value_of_scr_par->{seq_filterRisk}{ralpha}		;
	$href_value_of_scr_par->{menu_parameter}{upperlimit_d}	= $href_value_of_scr_par->{seq_filterDist}{tudist}		;
	$href_value_of_scr_par->{menu_parameter}{lowerlimit_d}	= $href_value_of_scr_par->{seq_filterDist}{tldist}		;
	$href_value_of_scr_par->{menu_parameter}{dscale}		= $href_value_of_scr_par->{seq_filterDist}{dscale}		;
	$href_value_of_scr_par->{menu_parameter}{dalpha}		= $href_value_of_scr_par->{seq_filterDist}{dalpha}		;
	$href_value_of_scr_par->{menu_parameter}{tdistlimit}	= $href_value_of_scr_par->{tax_filterDist}{tdist1st2nd} ;
	############################################################
}

sub assign_selected_methods{
	
	my	$href_value_of_scr_par		= $_[0] ; # key1: method_code; key2: method_option; value: float or integer		IN: defined; OUT: unchanged
	my	$href_algorithm_def_methods	= $_[1] ; # key1: methode type; key2: method code; value: integer				IN: undefined; OUT: defined
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	## assign true (==1) parameters of org. & filter methods to new hash #
	######################################################################
	#print "::&assign_selected_methods\n" ;
	
	for my $sm ( keys %{$href_value_of_scr_par->{seq_algorithm}} ){ 
		
		##### check main (original) method definitions, define method in %$href_algorithm_def_methods if true (1)....
		if ( $href_value_of_scr_par->{seq_algorithm}{$sm} == 1 ){ 
			
			$href_algorithm_def_methods->{original}{$sm} = 1	;	
		
			##### check filter risk methods and add them to expanded method $c of original method $sm if true (1)
			for my $fm_risk ( keys %{$href_value_of_scr_par->{seq_filterRisk}} ){ 
				
				if ( $href_value_of_scr_par->{seq_filterRisk}{$fm_risk} == 1 ){
					
					my $c = $sm."_".$fm_risk ;
					$href_algorithm_def_methods->{expanded}{$c} = 1
				}
			}
			
			##### check filter dist methods and add them to expanded method $c of original method $sm if true (1)
			for my $fm_dist ( keys %{$href_value_of_scr_par->{seq_filterDist}} ){ 
			
				if ( $href_value_of_scr_par->{seq_filterDist}{$fm_dist} == 1 ){
					
					my $c = $sm."_".$fm_dist ;
					$href_algorithm_def_methods->{expanded}{$c} = 1
				}
			}
		}
	}
	
	if ( $href_value_of_scr_par->{sys_parameter}{mean}			== 1 ){ push @{$href_algorithm_def_methods->{qsuppo}}	, 'mean' 		 }
	if ( $href_value_of_scr_par->{sys_parameter}{median}		== 1 ){ push @{$href_algorithm_def_methods->{qsuppo}}	, 'median' 		 }
	######################################################################
}

sub assign_split_pattern_to_code{

	my	$aref_ptaxa					=	$_[0]	; # array of quartet sequences																IN: defined; OUT: unchanged
	my	$href_clan_of_seq			=	$_[1]	; # hash with sequence name (key) assigned clan (value)										IN: defined; OUT: unchanged
	
	######################################################################
	## assigning split pattern to corresponding code                    ##
	######################################################################
	print "::&assign_split_pattern_to_code\n";
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	#########################################
	# assign split pattern signal to corresponding codes
	my	%splitcode	=	(
		
		'A'	=>	$href_clan_of_seq->{$aref_ptaxa->[0]}.",".$href_clan_of_seq->{$aref_ptaxa->[1]}."||".$href_clan_of_seq->{$aref_ptaxa->[2]}.",".$href_clan_of_seq->{$aref_ptaxa->[3]} ,
		'B'	=>	$href_clan_of_seq->{$aref_ptaxa->[0]}.",".$href_clan_of_seq->{$aref_ptaxa->[2]}."||".$href_clan_of_seq->{$aref_ptaxa->[1]}.",".$href_clan_of_seq->{$aref_ptaxa->[3]} ,
		'C'	=>	$href_clan_of_seq->{$aref_ptaxa->[0]}.",".$href_clan_of_seq->{$aref_ptaxa->[3]}."||".$href_clan_of_seq->{$aref_ptaxa->[1]}.",".$href_clan_of_seq->{$aref_ptaxa->[2]} ,
		'D'	=>	$href_clan_of_seq->{$aref_ptaxa->[0]}.",".$href_clan_of_seq->{$aref_ptaxa->[1]}."|-".$href_clan_of_seq->{$aref_ptaxa->[2]}.",".$href_clan_of_seq->{$aref_ptaxa->[3]} ,
		'E'	=>	$href_clan_of_seq->{$aref_ptaxa->[2]}.",".$href_clan_of_seq->{$aref_ptaxa->[3]}."|-".$href_clan_of_seq->{$aref_ptaxa->[0]}.",".$href_clan_of_seq->{$aref_ptaxa->[1]} ,
		'F'	=>	$href_clan_of_seq->{$aref_ptaxa->[1]}.",".$href_clan_of_seq->{$aref_ptaxa->[3]}."|-".$href_clan_of_seq->{$aref_ptaxa->[0]}.",".$href_clan_of_seq->{$aref_ptaxa->[2]} ,
		'G'	=>	$href_clan_of_seq->{$aref_ptaxa->[0]}.",".$href_clan_of_seq->{$aref_ptaxa->[2]}."|-".$href_clan_of_seq->{$aref_ptaxa->[1]}.",".$href_clan_of_seq->{$aref_ptaxa->[3]} ,
		'H'	=>	$href_clan_of_seq->{$aref_ptaxa->[0]}.",".$href_clan_of_seq->{$aref_ptaxa->[3]}."|-".$href_clan_of_seq->{$aref_ptaxa->[1]}.",".$href_clan_of_seq->{$aref_ptaxa->[2]} ,
		'J'	=>	$href_clan_of_seq->{$aref_ptaxa->[1]}.",".$href_clan_of_seq->{$aref_ptaxa->[2]}."|-".$href_clan_of_seq->{$aref_ptaxa->[0]}.",".$href_clan_of_seq->{$aref_ptaxa->[3]} ,
		'K'	=>	$href_clan_of_seq->{$aref_ptaxa->[0]}."|".$href_clan_of_seq->{$aref_ptaxa->[1]}.",".$href_clan_of_seq->{$aref_ptaxa->[2]}.",".$href_clan_of_seq->{$aref_ptaxa->[3]} ,
		'L'	=>	$href_clan_of_seq->{$aref_ptaxa->[1]}."|".$href_clan_of_seq->{$aref_ptaxa->[0]}.",".$href_clan_of_seq->{$aref_ptaxa->[2]}.",".$href_clan_of_seq->{$aref_ptaxa->[3]} ,
		'M'	=>	$href_clan_of_seq->{$aref_ptaxa->[2]}."|".$href_clan_of_seq->{$aref_ptaxa->[0]}.",".$href_clan_of_seq->{$aref_ptaxa->[1]}.",".$href_clan_of_seq->{$aref_ptaxa->[3]} ,
		'N'	=>	$href_clan_of_seq->{$aref_ptaxa->[3]}."|".$href_clan_of_seq->{$aref_ptaxa->[0]}.",".$href_clan_of_seq->{$aref_ptaxa->[1]}.",".$href_clan_of_seq->{$aref_ptaxa->[2]} ,
		'V'	=>	$href_clan_of_seq->{$aref_ptaxa->[0]}."|".$href_clan_of_seq->{$aref_ptaxa->[1]}."|".$href_clan_of_seq->{$aref_ptaxa->[2]}."|".$href_clan_of_seq->{$aref_ptaxa->[3]} ,
		'I'	=>	$href_clan_of_seq->{$aref_ptaxa->[0]}.",".$href_clan_of_seq->{$aref_ptaxa->[1]}.",".$href_clan_of_seq->{$aref_ptaxa->[2]}.",".$href_clan_of_seq->{$aref_ptaxa->[3]} ,
	) ;
	#########################################
	
	#########################################
	# print terminal info
	for my $t	( 			@$aref_ptaxa 	){ print "\n\t", $href_clan_of_seq->{$t}, "\t-> ", $t 	} print "\n";
	for my $pc	( sort keys %splitcode 		){ print "\n\tcode: ", $pc, " split: ", $splitcode{$pc} } print "\n";
	#########################################
	
	#########################################
	# retutn code assigment
	#return  %splitcode
	#########################################
	
	######################################################################
} 

sub assign_tree_order{
	
	my	$href_fsupport_of_type_of_qt			= $_[0] ;	# key1: support type (e.g. 'nap', 'napv'...); key2: rooted quartet tree; value: score									IN: defined; OUT: unchanged 
	my	$sref_method							= $_[1] ;	# actual support method																									IN: defined; OUT: unchanged 
	my	$tree									= $_[2] ;	# actual quartet tree																									IN: defined; OUT: unchanged 
	
	#print "::&assign_tree_order\n";
	
	######################################################################
	# called from subroutine(s)...
	# &print_tsv_Q3
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	## assign tree support ranking of given method and tree from a hash_##
	######################################################################
	my	$code1st	=	$$sref_method."_1stbestTree" ;
	my	$code2nd	=	$$sref_method."_2ndbestTree" ;
	my	$code3rd	=	$$sref_method."_3rdbestTree" ;
	
	if 		( $tree eq $href_fsupport_of_type_of_qt->{$code1st} ){ return "best"	}
	elsif 	( $tree eq $href_fsupport_of_type_of_qt->{$code2nd} ){ return "2nd" 	}
	elsif 	( $tree eq $href_fsupport_of_type_of_qt->{$code3rd} ){ return "3rd" 	}
	else	{ die "\nBUG-ERROR: Cannot assign quartet-tree ", $tree, " in subroutine &assign_tree_order!\n\tPlease, report BUG to software developer!\n\n" }
	######################################################################
}

sub build_multiclan_support_matrix{

	my	$href_hol_sup_of_type_of_mfilter	= $_[0]	; # key1: type of supportfilter ('nap_risk1', 'napr_risk1', or 'napv_risk1' ...); key2: mean|median; value: list of 'tree\tmean|medianScore'	IN: defined; OUT: unchanged
	my	$href_hoh_score_rooted_sister_pair	= $_[1]	; # key 1: rooted clan (y-axis e.g. RO), key 2: sister pair of rooted clan (x-axis e.g. AB), value: phyquart identified support score			IN: undefined; OUT: defined/changed
	my	$href_seen_sister_pair				= $_[2]	; # key : sister_pair combination (e.g. 'AB', value: counter of appearance in multiple clan analyses											IN: undefined; OUT: defined/changed
	my	$href_val_of_scr_par				= $_[3] ; # key1: sealion user available parameter options; value: optional parameter setting															IN: defined; OUT: unchanged
	my	$sref_stype							= $_[4] ; # mean|median																																	IN: defined; OUT: unchanged
	my	$sref_ftype							= $_[5] ; # supportfilter ('nap_risk1', 'napr_risk1', or 'napv_risk1' ...)																				IN: defined; OUT: unchanged
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# ...
	######################################################################
	
	######################################################################
	## Build multi-clan support matrix                                  ##
	######################################################################
	
	###############
	# sampling of filter and support type corresponding nwk trees and support (each joined by a tab sign)
	my	@topology_score	=	@{$href_hol_sup_of_type_of_mfilter->{$$sref_ftype}{$$sref_stype}} ;
	###############
	
	###############
	# terminal info print
	if ( $href_val_of_scr_par->{sys_parameter}{screen_prt} == 1 ){
		
		print "\n\t", $$sref_stype ," tree support of ", $$sref_ftype, ":" ;
		for my $str	( @topology_score ){ print "\n\t", $str }
		print "\n"
	}
	###############
	
	###############
	# sum up single 4clan tree support
	for my $top_and_score (	@topology_score ){
		
		# split 'topology\tscore'-string, code left clan group as rooted clans, right clan group as sister pair
		my	@top_and_score		=	split	"\t"	, $top_and_score	;
			$top_and_score[0]	=~	s/\(|\)|;//g						;
		my	@clan_taxa			=	split	","		, $top_and_score[0]	;
		my	$root_clan_1		=	'R|'.$clan_taxa[0]					;
		my	$root_clan_2		=	'R|'.$clan_taxa[1]					;
		my	$sister_pair		=	$clan_taxa[2]."|".$clan_taxa[3]		;
		
		$href_hoh_score_rooted_sister_pair->{$root_clan_1}{$sister_pair}	+=	$top_and_score[1] ;
		$href_hoh_score_rooted_sister_pair->{$root_clan_2}{$sister_pair}	+=	$top_and_score[1] ;
		$href_seen_sister_pair->{$sister_pair}++ ;
		
		###############
		# terminal info print
		if ( $href_val_of_scr_par->{sys_parameter}{screen_prt} == 1 ){
			
			print	"\n\t::&build_multiclan_result_matrix\n",
					"\n\t", $top_and_score	,
					"\n\troot-clan1\t"	, $root_clan_1 ,
					"\n\troot-clan2\t"	, $root_clan_2 ,
					"\n\tsister-clan1\t", $sister_pair , 
					"\n" ;
		}
		###############
	}#exit;
	######################################################################
}

sub build_quartet_alignment{

	my	$href_red_seq_of_taxa	= $_[0] ; # key1: taxon; value: gap&amb reduced sequence												IN: defined; OUT: unchanged
	my	$aref_qtaxa				= $_[1] ; # list of generated 4-taxon combination														IN: defined; OUT: unchanged
	my	$sref_alignment 		= $_[2] ; # new quartet alignment																		IN: undefined; OUT: defined
	
	######################################################################
	# called from subroutine(s)...
	# &quartet_building_and_p4_main
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	# Build all P4 quartet alignments of actual 4clan combination        #
	######################################################################
	my	$seqlength			=	length $href_red_seq_of_taxa->{$aref_qtaxa->[0]} ;
		$$sref_alignment	=	"p4_temp_qmsa.phy" ;
	my	$n_taxa				=	@$aref_qtaxa ;
	
	open  OUT, ">$$sref_alignment" or die "\nMSAFILE-ERROR: Cannot write to P4 quartet-alignment infile ", $$sref_alignment, "!\n\n" ;
	
	print OUT  "\t", $n_taxa, "\t", $seqlength, "\n" ,
				$aref_qtaxa->[0]."    ".$href_red_seq_of_taxa->{$aref_qtaxa->[0]} , "\n" ,
				$aref_qtaxa->[1]."    ".$href_red_seq_of_taxa->{$aref_qtaxa->[1]} , "\n" ,
				$aref_qtaxa->[2]."    ".$href_red_seq_of_taxa->{$aref_qtaxa->[2]} , "\n" ,
				$aref_qtaxa->[3]."    ".$href_red_seq_of_taxa->{$aref_qtaxa->[3]} ;
	
	close OUT ;
	######################################################################
}

sub build_quartet_trees{
	
	my	$aref_ptaxa				= $_[0] ; # array of quartet sequences											IN: defined; OUT: unchanged
	my	$href_clan_of_seq		= $_[1] ; # hash with sequence name (key) assigned clan (value)					IN: defined; OUT: unchanged
	my	$sref_outgroup			= $_[2] ; # outgroup-clan														IN: defined; OUT: unchanged
	my	$aref_newick_unrooted	= $_[3] ; # list of the three unrooted newick strings							IN: undefined; OUT: defined
	my	$aref_newick_rooted		= $_[4] ; # list of the three rooted newick strings								IN: undefined; OUT: defined
	my	$href_value_of_scr_par	= $_[5] ; # key1: sealion user available parameter options; value: setting		IN: defined; OUT: unchanged
	my	$href_clan4_of_tree		= $_[6] ; # key1: 4clan tree; value: 4clan code									IN: undefined; OUT: defined
	my	$sref_clan4				= $_[7] ; # Clan4 of given quartet combination									IN: defined; OUT: unchanged
	
	######################################################################
	## build quartet-trees in split code corresponding order            ##
	######################################################################
	#print "::&build_quartet_trees\n";
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	#########################################
	# build rooted quartet trees in split-pattern support corresponding order from alphabetically order list of taxa if outgroup appears as first taxon
	if		( $href_clan_of_seq->{$aref_ptaxa->[0]} eq $$sref_outgroup ){ 
		
		my	@sorted_ingroup_qt1	= sort ( $href_clan_of_seq->{$aref_ptaxa->[2]}, $href_clan_of_seq->{$aref_ptaxa->[3]}) ;
		my	@sorted_ingroup_qt2	= sort ( $href_clan_of_seq->{$aref_ptaxa->[1]}, $href_clan_of_seq->{$aref_ptaxa->[3]}) ;
		my	@sorted_ingroup_qt3	= sort ( $href_clan_of_seq->{$aref_ptaxa->[1]}, $href_clan_of_seq->{$aref_ptaxa->[2]}) ;
					
		$aref_newick_rooted->[0] = "(".$href_clan_of_seq->{$aref_ptaxa->[0]}.",(".$href_clan_of_seq->{$aref_ptaxa->[1]}.",(".$sorted_ingroup_qt1[0].",".$sorted_ingroup_qt1[1].")));" ;
		$aref_newick_rooted->[1] = "(".$href_clan_of_seq->{$aref_ptaxa->[0]}.",(".$href_clan_of_seq->{$aref_ptaxa->[2]}.",(".$sorted_ingroup_qt2[0].",".$sorted_ingroup_qt2[1].")));" ;
		$aref_newick_rooted->[2] = "(".$href_clan_of_seq->{$aref_ptaxa->[0]}.",(".$href_clan_of_seq->{$aref_ptaxa->[3]}.",(".$sorted_ingroup_qt3[0].",".$sorted_ingroup_qt3[1].")));" ;
	}
	#########################################
	
	#########################################
	# build rooted quartet trees in split-pattern support corresponding order from alphabetically order list of taxa if outgroup appears as second taxon
	elsif	( $href_clan_of_seq->{$aref_ptaxa->[1]} eq $$sref_outgroup ){ 
						
		my	@sorted_ingroup_qt1	= sort ( $href_clan_of_seq->{$aref_ptaxa->[2]}, $href_clan_of_seq->{$aref_ptaxa->[3]}) ;
		my	@sorted_ingroup_qt2	= sort ( $href_clan_of_seq->{$aref_ptaxa->[0]}, $href_clan_of_seq->{$aref_ptaxa->[2]}) ;
		my	@sorted_ingroup_qt3	= sort ( $href_clan_of_seq->{$aref_ptaxa->[0]}, $href_clan_of_seq->{$aref_ptaxa->[3]}) ;
		
		$aref_newick_rooted->[0] = "(".$href_clan_of_seq->{$aref_ptaxa->[1]}.",(".$href_clan_of_seq->{$aref_ptaxa->[0]}.",(".$sorted_ingroup_qt1[0].",".$sorted_ingroup_qt1[1].")));" ;
		$aref_newick_rooted->[1] = "(".$href_clan_of_seq->{$aref_ptaxa->[1]}.",(".$href_clan_of_seq->{$aref_ptaxa->[3]}.",(".$sorted_ingroup_qt2[0].",".$sorted_ingroup_qt2[1].")));" ;
		$aref_newick_rooted->[2] = "(".$href_clan_of_seq->{$aref_ptaxa->[1]}.",(".$href_clan_of_seq->{$aref_ptaxa->[2]}.",(".$sorted_ingroup_qt3[0].",".$sorted_ingroup_qt3[1].")));" ;
	}
	#########################################
	
	#########################################
	# build rooted quartet trees in split-pattern support corresponding order from alphabetically order list of taxa if outgroup appears as third taxon
	elsif	( $href_clan_of_seq->{$aref_ptaxa->[2]} eq $$sref_outgroup ){ 
						
		my	@sorted_ingroup_qt1	= sort ( $href_clan_of_seq->{$aref_ptaxa->[0]}, $href_clan_of_seq->{$aref_ptaxa->[1]}) ;
		my	@sorted_ingroup_qt2	= sort ( $href_clan_of_seq->{$aref_ptaxa->[1]}, $href_clan_of_seq->{$aref_ptaxa->[3]}) ;
		my	@sorted_ingroup_qt3	= sort ( $href_clan_of_seq->{$aref_ptaxa->[0]}, $href_clan_of_seq->{$aref_ptaxa->[3]}) ;
		
		$aref_newick_rooted->[0] = "(".$href_clan_of_seq->{$aref_ptaxa->[2]}.",(".$href_clan_of_seq->{$aref_ptaxa->[3]}.",(".$sorted_ingroup_qt1[0].",".$sorted_ingroup_qt1[1].")));" ;
		$aref_newick_rooted->[1] = "(".$href_clan_of_seq->{$aref_ptaxa->[2]}.",(".$href_clan_of_seq->{$aref_ptaxa->[0]}.",(".$sorted_ingroup_qt2[0].",".$sorted_ingroup_qt2[1].")));" ;
		$aref_newick_rooted->[2] = "(".$href_clan_of_seq->{$aref_ptaxa->[2]}.",(".$href_clan_of_seq->{$aref_ptaxa->[1]}.",(".$sorted_ingroup_qt3[0].",".$sorted_ingroup_qt3[1].")));" ;
	}
	
	#########################################
	# build rooted quartet trees in split-pattern support corresponding order from alphabetically order list of taxa if outgroup appears as fourth taxon
	elsif	( $href_clan_of_seq->{$aref_ptaxa->[3]} eq $$sref_outgroup ){ 
						
		my	@sorted_ingroup_qt1	= sort ( $href_clan_of_seq->{$aref_ptaxa->[0]}, $href_clan_of_seq->{$aref_ptaxa->[1]}) ;
		my	@sorted_ingroup_qt2	= sort ( $href_clan_of_seq->{$aref_ptaxa->[0]}, $href_clan_of_seq->{$aref_ptaxa->[2]}) ;
		my	@sorted_ingroup_qt3	= sort ( $href_clan_of_seq->{$aref_ptaxa->[1]}, $href_clan_of_seq->{$aref_ptaxa->[2]}) ;
		
		$aref_newick_rooted->[0] = "(".$href_clan_of_seq->{$aref_ptaxa->[3]}.",(".$href_clan_of_seq->{$aref_ptaxa->[2]}.",(".$sorted_ingroup_qt1[0].",".$sorted_ingroup_qt1[1].")));" ;
		$aref_newick_rooted->[1] = "(".$href_clan_of_seq->{$aref_ptaxa->[3]}.",(".$href_clan_of_seq->{$aref_ptaxa->[1]}.",(".$sorted_ingroup_qt2[0].",".$sorted_ingroup_qt2[1].")));" ;
		$aref_newick_rooted->[2] = "(".$href_clan_of_seq->{$aref_ptaxa->[3]}.",(".$href_clan_of_seq->{$aref_ptaxa->[0]}.",(".$sorted_ingroup_qt3[0].",".$sorted_ingroup_qt3[1].")));" ;
	}
	#########################################
	
	#########################################
	# assign rooted 4clan trees to corresponding 4clan code
	$href_clan4_of_tree->{$aref_newick_rooted->[0]} = $$sref_clan4 ;
	$href_clan4_of_tree->{$aref_newick_rooted->[1]} = $$sref_clan4 ;
	$href_clan4_of_tree->{$aref_newick_rooted->[2]} = $$sref_clan4 ;
	#########################################
	
	#########################################
	# build outgroup unrooted trees in split-pattern support corresponding order from alphabetically order list of taxa
	my	@sort_sister_qt1_l	=	sort ( $href_clan_of_seq->{$aref_ptaxa->[0]}, $href_clan_of_seq->{$aref_ptaxa->[1]} ) ;
	my	@sort_sister_qt1_r	=	sort ( $href_clan_of_seq->{$aref_ptaxa->[2]}, $href_clan_of_seq->{$aref_ptaxa->[3]} ) ;
	my	@sort_sister_qt2_l	=	sort ( $href_clan_of_seq->{$aref_ptaxa->[0]}, $href_clan_of_seq->{$aref_ptaxa->[2]} ) ;
	my	@sort_sister_qt2_r	=	sort ( $href_clan_of_seq->{$aref_ptaxa->[1]}, $href_clan_of_seq->{$aref_ptaxa->[3]} ) ;
	my	@sort_sister_qt3_l	=	sort ( $href_clan_of_seq->{$aref_ptaxa->[0]}, $href_clan_of_seq->{$aref_ptaxa->[3]} ) ;
	my	@sort_sister_qt3_r	=	sort ( $href_clan_of_seq->{$aref_ptaxa->[1]}, $href_clan_of_seq->{$aref_ptaxa->[2]} ) ;
	
	$aref_newick_unrooted->[0] = "(".$sort_sister_qt1_l[0].",".$sort_sister_qt1_l[1]."),(".$sort_sister_qt1_r[0].",".$sort_sister_qt1_r[1].");" ;
	$aref_newick_unrooted->[1] = "(".$sort_sister_qt2_l[0].",".$sort_sister_qt2_l[1]."),(".$sort_sister_qt2_r[0].",".$sort_sister_qt2_r[1].");" ;
	$aref_newick_unrooted->[2] = "(".$sort_sister_qt3_l[0].",".$sort_sister_qt3_l[1]."),(".$sort_sister_qt3_r[0].",".$sort_sister_qt3_r[1].");" ;
	#########################################
	
	#########################################
	# terminal info print
	if ( $href_value_of_scr_par->{sys_parameter}{screen_prt} == 1 ){
		
		for my $t ( @$aref_ptaxa ){ print "\n\t", $href_clan_of_seq->{$t}, "\t-> ", $t } #exit;
		
		print	"\n\n\trooted 4clan-trees:\n"					,
				"\tqt1: ",	$aref_newick_rooted->[0], "\n"		,
				"\tqt2: ",  $aref_newick_rooted->[1], "\n"		,
				"\tqt3: ",  $aref_newick_rooted->[2], "\n"		,
				"\n"											,
				"\n\tunrooted 4clan-trees:\n"					,
				"\tqt1: ",	$aref_newick_unrooted->[0], "\n"	,
				"\tqt2: ",  $aref_newick_unrooted->[1], "\n"	,
				"\tqt3: ",  $aref_newick_unrooted->[2], "\n\n"	
	} #exit;
	#########################################
	
	#########################################
	# test array list matching rooted and unrooted quartet-trees
	# matching is important for correct assignment of P4 tree expected split-pattern distributions
	my $prt_t1_l = $sort_sister_qt1_l[0].",".$sort_sister_qt1_l[1] ;	#print "\n\tt1_l\t", $prt_t1_l ;
	my $prt_t1_r = $sort_sister_qt1_r[0].",".$sort_sister_qt1_r[1] ;	#print "\n\tt1_r\t", $prt_t1_r ;
	unless ( $aref_newick_rooted->[0] =~ /$prt_t1_l|$prt_t1_r/  ){ die "\nBUG-ERROR: incorrect matching of rooted and unrooted tree 1 ", $aref_newick_rooted->[0], " in &build_quartet_trees!\n\tPlease, report BUG to software developer!\n\n" } 
	
	my $prt_t2_l = $sort_sister_qt2_l[0].",".$sort_sister_qt2_l[1] ;	#print "\n\tt2_l\t", $prt_t2_l ;
	my $prt_t2_r = $sort_sister_qt2_r[0].",".$sort_sister_qt2_r[1] ;	#print "\n\tt2_r\t", $prt_t2_r ;
	unless ( $aref_newick_rooted->[1] =~ /$prt_t2_l|$prt_t2_r/  ){ die "\nBUG-ERROR: incorrect matching of rooted and unrooted tree 2 ", $aref_newick_rooted->[1], " in &build_quartet_trees!\n\tPlease, report BUG to software developer!\n\n" }
	
	my $prt_t3_l = $sort_sister_qt3_l[0].",".$sort_sister_qt3_l[1] ;	#print "\n\tt3_l\t", $prt_t3_l ;
	my $prt_t3_r = $sort_sister_qt3_r[0].",".$sort_sister_qt3_r[1] ;	#print "\n\tt3_r\t", $prt_t3_r ;
	unless ( $aref_newick_rooted->[2] =~ /$prt_t3_l|$prt_t3_r/  ){ die "\nBUG-ERROR: incorrect matching of rooted and unrooted tree 3 ", $aref_newick_rooted->[2], " in &build_quartet_trees!\n\tPlease, report BUG to software developer!\n\n" }
	
	if ( $href_value_of_scr_par->{sys_parameter}{screen_prt} == 1 ){
		
		print	"::tree-match test between rooted and unrooted trees\n\n\tPassed!\n\n" ;
	}
	#exit;
	######################################################################
}

sub build_resultfolder{
	
	my	$href_value_of_scr_par	= $_[0] ;	# key1: method_code; key2: method_option; value: float or integer		IN: defined; OUT: unchanged
	my	$href_subfolder_q		= $_[1] ;	# key1: subfolder code ('pdf', 'tsv'... ); value: resultfolder path		IN: undefined; OUT: defined
	my	$href_subsubfolder_q	= $_[2] ;	# key1: subsubfolder code ('pdf_single', 'pdf_sumed' ); value: path		IN: undefined; OUT: defined
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	# Build result with subfolder for each output format                 #
	######################################################################
	
	#######################
	# MAIN
	my		$main_resultfolder = $href_value_of_scr_par->{sys_parameter}{o_main} ;
	mkdir	$main_resultfolder ;
	#######################
	
	#######################
	# SUBFOLDER
	$href_subfolder_q->{pdf} =  $main_resultfolder."/PDF" ;
	$href_subfolder_q->{tsv} =  $main_resultfolder."/TSV" ;
	$href_subfolder_q->{tex} =  $main_resultfolder."/TEX" ;
	$href_subfolder_q->{svg} =  $main_resultfolder."/SVG" ;
	$href_subfolder_q->{txt} =  $main_resultfolder."/TXT" ;
	$href_subfolder_q->{tre} =  $main_resultfolder."/TRE" ;
	$href_subfolder_q->{spd} =  $main_resultfolder."/SPD" ;
	$href_subfolder_q->{mcm} =  $main_resultfolder."/MCM" ;
	
	
	for my $subfolder ( sort keys %$href_subfolder_q ){ mkdir $href_subfolder_q->{$subfolder} }
	#######################
	
	#######################
	# SUBSUBFOLDER
	$href_subsubfolder_q->{pdf_single} = $main_resultfolder."/PDF/RPLOTS"  ;
	$href_subsubfolder_q->{pdf_sumed}  = $main_resultfolder."/PDF" 		;
	
	for my $subsubfolder ( sort keys %$href_subsubfolder_q ){ mkdir $href_subsubfolder_q->{$subsubfolder} }
	#######################
	
	#######################
	# DELETE possible .tsv files (in case if path /TSV/ is not new and contains old *.tsv files
	for my $tsv (<$main_resultfolder/TSV/*.tsv>){ unlink $tsv } ;
	for my $tex (<$main_resultfolder/TEX/*.tex>){ unlink $tex } ; # same with tex
	for my $svg (<$main_resultfolder/SVG/*.svg>){ unlink $svg } ; # same with clf
	for my $txt (<$main_resultfolder/TXT/*.txt>){ unlink $txt } ; # same with txt
	for my $tre (<$main_resultfolder/TRE/*.tre>){ unlink $tre } ; # same with tre
	for my $mcm (<$main_resultfolder/MCM/*.txt>){ unlink $mcm } ; # same with mcm
	for my $clf (<$main_resultfolder/CLF/*.txt>){ unlink $clf } ; # same with clf
	#######################
	
	#######################
	# if defined alignment infile, delete possible SPD files of previous analysis
	if ( $href_value_of_scr_par->{sys_parameter}{write_spd} == 1 ){ 
		
		my @sub_dirs = grep {-d} glob("$href_subfolder_q->{spd}/*");
		#for (@sub_dirs){ print $_, "\n"}
		for my $d (@sub_dirs){ rmtree $d }
	}
	#######################
	#exit;
	######################################################################
}

sub calc_4clan_tree_support{
	
	my	$sref_clan4											=	$_[0] ;	# clan4 combination of given quartet combination																										IN: defined; OUT: unchanged
	my	$sref_Nq_4clan_total								=	$_[1] ;	# total number of single 4clan quartet analyses (unfiltered)																							IN: defined; OUT: unchanged
	my	$href_hol_sup_of_clan4_of_calc_of_qtree				=	$_[2] ;	# multiple hash, listing for each of the three rooted quartet-trees of a given 4clan combination all single quartet inferred support (NAP or NAPV)		IN: defined; OUT: unchanged
	my	$href_hoh_sup_of_clan4_of_calc_of_qtree_of_fc		=	$_[3] ;	# multiple hash, storing the final support (nap or napv) related to each rooted quartet-tree															IN: defined; OUT: changed
	my	$href_hol_sup_of_type_of_mfilter					=	$_[4] ; # key1: type of supportfilter; key2: mean|median; value: list of 'tree\tmean|medianScore'																IN: undefined; OUT: defined
	my	$href_hoh_linedata_of_latexfilename_of_printcode	=	$_[5] ;	# key1: Latex table unique printcode (e.g. 'latexC'); key2: *.tex filename; value: @list of assigned table lines										IN: defined; OUT: changed
	my	$href_value_of_scr_par								=	$_[6] ;	# key1: sealion user available parameter options; value: optional parameter setting																	 	IN: defined; OUT: unchanged
	my	$href_subfolder_q									=	$_[7] ;	# hashpath to output-type specific subdirs of single quartet-pattern analyses; key: output-type (e.g. 'pdf'); value: pathdir							IN: defined; OUT: unchanged
	my	$href_textfilename_of_printcode						=	$_[8] ;	# key1: textfile unique printcode (e.g. 'textA'); value: filename for print																				IN: defined; OUT: unchanged
	
	#print "::&calc_4clan_tree_support\n";
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# &calc_mean
	# &calc_median
	######################################################################
	
	######################################################################
	# calc+sample mean and median support(s) of given 4clan ccombination # 
	######################################################################
	for my $type_support ( sort keys %{$href_hol_sup_of_clan4_of_calc_of_qtree->{$$sref_clan4}} ){
		
		######################################################################
		# calculate mean|median 4-clan tree support of single quartets given the actual support_type of the 4clan combination
		######################################################################
		my	( @trees, @median_of_trees );
		
		for my $tree ( sort keys %{$href_hol_sup_of_clan4_of_calc_of_qtree->{$$sref_clan4}{$type_support}} ){
			
			############################
			# calculate mean and median support for each suuport type of actual 4clan analysis
			my	$mean	= &calc_mean	( \@{$href_hol_sup_of_clan4_of_calc_of_qtree->{$$sref_clan4}{$type_support}{$tree}} ) ;
			my	$median	= &calc_median	( \@{$href_hol_sup_of_clan4_of_calc_of_qtree->{$$sref_clan4}{$type_support}{$tree}} ) ;
			############################
			
			############################
			push @trees, $tree ;
			push @median_of_trees, $median ;
			############################
			
			############################
			# store mean score for ternary plotting
			$href_hoh_sup_of_clan4_of_calc_of_qtree_of_fc->{$$sref_clan4}{$type_support}{$tree}{mean} = $mean	;
			############################
			
			############################
			# store tree and mean score for icebreaker matrix production
			push @{$href_hol_sup_of_type_of_mfilter->{$type_support}{mean}}, $tree."\t".$mean ;
			############################
		}
		############################
		
		############################
		# median values have to be normalized again, otherwise they don't have to sum up to 1
		my	@normed_median = &norm_values(\@median_of_trees) ;
		############################
		
		############################
		# normed list values stay in the sam eorder as given into subroutine &norm_values
		# thus, they can be easily tree assigned using the value corresponding tree-list order
		for ( 0 .. 2 ){ 
			
			$href_hoh_sup_of_clan4_of_calc_of_qtree_of_fc->{$$sref_clan4}{$type_support}{$trees[$_]}{median} = $normed_median[$_] ;
			
			############################
			# store tree and mean score for icebreaker matrix production
			push @{$href_hol_sup_of_type_of_mfilter->{$type_support}{median}}, $trees[$_]."\t".$normed_median[$_] ;
			############################
		}
		######################################################################
		
		
		######################################################################
		# data sampling for latex table LQ1                                 ##
		######################################################################
		
		####################################
		# single-line sampling for each of the three 4clan trees the rounded mean and median support
		my	$counter_tree = 0 		;	# counter for 'QT' print number
		my	%mean_support_of_tree	;	# key: 4clan-tree; value: corresponding mean support
		my	%median_support_of_tree ;	# key: 4clan-tree; value: corresponding median support
		my	%code_of_tree 			;	# key: 4clan-tree; value: 'QT' code, e.g. 'QT1'
		
		my	@sprts	= split "_", $type_support	;	# split mixture string of method and filter combination in single parts
		my	$m		= shift @sprts 				;	# take 1st list element as support-method code (e.g. 'nap'
		my	$f 									;	# remaining list elements can be only filter codes. No filter used if list is empty
		if ( @sprts >= 1 ){ $f = join "_", @sprts } else{ $f = 'none' }
		
		for my $tree ( sort keys %{$href_hol_sup_of_clan4_of_calc_of_qtree->{$$sref_clan4}{$type_support}} ){
			
			####################################
			# generate latex tableline and store tableline in latex table (latexQ1) associated %hol tablelines
			$counter_tree++ ;
			$code_of_tree{$tree} = 'QT'.$counter_tree ;
			
			my	$mean_rounded	 = sprintf( "%.2f", $href_hoh_sup_of_clan4_of_calc_of_qtree_of_fc->{$$sref_clan4}{$type_support}{$tree}{mean}	) ;
			my	$median_rounded	 = sprintf( "%.2f", $href_hoh_sup_of_clan4_of_calc_of_qtree_of_fc->{$$sref_clan4}{$type_support}{$tree}{median}	) ;
			
			$mean_support_of_tree{$tree}	= $mean_rounded ;
			$median_support_of_tree{$tree}	= $median_rounded ;
			
			if ( $href_value_of_scr_par->{sys_parameter}{mean} == 1 ){
				
				&latex_data_sample_LQ1(
					
					\%$href_hoh_linedata_of_latexfilename_of_printcode	, # key1: Latex table unique printcode (e.g. 'latexC'); key2: *.tex filename; value: @list of assigned table lines			IN: defined; OUT: changed
					\$$sref_clan4										, # clan4-combination																										IN: defined; OUT: unchanged
					\$m													, # calcultaion method (nap, napv, napr)																					IN: defined; OUT: unchanged
					\$f													, # filter method (risk1, risk2...)																							IN: defined; OUT: unchanged
					\$code_of_tree{$tree}								, # 4clan tree code (QT1, QT2, QT3)																							IN: defined; OUT: unchanged
					\$tree												, # 4clan tree																												IN: defined; OUT: unchanged
					\'mean'												, # mean|median																												IN: defined; OUT: unchanged
					\$mean_rounded										, # mean supported rounded																									IN: defined; OUT: unchanged
					\%$href_subfolder_q									, # hashpath to output-type specific subdirs of single quartet-pattern analyses; key: type (e.g. 'pdf'); value: pathdir		IN: defined; OUT: unchanged
					\%$href_textfilename_of_printcode					, # key1: textfile unique printcode (e.g. 'textA'); value: filename for print												IN: defined; OUT: unchanged
				) ;
			}
			
			if ( $href_value_of_scr_par->{sys_parameter}{median} == 1 ){
				
				&latex_data_sample_LQ1(
					
					\%$href_hoh_linedata_of_latexfilename_of_printcode	, # key1: Latex table unique printcode (e.g. 'latexC'); key2: *.tex filename; value: @list of assigned table lines			IN: defined; OUT: changed
					\$$sref_clan4										, # clan4-combination																										IN: defined; OUT: unchanged
					\$m													, # calcultaion method (nap, napv, napr)																					IN: defined; OUT: unchanged
					\$f													, # filter method (risk1, risk2...)																							IN: defined; OUT: unchanged
					\$code_of_tree{$tree}								, # 4clan tree code (QT1, QT2, QT3)																							IN: defined; OUT: unchanged
					\$tree												, # 4clan tree																												IN: defined; OUT: unchanged
					\'median'											, # mean|median																												IN: defined; OUT: unchanged
					\$median_rounded									, # median support rounded																									IN: defined; OUT: unchanged
					\%$href_subfolder_q									, # hashpath to output-type specific subdirs of single quartet-pattern analyses; key: type (e.g. 'pdf'); value: pathdir		IN: defined; OUT: unchanged
					\%$href_textfilename_of_printcode					, # key1: textfile unique printcode (e.g. 'textA'); value: filename for print												IN: defined; OUT: unchanged
				) ;
			}
			####################################
			
			####################################
			# test print
			if ( $href_value_of_scr_par->{sys_parameter}{screen_prt} == 1 ){
				
				print	"\t", $$sref_clan4, "\t", $type_support, "\t", $tree, 
						"\t", "mean: "	, $href_hoh_sup_of_clan4_of_calc_of_qtree_of_fc->{$$sref_clan4}{$type_support}{$tree}{mean}		, 
						"\t", "median: ", $href_hoh_sup_of_clan4_of_calc_of_qtree_of_fc->{$$sref_clan4}{$type_support}{$tree}{median}	, 
						"\n"
			}
			####################################
		}
		######################################################################
		
		
		######################################################################
		# data sampling for latex table LQ2                                 ##
		######################################################################
		
		####################################
		# sort trees by best support (in descending order)
		my @sorted_trees_mean	= sort { $mean_support_of_tree{$b}   <=> $mean_support_of_tree{$a}   or $b cmp $a } keys %mean_support_of_tree   ;
		my @sorted_trees_median = sort { $median_support_of_tree{$b} <=> $median_support_of_tree{$a} or $b cmp $a } keys %median_support_of_tree ;
		####################################
		
		####################################
		# calculate number of single quartet support after filtering
		my $Nq_remain = @{$href_hol_sup_of_clan4_of_calc_of_qtree->{$$sref_clan4}{$type_support}{$sorted_trees_mean[0]}} ;
        ####################################
        
        ####################################
        # sample single data lines for mean and median support
		if ( $href_value_of_scr_par->{sys_parameter}{mean} == 1 ){
			
			&latex_data_sample_LQ2(
				
				\%$href_hoh_linedata_of_latexfilename_of_printcode	, # key1: Latex table unique printcode (e.g. 'latexC'); key2: *.tex filename; value: @list of assigned table lines			IN: defined; OUT: changed
				\$$sref_clan4										, # clan4-combination																										IN: defined; OUT: unchanged
				\$m													, # calcultaion method (nap, napv, napr)																					IN: defined; OUT: unchanged
				\$f													, # filter method (risk1, risk2...)																							IN: defined; OUT: unchanged
				\$sorted_trees_mean[0]								, # best supported 4clan tree																								IN: defined; OUT: unchanged
				\$code_of_tree{$sorted_trees_mean[0]}				, # 4clan tree code (QT1, QT2, QT3)																							IN: defined; OUT: unchanged
				\'mean'												, # mean|median																												IN: defined; OUT: unchanged
				\$mean_support_of_tree{$sorted_trees_mean[0]}		, # best tree support																										IN: defined; OUT: unchanged
				\$$sref_Nq_4clan_total								, # total number of single 4clan quartet analyses (unfiltered)																IN: defined; OUT: unchanged
				\$Nq_remain											, # number of remaining quartets																							IN: defined; OUT: unchanged
				\%$href_subfolder_q									, # hashpath to output-type specific subdirs of single quartet-pattern analyses; key: type (e.g. 'pdf'); value: pathdir		IN: defined; OUT: unchanged
				\%$href_textfilename_of_printcode					, # key1: textfile unique printcode (e.g. 'textA'); value: filename for print												IN: defined; OUT: unchanged
			) ;
		}
        
        if ( $href_value_of_scr_par->{sys_parameter}{median} == 1 ){
			
			&latex_data_sample_LQ2(
				
				\%$href_hoh_linedata_of_latexfilename_of_printcode	, # key1: Latex table unique printcode (e.g. 'latexC'); key2: *.tex filename; value: @list of assigned table lines			IN: defined; OUT: changed
				\$$sref_clan4										, # clan4-combination																										IN: defined; OUT: unchanged
				\$m													, # calcultaion method (nap, napv, napr)																					IN: defined; OUT: unchanged
				\$f													, # filter method (risk1, risk2...)																							IN: defined; OUT: unchanged
				\$sorted_trees_median[0]							, # best supported 4clan tree																								IN: defined; OUT: unchanged
				\$code_of_tree{$sorted_trees_median[0]}				, # 4clan tree code (QT1, QT2, QT3)																							IN: defined; OUT: unchanged
				\'median'											, # mean|median																												IN: defined; OUT: unchanged
				\$median_support_of_tree{$sorted_trees_median[0]}	, # best tree support																										IN: defined; OUT: unchanged
				\$$sref_Nq_4clan_total								, # total number of single 4clan quartet analyses (unfiltered)																IN: defined; OUT: unchanged
				\$Nq_remain											, # number of remaining quartets																							IN: defined; OUT: unchanged
				\%$href_subfolder_q									, # hashpath to output-type specific subdirs of single quartet-pattern analyses; key: type (e.g. 'pdf'); value: pathdir		IN: defined; OUT: unchanged
				\%$href_textfilename_of_printcode					, # key1: textfile unique printcode (e.g. 'textA'); value: filename for print												IN: defined; OUT: unchanged
			) ;
		}
        ######################################################################
		
		######################################################################
		# data sampling for latex table LQT1                                ##
		######################################################################
		my	@trees_q = sort keys %{$href_hol_sup_of_clan4_of_calc_of_qtree->{$$sref_clan4}{$type_support}} ;
		
		&latex_data_sample_LQT1(
			
			\$type_support										, # support type (nap, napr...)																													IN: defined; OUT: unchanged
			\$type_support										, # support type (nap, napr...)																													IN: defined; OUT: unchanged
			\%$href_hoh_linedata_of_latexfilename_of_printcode	, # key1: Latex table unique printcode (e.g. 'latexC'); key2: *.tex filename; value: @list of assigned table lines								IN: defined; OUT: changed
			\$$sref_clan4										, # clan4-combination																															IN: defined; OUT: unchanged
			\@trees_q											, # list of 4clan-trees in same order as support in @normed_supp																				IN: defined; OUT: unchanged
			\%$href_subfolder_q									, # hashpath to output-type specific subdirs of single quartet-pattern analyses; key: output-type (e.g. 'pdf'); value: pathdir					IN: defined; OUT: unchanged
			\%$href_textfilename_of_printcode					, # key1: textfile unique printcode (e.g. 'textA'); value: filename for print																	IN: defined; OUT: unchanged
		) ;
		#push @{$href_hoh_linedata_of_latexfilename_of_printcode->{latexLQT1}{$type_support}}, $type_support."\t".$$sref_clan4.$latex_LQT1_dline."\n" ;
		######################################################################
	}
	######################################################################
}

sub calc_cf{
	
	my	$Nshortest_singleton	= $_[0] ;	# lowest number of observed|expected singleton pattern			IN: defined; OUT: unchanged
	my	$total_singleton		= $_[1] ;	# total number of the four observed|expected singleton pattern	IN: defined; OUT: unchanged
	
	######################################################################
	## calculate and return correction factor CF                        ##
	######################################################################
	#print "::&calc_cf\n";
	
	######################################################################
	# called from subroutine(s)...
	# &evaluation_quartet_singletons
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	#########################################
	# cf = 1 - ( ( 4*lowest singleton number ) / total N singletons )
	my	$cf	= sprintf ( "%.4f", 1 - ( ( $Nshortest_singleton * 4 ) / $total_singleton ) ) ;
	#########################################
	
	#########################################
	# return cf value
	return	$cf ;
	#########################################
	
	######################################################################
}

sub calc_mclan_besttree_support{
	
	my	$aref_support_scores	= $_[0]	;	# 
	
	######################################################################
	## Calculate final mean|median support of the best multi-clan tree  ##
	######################################################################
	
	######################################################################
	# called from subroutine(s)...
	# &treefishing_result_processing
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	####################################
	# if list-reference consists of three float values in descending order, calculate the final support
	if ( 	( @$aref_support_scores 	== 3 			)	&&
			( $aref_support_scores->[0]	=~ /\d+(.\d+)?/ )	&&
			( $aref_support_scores->[1]	=~ /\d+(.\d+)?/ )	&&
			( $aref_support_scores->[2]	=~ /\d+(.\d+)?/ )	){
		
		my	$final_best_tree_support = ( $aref_support_scores->[0] - $aref_support_scores->[1] ) / ( $aref_support_scores->[0] - $aref_support_scores->[2] ) ;
		
		return $final_best_tree_support ;
	}
	else{ die "\nBUG-ERROR: incorrect support values for final support calculation in &calc_mclan_besttree_support!\nPlease, report BUG to software developer!\n\n"  }
	####################################
}

sub calc_mean{

	my	$aref_values	= $_[0] 	;
	
	######################################################################
	## calculate and return mean from a list of values                  ##
	######################################################################
	#print "::&calc_mean\n";
	
	######################################################################
	# called from subroutine(s)...
	# &evaluation_quartet_singletons
	# &evaluation_quartet_inv_var
	# &calc_4clan_tree_support
	# &filter_quartets
	# &determining_mean_median_4clan_taxon_support
	# &taxon_clipping
	#
	############
	# using subroutine(s)...
	# &calc_sum
	######################################################################
	
	#########################################
	# determine number of list elements
	my	$N_values = @$aref_values	; 
	#for (@$aref_values){print "\nvalue\t", $_, "\t", $N_values }
	#########################################
	
	#########################################
	# calculate total of single values
	my	$sum_values	= &calc_sum( \@$aref_values ) ;
	#########################################
	
	#########################################
	# calculate mean value
	my	$mean = sprintf( "%.4f", ( $sum_values / $N_values ) );
	#########################################
	
	#########################################
	# return mean value
	return	$mean ;
	#########################################
	
	######################################################################
}

sub calc_median{

	my	$aref_vals	=	$_[0] ;	# list of single values		IN: defined; OUT: unchanged
	
	######################################################################
	## calculate and return median from a list of values                ##
	######################################################################
	#print "::&calc_median\n";
	
	######################################################################
	# called from subroutine(s)...
	# &calc_4clan_tree_support
	# &filter_quartets
	# &determining_mean_median_4clan_taxon_support
	# &taxon_clipping
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	#########################################
	# sort values in increasing order
	@$aref_vals	=	sort {$a <=> $b} @$aref_vals	;
	#########################################
	
	#########################################
	# determine number of list elements
	my $len =	@$aref_vals	;
	#########################################
	
	#########################################
	#odd? return median
	if	( $len%2 ){ return $aref_vals->[int($len/2)] }
	#########################################
	
	#########################################
	#even? return mean of lower and upper median
	else{ my $v_median = ( return $aref_vals->[int($len/2)-1] + $aref_vals->[int($len/2)])/2 }
	#########################################
	
	######################################################################
}

sub calc_percent{
	
	my	$value		= $_[0] ;	# single value		IN: defined; OUT: unchanged
	my	$total		= $_[1] ;	# total 			IN: defined; OUT: unchanged
	
	######################################################################
	## calculate percentage of a value based on $total value            ##
	######################################################################
	#print "::&calc_percent\n";
	
	######################################################################
	# called from subroutine(s)...
	# &main
	# &evaluation_quartet_singletons
	# &evaluation_quartet_NaNpNc_signal
	# &evaluation_quartet_corrected_NaNpNc
	# &evaluation_quartet_inv_var
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	#########################################
	# do calculation
	my	$percent 	= ( $value / $total ) * 100 ;
	#########################################
	
	#########################################
	# return calculate value
	return	$percent ;
	#########################################
	
	######################################################################
}

sub calc_scales{
	
	my	$sref_upper_limit 	= $_[0] ;	# upper scale threshold 		IN: defined; OUT: unchanged
	my	$sref_lower_limit 	= $_[1] ;	# lower scale threshold 		IN: defined; OUT: unchanged
	my	$sref_scale_factor	= $_[2] ;	# scale steps			 		IN: defined; OUT: unchanged
	my	$sref_filtertype	= $_[3] ;	# filtertype (risk or qdist)	IN: defined; OUT: unchanged
	
	######################################################################
	# called from subroutine(s)...
	# &assign_filter_scales
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	## generation of possible sclae values between defined limit-range  ##
	######################################################################
	
	#########################################
	# sample and return list of single scale values
	my	@scale_list	;
	my	$i = $$sref_lower_limit ;
	
	if ( $$sref_scale_factor >= 0 ){ 
		
		while ( $i <= $$sref_upper_limit ){ 
			
			push @scale_list, $i; 
			#print "\ni\t", $i, "\t", $scale, "\t", $upper_limit, "\t", ( $i+$scale) ;
			
			$i += $$sref_scale_factor; 
		}
		
		return @scale_list
	}
	
	else { return $$sref_lower_limit }
	#########################################
	
	######################################################################
}

sub calc_sum{
	
	my	$aref_list_values  = $_[0] ;	# list of value to sum up 	IN: defined; OUT: unchanged
	
	######################################################################
	# called from subroutine(s)...
	# &norm_values
	# &calc_mean
	# &evaluation_quartet_singletons
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	## calculate and return sum from list of single scores              ##
	######################################################################
	
	#########################################
	# calculate sum
	my	$sum ;
	for my $score (@$aref_list_values){ $sum += $score }
	#########################################
	
	#########################################
	# return sum
	return	$sum ;
	#########################################
	
	######################################################################
}

sub check_clan_infile{

	my	$href_setting_of_option		= $_[0]	;	# key: option command, value: option value						IN: defined; OUT: unchanged
	my	$sref_outgroup_clan			= $_[1]	;	# user specified outgroup clan									IN: defined; OUT: unchanged
	my	$href_hol_taxa_of_clan		= $_[2]	;	# key: defined clan name; value: list of clan assigned taxa		IN: undefined; OUT: defined
	my	$href_clan_of_seq			= $_[3] ;	# key: sequence-name; value: clan definition					IN: undefined, OUT: defined
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# &help
	######################################################################
	
	######################################################################
	# Read, check, and assign .txt formated clanfile                     #
	######################################################################
	
	####################################
	# READ IN txt file with defined clade subgroups
	# Allowed Format:
	#
	# clan-code1,seqname1,seqname2,seqname3,seqname4...\n
	# clan-code2,seqname5,seqname6,seqname7,seqname8...\n
	# clan-code3,seqname9,seqname10,seqname11,seqname12...\n
	# clan-code4,seqnameX,seqnameY,seqnameW,seqnameZ...\n
	# clan-codeX...
	#
	my $clan_infile	= $href_setting_of_option->{sys_parameter}{clfile} ;
	open IN, "<$clan_infile" or die "\nCLANFILE-ERROR: Cannot read ", $clan_infile, "!\n\n" ;
	####################################
	
	####################################
	# Defined seqnames must be unique within and between each clan!!
	# check eech clanfile-line
	my	%seen_clan_name 			;		# key1: cladename; value: N occurence
	my	%seen_tax					;		# key: taxon name; value: N occurence
	my	%seen_clan_startletter 		;
	
	while ( my $line = <IN> ){ chomp $line;
		
		$line =~ s/\s+//g ;
		
		if ( $line =~ /^\w+(,\w+)+$/ ){
			
			#############
			# split line at commata
			my @lineparts	=	split ",", $line	; 
			#for (@lineparts){print $_, "\n"}; #exit;
			#############
			
			#############
			# draw first element as clan-name
			my $clanname	=	$lineparts[0]	; 
			#print "clanname: ", $clanname, "\n"; #exit;
			#############
			
			#############
			# check first letter of defined clan to be unique
			my $startletter = substr $clanname, 0, 1 ;
			
			if ( $seen_clan_startletter{$startletter} ){
				
				print "\nCLANFILE-ERROR: Found start letter '".$startletter."' for different clans!\n\n" ;
				
				&help( \'p', \%$href_setting_of_option )
			}
			
			$seen_clan_startletter{$startletter}++ ;
			#############
			
			#############
			# go through remaining lineparts (taxa) of clanfile
			for my $tax ( 1 .. @lineparts-1 ){
				
				push @{$href_hol_taxa_of_clan->{$clanname}}, $lineparts[$tax] ;
					
				$href_clan_of_seq->{$lineparts[$tax]} = $clanname ; 
				
				#print "\nx", $lineparts[$tax], "\t", $href_clan_of_seq->{$lineparts[$tax]} ;
				
				#############
				# Check that each sequence name appears only once in given clan
				if	( $seen_tax{$lineparts[$tax]} ){	
					
					print "\nCLANFILE-ERROR: Defined sequence name ", $lineparts[$tax], " appears multiple times in clanfile!\n\n";
					
					&help( \'p', \%$href_setting_of_option )
				}
				
				else{ $seen_tax{$lineparts[$tax]}++ }
				#############
			}
			#############
			
			#############
			# Check that each defined clan name appears only once in the specified clan-infile
			if	( $seen_clan_name{$clanname} ){	
				
				print "\nCLANFILE-ERROR: Defined clan ", $clanname, " appears multiple times in clanfile!\n\n";
				
				&help( \'p', \%$href_setting_of_option )
			}
			
			else{ $seen_clan_name{$clanname}++ }
			#############
		}
	}
	close IN ; #exit;
	####################################
	
	####################################
	# terminal print of defined clans
	my $clan_counter = 0 ;
	
	print "\n\tDefined Clans:" ; 
	for ( sort keys %$href_hol_taxa_of_clan ){ $clan_counter++; print "\n\tClan ", $clan_counter, ": $_"  }
	print "\n" ;
	####################################
	
	####################################
	# Check if number of defined clans is at least 4
	my $N_clans = keys %seen_clan_name ; 
	
	print "\n\tN Clans: ", $N_clans, "\n" ;

	if ( $N_clans <= 3 ){	
		
		print "\nCLANFILE-ERROR: Number of defined clans < 4!\n\n";
		
		&help( \'p', \%$href_setting_of_option )
	}
	####################################
	
	####################################
	# Check if specified outgroup is part of the clanfile
	unless ( $href_hol_taxa_of_clan->{$$sref_outgroup_clan} ){
		
		print "\nCLANFILE-ERROR: Cannot finde defined outgroup clan: ", $$sref_outgroup_clan ,"!\n\n";
		
		&help( \'o', \%$href_setting_of_option )
	}
	####################################
}

sub check_msa_infile{

	my	$href_value_of_scr_par				= $_[0] ; # key1: sealion user available parameter options; value: optional parameter setting			IN: defined; OUT: unchanged
	my	$href_info_of_property				= $_[1]	; # key1: property; value : property value														IN: undefined; OUT: defined
	my	$href_unred_seq_of_taxon			= $_[2]	; # key1: taxonname; value sequence -> defined via ARGV (-q option)								IN: undefined; OUT: defined
	my	$href_hol_sequence_states_of_taxon	= $_[3]	; # key1: taxon name; key2: state position number; value: sequence state at that position		IN: undefined; OUT: defined
	my	$href_clan_of_seq					= $_[4]	; # key1: sequence-name; value: clan definition													IN: defined; OUT: changed
	my	$href_hol_taxa_of_clan				= $_[5]	; # key: defined clan name; value: list of clan assigned taxa									IN: defined; OUT: changed
	
	print	"\tcheck alignment file...\n";
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# &help
	# &readin_fasta
	# &readin_phylip
	######################################################################
	
	######################################################################
	# Generate latex data line, listing best final tree support data     #
	######################################################################
	
	#######################
	# read msa infile and store info in %unred_seq_of_tax
	# key: taxon name; value: taxon associated unreduced sequence
	# check of multiple taxon names and
	# non-alphanumeric signs in taxonnames
	if 		( $href_value_of_scr_par->{sys_parameter}{msa_nopath} =~ /\.fas$|\.fasta$/i	 ){ 
		
		%$href_unred_seq_of_taxon = &readin_fasta	( \%$href_value_of_scr_par ) 
	}
	
	elsif 	( $href_value_of_scr_par->{sys_parameter}{msa_nopath} =~ /\.phy$|\.phylip$/i ){ 
		
		%$href_unred_seq_of_taxon = &readin_phylip	( \%$href_value_of_scr_par ) 
	}
	
	else {  die "\nBUG-ERROR: Cannot assign 'msa_nopath' in subroutine &check_msa_infile!\nPlease, report BUG to software developer!\n\n" }
	###############
	
	###############
	# test print
	for my $taxon ( sort keys %$href_unred_seq_of_taxon ){ 
			
		#print "\n\tunred. sequence of OTU ", $taxon, ":\n", $href_unred_seq_of_taxon->{$taxon}, "\n" ;
		print "\n\tread input-sequence of OTU ", $taxon; 
	}   print "\n";
	#######################
	
	
	##############################################
	# START data check
	my %hol_sequence_states_of_taxon ;
	$href_info_of_property->{type} = 'nuc' ;
	my	%msa_unseen_clantaxa ;
	
	#######################
	# check seq lengths -> store msa length to $href_info_of_property->{length}
	# identify data type -> srore data type to $href_info_of_property->{type}
	# check of forbidden sequence signs
	for my $taxon ( sort keys %$href_clan_of_seq ){
		
		if ( $href_unred_seq_of_taxon->{$taxon} ){
			
			#######################
			# check equal sequence lengths
			unless ( $href_info_of_property->{length} ){ $href_info_of_property->{length} = length $href_unred_seq_of_taxon->{$taxon} }
			
			else{
			
				unless ( $href_info_of_property->{length} == length $href_unred_seq_of_taxon->{$taxon} ){
			
					print "\nMSAFILE-ERROR: Unequal sequence lengths in multiple-sequence alignment!\n\n";
					
					&help( \'i', \%$href_value_of_scr_par )
				}
			}
			#######################
			
			#######################
			# substitute lower cased site states to upper cased site states
			# substitute U (Uracil) to T (Thymine)
			$href_unred_seq_of_taxon->{$taxon} =~ s/(\w+)/\U$1/gi ;
			$href_unred_seq_of_taxon->{$taxon} =~ s/U/T/g ;
			# print "\nunreduced sequence after changing lower cases to upper cases and U to T taxon: $taxon\n", $unred_seq_of_tax{$taxon}, "\n";
			#######################
			
			#######################
			# store single site states in identic order as list in hol @hol_sequence_states_of_taxon
			# key: taxon names; values: list elements site states of corresponding taxon
			my @seq_states = split "", $href_unred_seq_of_taxon->{$taxon} ;
			for my $state_pos ( 0 .. @seq_states-1 ){
				
				#######################
				# Check of forbidden sequence signs
				unless ( $seq_states[$state_pos] =~ /A|C|G|T|N|Y|R|W|S|K|M|D|V|H|I|E|L|Q|F|P|B|X|-|\?/ ){
					
					print "\nMSAFILE-ERROR: Forbidden sign (", $seq_states[$state_pos], ") in sequence of OTU ", $taxon, "!\n\n";
					
					&help( \'i', \%$href_value_of_scr_par )
				}
				#######################
				
				#######################
				# Identify data  type
				# if sequence includes characters only defined for amino acid
				# change data type from nuc to aa
				# afterwards stop checking data type
				unless ( $href_info_of_property->{type} eq 'aa' ){
					
					if ( $seq_states[$state_pos] =~ /E|L|Q|F|P|X/i ){ $href_info_of_property->{type} = 'aa' }
				}
				#######################
				
				$href_hol_sequence_states_of_taxon->{$taxon}[$state_pos] = $seq_states[$state_pos]
			}
			#######################
		}
		
		else { 
			
			print "\n\t\t****","\n\t\t!NOTE! OTU ", $taxon, " of clan ", $href_clan_of_seq->{$taxon}, " missing in alignment file!\n\t\t****\n" ;
					
			$msa_unseen_clantaxa{$taxon} = $href_clan_of_seq->{$taxon} ;
			unlink $href_clan_of_seq->{$taxon} 
		}
	}
	#######################
	
	#######################
	# exclude clanfile defined taxa from further analyses if they are not present in the alignment
	for my $taxon ( sort keys %msa_unseen_clantaxa ){
		
		my @clantaxa = @{$href_hol_taxa_of_clan->{$href_clan_of_seq->{$taxon}}} ;
		
		@{$href_hol_taxa_of_clan->{$msa_unseen_clantaxa{$taxon}}} = () ;
		
		for my $taxon ( @clantaxa ){ 
			
			unless ( $msa_unseen_clantaxa{$taxon} ){ 
				
				push @{$href_hol_taxa_of_clan->{$href_clan_of_seq->{$taxon}}}, $taxon 
			}
		}
	}
	#######################
	
	#######################
	# terminal print out of choosen parameter settings
	print   "\n\n\tMSA Sequence Type:\t", $href_info_of_property->{type} ,
			"\n\tMSA Sequence Length:\t", $href_info_of_property->{length}, "bp\n" ;
	#######################

	##############################################
}

sub check_newick_chars{
	
	my	$href_value_of_scr_par	= $_[0] ;	# key1: sealion user available parameter options; value: optional parameter setting				IN: defined; OUT: changed
	my	$sref_newick_tree		= $_[1] ;	# newick-tree string																			IN: defined; OUT: unchanged
	my	$sref_cl4_code			= $_[2] ;	# actual 4clan code																				IN: defined; OUT: unchanged 
	
	######################################################################
	## check char-length of newick tree for R plot printing acceptance  ##
	######################################################################
	#print	"::&check_newick_chars\n";
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	if ( length $$sref_newick_tree > $href_value_of_scr_par->{sys_parameter}{t_qtrPrt} ){
		
		$href_value_of_scr_par->{sys_parameter}{flag_qtrPrt} = 1 ;
		
		if ( $href_value_of_scr_par->{sys_parameter}{screen_prt} == 1 ){
			
			print	"\n\tLength of rooted 4clan (".$$sref_cl4_code.") newick larger than defined for grafic print-out!" ,
					"\n\tAll 4clan tree are replaced in grafic prints by constant 'Qtree'coding...\n\n"
		}
	}
	else{
		
		$href_value_of_scr_par->{sys_parameter}{flag_qtrPrt} = 0 ;
		
		if ( $href_value_of_scr_par->{sys_parameter}{screen_prt} == 1 ){
		
			print	"\n\tLength of rooted 4clan (".$$sref_cl4_code.") newick below max N chars." ,
					"\n\tNewick string accepted in grafic prints\n\n";
		}
	}
	#exit;
}

sub check_taxon_quartet_number{

	my	$aref_taxon_quartets	= $_[0]	; # list of generated 4-taxon combinations														IN: defined; OUT: changed
	my	$href_set_opt			= $_[1]	; #	key1: sealion user available parameter options; value: optional parameter setting			IN: defined; OUT: unchanged
	my	$sref_clan4				= $_[2]	; # 4-clan string, separated by doublepoints (e.g.: c1:c2:c3:c4)								IN: defined; OUT: unchanged
	
	######################################################################
	# called from subroutine(s)...
	# &quartet_building_and_p4_main
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	# Screen request if Nq > Max. Nq -> Nq = Max. Nq by Random|Neglected #
	######################################################################
	
	##############################################
	# If Nq > Max. Limit Nq...
	if ( @$aref_taxon_quartets > $href_set_opt->{qlimit}{qnumber} ){
		
		#######################
		# user request (unless blocked by -u option)
		my $answer	=	0 ;
		my $nq 		=	 @$aref_taxon_quartets;
		
		unless ( $href_set_opt->{sys_parameter}{request} == 1 ){

			print	"\n\t---------------------------------------------------\n" ,
					"\n\tNOTE: N quartets > max. limit of allowed quartets\n\t" ,
					"for 4clan combination ", $sref_clan4, "!\n\n\t"			,
					"N quartets: ", $nq, "\n\t" 								,
					"Max. limit: ", $href_set_opt->{qlimit}{qnumber}, "\n\n\t"	,
					"Press <q>	<enter>\tto Quit\n\t"							,
					"Press <r>	<enter>\tto Choose ", $href_set_opt->{sys_parameter}{request} ," Quartet(s) Randomly\n\t",
					"Press		<enter>\tto Proceed\n\n\tCommand: "				,
					"\n\t---------------------------------------------------\n" ;

			chomp ( $answer = <STDIN> ) ;
			
			#######################
			# if user choose q (quit)
			if		( $answer =~ /q/i )	{ print "\n\tQuit ", $href_set_opt->{sys_parameter}{pen_bin} ,"!\n\n"; exit }
			#######################
		}
		#######################
		
		#######################
		# if user choose random sampling up to Max. Limit Nq
		if	(	(	$answer									=~	/r/i	)	||
				(	$href_set_opt->{sys_parameter}{request}	==	1		)	){

			print "\n\t", $href_set_opt->{qlimit}{qnumber} ," Quartet(s) drawn by random!\n";

			my $quartet_counter	= 1	;
			my @new_quartets	= ();
			
			#######################
			# random quartet sampling
			while ( $quartet_counter <= $href_set_opt->{qlimit}{qnumber} ){
				
				my		$random_number = int(rand(@$aref_taxon_quartets)) - 1 ;
				push 	@new_quartets, $aref_taxon_quartets->[$random_number] ;
				splice	@$aref_taxon_quartets, $random_number, 1 ;
				
				$quartet_counter++
			}
			#######################
			
			@$aref_taxon_quartets = @new_quartets ;
			@new_quartets = () ;
		}
		#######################
		
		#######################
		# proceed with original Nq
		else	{ print "\n\tProceeding with all quartets!\n" }
		#######################
	}
	##############################################

	##############################################
	# Screen print
	my	$quartet_counter = 1 ;
	
	open OUTlist, ">>SeaLion_list_of_analysed_quartets.txt" ;
	for ( @$aref_taxon_quartets ){  print OUTlist $$sref_clan4, "\tq", $quartet_counter, "\t", $_, "\n"; $quartet_counter++ }
	close OUTlist ; 
	
	$quartet_counter = () ; #exit;
	##############################################
	
	######################################################################
}

sub config_menu{
	
	my	$href_value_of_scr_par		= $_[0] ; # key1: available options; value: parameter setting		IN: defined; OUT: unchanged
	my	$sref_config_line_new		= $_[1] ; # via config menu specified parameter line				IN: undefined; OUT: defined
	
	######################################################################
	# called from subroutine(s)...
	# &argv_handling
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	##############################
	## Print parameter setting and available command options
	print	"\n",
			"\t------------------------------------------------------------\n",
			"\tSealion Parameter Configuration Menu\n",
			"\t------------------------------------------------------------\n",
			"\tSTART SeaLion:\t\t type\t-s\t\t<enter>\t\t|(Setup)",
			"\n\t\t\t\t\t\t\t\t\t--------------------"																	,
			"\n\tMSA infile:\t\t type\t-i filename\t<enter>\t\t('"					, $href_value_of_scr_par->{menu_parameter}{msa_infile}				,"')",
			"\n\tSPD infile Analysis:\t type\t-restart\t<enter>\t\t('"				, $href_value_of_scr_par->{menu_parameter}{restart_setting}			,"')",
			"\n\tSPD infile DIR:\t\t type\t-imain SPDDIR\t<enter>\t\t('"			, $href_value_of_scr_par->{menu_parameter}{spd_folder}				,"/*')",
			"\n"																	,
			"\n\tClade infile:\t\t type\t-p filename\t<enter>\t\t('"				, $href_value_of_scr_par->{menu_parameter}{clfile}					,"')",
			"\n\tOutgroup-clade:\t\t type\t-o clade-code\t<enter>\t\t('"			, $href_value_of_scr_par->{menu_parameter}{outgroup}				,"')",
			"\n\tNmin quartet sites:\t type\t-l integer\t<enter>\t\t('"				, $href_value_of_scr_par->{menu_parameter}{minlengths}				,"')",
			"\n\tNmax quartets:\t\t type\t-M integer\t<enter>\t\t('"					, $href_value_of_scr_par->{menu_parameter}{maxquartets}				,"')",
			"\n\tSupport Type:\t\t type\t-average\t<enter>\t\t('"					, $href_value_of_scr_par->{menu_parameter}{average}					,"')",
			"\n"																	,
			"\n\tQuartet Filter 'RISK':\t type\t-r integer\t<enter>\t\t('"			, $href_value_of_scr_par->{menu_parameter}{risk_default}			,"')",
			"\n\t'RISK' upper limit:\t type\t-turisk float\t<enter>\t\t('"			, $href_value_of_scr_par->{menu_parameter}{upperlimit}				,"')",
			"\n\t'RISK' lower limit:\t type\t-tlrisk float\t<enter>\t\t('"			, $href_value_of_scr_par->{menu_parameter}{lowerlimit}				,"')",
			"\n\t'RISK' scaling:\t\t type\t-tsrisk float\t<enter>\t\t('"			, $href_value_of_scr_par->{menu_parameter}{rscale}					,"')",
			"\n"																	,
			"\n\tQuartet Filter 'DIST':\t type\t-d\t\t<enter>\t\t('"				, $href_value_of_scr_par->{menu_parameter}{qdist_default}			,"')",
			"\n\t'DIST' upper limit:\t type\t-tudist float\t<enter>\t\t('"			, $href_value_of_scr_par->{menu_parameter}{upperlimit_d}			,"')",
			"\n\t'DIST' lower limit:\t type\t-tldist float\t<enter>\t\t('"			, $href_value_of_scr_par->{menu_parameter}{lowerlimit_d}			,"')",
			"\n\t'DIST' scaling:\t\t type\t-tsdist float\t<enter>\t\t('"			, $href_value_of_scr_par->{menu_parameter}{dscale}					,"')",
			"\n"																	,
			"\n\tOutput DIR:\t\t type\t-omain OutDIR\t<enter>\t\t('"				, $href_value_of_scr_par->{menu_parameter}{output_folder}			,"/*')",
			"\n\tTable/Graphic Plotting:\t type\t-prt integer\t<enter>\t\t('"		, $href_value_of_scr_par->{menu_parameter}{prt_setting}				,"')",
			"\n"																	,
			"\n\tP4 Model (AA):\t\t type\t-m integer\t<enter>\t\t('"				, $href_value_of_scr_par->{menu_parameter}{model_aa}				,"')",
			"\n\tP4 Model (NUC):\t\t type\t-m integer\t<enter>\t\t('"				, $href_value_of_scr_par->{menu_parameter}{model_nu}				,"')",
			"\n\tP4 Alpha:\t\t type\t-a float\t<enter>\t\t('"						, $href_value_of_scr_par->{menu_parameter}{start_alpha}				,"')",
			"\n\tP4 pINV:\t\t type\t-I float\t<enter>\t\t('"						, $href_value_of_scr_par->{menu_parameter}{start_pinv}				,"')",
			"\n"																	,
			"\n\tRY coding:\t\t type\t-ry\t\t<enter>\t\t('"							, $href_value_of_scr_par->{menu_parameter}{ry_setting}				,"')",
			"\n\tScript Query:\t\t type\t-u\t\t<enter>\t\t('"						, $href_value_of_scr_par->{menu_parameter}{query_setting}			,"')",
			"\n"																	,
			"\n\tHELP General:\t\t type\t-h\t\t<enter>"								,
			"\n\tHELP Specific:\t\t type\t-h command\t<enter>"						,
			"\n\tQUIT:\t\t\t type\t-q\t\t<enter>"									,
			"\n\tPREFACE:\t\t type\t-P\t\t<enter>"									,
			"\n\t------------------------------------------------------------\n",
			"\n\tCOMMAND:\t "
	;
	##############################
	#exit;
	##############################
	## READ IN and return user command
	chomp ( $$sref_config_line_new = <STDIN> );
	#print "\n\t", $$sref_config_line_new, "\n"; #exit;
	##############################
}

sub check_user_parameter{
	
	my	$href_value_of_scr_par	= $_[0] ;	# key1: sealion user available parameter options; value: optional parameter setting			IN: defined; OUT: changed
	my	$sref_commandline		= $_[1] ;	# command line or config menu specied string: 'command parameter'							IN: defined; OUT: unchanged
	
	######################################################################
	# called from subroutine(s)...
	# &argv_handling
	#
	############
	# using subroutine(s)...
	# &help
	# &preface
	######################################################################
	
	##############################################
	# parameter options:
	# -r (visible)				: quartet-filter method (risk)		; integer: 0 to 15 ;
	# -d (visible)				: quartet-filter method (dist)		; no extra parameter needed ;
	# -i (visible)				: MSA-infile						; string with allowed path extension ;
	# -p (visible)				: Clan-infile						; string with allowed path extension ;
	# -a (visible)				: P4 alpha start parameter			; float: 0 to 100 ;
	# -I (visible)				: P4 pinv start parameter			; float: 0 to 1 ;
	# -m (visible)				: P4 model							; string ;
	# -M (visible)				: Min. quartet sequence-length  	; integer ;
	# -l (visible)				: Max. quartet number of 4clan		; integer ;
	# -o (visible)				: Outgroup-clan						; string ;
	# -u (visible)				: Reject user requests				; no extra parameter needed ;
	# -tlrisk (visible)			; lower risk-filter bound			; float 0 to 1 ;
	# -turisk (visible)			; upper risk-filter bound			; float 0 to 1 ;
	# -tsrisk (visible)			; scale risk-filter bound			; float 0 to 1 ;
	# -tldist (visible)			; lower qdist-filter bound			; float 0 to 1 ;
	# -tudist (visible)			; upper qdist-filter bound			; float 0 to 1 ;
	# -tsdist (visible)			; scale qdist-filter bound			; float 0 to 1 ;
	# -prt (visible)			: latex & R printing				; integer 0 to 3 ;
	# -average (visible)		: measure of average support		; integer 1 to 3 ;
	# -imain (visible)			: main SPD-infile pathDIR			: string ;
	# -omain (visible)			: Main result pathDIR				: string ;
	# -restart (visible)		: activate|deactivate SPD readIN	; no extra parameter needed ;
	# -ry (visible)				: activate|deactivate RY coding		; no extra parameter needed ;
	# -s						; process start-command				; no extra parameter needed ;
	# -P						; print script-preface to terminal	; no extra parameter needed ;
	# -h						; help menu							; no extra parameter needed, but possible ;
	# -q						; quit script						; no extra parameter needed ;
	##############################################

	##############################################
	# help menu options
	# -h					->	general help menu
	# -h -parameter option	->	detailed help menu of parameter option
	##############################################
	
	##############################################
	# Check commandline options and parameter
	##############################################
	my	$cmd		=	$$sref_commandline		;
		$cmd		=~	s/\s+|\t+/ /g 			;
		$cmd		=~	s/^ | $//g 				;	# substitute possible blanks at ethe beginning or end of the string
		
	my	@commands	=	split "-", $cmd ; 
	#for (@commands){print "\n\tcmd\t", $_, "\n"} #exit;
	
	shift	@commands ;
	
	for my $cmdline ( @commands ){
		
		my	@command_signs	= split	" ", $cmdline ;
		my	$command		= $command_signs[0] ;
		my	$parameter		= 0 ;
		if 	( $command_signs[1] ){ $parameter = $command_signs[1] } 
		#print "\nfound cmd\t", $command, "\nfound par\t", $parameter,"\n"; exit;
		
		if 		( ( $command =~ /^m$/ 			 ) && ( $parameter ) ){ 
			
			if 	  ( $parameter == 1  	){ $href_value_of_scr_par->{p4_parameter}{model_nu} = 'JC'  		}
			elsif ( $parameter == 2  	){ $href_value_of_scr_par->{p4_parameter}{model_nu} = 'F81' 		}
			elsif ( $parameter == 3  	){ $href_value_of_scr_par->{p4_parameter}{model_nu} = 'K2P' 		}
			elsif ( $parameter == 4  	){ $href_value_of_scr_par->{p4_parameter}{model_nu} = 'HKY' 		}
			elsif ( $parameter == 5  	){ $href_value_of_scr_par->{p4_parameter}{model_nu} = 'GTR' 		}
			elsif ( $parameter == 6  	){ $href_value_of_scr_par->{p4_parameter}{model_aa} = 'd78' 		}
			elsif ( $parameter == 7  	){ $href_value_of_scr_par->{p4_parameter}{model_aa} = 'jtt' 		}
			elsif ( $parameter == 8  	){ $href_value_of_scr_par->{p4_parameter}{model_aa} = 'wag' 		}
			elsif ( $parameter == 9  	){ $href_value_of_scr_par->{p4_parameter}{model_aa} = 'mtrev24'		}
			elsif ( $parameter == 10 	){ $href_value_of_scr_par->{p4_parameter}{model_aa} = 'lg' 			}
			elsif ( $parameter == 11 	){ $href_value_of_scr_par->{p4_parameter}{model_aa} = 'blosum62'	}
			
			else  { print "\n\tPARAMETER-ERROR: Unallowed value: ", $parameter, " for command ", $command, "!\n\n"; &help( \'m', \%$href_value_of_scr_par ) }
		}
		elsif 	( 	$command =~ /^m$/	 	 	 ) { print "\n\tPARAMETER-ERROR: Missing parameter for command '-", $command, "'!\n\n"; &help( \$command, \%$href_value_of_scr_par 			) }
		
		
		elsif 	( ( $command =~ /^r$/ 			 ) && ( $parameter ) ){ 
			
			if 	  ( $parameter =~ /x/ 	){ $href_value_of_scr_par->{seq_filterRisk}{risk1} = 0 ; $href_value_of_scr_par->{seq_filterRisk}{risk2} = 0 }
			elsif ( $parameter == 1 	){ $href_value_of_scr_par->{seq_filterRisk}{risk1} = 1 ; $href_value_of_scr_par->{seq_filterRisk}{risk2} = 0 }
			elsif ( $parameter == 2 	){ $href_value_of_scr_par->{seq_filterRisk}{risk1} = 0 ; $href_value_of_scr_par->{seq_filterRisk}{risk2} = 1 }
			elsif ( $parameter == 3 	){ $href_value_of_scr_par->{seq_filterRisk}{risk1} = 1 ; $href_value_of_scr_par->{seq_filterRisk}{risk2} = 1 }
				
			else  { print "\n\tPARAMETER-ERROR: Unallowed value: ", $parameter, " for command ", $command, "!\n\n"; &help( \'r', \%$href_value_of_scr_par ) }
		}
		elsif 	( 	$command =~ /^r$/	 	 	 ) { print "\n\tPARAMETER-ERROR: Missing parameter for command '-", $command, "'!\n\n"; &help( \$command, \%$href_value_of_scr_par 			) }
		
		
		elsif 	( ( $command =~ /^prt$/ 		 ) && ( $parameter ) ){ 
			
			if 	  ( $parameter =~ /x/ 	){ $href_value_of_scr_par->{sys_parameter}{r_prt} = 0 ; $href_value_of_scr_par->{sys_parameter}{l_prt} = 0 }
			elsif ( $parameter == 1 	){ $href_value_of_scr_par->{sys_parameter}{r_prt} = 1 ; $href_value_of_scr_par->{sys_parameter}{l_prt} = 0 }
			elsif ( $parameter == 2 	){ $href_value_of_scr_par->{sys_parameter}{r_prt} = 0 ; $href_value_of_scr_par->{sys_parameter}{l_prt} = 1 }
			elsif ( $parameter == 3 	){ $href_value_of_scr_par->{sys_parameter}{r_prt} = 1 ; $href_value_of_scr_par->{sys_parameter}{l_prt} = 1 }
			
			else  { print "\n\tPARAMETER-ERROR: Unallowed value: ", $parameter, " for command ", $command, "!\n\n"; &help( \'prt', \%$href_value_of_scr_par ) }
		}
		elsif 	( 	$command =~ /^prt$/ 	 	 ) { print "\n\tPARAMETER-ERROR: Missing parameter for command '-", $command, "'!\n\n"; &help( \$command, \%$href_value_of_scr_par 			) }
		
		
		elsif 	( ( $command =~ /^average$/ 	 ) && ( $parameter ) ){ 
			
			if 	  ( $parameter == 1 	){ $href_value_of_scr_par->{sys_parameter}{mean} = 1 ; $href_value_of_scr_par->{sys_parameter}{median} = 0 }
			elsif ( $parameter == 2 	){ $href_value_of_scr_par->{sys_parameter}{mean} = 0 ; $href_value_of_scr_par->{sys_parameter}{median} = 1 }
			elsif ( $parameter == 3 	){ $href_value_of_scr_par->{sys_parameter}{mean} = 1 ; $href_value_of_scr_par->{sys_parameter}{median} = 1 }
			
			else  { print "\n\tPARAMETER-ERROR: Unallowed value: ", $parameter, " for command ", $command, "!\n\n"; &help( \'average', \%$href_value_of_scr_par ) }
		}
		elsif 	( 	$command =~ /^average$/ 	 ) { print "\n\tPARAMETER-ERROR: Missing parameter for command '-", $command, "'!\n\n"; &help( \$command, \%$href_value_of_scr_par 			) }
		
		
		elsif 	( ( $command =~ /^p$/ 		 	 ) && ( $parameter ) ){ 
			
			if 	( 	$parameter =~ /^\w+\.txt$/ 	 ) { $href_value_of_scr_par->{sys_parameter}{clfile}		= $parameter }
			else								   { print "\n\tPARAMETER-ERROR: Unallowed '-", $command, "' parameter: ", $parameter, "!\n\n"; &help( \$command, \%$href_value_of_scr_par 	) }
		}
		elsif 	( 	$command =~ /^p$/ 	  		 ) { print "\n\tPARAMETER-ERROR: Missing parameter for command '-", $command, "'!\n\n"; &help( \$command, \%$href_value_of_scr_par 			) }
		
		
		elsif 	( ( $command =~ /^i$/ 		 	 ) && ( $parameter ) ){ 
			
			if 	( 	$parameter =~ /\w+\.fas$|\w+\.phy$|\w+\.fasta$/i ){ $href_value_of_scr_par->{sys_parameter}{msafile}	= $parameter }
			else								   { print "\n\tPARAMETER-ERROR: Unallowed '-", $command, "' parameter: ", $parameter, "!\n\n"; &help( \$command, \%$href_value_of_scr_par 	) }
		}
		elsif 	( 	$command =~ /^i$/ 	  		 ) { print "\n\tPARAMETER-ERROR: Missing parameter for command '-", $command, "'!\n\n"; &help( \$command, \%$href_value_of_scr_par 			) }
		
		
		elsif 	( ( $command =~ /^o$/ 		 	 ) && ( $parameter ) ){ 
			
			if 	( 	$parameter =~ /^\w+$/ 		 ) { $href_value_of_scr_par->{sys_parameter}{outgroup}	 = $parameter }
			else								   { print "\n\tPARAMETER-ERROR: Unallowed '-", $command, "' parameter: ", $parameter, "!\n\n"; &help( \$command, \%$href_value_of_scr_par 	) }
		}
		elsif 	( 	$command =~ /^o$/ 	  		 ) { print "\n\tPARAMETER-ERROR: Missing parameter for command '-", $command, "'!\n\n"; &help( \$command, \%$href_value_of_scr_par 			) }
		
		
		elsif 	( ( $command =~ /^imain$/ 		 ) && ( $parameter ) ){ 
			
			if 	( 	$parameter =~ /^.?(\/?\w+)+$/ ) { $href_value_of_scr_par->{sys_parameter}{i_main}	 = $parameter }
			else								   { print "\n\tPARAMETER-ERROR: Unallowed '-", $command, "' parameter: ", $parameter, "!\n\n"; &help( \$command, \%$href_value_of_scr_par 	) }
		}
		elsif 	( 	$command =~ /^imain$/ 	  	 ) { print "\n\tPARAMETER-ERROR: Missing parameter for command '-", $command, "'!\n\n"; &help( \$command, \%$href_value_of_scr_par 			) }
		
		
		elsif 	( ( $command =~ /^omain$/ 		 ) && ( $parameter ) ){ 
			
			if 	( 	$parameter =~ /^.?(\/?\w+)+$/ ) { $href_value_of_scr_par->{sys_parameter}{o_main}	 = $parameter }
			else								   { print "\n\tPARAMETER-ERROR: Unallowed '-", $command, "' parameter: ", $parameter, "!\n\n"; &help( \$command, \%$href_value_of_scr_par 	) }
		}
		elsif 	( 	$command =~ /^omain$/ 	  	 ) { print "\n\tPARAMETER-ERROR: Missing parameter for command '-", $command, "'!\n\n"; &help( \$command, \%$href_value_of_scr_par 			) }
		
		
		elsif 	( ( $command =~ /^M$/ 		 	 ) && ( $parameter ) ){ 
			
			if 		( ( $parameter =~ /^\d+(\.\d+)?$/ ) && ( $parameter <  100 ) ) { print "\n\tPARAMETER-ERROR: Unallowed '-", $command, "' parameter: ", $parameter, "!\n\n"; &help( \$command, \%$href_value_of_scr_par ) }
			elsif 	( ( $parameter =~ /^\d+(\.\d+)?$/ ) && ( $parameter >= 100 ) ) { $href_value_of_scr_par->{qlimit}{qlength}		 = $parameter }
			else																   { print "\n\tPARAMETER-ERROR: Unallowed '-", $command, "' parameter: ", $parameter, "!\n\n"; &help( \$command, \%$href_value_of_scr_par ) }
		}
		elsif 	( 	$command =~ /^M$/ 		 	 ) { print "\n\tPARAMETER-ERROR: Missing parameter for command '-", $command, "'!\n\n"; &help( \$command, \%$href_value_of_scr_par 		   ) }
			
		
		elsif 	( ( $command =~ /^l$/ 		 	 ) && ( $parameter ) ){ 
			
			if 		( ( $parameter =~ /^\d+(\.\d+)?$/ ) && ( $parameter <  1   ) ) { print "\n\tPARAMETER-ERROR: Unallowed '", $command, "' parameter: ", $parameter, "!\n\n"; &help( \$command, \%$href_value_of_scr_par  ) }
			elsif 	( ( $parameter =~ /^\d+(\.\d+)?$/ ) && ( $parameter >= 1   ) ) { $href_value_of_scr_par->{qlimit}{qnumber}		 = $parameter }
			else																   { print "\n\tPARAMETER-ERROR: Unallowed '", $command, "' parameter: ", $parameter, "!\n\n"; &help( \$command, \%$href_value_of_scr_par  ) }
		}
		elsif 	( 	$command =~ /^l$/ 		 	 ) { print "\n\tPARAMETER-ERROR: Missing parameter for command '-", $command, "'!\n\n"; &help( \$command, \%$href_value_of_scr_par 		   ) }
		
		
		elsif 	( ( $command =~ /^tlrisk$/ 		 ) && ( $parameter ) ){ 
			
			if 		( ( $parameter =~ /^\d+(\.\d+)?$/ ) && ( $parameter <  0   ) ) { print "\n\tPARAMETER-ERROR: Unallowed '", $command, "' parameter: ", $parameter, "!\n\n"; &help( \$command, \%$href_value_of_scr_par  ) }
			elsif 	( ( $parameter =~ /^\d+(\.\d+)?$/ ) && ( $parameter >  1   ) ) { print "\n\tPARAMETER-ERROR: Unallowed '", $command, "' parameter: ", $parameter, "!\n\n"; &help( \$command, \%$href_value_of_scr_par  ) }
			elsif 	( ( $parameter =~ /^\d+(\.\d+)?$/ ) && ( $parameter <= 1   ) ) { $href_value_of_scr_par->{seq_filterRisk}{tlrisk}		 = $parameter }
			else																   { print "\n\tPARAMETER-ERROR: Unallowed '", $command, "' parameter: ", $parameter, "!\n\n"; &help( \$command, \%$href_value_of_scr_par  ) }
		}
		elsif 	( 	$command =~ /^tlrisk$/ 		 ) { print "\n\tPARAMETER-ERROR: Missing parameter for command '-", $command, "'!\n\n"; &help( \$command, \%$href_value_of_scr_par 		   ) }
		
		
		elsif 	( ( $command =~ /^turisk$/ 		 ) && ( $parameter ) ){ 
			
			if 		( ( $parameter =~ /^\d+(\.\d+)?$/ ) && ( $parameter <  0   ) ) { print "\n\tPARAMETER-ERROR: Unallowed '", $command, "' parameter: ", $parameter, "!\n\n"; &help( \$command, \%$href_value_of_scr_par  ) }
			elsif 	( ( $parameter =~ /^\d+(\.\d+)?$/ ) && ( $parameter >  2   ) ) { print "\n\tPARAMETER-ERROR: Unallowed '", $command, "' parameter: ", $parameter, "!\n\n"; &help( \$command, \%$href_value_of_scr_par  ) }
			elsif 	( ( $parameter =~ /^\d+(\.\d+)?$/ ) && ( $parameter <= 2   ) ) { $href_value_of_scr_par->{seq_filterRisk}{turisk}		 = $parameter }
			else																   { print "\n\tPARAMETER-ERROR: Unallowed '", $command, "' parameter: ", $parameter, "!\n\n"; &help( \$command, \%$href_value_of_scr_par  ) }
		}
		elsif 	( 	$command =~ /^turisk$/ 		 ) { print "\n\tPARAMETER-ERROR: Missing parameter for command '-", $command, "'!\n\n"; &help( \$command, \%$href_value_of_scr_par 		   ) }
		
		
		elsif 	( ( $command =~ /^tsrisk$/ 		 ) && ( $parameter ) ){ 
			
			if 		( ( $parameter =~ /^\d+(\.\d+)?$/ ) && ( $parameter <  0   ) ) { print "\n\tPARAMETER-ERROR: Unallowed '-", $command, "' parameter: ", $parameter, "!\n\n"; &help( \$command, \%$href_value_of_scr_par ) }
			elsif 	( ( $parameter =~ /^\d+(\.\d+)?$/ ) && ( $parameter >  1   ) ) { print "\n\tPARAMETER-ERROR: Unallowed '-", $command, "' parameter: ", $parameter, "!\n\n"; &help( \$command, \%$href_value_of_scr_par ) }
			elsif 	( ( $parameter =~ /^\d+(\.\d+)?$/ ) && ( $parameter <= 1   ) ) { $href_value_of_scr_par->{seq_filterRisk}{rscale}		 = $parameter }
			else																   { print "\n\tPARAMETER-ERROR: Unallowed '-", $command, "' parameter: ", $parameter, "!\n\n"; &help( \$command, \%$href_value_of_scr_par ) }
		}
		elsif 	( 	$command =~ /^tsrisk$/ 		 ) { print "\n\tPARAMETER-ERROR: Missing parameter for command '-", $command, "'!\n\n"; &help( \$command, \%$href_value_of_scr_par 		   ) }
		
		
		elsif 	( ( $command =~ /^tudist$/ 		 ) && ( $parameter ) ){ 
			
			if 		( ( $parameter =~ /^\d+(\.\d+)?$/ ) && ( $parameter <  0   ) ) { print "\n\tPARAMETER-ERROR: Unallowed '-", $command, "' parameter: ", $parameter, "!\n\n"; &help( \$command, \%$href_value_of_scr_par ) }
			elsif 	( ( $parameter =~ /^\d+(\.\d+)?$/ ) && ( $parameter >  1   ) ) { print "\n\tPARAMETER-ERROR: Unallowed '-", $command, "' parameter: ", $parameter, "!\n\n"; &help( \$command, \%$href_value_of_scr_par ) }
			elsif 	( ( $parameter =~ /^\d+(\.\d+)?$/ ) && ( $parameter <= 1   ) ) { $href_value_of_scr_par->{seq_filterDist}{tudist}		 = $parameter }
			else																   { print "\n\tPARAMETER-ERROR: Unallowed '-", $command, "' parameter: ", $parameter, "!\n\n"; &help( \$command, \%$href_value_of_scr_par ) }
		}
		elsif 	( 	$command =~ /^tudist$/ 		 ) { print "\n\tPARAMETER-ERROR: Missing parameter for command '-", $command, "'!\n\n"; &help( \$command, \%$href_value_of_scr_par 		   ) }
		
		
		elsif 	( ( $command =~ /^tldist$/ 		 ) && ( $parameter ) ){ 
			
			if 		( ( $parameter =~ /^\d+(\.\d+)?$/ ) && ( $parameter <  0   ) ) { print "\n\tPARAMETER-ERROR: Unallowed '-", $command, "' parameter: ", $parameter, "!\n\n"; &help( \$command, \%$href_value_of_scr_par ) }
			elsif 	( ( $parameter =~ /^\d+(\.\d+)?$/ ) && ( $parameter >  1   ) ) { print "\n\tPARAMETER-ERROR: Unallowed '-", $command, "' parameter: ", $parameter, "!\n\n"; &help( \$command, \%$href_value_of_scr_par ) }
			elsif 	( ( $parameter =~ /^\d+(\.\d+)?$/ ) && ( $parameter <= 1   ) ) { $href_value_of_scr_par->{seq_filterDist}{tldist}		 = $parameter }
			else																   { print "\n\tPARAMETER-ERROR: Unallowed '-", $command, "' parameter: ", $parameter, "!\n\n"; &help( \$command, \%$href_value_of_scr_par ) }
		}
		elsif 	( 	$command =~ /^tldist$/ 		 ) { print "\n\tPARAMETER-ERROR: Missing parameter for command '-", $command, "'!\n\n"; &help( \$command, \%$href_value_of_scr_par 		   ) }
		
		
		elsif 	( ( $command =~ /^tsdist$/ 		 ) && ( $parameter ) ){ 
			
			if 		( ( $parameter =~ /^\d+(\.\d+)?$/ ) && ( $parameter <  0   ) ) { print "\n\tPARAMETER-ERROR: Unallowed '-", $command, "' parameter: ", $parameter, "!\n\n"; &help( \$command, \%$href_value_of_scr_par ) }
			elsif 	( ( $parameter =~ /^\d+(\.\d+)?$/ ) && ( $parameter >  1   ) ) { print "\n\tPARAMETER-ERROR: Unallowed '-", $command, "' parameter: ", $parameter, "!\n\n"; &help( \$command, \%$href_value_of_scr_par ) }
			elsif 	( ( $parameter =~ /^\d+(\.\d+)?$/ ) && ( $parameter <= 1   ) ) { $href_value_of_scr_par->{seq_filterDist}{dscale}		 = $parameter }
			else																   { print "\n\tPARAMETER-ERROR: Unallowed '-", $command, "' parameter: ", $parameter, "!\n\n"; &help( \$command, \%$href_value_of_scr_par ) }
		}
		elsif 	( 	$command =~ /^tsdist$/ 		 ) { print "\n\tPARAMETER-ERROR: Missing parameter for command '-", $command, "'!\n\n"; &help( \$command, \%$href_value_of_scr_par 		   ) }
		
		
		elsif 	( ( $command =~ /^I$/ 			 ) && ( $parameter ) ){ 
			
			if 		( ( $parameter =~ /^\d+(\.\d+)?$/ ) && ( $parameter <  0   ) ) { print "\n\tPARAMETER-ERROR: Unallowed '-", $command, "' parameter: ", $parameter, "!\n\n"; &help( \$command, \%$href_value_of_scr_par ) }
			elsif 	( ( $parameter =~ /^\d+(\.\d+)?$/ ) && ( $parameter >  1   ) ) { print "\n\tPARAMETER-ERROR: Unallowed '-", $command, "' parameter: ", $parameter, "!\n\n"; &help( \$command, \%$href_value_of_scr_par ) }
			elsif 	( ( $parameter =~ /^\d+(\.\d+)?$/ ) && ( $parameter <= 1   ) ) { $href_value_of_scr_par->{p4_parameter}{pinv} 		= $parameter }
			else																   { print "\n\tPARAMETER-ERROR: Unallowed '-", $command, "' parameter: ", $parameter, "!\n\n"; &help( \$command, \%$href_value_of_scr_par ) }
		}
		elsif 	( 	$command =~ /^I$/ 			 ) { print "\n\tPARAMETER-ERROR: Missing parameter for command '-", $command, "'!\n\n"; &help( \$command, \%$href_value_of_scr_par 		   ) }
		
		
		elsif 	( ( $command =~ /^a$/ 			 ) && ( $parameter ) ){ 
			
			if 		( ( $parameter =~ /^\d+(\.\d+)?$/ ) && ( $parameter <  0   ) ) { print "\n\tPARAMETER-ERROR: Unallowed '-", $command, "' parameter: ", $parameter, "!\n\n"; &help( \$command, \%$href_value_of_scr_par ) }
			elsif 	( ( $parameter =~ /^\d+(\.\d+)?$/ ) && ( $parameter >  100 ) ) { print "\n\tPARAMETER-ERROR: Unallowed '-", $command, "' parameter: ", $parameter, "!\n\n"; &help( \$command, \%$href_value_of_scr_par ) }
			elsif 	( ( $parameter =~ /^\d+(\.\d+)?$/ ) && ( $parameter <= 100 ) ) { $href_value_of_scr_par->{p4_parameter}{alpha} 		= $parameter }
			else																   { print "\n\tPARAMETER-ERROR: Unallowed '-", $command, "' parameter: ", $parameter, "!\n\n"; &help( \$command, \%$href_value_of_scr_par ) }
		}
		elsif 	( 	$command =~ /^a$/ 			 ) { print "\n\tPARAMETER-ERROR: Missing parameter for command '-", $command, "'!\n\n"; &help( \$command, \%$href_value_of_scr_par 		   ) }
		
			
		elsif 	( 	$command =~ /^d$/			 ){
			
			if ( $href_value_of_scr_par->{seq_filterDist}{qdist1} ) { $href_value_of_scr_par->{seq_filterDist}{qdist1} 	= 0 } # de-activation of quartet-filter 'DIST' if preset is defined ('1')	
			else													{ $href_value_of_scr_par->{seq_filterDist}{qdist1} 	= 1 } # activation of quartet-filter 'DIST' if preset is undefined ('0')
		}
		
		
		elsif	( 	$command =~ /^restart$/		 ){
			
			if ( $href_value_of_scr_par->{sys_parameter}{restart} ) { $href_value_of_scr_par->{sys_parameter}{restart}	= 0 }	# enables P4 quartet calculations via config menu if restart has been previously defined
			else													{ $href_value_of_scr_par->{sys_parameter}{restart}	= 1 }	# disables P4 quartet calculations, needs SPD-infile path (via 'imain' option)
		}
		
		
		elsif 	( 	$command =~ /^u$/			 ){ 
			
			if ( $href_value_of_scr_par->{sys_parameter}{request} ){ $href_value_of_scr_par->{sys_parameter}{request}	= 0 }	# enables terminal user requests via config menu if defined earlier via config meu
			else												   { $href_value_of_scr_par->{sys_parameter}{request}	= 1 }	# disbales terminal user requests if N quartets in a 4clan combination is too large (maximum number of allowed quartets are drawn randomly instead)
		}
		
		
		elsif 	( 	$command =~ /^ry$/			 ){ 
			
			if ( $href_value_of_scr_par->{sys_parameter}{rycoding} ){ $href_value_of_scr_par->{sys_parameter}{rycoding}	= 0 }	# disables RY-coding of site characters via config menu if ry coding has been previously defined
			else													{ $href_value_of_scr_par->{sys_parameter}{rycoding}	= 1 }	# RY-coding of site characters in advance of the analysis
		}
		
		
		elsif 	( 	$command =~ /^s$/			 ){ $href_value_of_scr_par->{sys_parameter}{sysstart}	 = 1 }
		
		
		elsif	( ( $command =~ /^h$/i			 ) && ( $parameter == 0 ) ){ &help( \'h', \%$href_value_of_scr_par ) }
		
		elsif 	( (	$command =~ /^h$/i			 ) && ( $parameter ) ){
			
			if (  $parameter =~ /^s$|^r$|^d$|^i$|^p$|^a$|^I$|^m$|^M$|^l$|^o$|^u$|^tlrisk$|^turisk$|^tsrisk$|^rscale$|^ralpha$|^tldist$|^tudist$|^tsdist$|^dscale$|^dalpha$|^prt$|^average$|^imain$|^omain$|^restart$|^ry$/ ){ 
				
				&help( \$parameter, \%$href_value_of_scr_par ) 
			}
			
			else { print "\n\tPARAMETER-ERROR: Unallowed command or value: -", $command, " ",$parameter, "!\n\n"; &help( \'h', \%$href_value_of_scr_par ) }
		}
		
		elsif 	(	$command =~ /^P$/ 			 ){ &preface( \%$href_value_of_scr_par 	  ) }
		
		
		elsif 	(	$command =~ /^q$/i 			 ){ &print_endscreen( \%$href_value_of_scr_par ); exit	}
		
		
		else 	{ print "\n\tPARAMETER-ERROR: Unallowed command-line: '-", $cmdline, "'!\n\n"; &help( \'h', \%$href_value_of_scr_par ) }
	}
	##############################################
	
	
	##############################################
	#  MSA + SPD without|with RESTART
	##############################################
	
	#######################
	# Read MSA infile (with opt.path) if RESTART is undefined
	if		( ( $href_value_of_scr_par->{sys_parameter}{msafile}  ) && ( $href_value_of_scr_par->{sys_parameter}{restart} == 0 ) ){ 
		
		my	@infile_path	= 	split( /\/|\\/, $href_value_of_scr_par->{sys_parameter}{msafile} ) ;
		my	$msa_nopath		= 	$infile_path[-1] ;
		
		if ( $msa_nopath =~ /^\w+.\w+$/ ){
			
			$href_value_of_scr_par->{sys_parameter}{msa_nopath}	= $infile_path[-1] ;
			
			#######################
			# set spd-print to 1 (active)
			$href_value_of_scr_par->{sys_parameter}{write_spd}	= 1 ;
			#######################
			
			#######################
			# define main penguin folder by adding alignment name
			$href_value_of_scr_par->{sys_parameter}{o_main} = $href_value_of_scr_par->{sys_parameter}{o_main}."_".$href_value_of_scr_par->{sys_parameter}{msa_nopath} ;
			$href_value_of_scr_par->{sys_parameter}{o_main} =~	s/.phy$|.fas$// ;
			#######################
		}
		
		else{  print "\n\tPARAMETER-ERROR: alignment-filename includes non-alphanumeric signs!\n\n" ; &help( \'p', \%$href_value_of_scr_par ) }
	}
	#######################
	
	#######################
	# SPD-file handling if RESTART is defined
	if	( $href_value_of_scr_par->{sys_parameter}{restart} == 1 ){
		
		# delete possible alignment definition
		$href_value_of_scr_par->{sys_parameter}{msafile} = () ;
		$href_value_of_scr_par->{sys_parameter}{msa_nopath}	= () ;
		
		# if main (spd) infile and outfile pathdir are identic: no spd-file print 
		if		( $href_value_of_scr_par->{sys_parameter}{o_main} eq $href_value_of_scr_par->{sys_parameter}{i_main} ){
			
			$href_value_of_scr_par->{sys_parameter}{write_spd}	= 0 ;
		}
		
		# if main (spd) infile and outfile pathdir are not identic: spd-file print to new outfilepath
		elsif 	( $href_value_of_scr_par->{sys_parameter}{o_main} ne $href_value_of_scr_par->{sys_parameter}{i_main} ){
			
			$href_value_of_scr_par->{sys_parameter}{write_spd}	= 1 ;
		}
	}
	#######################
	
	
	##############################################
	#  CLANFILE
	##############################################
	
	#######################
	# Check if clade definition file has been defined
	if		( $href_value_of_scr_par->{sys_parameter}{clfile} ){
		
		my	@infile_path 	= split( /\/|\\/, $href_value_of_scr_par->{sys_parameter}{clfile} ) ;
		my	$clf_nopath		= $infile_path[-1] ; 
			@infile_path 	= () ;
			
		if ( $clf_nopath =~ /^\w+\.\w+$/ ){
				
			#######################
			# store infilename without path extension in {sys_parameter}{clfile}->{clf_nopath}
			$href_value_of_scr_par->{sys_parameter}{clf_nopath}	= $clf_nopath ;
			#######################
		}
		
		else{  print "\nPARAMETER-ERROR: clan-filename includes non-alphanumeric signs!\n\n" ; &help( \'p', \%$href_value_of_scr_par ) }
	}
	#######################
	
	#######################
	# if clanfile is missing go to help
	else	{ print "\nPARAMETER-ERROR: Missing clanfile!\n\n" ; &help( \'p', \%$href_value_of_scr_par ) }
	##############################################
	
	
	##############################################
	#  OUTGROUP		
	##############################################
		
	#######################
	# Check for defined outgroup clan
	unless ( $href_value_of_scr_par->{sys_parameter}{outgroup} =~ /\w+/ ){
		
		print "\nPARAMETER-ERROR: Missing or incorrectly defined outgroup-clan ", $href_value_of_scr_par->{sys_parameter}{outgroup}, "!\n\n" ; &help( \'o', \%$href_value_of_scr_par )
	}
	##############################################
	
	
	##############################################
	#  'RISK' FILTER SETTINGS                    #
	##############################################
	
	#######################
	# check the upper bound limit of risk filtering is not lower than the lower bound limit
	if	( $href_value_of_scr_par->{seq_filterRisk}{turisk} < 0 ){ 
		
		print "\nPARAMETER-ERROR: upper threshold-limit of 'risk' filter < 0!\n\n" ; &help( \'turisk', \%$href_value_of_scr_par )
	}
	
	if	( $href_value_of_scr_par->{seq_filterRisk}{tlrisk} < 0 ){ 
		
		print "\nPARAMETER-ERROR: lower threshold-limit of 'risk' filter < 0!\n\n" ; &help( \'tlrisk', \%$href_value_of_scr_par )
	}
	
	if	( $href_value_of_scr_par->{seq_filterRisk}{turisk} < $href_value_of_scr_par->{seq_filterRisk}{tlrisk} ){ 
		
		print "\nPARAMETER-ERROR: upper bound limit of 'risk' filter < lower bound limit!\n\n" ; &help( \'tlrisk', \%$href_value_of_scr_par ) 
	}
	
	if	( $href_value_of_scr_par->{seq_filterRisk}{rscale} < 0 ){ 
		
		print "\nPARAMETER-ERROR: threshold-scaling of 'risk' filter < 0!\n\n" ; &help( \'tsrisk', \%$href_value_of_scr_par )
	}
	
	if	( $href_value_of_scr_par->{seq_filterRisk}{ralpha} < 0 ){ 
		
		print "\nPARAMETER-ERROR: alpha-score of 'risk' filter optimization step < 0!\n\n" ; &help( \'tarisk', \%$href_value_of_scr_par )
	}
	##############################################
	
	
	##############################################
	#  'DIST' FILTER SETTINGS                    #
	##############################################
	
	#######################
	# check the upper bound limit of risk filtering is not lower than the lower bound limit
	if	( $href_value_of_scr_par->{seq_filterDist}{tudist} < 0 ){ 
		
		print "\nPARAMETER-ERROR: upper threshold-limit of 'qdist' filter < 0!\n\n" ; &help( \'tudist', \%$href_value_of_scr_par )
	}
	
	if	( $href_value_of_scr_par->{seq_filterDist}{tldist} < 0 ){ 
		
		print "\nPARAMETER-ERROR: lower threshold-limit of 'qdist' filter < 0!\n\n" ; &help( \'tldist', \%$href_value_of_scr_par )
	}
	
	if	( $href_value_of_scr_par->{seq_filterDist}{tudist} < $href_value_of_scr_par->{seq_filterDist}{tldist} ){ 
		
		print "\nPARAMETER-ERROR: upper bound limit of 'qdist' filter < lower bound limit!\n\n" ; &help( \'tldist', \%$href_value_of_scr_par ) 
	}
	
	if	( $href_value_of_scr_par->{seq_filterDist}{dscale} < 0 ){ 
		
		print "\nPARAMETER-ERROR: threshold-scaling of 'qdist' filter < 0!\n\n" ; &help( \'tsdist', \%$href_value_of_scr_par )
	}
	
	if	( $href_value_of_scr_par->{seq_filterDist}{dalpha} < 0 ){ 
		
		print "\nPARAMETER-ERROR: alpha-score of 'qdist' filter optimization step < 0!\n\n" ; &help( \'tadist', \%$href_value_of_scr_par )
	}
	##############################################
	
	######################################################################
}

sub clan4_extraction_from_filepaths{
	
	my	$aref_qcalc_filepaths			=	$_[0] ;	# list of quartet individual split-claculation filepaths 									IN: defined; OUT: unchanged
	my	$sref_main_resultfolder_code	=	$_[1] ;	# constant string of the main penguin resultfolder name										IN: defined; OUT: unchanged
	my	$href_hol_fp_qcal_of_4clan_code	=	$_[2] ;	# key: 4clan combination code; value: list of qcalc filepaths								IN: undefined; OUT: defined
	my	$href_value_of_scr_param		=	$_[3] ;	# key1: sealion user available parameter options; value: optional parameter setting			IN: defined; OUT: unchanged
	my	$href_code_of_clan				=	$_[4] ;	# key1: clanname; value assigned code														IN: undefined; OUT: defined
	my	$href_code_of_4clan				=	$_[5] ;	# key1: 4clan assigned clan names; value assigned code										IN: undefined; OUT: defined
	my	$href_subfolder_q				=	$_[6] ;	# key1: output-type (e.g. 'pdf'); value: pathdir											IN: defined; OUT: unchanged
	
	#print "::&clan4_extraction_from_filepaths\n" ;
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# &help
	######################################################################
	
	######################################################################
	## extract 4clan combination-code from filepath name                ##
	######################################################################
	for my $fp ( sort @$aref_qcalc_filepaths ){
		
		####################################
		# substitute all pathinfo in front of the 4clan combination code in each quartet calculation filepath
		# allows multiple paths in front of the main penguin resultfolder $sref_main_resultfolder_code
		( my $fp_red = $fp ) =~ s/^.*$$sref_main_resultfolder_code\/(SPD\/)?// ;	
		####################################
		
		####################################
		# substitute all pathinfo after the 4clan combination code in each quartet calculation filepath
		( my $clan4_code_complete = $fp_red ) =~ s/\/.*$// ;
		####################################
		
		####################################
		# for each of the four clans, take the first start letter as clan code
		# add integer value+1 if start letter appears multiple times (for different clans)
		# combine them to $clan4_code_reduced
		my $clan4_code_reduced	; 
		my @clans = split "_", $clan4_code_complete ;
		my %seen_reduced_code ;
		
		for my $clan ( @clans ){
			
			$href_code_of_clan->{$clan} = substr $clan, 0, 1 ;
			
			if ( $seen_reduced_code{$href_code_of_clan->{$clan}} ){
				
				print "\nSPDFILE-ERROR: Found identic start letter '".$href_code_of_clan->{$clan}."' for different clans in SPD file ", $fp ,"!\n\n" ;
				&help( \'p', \%$href_value_of_scr_param ) ;
			}
			
			$seen_reduced_code{$href_code_of_clan->{$clan}}++ ;
			
			$clan4_code_reduced .= $href_code_of_clan->{$clan} ;
		}
		####################################
		
		####################################
		# assign complete qcalc filepath to corresponding reduced 4clan code list
		push @{$href_hol_fp_qcal_of_4clan_code->{$clan4_code_reduced}},  $fp ;
		####################################
		
		####################################
		# assign original 4clan combi with corresponding code
		my	$cl4 = join ":", @clans ; 
		$href_code_of_4clan->{$cl4} = $clan4_code_reduced ;
		####################################
		
		####################################
		# test print
		if ( $href_value_of_scr_param->{sys_parameter}{screen_prt} == 1 ){
			
			print	"\nfp: ", $fp,"\nred fp: ", $fp_red 		,
					"\n4clan full: "	, $clan4_code_complete	,
					"\n4clan reduced: "	, $clan4_code_reduced	,
					"\n";
		}#exit;
		####################################
	}
	######################################################################
}

sub determining_4clan_conflicted_taxon_support{
	
	my	$href_value_of_scr_param								= $_[0]	 ;	# key1: sealion user available parameter options; value: optional parameter setting																 							IN: defined; OUT: unchanged
	my	$sref_cmtf												= $_[1]	 ;	# code of quartet-method (e.g. nap_risk1) and taxon-filter (strict_in, strict_out, strict_all)																				IN: defined; OUT: unchanged
	my	$sref_clan4												= $_[2]	 ;	# actual 4clan combi code																																					IN: defined; OUT: unchanged
	my	$href_algorithm_def_met									= $_[3]  ;	# key1: 'ftaxon'; value: list of defined taxon filter																														IN: defined; OUT: unchanged
	my	$href_dist_sup_best2nd_of_stype							= $_[4]	 ;	# key1: mean|median; value: distance best to 2nd best supported 4clan-tree																									IN: undefined; OUT: defined
	my	$href_clan_of_sequence									= $_[5]	 ;	# key1: taxon; value: assigned clan																																			IN: defined; OUT: unchanged 
	my	$href_hoh_seen_confltype_of_clan4_of_taxon				= $_[6]	 ;	# key1: combined filtertype; key2: clan4; key3: taxon; value: counter++																										IN: undefined; OUT: defined  
	my	$href_hoh_seen_confl_taxon_of_suptype					= $_[7]	 ;	# key1: taxon_filter (stric_out, strict_in, strict_all); key2: sequence-name; key3: type of support (mean|median); key4: 'counter'; value: N occurence						IN: undefined; OUT: defined  
	my	$href_prt_path											= $_[8]  ;	# key1: R grafic unique printcode (e.g. 'metafileQ1'); key2: *.tsv filename; value: *tsv filename data addings 																IN: defined; OUT: unchanged
	my	$href_hoh_dat_of_tsvfile_of_prtcode						= $_[9]  ;	# key: output-type (e.g. 'pdf'); value: pathdir																																IN: defined; OUT: unchanged
	my	$href_hoh_linedata_of_latex_of_prtcode					= $_[10] ;	# key1: Latex table unique printcode (e.g. 'latexC'); key2: *.tex filename; value: @list of assigned table lines															IN: defined; OUT: changed
	my	$href_seen_tax_of_4clan									= $_[11] ;	# key1: 4clan code; key2: remaining taxon after OTU-filtering; value: counter number																						IN: defined; OUT: unchanged
	my	$href_hoh_seen_confl_taxon_of_type						= $_[12] ;	# key1: taxon_filter (stric_out, strict_in, strict_all); key2: type of support (mean|median); key3: sequence-name; 															IN: undefined; OUT: defined  
	my	$href_subfolder_q										= $_[13] ;	# hashpath to output-type specific subdirs of single quartet-pattern analyses; key: output-type (e.g. 'pdf'); value: pathdir												IN: defined; OUT: unchanged
	my	$href_textfilename_of_prtcode							= $_[14] ;	# key1: textfile unique printcode (e.g. 'textA'); value: filename for print																									IN: defined; OUT: unchanged
	
	#print "::&determining_4clan_conflicted_taxon_support\n";
	
	######################################################################
	# called from subroutine(s)...
	# &evaluation_otu_support_main
	#
	############
	# using subroutine(s)...
	# &print_tsv_T4
	# &latex_data_sample_LT1
	# &latex_data_sample_LT3
	######################################################################
	
	######################################################################
	# finding taxa with mean|median supp dist below defined threshold    #
	######################################################################
	for my $taxon ( sort keys %{$href_seen_tax_of_4clan->{$$sref_clan4}} ){
		
		my	$threshold_conflict = $href_value_of_scr_param->{tax_filterDist}{tdist1st2nd} ;
		my	%seen_stype_of_taxon ;
			
		for	my $taxon_filter ( @{$href_algorithm_def_met->{ftaxon}} ){
			
			my	$new_ftype = $$sref_cmtf."_".$taxon_filter ;
				
			for my $type_support ( sort @{$href_algorithm_def_met->{tsuppo}} ){
				
				my	$score_dist = $href_dist_sup_best2nd_of_stype->{$taxon}{$type_support} ;
				
				my	$specified_clan_conflict_handling;
				if 	(	(	(	$taxon_filter	=~ /strict_out|strict_all/i	) 	&&	( $href_clan_of_sequence->{$taxon} eq $href_value_of_scr_param->{sys_parameter}{outgroup}	) && ( $score_dist	<= $threshold_conflict	)	)	||
					(	(		$taxon_filter	=~ /strict_in|strict_all/i	) 	&&	( $href_clan_of_sequence->{$taxon} ne $href_value_of_scr_param->{sys_parameter}{outgroup}	) && ( $score_dist	<= $threshold_conflict	)	)	){
					
					####################################
					# (a) sample taxon as conflicted for given 4-clan combination and taxon filter
					$href_hoh_seen_confltype_of_clan4_of_taxon->{$new_ftype}{$$sref_clan4}{$type_support}{$taxon}++ ;
					$href_hoh_seen_confl_taxon_of_type->{$new_ftype}{$type_support}{$taxon}++ ;
					####################################
					
					####################################
					# (b) mark taxon as 'highly conflicted' if mean|median support difference of best and second best quartet tree are smaller as specified threshold ($threshold_conflict)
					# (c) store in %hoh_seen_conflicted_taxon_of_clan4_of_ctype_of_stype all neccessary info about conflicted taxon as hashkeys
					my		  $tax_status			=	'reject' ;
							  $href_hoh_seen_confl_taxon_of_suptype->{$new_ftype}{$taxon}{$type_support}{counter}++ ;
					push	@{$href_hoh_seen_confl_taxon_of_suptype->{$new_ftype}{$taxon}{$type_support}{clan4}}, $$sref_clan4  ;
					####################################
					
					############################
					# test print
					if ( $href_value_of_scr_param->{sys_parameter}{screen_prt} == 1 ){
						
						print "\n\tfailed SD:\t", $$sref_clan4, "\t", $new_ftype, "\t", $type_support, "\tdBest:2nd:\t", $score_dist, "\t", $taxon; 
					}
					############################
					
					#################################### LATEX LT1
					# sample 'highly conflicted' taxon with additional info as latex info line ($latex_info_line) in %hoh_linedata_of_latexfilename_of_printcode
					&latex_data_sample_LT1(
						
						\%$href_hoh_linedata_of_latex_of_prtcode	,	# key1: Latex table unique printcode (e.g. 'latexC'); key2: *.tex filename; value: @list of assigned table lines								IN: defined; OUT: changed
						\$$sref_cmtf								,	# code of quartet-method (e.g. nap_risk1) and taxon-filter (strict_in, strict_out, strict_all)													IN: defined; OUT: unchanged
						\$taxon_filter								,	# actual taxonfilter code																														IN: defined; OUT: unchanged
						\$tax_status								,	# 'reject' string																																IN: defined; OUT: unchanged
						\$$sref_clan4								,	# actual 4clan combi code																														IN: defined; OUT: unchanged
						\$taxon										,	# sequence name																																	IN: defined; OUT: unchanged
						\$href_clan_of_sequence->{$taxon}			,	# sequence assigned clan																														IN: defined; OUT: unchanged	
						\$type_support								,	# type of support calculation (nap, napv, napr)																									IN: defined; OUT: unchanged
						\$score_dist								,	# distance best OTU tree and 2nd best OTU tree																									IN: defined; OUT: unchanged
						\%$href_subfolder_q							,	# hashpath to output-type specific subdirs of single quartet-pattern analyses; key: output-type (e.g. 'pdf'); value: pathdir					IN: defined; OUT: unchanged
						\%$href_textfilename_of_prtcode				,	# key1: textfile unique printcode (e.g. 'textA'); value: filename for print																		IN: defined; OUT: unchanged
	
					) ;
					####################################
					
					unless ( $seen_stype_of_taxon{$taxon}{$type_support} ){
						
						######################################################################
						# write data to tsv. table T4                                        #
						######################################################################
						#print	"::printing data to .tsv metatable T4\n"		;
						
						&print_tsv_T4(
							
							\%$href_prt_path								,	# key1: R grafic unique printcode (e.g. 'metafileQ1'); key2: *.tsv filename; value: *tsv filename data addings 							IN: defined; OUT: unchanged
							\%$href_hoh_dat_of_tsvfile_of_prtcode			,	# key: output-type (e.g. 'pdf'); value: pathdir																							IN: defined; OUT: unchanged
							\$$sref_cmtf									,	# type of calculation (nap, napv...) and or filter (rsik_1, risk1_dist1...) method														IN: defined; OUT: unchanged 
							\$$sref_clan4									,	# clan4 code																															IN: defined; OUT: unchanged 
							\$taxon											,	# taxon name																															IN: defined; OUT: unchanged
							\$href_clan_of_sequence->{$taxon}				,	# clan name of given taxon																												IN: defined; OUT: unchanged
							\$score_dist									,	# distance best OTU tree and 2nd best OTU tree																							IN: defined; OUT: unchanged
							\$type_support									,	# support type (mean || median)																											IN: defined; OUT: unchanged
						) ;
						######################################################################
					}
					
					$seen_stype_of_taxon{$taxon}{$type_support}++
				}
				
				else{
					
					my	$tax_status	= 'remain' ;
					
					#################################### LATEX LT3
					# sample 'unconflicted' taxon with additional info as latex info line ($latex_info_line) in %hoh_linedata_of_latexfilename_of_printcode
					&latex_data_sample_LT3(
						
						\%$href_hoh_linedata_of_latex_of_prtcode	,	# key1: Latex table unique printcode (e.g. 'latexC'); key2: *.tex filename; value: @list of assigned table lines								IN: defined; OUT: changed
						\$$sref_cmtf								,	# code of quartet-method (e.g. nap_risk1) and taxon-filter (strict_in, strict_out, strict_all)													IN: defined; OUT: unchanged
						\$taxon_filter								,	# actual taxonfilter code																														IN: defined; OUT: unchanged
						\$tax_status								,	# 'reject' string																																IN: defined; OUT: unchanged
						\$$sref_clan4								,	# actual 4clan combi code																														IN: defined; OUT: unchanged
						\$taxon										,	# sequence name																																	IN: defined; OUT: unchanged
						\$href_clan_of_sequence->{$taxon}			,	# sequence assigned clan																														IN: defined; OUT: unchanged	
						\$type_support								,	# type of support calculation (nap, napv, napr)																									IN: defined; OUT: unchanged
						\$score_dist								,	# distance best OTU tree and 2nd best OTU tree																									IN: defined; OUT: unchanged
						\%$href_subfolder_q							,	# hashpath to output-type specific subdirs of single quartet-pattern analyses; key: output-type (e.g. 'pdf'); value: pathdir					IN: defined; OUT: unchanged
						\%$href_textfilename_of_prtcode				,	# key1: textfile unique printcode (e.g. 'textA'); value: filename for print																		IN: defined; OUT: unchanged
					) ;
					####################################
				}
			}
		}
	} #exit;
	######################################################################
	######################################################################
}

sub determining_4clan_taxon_treesupport_distances{
	
	my	$href_value_of_scr_param								= $_[0]	;	# key1: sealion user available parameter options; value: optional parameter setting																IN: defined; OUT: unchanged
	my	$sref_clan4												= $_[1]	;	# actual 4clan combi code																														IN: defined; OUT: unchanged
	my	$href_seen_tax_of_4clan									= $_[2]	;	# key1: 4clan code; key2: remaining taxon after OTU-filtering; value: counter number															IN: defined; OUT: unchanged
	my	$href_supp_of_4clan_of_tax_of_type_of_tree				= $_[3]	;	# key1: 4clan code; key2: 4taxon-string; key3: mean|median; key4: quartet-tree; value: support value 											IN: defined; OUT: unchanged
	my	$href_dist_sup_best2nd_of_stype							= $_[4]	;	# key1: taxon; key2: mean|median; value: distance best to 2nd best supported 4clan-tree															IN: undefined; OUT: defined
	my	$href_dist_sup_best3rd_of_stype							= $_[5]	;	# key1: taxon; key2: mean|median; value: distance best to 3rd best supported 4clan-tree															IN: undefined; OUT: defined
	
	######################################################################
	# called from subroutine(s)...
	# &evaluation_otu_support_main
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	# calculation of taxon specific tree-support differences             #
	######################################################################
	for my $taxon ( sort keys %{$href_seen_tax_of_4clan->{$$sref_clan4}} ){
		
		my	@qtree_sorted_by_mean_support						=	(sort { $href_supp_of_4clan_of_tax_of_type_of_tree->{$$sref_clan4}{$taxon}{mean}{$b} <=> $href_supp_of_4clan_of_tax_of_type_of_tree->{$$sref_clan4}{$taxon}{mean}{$a} } keys %{$href_supp_of_4clan_of_tax_of_type_of_tree->{$$sref_clan4}{$taxon}{mean}} ) ; 
			$href_dist_sup_best2nd_of_stype->{$taxon}{mean}		=	$href_supp_of_4clan_of_tax_of_type_of_tree->{$$sref_clan4}{$taxon}{mean}{$qtree_sorted_by_mean_support[0]} - $href_supp_of_4clan_of_tax_of_type_of_tree->{$$sref_clan4}{$taxon}{mean}{$qtree_sorted_by_mean_support[1]} ;
			$href_dist_sup_best3rd_of_stype->{$taxon}{mean}		=	$href_supp_of_4clan_of_tax_of_type_of_tree->{$$sref_clan4}{$taxon}{mean}{$qtree_sorted_by_mean_support[0]} - $href_supp_of_4clan_of_tax_of_type_of_tree->{$$sref_clan4}{$taxon}{mean}{$qtree_sorted_by_mean_support[2]} ;
			
		my	@qtree_sorted_by_median_support						=	(sort { $href_supp_of_4clan_of_tax_of_type_of_tree->{$$sref_clan4}{$taxon}{median}{$b} <=> $href_supp_of_4clan_of_tax_of_type_of_tree->{$$sref_clan4}{$taxon}{median}{$a} } keys %{$href_supp_of_4clan_of_tax_of_type_of_tree->{$$sref_clan4}{$taxon}{median}} ) ; 
			$href_dist_sup_best2nd_of_stype->{$taxon}{median}	=	$href_supp_of_4clan_of_tax_of_type_of_tree->{$$sref_clan4}{$taxon}{median}{$qtree_sorted_by_median_support[0]}	- $href_supp_of_4clan_of_tax_of_type_of_tree->{$$sref_clan4}{$taxon}{median}{$qtree_sorted_by_median_support[1]} ;
			$href_dist_sup_best3rd_of_stype->{$taxon}{median}	=	$href_supp_of_4clan_of_tax_of_type_of_tree->{$$sref_clan4}{$taxon}{median}{$qtree_sorted_by_median_support[0]}	- $href_supp_of_4clan_of_tax_of_type_of_tree->{$$sref_clan4}{$taxon}{median}{$qtree_sorted_by_median_support[2]} ;
		
		####################################
		## test print
		if ( $href_value_of_scr_param->{sys_parameter}{screen_prt} == 1 ){
			
			print 	"::&determining_4clan_taxon_treesupport_distances\n",
					"\n\t", $taxon, "\t1st Best mean Tree:\t"		, $qtree_sorted_by_mean_support[0]					,
					"\n\t", $taxon, "\t2nd Best mean Tree:\t"		, $qtree_sorted_by_mean_support[1]					,
					"\n\t", $taxon, "\t3rd Best mean Tree:\t"		, $qtree_sorted_by_mean_support[2]					,
					"\n\t", $taxon, "\tDis Best mean to 2nd\t"		, $href_dist_sup_best2nd_of_stype->{$taxon}{mean}   ,
					"\n\t", $taxon, "\tDis Best mean to 3rd\t"		, $href_dist_sup_best3rd_of_stype->{$taxon}{mean}	,
					"\n"  ,
					"\n\t", $taxon, "\t1st Best median Tree:\t"		, $qtree_sorted_by_median_support[0]				,
					"\n\t", $taxon, "\t2nd Best median Tree:\t"		, $qtree_sorted_by_median_support[1]				,
					"\n\t", $taxon, "\t3rd Best median Tree:\t"		, $qtree_sorted_by_median_support[2]				,
					"\n\t", $taxon, "\tDis Best median to 2nd\t"	, $href_dist_sup_best2nd_of_stype->{$taxon}{median}	,
					"\n\t", $taxon, "\tDis Best median to 3rd\t"	, $href_dist_sup_best3rd_of_stype->{$taxon}{median} ,
					"\n"  ;
		} #exit;
		####################################
		
		####################################
		## control check
		if ( ($href_dist_sup_best2nd_of_stype->{$taxon}{mean}   < 0  )	||
			 ($href_dist_sup_best3rd_of_stype->{$taxon}{mean}   < 0  )	||
			 ($href_dist_sup_best2nd_of_stype->{$taxon}{median} < 0  )	||
			 ($href_dist_sup_best3rd_of_stype->{$taxon}{median} < 0  )	){
		
			die "\nBUG-ERROR: Negative support-distance score between best and lower supported trees for OTU ", $taxon, " in &eval_otu_4clan_support!\nPlease, report BUG to software developer!\n\n"
		}
		####################################
	}
	######################################################################
	######################################################################
}

sub determining_mean_median_4clan_taxon_support{
	
	my	$href_value_of_scr_param								= $_[0]	;	# key1: sealion user available parameter options; value: optional parameter setting																IN: defined; OUT: unchanged
	my	$sref_clan4												= $_[1]	;	# actual 4clan combi code																														IN: defined; OUT: unchanged
	my	$sref_filttype											= $_[2]	;	# type of calculation (nap, napv...) and or filter (rsik_1, risk1_dist1...) method																IN: defined; OUT: unchanged 
	my	$href_hol_sup_of_4clan_of_tax_of_tree					= $_[3]	;	# key1: 4clan code; key2: 4taxon-string; key3: quartet-tree; value: list of taxon																IN: defined; OUT: unchanged
	my	$href_seen_tax_of_4clan									= $_[4]	;	# key1: 4clan code; key2: remaining taxon after OTU-filtering; value: counter number															IN: defined; OUT: unchanged
	my	$href_supp_of_4clan_of_tax_of_type_of_tree				= $_[5]	;	# key1: 4clan code; key2: 4taxon-string; key3: mean|median; key4: quartet-tree; value: support value 											IN: undefined; OUT: defined
	
	######################################################################
	# called from subroutine(s)...
	# &evaluation_otu_support_main
	#
	############
	# using subroutine(s)...
	# &calc_mean
	# &calc_median
	######################################################################
	
	######################################################################
	# calculation of taxon specific median (normalised) tree-support     #
	######################################################################
	for my $taxon ( sort keys %{$href_seen_tax_of_4clan->{$$sref_clan4}} ){
		
		my  ( @trees, @median, @mean ) ;
		
		for my $tree ( sort keys %{$href_hol_sup_of_4clan_of_tax_of_tree->{$$sref_clan4}{$taxon}} ){
			
			my $mean	= &calc_mean	( \@{$href_hol_sup_of_4clan_of_tax_of_tree->{$$sref_clan4}{$taxon}{$tree}} ) ; 
			my $median	= &calc_median	( \@{$href_hol_sup_of_4clan_of_tax_of_tree->{$$sref_clan4}{$taxon}{$tree}} ) ;
			
			push @mean	, $mean 	;
			push @median, $median	;
			push @trees , $tree		;
		}
		
		my @normed_mean		= &norm_values ( \@mean   ) ;
		my @normed_median	= &norm_values ( \@median ) ;
		
		$href_supp_of_4clan_of_tax_of_type_of_tree->{$$sref_clan4}{$taxon}{mean}{$trees[0]}		= $normed_mean[0] ;
		$href_supp_of_4clan_of_tax_of_type_of_tree->{$$sref_clan4}{$taxon}{mean}{$trees[1]}		= $normed_mean[1] ;
		$href_supp_of_4clan_of_tax_of_type_of_tree->{$$sref_clan4}{$taxon}{mean}{$trees[2]}		= $normed_mean[2] ;
		
		$href_supp_of_4clan_of_tax_of_type_of_tree->{$$sref_clan4}{$taxon}{median}{$trees[0]}	= $normed_median[0] ;
		$href_supp_of_4clan_of_tax_of_type_of_tree->{$$sref_clan4}{$taxon}{median}{$trees[1]}	= $normed_median[1] ;
		$href_supp_of_4clan_of_tax_of_type_of_tree->{$$sref_clan4}{$taxon}{median}{$trees[2]}	= $normed_median[2] ;
		
		####################################
		## test print
		if ( $href_value_of_scr_param->{sys_parameter}{screen_prt} == 1 ){
			
			print 	"::&determining_mean_median_4clan_taxon_support\n",
					"\n\t", $$sref_filttype, "\t", $$sref_clan4, "\t", $taxon, "\t", $trees[0], "\tmean\t", $href_supp_of_4clan_of_tax_of_type_of_tree->{$$sref_clan4}{$taxon}{mean}{$trees[0]} ,
					"\n\t", $$sref_filttype, "\t", $$sref_clan4, "\t", $taxon, "\t", $trees[1], "\tmean\t", $href_supp_of_4clan_of_tax_of_type_of_tree->{$$sref_clan4}{$taxon}{mean}{$trees[1]} ,
					"\n\t", $$sref_filttype, "\t", $$sref_clan4, "\t", $taxon, "\t", $trees[2], "\tmean\t", $href_supp_of_4clan_of_tax_of_type_of_tree->{$$sref_clan4}{$taxon}{mean}{$trees[2]} ,
					"\n"  ,
					"\n\t", $$sref_filttype, "\t", $$sref_clan4, "\t", $taxon, "\t", $trees[0], "\tmedian\t", $href_supp_of_4clan_of_tax_of_type_of_tree->{$$sref_clan4}{$taxon}{median}{$trees[0]} ,
					"\n\t", $$sref_filttype, "\t", $$sref_clan4, "\t", $taxon, "\t", $trees[1], "\tmedian\t", $href_supp_of_4clan_of_tax_of_type_of_tree->{$$sref_clan4}{$taxon}{median}{$trees[1]} ,
					"\n\t", $$sref_filttype, "\t", $$sref_clan4, "\t", $taxon, "\t", $trees[2], "\tmedian\t", $href_supp_of_4clan_of_tax_of_type_of_tree->{$$sref_clan4}{$taxon}{median}{$trees[2]} ,
					"\n"  ;	
		}
		####################################
	}#exit;
	######################################################################
	######################################################################
}

sub evaluation_otu_support_main{
	
	my	$href_printpath											= $_[0] ;	# key1: R grafic unique printcode (e.g. 'metafileQ1'); key2: *.tsv filename; value: *tsv filename data addings 																IN: defined; OUT: unchanged
	my	$href_hoh_data_of_tsvfile_of_prcode						= $_[1] ;	# key: output-type (e.g. 'pdf'); value: pathdir																																IN: defined; OUT: unchanged
	my	$href_hoh_sup_of_clan4_of_calc_of_qtree_of_fc			= $_[2]	;	# multiple hash, storing the final support (nap or napv) related to each rooted quartet-tree																				IN: defined; OUT: unchanged
	my	$href_hoh_seen_rejected_of_4clan_of_filtertype_of_fn	= $_[3]	;	# key1: 4clan code; key2: type of supportfilter ('nap_risk1', 'napr_risk1', or 'napv_risk1'); key3: actual sealion quartet number; value: filter corresponding score-value	IN: defined; OUT: unchanged 
	my	$href_taxa4_of_fc										= $_[4]	;	# key1: quartet-number ($fc); value : quartet assigned taxa (alphabetically ordered, joined by ':'																			IN: defined; OUT: unchanged 
	my	$sref_ftype												= $_[5]	;	# type of calculation (nap, napv...) and or filter (rsik_1, risk1_dist1...) method																							IN: defined; OUT: unchanged 
	my	$href_clan_of_seq										= $_[6]	;	# key1: taxon; value: assigned clan																																			IN: defined; OUT: unchanged 
	my	$href_hol_sup_of_type_of_mfilter						= $_[7]	;	# key1: type of supportfilter ('nap_risk1', 'napr_risk1', or 'napv_risk1' ...); key2: mean|median; value : list of 'tree\tmean|medianScore'									IN: defined; OUT: changed
	my	$href_algorithm_def_methods								= $_[8] ;	# key1: 'ftaxon'; value: list of defined taxon filter																														IN: defined; OUT: unchanged
	my	$href_value_of_scr_par									= $_[9] ;	# key1: sealion user available parameter options; value: optional parameter setting																 							IN: defined; OUT: unchanged
	my	$href_hoh_linedata_of_latexfilename_of_printcode		= $_[10] ;	# key1: Latex table unique printcode (e.g. 'latexC'); key2: *.tex filename; value: @list of assigned table lines															IN: defined; OUT: changed
	my	$href_seen_remain_of_cmfm_of_fc							= $_[11] ;	# key1: calculation method + filter method (nap_risk1, napv_dist2, napr_dist3); key2: quartet-number ($fc); value: seen +1													IN: defined; OUT: changed
	my	$href_seen_filter_method								= $_[12] ;	# key1: filter method incl OTU clipping approaches; value: seen +1																											IN: defined; OUT: changed
	my	$href_textfilename_of_printcode							= $_[13] ;	# key1: textfile unique printcode (e.g. 'textA'); value: filename for print																									IN: defined; OUT: unchanged
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# &sampling_single_4clan_taxon_support
	# &determining_mean_median_4clan_taxon_support
	# &determining_4clan_taxon_treesupport_distances
	# &determining_4clan_conflicted_taxon_support
	# &print_tsv_T1
	# &print_tsv_T2
	# &print_tsv_T3
	# &latex_data_sample_LT2
	# &print_taxon_reduced_clanfile
	# &taxon_clipping
	######################################################################
	
	
	######################################################################
	# global variables within subroutine
	my	(
			%hoh_seen_confltype_of_clan4_of_taxon	,	# key1: combined filtertype; key2: clan4; key3: taxon; value: counter++
			%hoh_seen_confl_taxon_of_suptype		,	# key1: taxon_filter (stric_out, strict_in, strict_all); key2: sequence-name; key3: type of support (mean|median); key4: 'counter'; value: N occurence
			%hoh_seen_confl_taxon_of_type			,	# key1: taxon_filter (stric_out, strict_in, strict_all); key2: type of support (mean|median); key3: sequence-name; 
			%seen_clan_of_taxon						,	# key1: sequence-name; value: associated clan of sequence-name
			%seen_qn_of_4clan						,	# key1: 4clan code; key2: quartet-number; value: list of 4clan assigned quartet number
			%seen_tax_of_4clan						,	# key1: 4clan code; key2: remaining taxon after OTU-filtering; value: counter number
			%seen_alltax_of_4clan					,	# key1: 4clan code; key2: each 4clan assigned taxon; value: counter number
			%hol_nwk_and_support_of_clan4_of_taxa4	,	# key1: 4clan code; key2: 4taxon-string; value: list of 'tree::supp:' of the three 4clan relationships
			%hol_sup_of_4clan_of_tax_of_tree		,	# key1: 4clan code; key2: 4taxon-string; key3: quartet-tree; value: list of taxon
			%supp_of_4clan_of_tax_of_type_of_tree	,	# key1: 4clan code; key2: 4taxon-string; key3: mean|median; key4: quartet-tree; value: support value
	) ;
	######################################################################
	
	
	######################################################################
	# extract original quartet calculation method ($cm) from possible 
	# quartet-filter method extension in $$sref_ftype
	my	@cm_fm		= split "_", $$sref_ftype ;
	my	$cm			= shift @cm_fm ;
	######################################################################
	
	
	######################################################################
	# for each 4clan combination...
	# ...sampling 4clan-tree support for each taxon from related quartets
	# ...calc. mean|median 4clan-tree (unfilt) support of each taxon
	# ...calc. mean|median 4clan-tree (unfilt) support diff. of each taxon
	# ...identification of pot. support confl. taxa (dist supp < threshold)
	# ...write data to tsv. table T1
	# ...write data to tsv. table T2
	# ...write data to tsv. table T3
	print 	"\n\t----------------------------------\n"	,
			"\tanalysing single OTU support for 4clan trees of...\n\n" ;

	for	my $clan4 ( sort keys %{$href_hoh_sup_of_clan4_of_calc_of_qtree_of_fc} ){
		
		######################################################################
		# sampling 4clan-tree support for each taxon from related quartets ###
		######################################################################
		
		print	"\t\t...", $clan4 , " (".$$sref_ftype.")\n" ;
		
		####################################
		# test for (unallowed) multiple 4-taxon combinations sharing the same set of taxa
		my	%seen_taxa ;
		for my $no ( sort {$a<=>$b} keys %$href_taxa4_of_fc ){ push @{$seen_taxa{$href_taxa4_of_fc->{$no}}}, $no  }
		for my $no ( sort {$a<=>$b} keys %$href_taxa4_of_fc ){ 
			
			if ( @{$seen_taxa{$href_taxa4_of_fc->{$no}}} > 1 ){
				
				my $l = join "_", @{$seen_taxa{$href_taxa4_of_fc->{$no}}} ;
				die "\nQUARTET-ERROR: sequence-quartet found multiple times:", $href_taxa4_of_fc->{$no}, "\tquartet ", $l, "\n\n"
			}
		}
		%seen_taxa = () ;
		####################################
		
		####################################
		# start sampling taxon support
		&sampling_single_4clan_taxon_support(
			
			\%$href_value_of_scr_par								,	# key1: sealion user available parameter options; value: optional parameter setting										IN: defined; OUT: unchanged
			\%hol_nwk_and_support_of_clan4_of_taxa4					,	# key1: 4clan code; key2: 4taxon-string; value: list of 'tree::supp:' of the three 4clan relationships					IN: undefined; OUT: defined
			\$cm													,	# original calc. method (nap, napr, napv)																				IN: defined; OUT: unchanged
			\$clan4													,	# actual 4clan combi code																								IN: defined; OUT: unchanged
			\%$href_hoh_sup_of_clan4_of_calc_of_qtree_of_fc			,	# multiple hash, storing the final support (nap or napv) related to each rooted quartet-tree							IN: defined; OUT: unchanged
			\%$href_hoh_seen_rejected_of_4clan_of_filtertype_of_fn	,	# key1: 4clan code; key2: type of supportfilter (e.g. 'nap_risk1'); key3: sealion quartet number; value: score			IN: defined; OUT: unchanged 
			\%$href_taxa4_of_fc										,	# key1: quartet-number ($fc); value : quartet assigned taxa (alphabetically ordered, joined by ':'						IN: defined; OUT: unchanged 
			\$$sref_ftype											,	# type of calculation (nap, napv...) and or filter (rsik_1, risk1_dist1...) method										IN: defined; OUT: unchanged  
			\%hol_sup_of_4clan_of_tax_of_tree						,	# key1: 4clan code; key2: 4taxon-string; key3: quartet-tree; value: list of taxon										IN: undefined; OUT: defined
			\%seen_qn_of_4clan										,	# key1: 4clan code; key2: quartet-number; value: list of 4clan assigned quartet number									IN: undefined; OUT: defined
			\%seen_tax_of_4clan										,	# key1: 4clan code; key2: remaining taxon after OTU-filtering; value: counter number									IN: undefined; OUT: defined
			\%seen_alltax_of_4clan									,	# key1: 4clan code; key2: each 4clan assigned taxon; value: counter number												IN: undefined; OUT: defined
		) ;
		######################################################################
		
		
		######################################################################
		# calc. mean|median 4clan-tree (unfilt) support of each taxon      ###
		######################################################################
		&determining_mean_median_4clan_taxon_support(
			
			\%$href_value_of_scr_par								,	# key1: sealion user available parameter options; value: optional parameter setting										IN: defined; OUT: unchanged
			\$clan4													,	# actual 4clan combi code																								IN: defined; OUT: unchanged
			\$$sref_ftype											,	# type of calculation (nap, napv...) and or filter (rsik_1, risk1_dist1...) method										IN: defined; OUT: unchanged  
			\%hol_sup_of_4clan_of_tax_of_tree						,	# key1: 4clan code; key2: 4taxon-string; key3: quartet-tree; value: list of taxon										IN: defined; OUT: unchanged
			\%seen_tax_of_4clan										,	# key1: 4clan code; key2: remaining taxon after OTU-filtering; value: counter number									IN: defined; OUT: unchanged
			\%supp_of_4clan_of_tax_of_type_of_tree					,	# key1: 4clan code; key2: 4taxon-string; key3: mean|median; key4: quartet-tree; value: support value 					IN: undefined; OUT: defined
		);
		######################################################################
		

		######################################################################
		# calc. mean|median 4clan-tree (unfilt) support diff. of each taxon ##
		######################################################################
		my	(
				
			%dist_sup_best2nd_of_stype	,	# key1: taxon; key2: mean|median; value: distance best to 2nd best supported 4clan-tree
			%dist_sup_best3rd_of_stype	,	# key1: taxon; key2: mean|median; value: distance best to 3rd best supported 4clan-tree
		) ;

		&determining_4clan_taxon_treesupport_distances(
			
			\%$href_value_of_scr_par								,	# key1: sealion user available parameter options; value: optional parameter setting										IN: defined; OUT: unchanged
			\$clan4													,	# actual 4clan combi code																								IN: defined; OUT: unchanged
			\%seen_tax_of_4clan										,	# key1: 4clan code; key2: remaining taxon after OTU-filtering; value: counter number									IN: defined; OUT: unchanged
			\%supp_of_4clan_of_tax_of_type_of_tree					,	# key1: 4clan code; key2: 4taxon-string; key3: mean|median; key4: quartet-tree; value: support value 					IN: defined; OUT: unchanged
			\%dist_sup_best2nd_of_stype								,	# key1: taxon; key2: mean|median; value: distance best to 2nd best supported 4clan-tree									IN: undefined; OUT: defined
			\%dist_sup_best3rd_of_stype								,	# key1: taxon; key2: mean|median; value: distance best to 3rd best supported 4clan-tree									IN: undefined; OUT: defined
		) ;
		######################################################################
		
		
		######################################################################
		# identification of pot. support confl. taxa (dist supp < threshold) #
		######################################################################
		&determining_4clan_conflicted_taxon_support(
			
			\%$href_value_of_scr_par								,	# key1: sealion user available parameter options; value: optional parameter setting																 							IN: defined; OUT: unchanged
			\$$sref_ftype											,	# type of calculation (nap, napv...) and or filter (rsik_1, risk1_dist1...) method																							IN: defined; OUT: unchanged  
			\$clan4													,	# actual 4clan combi code																																					IN: defined; OUT: unchanged
			\%$href_algorithm_def_methods							,	# key1: 'ftaxon'; value: list of defined taxon filter																														IN: defined; OUT: unchanged
			\%dist_sup_best2nd_of_stype								,	# key1: taxon; key2: mean|median; value: distance best to 2nd best supported 4clan-tree																						IN: defined; OUT: unchanged
			\%$href_clan_of_seq										,	# key1: taxon; value: assigned clan																																			IN: defined; OUT: unchanged 
			\%hoh_seen_confltype_of_clan4_of_taxon					,	# key1: combined filtertype; key2: clan4; key3: taxon; value: counter++																										IN: undefined; OUT: defined  
			\%hoh_seen_confl_taxon_of_suptype						,	# key1: taxon_filter (stric_out, strict_in, strict_all); key2: sequence-name; key3: type of support (mean|median); key4: 'counter'; value: N occurence						IN: undefined; OUT: defined
			\%$href_printpath										,	# key1: R grafic unique printcode (e.g. 'metafileQ1'); key2: *.tsv filename; value: *tsv filename data addings 																IN: defined; OUT: unchanged
			\%$href_hoh_data_of_tsvfile_of_prcode					,	# key: output-type (e.g. 'pdf'); value: pathdir																																IN: defined; OUT: unchanged
			\%$href_hoh_linedata_of_latexfilename_of_printcode		,	# key1: Latex table unique printcode (e.g. 'latexC'); key2: *.tex filename; value: @list of assigned table lines															IN: defined; OUT: changed
			\%seen_tax_of_4clan										,	# key1: 4clan code; key2: remaining taxon after OTU-filtering; value: counter number																						IN: defined; OUT: unchanged
			\%hoh_seen_confl_taxon_of_type							,	# key1: taxon_filter (stric_out, strict_in, strict_all); key2: type of support (mean|median); key3: sequence-name; 															IN: undefined; OUT: defined  
			\%$href_printpath										,	# key1: R grafic unique printcode (e.g. 'metafileQ1'); key2: *.tsv filename; value: *tsv filename data addings 																IN: defined; OUT: unchanged
			\%$href_textfilename_of_printcode						,	# key1: textfile unique printcode (e.g. 'textA'); value: filename for print																									IN: defined; OUT: unchanged
		) ;
		######################################################################
		
		
		######################################################################
		# write data to tsv. table T1                                        #
		######################################################################
		&print_tsv_T1(
			
			\%$href_printpath										,	# key1: R grafic unique printcode (e.g. 'metafileQ1'); key2: *.tsv filename; value: *tsv filename data addings 			IN: defined; OUT: unchanged
			\%$href_hoh_data_of_tsvfile_of_prcode					,	# key: output-type (e.g. 'pdf'); value: pathdir																			IN: defined; OUT: unchanged
			\%supp_of_4clan_of_tax_of_type_of_tree					,	# key1: 4clan code; key2: 4taxon-string; key3: mean|median; key4: quartet-tree; value: support value 					IN: defined; OUT: unchanged 
			\$$sref_ftype											,	# type of calculation (nap, napv...) and or filter (rsik_1, risk1_dist1...) method										IN: defined; OUT: unchanged 
			\$clan4													,	# clan4 code																											IN: defined; OUT: unchanged 
			\%$href_clan_of_seq										,	# key1: taxon; value: assigned clan																						IN: defined; OUT: unchanged 
			\%$href_algorithm_def_methods							,	# key1: 'ftaxon'; value: list of defined taxon filter																	IN: defined; OUT: unchanged
			\%seen_alltax_of_4clan									,	# key1: 4clan code; key2: each 4clan assigned taxon; value: counter number												IN: defined; OUT: unchanged 
		) ;
		######################################################################
		
		
		######################################################################
		# write data to tsv. table T2                                        #
		######################################################################
		&print_tsv_T2(
			
			\%$href_printpath										,	# key1: R grafic unique printcode (e.g. 'metafileQ1'); key2: *.tsv filename; value: *tsv filename data addings 			IN: defined; OUT: unchanged
			\%$href_hoh_data_of_tsvfile_of_prcode					,	# key: output-type (e.g. 'pdf'); value: pathdir																			IN: defined; OUT: unchanged
			\%supp_of_4clan_of_tax_of_type_of_tree					,	# key1: 4clan code; key2: 4taxon-string; key3: mean|median; key4: quartet-tree; value: support value 					IN: defined; OUT: unchanged
			\$$sref_ftype											,	# type of calculation (nap, napv...) and or filter (rsik_1, risk1_dist1...) method										IN: defined; OUT: unchanged 
			\$clan4													,	# clan4 code																											IN: defined; OUT: unchanged 
			\%$href_clan_of_seq										,	# key1: taxon; value: assigned clan																						IN: defined; OUT: unchanged 
			\%$href_algorithm_def_methods							,	# key1: 'ftaxon'; value: list of defined taxon filter																	IN: defined; OUT: unchanged
			\%seen_alltax_of_4clan									,	# key1: 4clan code; key2: each 4clan assigned taxon; value: counter number												IN: defined; OUT: unchanged 
		) ;
		######################################################################
	}
	######################################################################
}

sub evaluation_quartet_corrected_NaNpNc{
	
	my	$aref_rooted_qtrees					= $_[0] ; # list of rooted quartet trees 																								IN: defined; OUT: unchanged
	my	$href_cf_of_type_of_qt				= $_[1] ; # key1: info type (e.g. ('cf_obs'); key2: rooted quartet-tree; value: correction factor for site patterns of actual tree		IN:	defined; OUT: unchanged
	my	$href_NaNpNc_inf_of_type_of_tree	= $_[2] ; # key1: tree dependent info types (.e.g 'Na', 'Np'...); key2: rooted tree, value_ type value									IN: defined, OUT: unchanged
	my	$sref_total_Npattern				= $_[3] ; # total number of observed quartet site positions																				IN: defined, OUT: unchanged
	my	$href_value_of_scr_par				= $_[4] ; # key1: sealion user available parameter options; value: optional parameter setting											IN: defined; OUT: unchanged
	
	
	######################################################################
	# called from...
	# %evaluation_quartet_pattern_main
	#
	############
	# using 
	# &calc_percent
	######################################################################
	
	######################################################################
	## for each quartet-tree red. Na, Np, Nc due to CF; do red. rNc/Na  ##
	######################################################################
	
	###################################
	# for each of the three rooted trees (qt1, qt2, qt3)...
	# (1) calculate cf corrected number of reduced pattern support (Nb_obs_red, Na_obs_red, Np_obs_red, Nc_exp_red)
	# ...using observed calculated cf_obs value for Na_obs and Np_obs
	# ...using tree related cf_exp value for $mean_Nc_exp
	# (2) calculate proportion of reduced pattern support (pNa_obs_red, pNp_obs_red, pNc_exp_red)
	# (3) ratio rNc/Na (corrected)
	for my $tree ( @$aref_rooted_qtrees ){
		
		###################################
		# (1) calculate cf corrected number of reduced pattern support (Nb_obs_red, Na_obs_red, Np_obs_red, Nc_exp_red)
		$href_NaNpNc_inf_of_type_of_tree->{Na_nap_red_obs}{$tree}		=	$href_NaNpNc_inf_of_type_of_tree->{Na_obs}{$tree}		- ( $href_NaNpNc_inf_of_type_of_tree->{Na_obs}{$tree} * $href_cf_of_type_of_qt->{cf_obs}{$tree} ) ;
		$href_NaNpNc_inf_of_type_of_tree->{Np_nap_red_obs}{$tree}		=	$href_NaNpNc_inf_of_type_of_tree->{Np_obs}{$tree}		- ( $href_NaNpNc_inf_of_type_of_tree->{Np_obs}{$tree} * $href_cf_of_type_of_qt->{cf_obs}{$tree} ) ;
		$href_NaNpNc_inf_of_type_of_tree->{Nc_nap_red_exp}{$tree}		=	$href_NaNpNc_inf_of_type_of_tree->{Nc_exp}{$tree}		- ( $href_NaNpNc_inf_of_type_of_tree->{Nc_exp}{$tree} * $href_cf_of_type_of_qt->{cf_exp}{$tree} ) ;
		
		$href_NaNpNc_inf_of_type_of_tree->{Na_nap_red_obs_rev}{$tree}	=	$href_NaNpNc_inf_of_type_of_tree->{Na_obs_rev}{$tree}	- ( $href_NaNpNc_inf_of_type_of_tree->{Na_obs_rev}{$tree} * $href_cf_of_type_of_qt->{cf_obs}{$tree} ) ;
		$href_NaNpNc_inf_of_type_of_tree->{Nc_nap_red_exp_rev}{$tree}	=	$href_NaNpNc_inf_of_type_of_tree->{Nc_exp_rev}{$tree}	- ( $href_NaNpNc_inf_of_type_of_tree->{Nc_exp_rev}{$tree} * $href_cf_of_type_of_qt->{cf_exp}{$tree} ) ;
		###################################
		
		###################################
		# (2) calculate proportion of reduced pattern support (pNa_obs_red, pNp_obs_red, pNc_exp_red)
		$href_NaNpNc_inf_of_type_of_tree->{pNa_nap_red_obs}{$tree}		=	&calc_percent( $href_NaNpNc_inf_of_type_of_tree->{Na_nap_red_obs}{$tree}, $$sref_total_Npattern ) ;
		$href_NaNpNc_inf_of_type_of_tree->{pNp_nap_red_obs}{$tree}		=	&calc_percent( $href_NaNpNc_inf_of_type_of_tree->{Np_nap_red_obs}{$tree}, $$sref_total_Npattern ) ;
		$href_NaNpNc_inf_of_type_of_tree->{pNc_nap_red_exp}{$tree}		=	&calc_percent( $href_NaNpNc_inf_of_type_of_tree->{Nc_nap_red_exp}{$tree}, $$sref_total_Npattern ) ;
		###################################
		
		###################################
		# (3) ratio Nc/Na (corrected)
		$href_NaNpNc_inf_of_type_of_tree->{rNcNa_nap_red}{$tree}		=	$href_NaNpNc_inf_of_type_of_tree->{Nc_nap_red_exp}{$tree} / $href_NaNpNc_inf_of_type_of_tree->{Na_nap_red_obs}{$tree} ;
		###################################
		
		###################################
		# test print
		if ( $href_value_of_scr_par->{sys_parameter}{screen_prt} == 1 ){
			
			print	"::&evaluation_quartet_corrected_NaNpNc\n",
					"\n\tNa_obs\t\t"	, $tree, "\t", $href_NaNpNc_inf_of_type_of_tree->{Na_obs}{$tree}			,
					"\n\tNp_obs\t\t"	, $tree, "\t", $href_NaNpNc_inf_of_type_of_tree->{Np_obs}{$tree}			,
					"\n\tNc exp\t\t"	, $tree, "\t", $href_NaNpNc_inf_of_type_of_tree->{Nc_exp}{$tree}			,
					"\n"																							,
					"\n\tcf_obs\t\t"	, $tree, "\t", $href_cf_of_type_of_qt->{cf_obs}{$tree}						,
					"\n\tcf_exp\t\t"	, $tree, "\t", $href_cf_of_type_of_qt->{cf_exp}{$tree}						,
					"\n"																							,
					"\n\tNa_obs_red\t"	, $tree, "\t", $href_NaNpNc_inf_of_type_of_tree->{Na_nap_red_obs}{$tree} 	, 
					"\n\tNp_obs red\t"	, $tree, "\t", $href_NaNpNc_inf_of_type_of_tree->{Np_nap_red_obs}{$tree}	, 
					"\n\tNc exp_red\t"	, $tree, "\t", $href_NaNpNc_inf_of_type_of_tree->{Nc_nap_red_exp}{$tree}	,
					"\n"																							,
					"\n\trNc_Na_red\t"	, $tree, "\t", $href_NaNpNc_inf_of_type_of_tree->{rNcNa_nap_red}{$tree}		, 
					"\n" 				;
		}
	}
	######################################################################
}

sub evaluation_quartet_inv_var{
	
	my	$aref_rooted_qtrees					= $_[0]  ; # list of rooted quartet trees 																		IN: defined; OUT: unchanged
	my	$href_Nsites_obs_of_pattern			= $_[1]  ; # key: split-pattern code (A to V); value N observed pattern											IN: defined; OUT: unchanged
	my	$href_Nsites_exp_qt1_of_pattern		= $_[2]  ; # key: split-pattern code (A to V); value N expected given quartet-tree 1							IN: defined; OUT: unchanged
	my	$href_Nsites_exp_qt2_of_pattern		= $_[3]  ; # key: split-pattern code (A to V); value N expected given quartet-tree 2							IN: defined; OUT: unchanged
	my	$href_Nsites_exp_qt3_of_pattern		= $_[4]  ; # key: split-pattern code (A to V); value N expected given quartet-tree 3							IN: defined; OUT: unchanged
	my	$sref_tot_N_pattern_obs				= $_[5]  ; # total number of (gap excluded) site positions														IN: defined; OUT: unchanged
	my	$href_sing_inf_of_type				= $_[6]  ; # key: info type (.e.g 'sNs'); value_ type value														IN: defined, OUT: unchanged
	my	$href_sing_inf_of_type_of_tree		= $_[7]  ; # key1: tree dependent info types (.e.g 'sNs_exp_from_qt1'); key2: rooted tree, value_ type value	IN: defined, OUT: unchanged
	my	$href_invvar_inf_of_type			= $_[8]  ; # key: info type (.e.g 'sNs'); value_ type value														IN: undefined; OUT: defined
	my	$href_invvar_inf_of_type_of_tree	= $_[9]  ; # key1: tree dependent info types (.e.g 'sNs_exp_from_qt1'); key2: rooted tree, value_ type value	IN: undefined; OUT: defined
	my	$href_val_of_scr_par				= $_[10] ; # key1: sealion user available parameter options; value: optional parameter setting					IN: defined; OUT: unchanged
	
	######################################################################
	## analysis of complete variable and invariable pattern frequencies ##
	######################################################################
	
	######################################################################
	# called from...
	# &evaluation_quartet_pattern_main
	#
	############
	# using 
	# &calc_percent
	# &calc_mean
	######################################################################
	
	###################################
	# assignment of alignment observed invariable I and complete variable V pattern calculations
	# total number of observed invariable sites I 'tNi'
	$href_invvar_inf_of_type->{tNi_obs}		=	$href_Nsites_obs_of_pattern->{I} ;
	
	# total number of observed variable sites v 'tNv'
	$href_invvar_inf_of_type->{tNv_obs}		=	$href_Nsites_obs_of_pattern->{V} ;
	
	# proportion total number of observed invariable sites I related to total number of alignment observed quartet sites 'tPi'
	$href_invvar_inf_of_type->{tPi_obs}		=	&calc_percent( $href_invvar_inf_of_type->{tNi_obs}	, $$sref_tot_N_pattern_obs ) ;
	
	# proportion total number of observed variable sites V related to total number of alignment observed quartet sites 'tPv'
	$href_invvar_inf_of_type->{tPv_obs}		=	&calc_percent( $href_invvar_inf_of_type->{tNv_obs}	, $$sref_tot_N_pattern_obs ) ;
	 
	# proportion of observed invariable to variable site positions
	$href_invvar_inf_of_type->{rNi_Nv_obs}	=	$href_invvar_inf_of_type->{tNi_obs} / $href_invvar_inf_of_type->{tNv_obs} ;
	
	# proportion of lowest observed singleton proportion  + 1/4 of the variable site proportion
	$href_invvar_inf_of_type->{sPsv_obs}	=	$href_sing_inf_of_type->{sPs_obs} + ( $href_invvar_inf_of_type->{tPv_obs} / 4 ) ;
	###################################
	
	###################################
	# test print
	if ( $href_val_of_scr_par->{sys_parameter}{screen_prt} == 1 ){
		
		print	"::&evaluation_quartet_inv_var\n"								,
				"\n\ttNi_obs:\t"	, $href_invvar_inf_of_type->{tNi_obs} 		,
				"\n\ttNv_obs:\t"	, $href_invvar_inf_of_type->{tNv_obs} 		,
				"\n\ttPi_obs:\t"	, $href_invvar_inf_of_type->{tPi_obs} 		,
				"\n\ttPv_obs:\t"	, $href_invvar_inf_of_type->{tPv_obs} 		,
				"\n\trNi_Nv_obs:\t"	, $href_invvar_inf_of_type->{rNi_Nv_obs}	,
				"\n\tsPsv_obs:\t"	, $href_invvar_inf_of_type->{sPsv_obs}		; #exit ;
	}
	###################################
	
	###################################
	# assignment of quartet-tree depending INV (i) and VAR (v) calculations 'tNi_exp' and 'tNv_exp' (needed for following for-loop calculations)
	for my $tr_no ( 0 .. 2 ){
		
		# for given rooted tree assign total number of tree expected number of invariable site pattern 'tNi_exp'
		if 		( $tr_no == 0 ){ $href_invvar_inf_of_type_of_tree->{tNi_exp}{$aref_rooted_qtrees->[$tr_no]}	=	$href_Nsites_exp_qt1_of_pattern->{I} }
		elsif 	( $tr_no == 1 ){ $href_invvar_inf_of_type_of_tree->{tNi_exp}{$aref_rooted_qtrees->[$tr_no]}	=	$href_Nsites_exp_qt2_of_pattern->{I} }
		elsif 	( $tr_no == 2 ){ $href_invvar_inf_of_type_of_tree->{tNi_exp}{$aref_rooted_qtrees->[$tr_no]}	=	$href_Nsites_exp_qt3_of_pattern->{I} }
		
		# for given rooted tree assign total number of tree expected number of variable site pattern 'tNv_exp'
		if 		( $tr_no == 0 ){ $href_invvar_inf_of_type_of_tree->{tNv_exp}{$aref_rooted_qtrees->[$tr_no]}	=	$href_Nsites_exp_qt1_of_pattern->{V} }
		elsif 	( $tr_no == 1 ){ $href_invvar_inf_of_type_of_tree->{tNv_exp}{$aref_rooted_qtrees->[$tr_no]}	=	$href_Nsites_exp_qt2_of_pattern->{V} }
		elsif 	( $tr_no == 2 ){ $href_invvar_inf_of_type_of_tree->{tNv_exp}{$aref_rooted_qtrees->[$tr_no]}	=	$href_Nsites_exp_qt3_of_pattern->{V} }
		
		# for given rooted tree calculate the proportion of tree expected invariable sites in relation to the total number of observed quartet sites 'tPi_exp'
		$href_invvar_inf_of_type_of_tree->{tPi_exp}{$aref_rooted_qtrees->[$tr_no]}	=	&calc_percent( $href_invvar_inf_of_type_of_tree->{tNi_exp}{$aref_rooted_qtrees->[$tr_no]} , $$sref_tot_N_pattern_obs ) ;
		
		# for given rooted tree calculate the proportion of tree expected variable sites in relation to the total number of observed quartet sites 'tPv_exp'
		$href_invvar_inf_of_type_of_tree->{tPv_exp}{$aref_rooted_qtrees->[$tr_no]}	=	&calc_percent( $href_invvar_inf_of_type_of_tree->{tNv_exp}{$aref_rooted_qtrees->[$tr_no]} , $$sref_tot_N_pattern_obs ) ;
		
		###################################
		# test print
		if ( $href_val_of_scr_par->{sys_parameter}{screen_prt} == 1 ){
			
			print	"\n\ttNi_exp ", $aref_rooted_qtrees->[$tr_no] , ":\t"	, $href_invvar_inf_of_type_of_tree->{tNi_exp}{$aref_rooted_qtrees->[$tr_no]}	,
					"\n\ttNv_exp ", $aref_rooted_qtrees->[$tr_no] , ":\t"	, $href_invvar_inf_of_type_of_tree->{tNv_exp}{$aref_rooted_qtrees->[$tr_no]}	,
					"\n\ttPi_exp ", $aref_rooted_qtrees->[$tr_no] , ":\t"	, $href_invvar_inf_of_type_of_tree->{tPi_exp}{$aref_rooted_qtrees->[$tr_no]}	,
					"\n\ttPv_exp ", $aref_rooted_qtrees->[$tr_no] , ":\t"	, $href_invvar_inf_of_type_of_tree->{tPv_exp}{$aref_rooted_qtrees->[$tr_no]}	;
		}
	}#exit;
	###################################
	
	###################################
	# singleton calculations for each quartet-tree based on 'sNs_exp' and 'tNs_exp'
	for my $tr_no ( 0 .. 2 ){
		
		# for given rooted tree calculate the mean total number of convergent expected number of invariable sites I (expected from the other two trees) 'tNi_exp_for'
		my @scores ;
		if 		( $tr_no == 0 ){ @scores = ( $href_invvar_inf_of_type_of_tree->{tNi_exp}{$aref_rooted_qtrees->[1]}, $href_invvar_inf_of_type_of_tree->{tNi_exp}{$aref_rooted_qtrees->[2]} ) }
		elsif 	( $tr_no == 1 ){ @scores = ( $href_invvar_inf_of_type_of_tree->{tNi_exp}{$aref_rooted_qtrees->[0]}, $href_invvar_inf_of_type_of_tree->{tNi_exp}{$aref_rooted_qtrees->[2]} ) }
		elsif 	( $tr_no == 2 ){ @scores = ( $href_invvar_inf_of_type_of_tree->{tNi_exp}{$aref_rooted_qtrees->[0]}, $href_invvar_inf_of_type_of_tree->{tNi_exp}{$aref_rooted_qtrees->[1]} ) }
		$href_invvar_inf_of_type_of_tree->{tNi_exp_for}{$aref_rooted_qtrees->[$tr_no]}		=	&calc_mean ( \@scores ) ;
		
		# for given rooted tree calculate the mean total number of convergent expected number of variable sites V (expected from the other two trees) 'tNv_exp_for'
		if 		( $tr_no == 0 ){ @scores = ( $href_invvar_inf_of_type_of_tree->{tNv_exp}{$aref_rooted_qtrees->[1]}, $href_invvar_inf_of_type_of_tree->{tNv_exp}{$aref_rooted_qtrees->[2]} ) }
		elsif 	( $tr_no == 1 ){ @scores = ( $href_invvar_inf_of_type_of_tree->{tNv_exp}{$aref_rooted_qtrees->[0]}, $href_invvar_inf_of_type_of_tree->{tNv_exp}{$aref_rooted_qtrees->[2]} ) }
		elsif 	( $tr_no == 2 ){ @scores = ( $href_invvar_inf_of_type_of_tree->{tNv_exp}{$aref_rooted_qtrees->[0]}, $href_invvar_inf_of_type_of_tree->{tNv_exp}{$aref_rooted_qtrees->[1]} ) }
		$href_invvar_inf_of_type_of_tree->{tNv_exp_for}{$aref_rooted_qtrees->[$tr_no]}		=	&calc_mean	 ( \@scores ) ;
		
		# proportion of observed invariable to variable site positions
		$href_invvar_inf_of_type->{rNi_Nv_exp_for}{$aref_rooted_qtrees->[$tr_no]}			=	$href_invvar_inf_of_type_of_tree->{tNi_exp_for}{$aref_rooted_qtrees->[$tr_no]} / $href_invvar_inf_of_type_of_tree->{tNi_exp_for}{$aref_rooted_qtrees->[$tr_no]} ;
		
		# for given rooted tree calculate the proportion of convergent expected number of invariable sites I to the overall number of observed quartet sites 'tPi_exp_for'
		$href_invvar_inf_of_type_of_tree->{tPi_exp_for}{$aref_rooted_qtrees->[$tr_no]}		=	&calc_percent ( $href_invvar_inf_of_type_of_tree->{tNi_exp_for}{$aref_rooted_qtrees->[$tr_no]}, $$sref_tot_N_pattern_obs ) ;
		
		# for given rooted tree calculate the proportion of convergent expected number of variable sites V to the overall number of observed quartet sites 'tPv_exp_for'
		$href_invvar_inf_of_type_of_tree->{tPv_exp_for}{$aref_rooted_qtrees->[$tr_no]}		=	&calc_percent ( $href_invvar_inf_of_type_of_tree->{tNv_exp_for}{$aref_rooted_qtrees->[$tr_no]}, $$sref_tot_N_pattern_obs ) ;
		
		# proportion of invariable + 1/4 variable site positions expected from he other two quartet-trees 'sPsv_exp_for'
		$href_invvar_inf_of_type_of_tree->{sPsv_exp_for}{$aref_rooted_qtrees->[$tr_no]}		=	$href_sing_inf_of_type_of_tree	->{sPs_exp_for}	{$aref_rooted_qtrees->[$tr_no]} + ($href_invvar_inf_of_type_of_tree->{tPv_exp_for}{$aref_rooted_qtrees->[$tr_no]} / 4 ) ;
		
		# for given rooted tree calculate the difference between lowest P of singletons+variable sites (v) expected for actual tree to lowest P of observed singletons+v 'dsPsv_ExpObs_for'
		$href_invvar_inf_of_type_of_tree->{dsPsv_ExpObs_for}{$aref_rooted_qtrees->[$tr_no]}	=	$href_invvar_inf_of_type_of_tree->{sPsv_exp_for}{$aref_rooted_qtrees->[$tr_no]} - $href_invvar_inf_of_type->{sPsv_obs} ;
	}
	###################################
	
	######################################################################
}

sub evaluation_quartet_NaNpNc_signal{
	
	my	$aref_rooted_qtrees					= $_[0] ; # list of rooted quartet trees 																	IN: defined; OUT: unchanged
	my	$aref_quartet_taxa					= $_[1] ; # list of quartet taxon-name in seq-quartet order													IN: defined; OUT: unchanged
	my	$href_clan_of_tax					= $_[2] ; # key: clan-code; value: clan assigned taxon name													IN: defined; OUT: unchanged
	my	$href_Nobs_of_pattern				= $_[3] ; # key: split pattern code ( A to V ); value: N observed in given sequnce quartet					IN: defined; OUT: unchanged
	my	$href_Nexp_of_tree1_pattern			= $_[4] ; # pattern frequencies expected given tree1 (key: pattern code; value: frequency number			IN: defined, OUT: unchanged
	my	$href_Nexp_of_tree2_pattern			= $_[5] ; # pattern frequencies expected given tree2 (key: pattern code; value: frequency number			IN: defined, OUT: unchanged
	my	$href_Nexp_of_tree3_pattern			= $_[6] ; # pattern frequencies expected given tree3 (key: pattern code; value: frequency number			IN: defined, OUT: unchanged
	my	$sref_N_quartet_observed_sites		= $_[7] ; # total number of observed quartet site positions													IN: defined, OUT: unchanged
	my	$href_NaNpNc_inf_of_type_of_tree	= $_[8] ; # key1: tree dependent info types (.e.g 'Na', 'Np'...); key2: rooted tree, value_ type value		IN: undefined; OUT: defined
	my	$href_val_of_scr_par				= $_[9] ; # key1: sealion user available parameter options; value: optional parameter setting				IN: defined; OUT: unchanged
	
	######################################################################
	# called from...
	# &evaluation_quartet_pattern_main
	#
	############
	# using 
	# &calc_percent
	######################################################################
	
	######################################################################
	## evaluate Na, Np, Nc for each quartet tree; do rNc/Na             ##
	######################################################################
	
	######################################################################
	# test print
	if ( $href_val_of_scr_par->{sys_parameter}{screen_prt} == 1 ){
		
		print "::&evaluation_quartet_NaNpNc_signal\n";
		for	my $n ( 0 .. @$aref_quartet_taxa-1 ){ print "\n\tT", $n,":\t", $href_clan_of_tax->{$aref_quartet_taxa->[$n]} } print "\n"
	}
	######################################################################
	
	######################################################################
	# for each rooted quartet-tree...
	# (1) identification of Na, Np, and Nc, based on corresponding ingroup relationship
	# (2) calculation of proportion pNa, pNp, and pNc in respect of the total number of quartet site positions
	# (3) calculate ratio Nc/Na
	for	my $qtree ( @$aref_rooted_qtrees ){
		
		
		###################################
		# (1) calculate rooted quartet-tree depending Na, Np, and Nc split pattern signal
		# based on ingroup identification of (support corresponding) rooted quartet-tree
		# ingroup order can be skiped due to clan sorting through generation of rooted quartet-trees
		if			( 	( $qtree =~ /$href_clan_of_tax->{$aref_quartet_taxa->[2]},$href_clan_of_tax->{$aref_quartet_taxa->[3]}/ ) ||
						( $qtree =~ /$href_clan_of_tax->{$aref_quartet_taxa->[3]},$href_clan_of_tax->{$aref_quartet_taxa->[2]}/ ) ){ 
			
			if ( $href_val_of_scr_par->{sys_parameter}{screen_prt} == 1 ){ print	"\n\t", $qtree, "\t(T2,T3)", "\tNa = A + E" }
			
			$href_NaNpNc_inf_of_type_of_tree->{Na_obs}{$qtree} 		= 	$href_Nobs_of_pattern->{A} + $href_Nobs_of_pattern->{E} ; 
			$href_NaNpNc_inf_of_type_of_tree->{Np_obs}{$qtree}		= 	$href_Nobs_of_pattern->{D} ; 
			$href_NaNpNc_inf_of_type_of_tree->{Nc_exp}{$qtree}		= ( $href_Nexp_of_tree2_pattern->{A} + $href_Nexp_of_tree2_pattern->{E} + $href_Nexp_of_tree3_pattern->{A} + $href_Nexp_of_tree3_pattern->{E} ) / 2	;
			
			$href_NaNpNc_inf_of_type_of_tree->{Na_obs_rev}{$qtree} 	= 	$href_Nobs_of_pattern->{A} + $href_Nobs_of_pattern->{D} ; 
			$href_NaNpNc_inf_of_type_of_tree->{Nc_exp_rev}{$qtree}	= ( $href_Nexp_of_tree2_pattern->{A} + $href_Nexp_of_tree2_pattern->{D} + $href_Nexp_of_tree3_pattern->{A} + $href_Nexp_of_tree3_pattern->{D} ) / 2	;
		}
			
		elsif		( 	( $qtree =~ /$href_clan_of_tax->{$aref_quartet_taxa->[0]},$href_clan_of_tax->{$aref_quartet_taxa->[1]}/ ) ||
						( $qtree =~ /$href_clan_of_tax->{$aref_quartet_taxa->[1]},$href_clan_of_tax->{$aref_quartet_taxa->[0]}/ ) ){
			
			if ( $href_val_of_scr_par->{sys_parameter}{screen_prt} == 1 ){ print	"\n\t", $qtree, "\t(T0,T1)", "\tNa = A + D" }
			
			$href_NaNpNc_inf_of_type_of_tree->{Na_obs}{$qtree}		= 	$href_Nobs_of_pattern->{A} + $href_Nobs_of_pattern->{D} ; 
			$href_NaNpNc_inf_of_type_of_tree->{Np_obs}{$qtree}		= 	$href_Nobs_of_pattern->{E} ; 
			$href_NaNpNc_inf_of_type_of_tree->{Nc_exp}{$qtree}		= ( $href_Nexp_of_tree2_pattern->{A} + $href_Nexp_of_tree2_pattern->{D} + $href_Nexp_of_tree3_pattern->{A} + $href_Nexp_of_tree3_pattern->{D} ) / 2	;
			
			$href_NaNpNc_inf_of_type_of_tree->{Na_obs_rev}{$qtree} 	= 	$href_Nobs_of_pattern->{A} + $href_Nobs_of_pattern->{E} ; 
			$href_NaNpNc_inf_of_type_of_tree->{Nc_exp_rev}{$qtree}	= ( $href_Nexp_of_tree2_pattern->{A} + $href_Nexp_of_tree2_pattern->{E} + $href_Nexp_of_tree3_pattern->{A} + $href_Nexp_of_tree3_pattern->{E} ) / 2	;
		}
		
		elsif		( 	( $qtree =~ /$href_clan_of_tax->{$aref_quartet_taxa->[1]},$href_clan_of_tax->{$aref_quartet_taxa->[3]}/ ) ||
						( $qtree =~ /$href_clan_of_tax->{$aref_quartet_taxa->[3]},$href_clan_of_tax->{$aref_quartet_taxa->[1]}/ ) ){ 
			
			if ( $href_val_of_scr_par->{sys_parameter}{screen_prt} == 1 ){ print	"\n\t", $qtree, "\t(T1,T3)", "\tNa = B + F" }
			
			$href_NaNpNc_inf_of_type_of_tree->{Na_obs}{$qtree}		= 	$href_Nobs_of_pattern->{B} + $href_Nobs_of_pattern->{F} ; 
			$href_NaNpNc_inf_of_type_of_tree->{Np_obs}{$qtree}		= 	$href_Nobs_of_pattern->{G} ; 
			$href_NaNpNc_inf_of_type_of_tree->{Nc_exp}{$qtree}		= ( $href_Nexp_of_tree1_pattern->{B} + $href_Nexp_of_tree1_pattern->{F} + $href_Nexp_of_tree3_pattern->{B} + $href_Nexp_of_tree3_pattern->{F} ) / 2 ;
			
			$href_NaNpNc_inf_of_type_of_tree->{Na_obs_rev}{$qtree} 	= 	$href_Nobs_of_pattern->{B} + $href_Nobs_of_pattern->{G} ; 
			$href_NaNpNc_inf_of_type_of_tree->{Nc_exp_rev}{$qtree}	= ( $href_Nexp_of_tree1_pattern->{B} + $href_Nexp_of_tree1_pattern->{G} + $href_Nexp_of_tree3_pattern->{B} + $href_Nexp_of_tree3_pattern->{G} ) / 2	;
		}
		
		elsif		( 	( $qtree =~ /$href_clan_of_tax->{$aref_quartet_taxa->[0]},$href_clan_of_tax->{$aref_quartet_taxa->[2]}/ ) ||
						( $qtree =~ /$href_clan_of_tax->{$aref_quartet_taxa->[2]},$href_clan_of_tax->{$aref_quartet_taxa->[0]}/ ) ){ 
			
			if ( $href_val_of_scr_par->{sys_parameter}{screen_prt} == 1 ){ print	"\n\t", $qtree, "\t(T0,T2)", "\tNa = B + G" }
			
			$href_NaNpNc_inf_of_type_of_tree->{Na_obs}{$qtree}		=	$href_Nobs_of_pattern->{B} + $href_Nobs_of_pattern->{G} ; 
			$href_NaNpNc_inf_of_type_of_tree->{Np_obs}{$qtree}		=	$href_Nobs_of_pattern->{F} ; 
			$href_NaNpNc_inf_of_type_of_tree->{Nc_exp}{$qtree}		= ( $href_Nexp_of_tree1_pattern->{B} + $href_Nexp_of_tree1_pattern->{G} + $href_Nexp_of_tree3_pattern->{B} + $href_Nexp_of_tree3_pattern->{G} ) / 2 ;
			
			$href_NaNpNc_inf_of_type_of_tree->{Na_obs_rev}{$qtree} 	= 	$href_Nobs_of_pattern->{B} + $href_Nobs_of_pattern->{F} ; 
			$href_NaNpNc_inf_of_type_of_tree->{Nc_exp_rev}{$qtree}	= ( $href_Nexp_of_tree1_pattern->{B} + $href_Nexp_of_tree1_pattern->{F} + $href_Nexp_of_tree3_pattern->{B} + $href_Nexp_of_tree3_pattern->{F} ) / 2	;
		}
		
		elsif		( 	( $qtree =~ /$href_clan_of_tax->{$aref_quartet_taxa->[1]},$href_clan_of_tax->{$aref_quartet_taxa->[2]}/ ) ||
						( $qtree =~ /$href_clan_of_tax->{$aref_quartet_taxa->[2]},$href_clan_of_tax->{$aref_quartet_taxa->[1]}/ ) ){ 
			
			if ( $href_val_of_scr_par->{sys_parameter}{screen_prt} == 1 ){ print	"\n\t", $qtree, "\t(T1,T2)", "\tNa = C + H" }
			
			$href_NaNpNc_inf_of_type_of_tree->{Na_obs}{$qtree}		=	$href_Nobs_of_pattern->{C} + $href_Nobs_of_pattern->{H} ; 
			$href_NaNpNc_inf_of_type_of_tree->{Np_obs}{$qtree}		=	$href_Nobs_of_pattern->{J} ; 
			$href_NaNpNc_inf_of_type_of_tree->{Nc_exp}{$qtree}		= ( $href_Nexp_of_tree1_pattern->{C} + $href_Nexp_of_tree1_pattern->{H} + $href_Nexp_of_tree2_pattern->{C} + $href_Nexp_of_tree2_pattern->{H} ) / 2 ;
			
			$href_NaNpNc_inf_of_type_of_tree->{Na_obs_rev}{$qtree} 	= 	$href_Nobs_of_pattern->{C} + $href_Nobs_of_pattern->{J} ; 
			$href_NaNpNc_inf_of_type_of_tree->{Nc_exp_rev}{$qtree}	= ( $href_Nexp_of_tree1_pattern->{C} + $href_Nexp_of_tree1_pattern->{J} + $href_Nexp_of_tree2_pattern->{C} + $href_Nexp_of_tree2_pattern->{J} ) / 2	;
		}
		elsif		( 	( $qtree =~ /$href_clan_of_tax->{$aref_quartet_taxa->[0]},$href_clan_of_tax->{$aref_quartet_taxa->[3]}/ ) ||
						( $qtree =~ /$href_clan_of_tax->{$aref_quartet_taxa->[3]},$href_clan_of_tax->{$aref_quartet_taxa->[0]}/ ) ){ 
			
			if ( $href_val_of_scr_par->{sys_parameter}{screen_prt} == 1 ){ print	"\n\t", $qtree, "\t(T0,T3)", "\tNa = C + J" }
						
			$href_NaNpNc_inf_of_type_of_tree->{Na_obs}{$qtree}		=	$href_Nobs_of_pattern->{C} + $href_Nobs_of_pattern->{J} ; 
			$href_NaNpNc_inf_of_type_of_tree->{Np_obs}{$qtree}		=	$href_Nobs_of_pattern->{H} ; 
			$href_NaNpNc_inf_of_type_of_tree->{Nc_exp}{$qtree}		= ( $href_Nexp_of_tree1_pattern->{C} + $href_Nexp_of_tree1_pattern->{J} + $href_Nexp_of_tree2_pattern->{C} + $href_Nexp_of_tree2_pattern->{J} ) / 2 ;
			
			$href_NaNpNc_inf_of_type_of_tree->{Na_obs_rev}{$qtree} 	= 	$href_Nobs_of_pattern->{C} + $href_Nobs_of_pattern->{H} ; 
			$href_NaNpNc_inf_of_type_of_tree->{Nc_exp_rev}{$qtree}	= ( $href_Nexp_of_tree1_pattern->{C} + $href_Nexp_of_tree1_pattern->{H} + $href_Nexp_of_tree2_pattern->{C} + $href_Nexp_of_tree2_pattern->{H} ) / 2	;
		}
		###################################
		
		
		###################################
		# (2) calculation of proportion pNa, pNp, and pNc in respect of the total number of quartet site positions
		$href_NaNpNc_inf_of_type_of_tree->{pNa_obs}{$qtree}	=	&calc_percent( $href_NaNpNc_inf_of_type_of_tree->{Na_obs}{$qtree}, $$sref_N_quartet_observed_sites ) ;
		$href_NaNpNc_inf_of_type_of_tree->{pNp_obs}{$qtree}	=	&calc_percent( $href_NaNpNc_inf_of_type_of_tree->{Np_obs}{$qtree}, $$sref_N_quartet_observed_sites ) ;
		$href_NaNpNc_inf_of_type_of_tree->{pNc_exp}{$qtree}	=	&calc_percent( $href_NaNpNc_inf_of_type_of_tree->{Nc_exp}{$qtree}, $$sref_N_quartet_observed_sites ) ;
		###################################
		
		###################################
		# (3) calculate ratio Nc/Na of actual rooted quartet-tree
		$href_NaNpNc_inf_of_type_of_tree->{rNcNa_org}{$qtree}	=	$href_NaNpNc_inf_of_type_of_tree->{Nc_exp}{$qtree} / $href_NaNpNc_inf_of_type_of_tree->{Na_obs}{$qtree} ;
		###################################
		
		###################################
		# test print
		if ( $href_val_of_scr_par->{sys_parameter}{screen_prt} == 1 ){
			
			print	"\n\tNa_obs\t\t"	, $qtree , "\t",	$href_NaNpNc_inf_of_type_of_tree->{Na_obs}{$qtree} 		,
					"\n\tNp_obs\t\t"	, $qtree , "\t",	$href_NaNpNc_inf_of_type_of_tree->{Np_obs}{$qtree} 		,
					"\n\tNc_exp\t\t"	, $qtree , "\t",	$href_NaNpNc_inf_of_type_of_tree->{Nc_exp}{$qtree} 		, 
					"\n\tNa_obs_rev\t"	, $qtree , "\t",	$href_NaNpNc_inf_of_type_of_tree->{Na_obs_rev}{$qtree} 	,
					"\n\tNc_exp_rev\t"	, $qtree , "\t",	$href_NaNpNc_inf_of_type_of_tree->{Nc_exp_rev}{$qtree} 	, "\n" ;
		}
		###################################
	}
	######################################################################
}

sub	evaluation_quartet_napv_pattern_info{
	
	my	$aref_rooted_qtrees					= $_[0] ;	# list of rooted quartet trees 																									IN: defined; OUT: unchanged
	my	$href_sing_inf_of_type				= $_[1] ;	# key: info type (.e.g 'sNs'); value_ type value																				IN: undefined; OUT: defined
	my	$href_sing_inf_of_type_of_tree		= $_[2] ;	# key1: tree dependent info types (.e.g 'sNs_exp_from_qt1'); key2: rooted tree, value_ type value								IN: undefined; OUT: defined
	my	$href_invvar_inf_of_type			= $_[3] ;	# key: info type (.e.g 'sNs'); value_ type value																				IN: undefined; OUT: defined
	my	$href_invvar_inf_of_type_of_tree	= $_[4] ;	# key1: tree dependent info types (.e.g 'sNs_exp_from_qt1'); key2: rooted tree, value_ type value								IN: undefined; OUT: defined
	my	$href_NaNpNc_inf_of_type_of_tree	= $_[5] ;	# key1: tree dependent info types (.e.g 'Na', 'Np'...); key2: rooted tree, value_ type value									IN: defined; OUT: unchanged
	my	$href_cf_of_type_of_qt				= $_[6] ;	# key1: info type (e.g. ('cf_obs'); key2: rooted quartet-tree; value: correction factor for site patterns of actual tree		IN: undefined; OUT: defined	
	my	$href_val_of_scr_par				= $_[7] ;	# key1: sealion user available parameter options; value: optional parameter setting							IN: defined; OUT: unchanged
	
	######################################################################
	# called from...
	# &evaluation_quartet_pattern_main
	#
	############
	# using 
	# &calc_cf
	# &calc_mean
	######################################################################
	
	######################################################################
	## analysis of Nav, Npv, Ncv and singleton data in sense of napv    ##
	######################################################################
	
	# shortest number of observed singletons plus 1/4 of the total number of observed variable sites v 'sNsv_obs'
	$href_sing_inf_of_type->{sNsv_obs}	=	$href_sing_inf_of_type->{sNs_obs} + ( $href_invvar_inf_of_type->{tNv_obs} / 4 ) ;
	
	# total number of observed singletons plus total number of observed variable sites v 'tNsv_obs'
	$href_sing_inf_of_type->{tNsv_obs}	=	$href_sing_inf_of_type->{tNs_obs} +   $href_invvar_inf_of_type->{tNv_obs} ;
	
	# lowest P of observed singletons (sPs_obs) + 1/4 of the total site proportion of observed variable sites v 'sNsv_obs'
	$href_sing_inf_of_type->{sPsv_obs}	=	$href_sing_inf_of_type->{sPs_obs} + ( $href_invvar_inf_of_type->{tPv_obs} / 4 ) ;
	
	# total P of observed singletons + observed variable sites
	$href_sing_inf_of_type->{tPsv_obs}	=	$href_sing_inf_of_type->{tPs_obs} +   $href_invvar_inf_of_type->{tPv_obs} ;
	
	for my $tr_no ( 0 .. 2 ){
		
		# shortest number of tree expected singletons plus 1/4 of the total number of expected variable sites v 'sNsv_exp'
		$href_sing_inf_of_type_of_tree->{sNsv_exp}{$aref_rooted_qtrees->[$tr_no]}		=	$href_sing_inf_of_type_of_tree->{sNs_exp}{$aref_rooted_qtrees->[$tr_no]} + ( $href_invvar_inf_of_type_of_tree->{tNv_exp}{$aref_rooted_qtrees->[$tr_no]} / 4 );
		
		# total number of tree expected singletons plus total number of expected variable sites v 'tNsv_exp'
		$href_sing_inf_of_type_of_tree->{tNsv_exp}{$aref_rooted_qtrees->[$tr_no]}		=	$href_sing_inf_of_type_of_tree->{tNs_exp}{$aref_rooted_qtrees->[$tr_no]} + 	 $href_invvar_inf_of_type_of_tree->{tNv_exp}{$aref_rooted_qtrees->[$tr_no]} ;
		
		# shortest percentage of tree expected singletons plus 1/4 of the total percentage of expected variable sites v 'sNsv_exp'
		$href_sing_inf_of_type_of_tree->{sPsv_exp}{$aref_rooted_qtrees->[$tr_no]}		=	$href_sing_inf_of_type_of_tree->{sPs_exp}{$aref_rooted_qtrees->[$tr_no]} + ( $href_invvar_inf_of_type_of_tree->{tPv_exp}{$aref_rooted_qtrees->[$tr_no]} / 4 );
		
		# total percentage of tree expected singletons plus total percentage of expected variable sites v 'tPsv_exp'
		$href_sing_inf_of_type_of_tree->{tPsv_exp}{$aref_rooted_qtrees->[$tr_no]}		=	$href_sing_inf_of_type_of_tree->{tPs_exp}{$aref_rooted_qtrees->[$tr_no]} + 	 $href_invvar_inf_of_type_of_tree->{tPv_exp}{$aref_rooted_qtrees->[$tr_no]} ;
	}
	
	for my $tr_no ( 0 .. 2 ){
		
		# distance percentage of tree expected singletons+variable sites - percentage of observed singletons+variable sites 'dsPsv_ExpObs'
		$href_sing_inf_of_type_of_tree	->{dsPsv_ExpObs}{$aref_rooted_qtrees->[$tr_no]}	=	$href_sing_inf_of_type_of_tree->{sPsv_exp}{$aref_rooted_qtrees->[$tr_no]} - $href_sing_inf_of_type->{sPsv_obs} ;
		
		# calculation of cf correction factor for observed split pattern 'cfv_obs'
		$href_cf_of_type_of_qt			->{cfv_obs}		{$aref_rooted_qtrees->[$tr_no]} = 	&calc_cf( $href_sing_inf_of_type->{sNsv_obs}, $href_sing_inf_of_type->{tNsv_obs} 	) ;
		
		# calculate correction factor expected for each tree, including variable sites in singletons 'cfv_exp'
		my	@scores ;
		if 		( $tr_no == 0 ){ 
			
			$scores[0] = &calc_cf( $href_sing_inf_of_type_of_tree->{sNsv_exp}{$aref_rooted_qtrees->[1]}, $href_sing_inf_of_type_of_tree->{tNsv_exp}{$aref_rooted_qtrees->[1]} ) ;
			$scores[1] = &calc_cf( $href_sing_inf_of_type_of_tree->{sNsv_exp}{$aref_rooted_qtrees->[2]}, $href_sing_inf_of_type_of_tree->{tNsv_exp}{$aref_rooted_qtrees->[2]} ) ;
		}
		elsif 	( $tr_no == 1 ){ 
			
			$scores[0] = &calc_cf( $href_sing_inf_of_type_of_tree->{sNsv_exp}{$aref_rooted_qtrees->[0]}, $href_sing_inf_of_type_of_tree->{tNsv_exp}{$aref_rooted_qtrees->[0]} ) ;
			$scores[1] = &calc_cf( $href_sing_inf_of_type_of_tree->{sNsv_exp}{$aref_rooted_qtrees->[2]}, $href_sing_inf_of_type_of_tree->{tNsv_exp}{$aref_rooted_qtrees->[2]} ) ;
		}	
		elsif 	( $tr_no == 2 ){ 
			
			$scores[0] = &calc_cf( $href_sing_inf_of_type_of_tree->{sNsv_exp}{$aref_rooted_qtrees->[0]}, $href_sing_inf_of_type_of_tree->{tNsv_exp}{$aref_rooted_qtrees->[0]} ) ;
			$scores[1] = &calc_cf( $href_sing_inf_of_type_of_tree->{sNsv_exp}{$aref_rooted_qtrees->[1]}, $href_sing_inf_of_type_of_tree->{tNsv_exp}{$aref_rooted_qtrees->[1]} ) ;
		}
		$href_cf_of_type_of_qt->{cfv_exp}{$aref_rooted_qtrees->[$tr_no]}						=	&calc_mean( \@scores ) ;
		
		# calculate number of potential, alignment observed, apomorph split-signal Napv, corrected by the variable site influenced (observed) correction factor 'Nav_obs_red'
		$href_NaNpNc_inf_of_type_of_tree->{Nav_obs_red}		{$aref_rooted_qtrees->[$tr_no]}		=	$href_NaNpNc_inf_of_type_of_tree->{Na_obs}		{$aref_rooted_qtrees->[$tr_no]} - ( $href_NaNpNc_inf_of_type_of_tree->{Na_obs}		{$aref_rooted_qtrees->[$tr_no]} * $href_cf_of_type_of_qt->{cfv_obs}{$aref_rooted_qtrees->[$tr_no]} ) ;
		
		# calculate number of alignment observed, plesiomorph split-signal Nav, corrected by the variable site influenced (observed) correction factor 'Npv_obs_red'
		$href_NaNpNc_inf_of_type_of_tree->{Npv_obs_red}		{$aref_rooted_qtrees->[$tr_no]}		=	$href_NaNpNc_inf_of_type_of_tree->{Np_obs}		{$aref_rooted_qtrees->[$tr_no]} - ( $href_NaNpNc_inf_of_type_of_tree->{Np_obs}		{$aref_rooted_qtrees->[$tr_no]} * $href_cf_of_type_of_qt->{cfv_obs}{$aref_rooted_qtrees->[$tr_no]} ) ;
		
		# calculate number of expected convergent split-signal Nc, corrected by the variable site influenced (expected) correction factor 'Ncv_exp_red'
		$href_NaNpNc_inf_of_type_of_tree->{Ncv_exp_red}		{$aref_rooted_qtrees->[$tr_no]}		=	$href_NaNpNc_inf_of_type_of_tree->{Nc_exp}		{$aref_rooted_qtrees->[$tr_no]} - ( $href_NaNpNc_inf_of_type_of_tree->{Nc_exp}		{$aref_rooted_qtrees->[$tr_no]} * $href_cf_of_type_of_qt->{cfv_exp}{$aref_rooted_qtrees->[$tr_no]} ) ;
		
		# calculate ratio of expected convergent split-signal Ncv and potential apomorph split-signal Nav 'rNcv_Nav_red'
		$href_NaNpNc_inf_of_type_of_tree->{rNcNa_napv_red}	{$aref_rooted_qtrees->[$tr_no]}		=	$href_NaNpNc_inf_of_type_of_tree->{Ncv_exp_red}	{$aref_rooted_qtrees->[$tr_no]} / 	$href_NaNpNc_inf_of_type_of_tree->{Nav_obs_red}	{$aref_rooted_qtrees->[$tr_no]} ;
		
		# calculate ratio correction factor expected for each tree related to correction factor for observed site-pattern 'rcfv'
		$href_cf_of_type_of_qt			->{rcfv_ExpObs}		{$aref_rooted_qtrees->[$tr_no]}		=	$href_cf_of_type_of_qt			->{cfv_exp}		{$aref_rooted_qtrees->[$tr_no]} / 	$href_cf_of_type_of_qt			->{cfv_obs}		{$aref_rooted_qtrees->[$tr_no]} ;
		
		# calculate number of potential, alignment observed, apomorph split-signal Navof contrary polarized tree, corrected by the variable site influenced (observed) correction factor 'Nav_obs_red'
		$href_NaNpNc_inf_of_type_of_tree->{Nav_obs_red_rev}	{$aref_rooted_qtrees->[$tr_no]}		=	$href_NaNpNc_inf_of_type_of_tree->{Na_obs_rev}	{$aref_rooted_qtrees->[$tr_no]} - ( $href_NaNpNc_inf_of_type_of_tree->{Na_obs_rev}	{$aref_rooted_qtrees->[$tr_no]} * $href_cf_of_type_of_qt->{cfv_obs}{$aref_rooted_qtrees->[$tr_no]} ) ;
		
		# calculate number of expected convergent split-signal Nc of contrary polarized tree, corrected by the variable site influenced (expected) correction factor 'Ncv_exp_red_rev'
		$href_NaNpNc_inf_of_type_of_tree->{Ncv_exp_red_rev}	{$aref_rooted_qtrees->[$tr_no]}		=	$href_NaNpNc_inf_of_type_of_tree->{Nc_exp_rev}	{$aref_rooted_qtrees->[$tr_no]} - ( $href_NaNpNc_inf_of_type_of_tree->{Nc_exp_rev}	{$aref_rooted_qtrees->[$tr_no]} * $href_cf_of_type_of_qt->{cfv_exp}{$aref_rooted_qtrees->[$tr_no]} ) ;
		
		###################################
		# test print
		if ( $href_val_of_scr_par->{sys_parameter}{screen_prt} == 1 ){
			
			print	"::&evaluation_quartet_napv_pattern_info\n" ,
					"\n\tNav_obs\t\t"	, $aref_rooted_qtrees->[$tr_no], "\t", $href_NaNpNc_inf_of_type_of_tree	->{Na_obs}		{$aref_rooted_qtrees->[$tr_no]}	,
					"\n\tNpv_obs\t\t"	, $aref_rooted_qtrees->[$tr_no], "\t", $href_NaNpNc_inf_of_type_of_tree	->{Np_obs}		{$aref_rooted_qtrees->[$tr_no]}	,
					"\n\tNcv exp\t\t"	, $aref_rooted_qtrees->[$tr_no], "\t", $href_NaNpNc_inf_of_type_of_tree	->{Nc_exp}		{$aref_rooted_qtrees->[$tr_no]}	,
					"\n"																																		,
					"\n\tcfv_obs\t\t"	, $aref_rooted_qtrees->[$tr_no], "\t", $href_cf_of_type_of_qt			->{cfv_obs}		{$aref_rooted_qtrees->[$tr_no]}	,
					"\n\tcfv_exp\t\t"	, $aref_rooted_qtrees->[$tr_no], "\t", $href_cf_of_type_of_qt			->{cfv_exp}		{$aref_rooted_qtrees->[$tr_no]}	,
					"\n"																																		,
					"\n\tNav_obs_red\t"	, $aref_rooted_qtrees->[$tr_no], "\t", $href_NaNpNc_inf_of_type_of_tree	->{Nav_obs_red}	{$aref_rooted_qtrees->[$tr_no]} , 
					"\n\tNpv_obs red\t"	, $aref_rooted_qtrees->[$tr_no], "\t", $href_NaNpNc_inf_of_type_of_tree	->{Npv_obs_red}	{$aref_rooted_qtrees->[$tr_no]}	, 
					"\n\tNcv exp_red\t"	, $aref_rooted_qtrees->[$tr_no], "\t", $href_NaNpNc_inf_of_type_of_tree	->{Ncv_exp_red}	{$aref_rooted_qtrees->[$tr_no]}	,
					"\n\trNcv_Nav_red\t", $aref_rooted_qtrees->[$tr_no], "\t", $href_NaNpNc_inf_of_type_of_tree	->{rNcNa_napv_red}{$aref_rooted_qtrees->[$tr_no]}	, "\n" ; #exit;
		}
	}
}

sub evaluation_quartet_pattern_main{				
	
	my	$aref_rooted_qtrees										=	$_[0] ;		# rooted 4-clan tree 3																																								IN: defined; OUT: unchanged
	my	$aref_ptaxa												=	$_[1] ;		# array of quartet sequences																																						IN: defined; OUT: unchanged
	my	$href_obs_n_of_pattern									=	$_[2] ;		# observed pattern frequencies (key: pattern code; value: frequency number																											IN: defined; OUT: unchanged
	my	$href_exp_n_of_tree1_pattern							=	$_[3] ;		# pattern frequencies expected given tree1 (key: pattern code; value: frequency number																								IN: defined; OUT: unchanged
	my	$href_exp_n_of_tree2_pattern							=	$_[4] ;		# pattern frequencies expected given tree2 (key: pattern code; value: frequency number																								IN: defined; OUT: unchanged
	my	$href_exp_n_of_tree3_pattern							=	$_[5] ;		# pattern frequencies expected given tree3 (key: pattern code; value: frequency number																								IN: defined; OUT: unchanged
	my	$href_obs_p_of_pattern									=	$_[6] ;		# observed pattern percentages (key: pattern code; value: frequency number																											IN: defined; OUT: unchanged
	my	$href_exp_p_of_tree1_pattern							=	$_[7] ;		# pattern percentages expected given tree1 (key: pattern code; value: frequency number																								IN: defined; OUT: unchanged
	my	$href_exp_p_of_tree2_pattern							=	$_[8] ;		# pattern percentages expected given tree2 (key: pattern code; value: frequency number																								IN: defined; OUT: unchanged
	my	$href_exp_p_of_tree3_pattern							=	$_[9] ;		# pattern percentages expected given tree3 (key: pattern code; value: frequency number																								IN: defined; OUT: unchanged
	my	$sref_Npat_total_obs									=	$_[10] ;	# total number of (gap excluded) site positions																																		IN: defined; OUT: unchanged
	my	$href_clan_of_seq										=	$_[11] ;	# hash with sequence name (key) assigned clan (value)																																IN: defined; OUT: unchanged
	my	$sref_outgroup											=	$_[12] ;	# clan name of the outgroup																																							IN: defined; OUT: unchanged
	my	$sref_clan4												=	$_[13] ;	# clan4 combination of given quartet combination																																	IN: defined; OUT: unchanged
	my	$sref_fc												=	$_[14] ;	# actual sealion quartet number																																						IN: defined; OUT: unchanged
	my	$href_hoh_sup_of_clan4_of_calc_of_qtree_of_fc			=	$_[15] ;	# multiple hash, storing the final support (nap or napv) related to each rooted quartet-tree																						IN: defined; OUT: changed
	my	$href_hol_sup_of_clan4_of_calc_of_qtree					=	$_[16] ;	# multiple hash, listing for each of the three rooted quartet-trees of a given 4clan combination all single quartet inferred support (NAP or NAPV)									IN: defined; OUT: changed
	my	$href_hoh_tsvdata_of_tsvfilename_of_printcode			=	$_[17] ;	# key1: R grafic unique printcode (e.g. 'metafileT'); key2: *.tsv filename; value: *tsv filename data addings (e.g. clanfilename, type of support...)								IN: defined; OUT: changed
	my	$href_subfolder_q										=	$_[18] ;	# hashpath to output-type specific subdirs of single quartet-pattern analyses; key: output-type (e.g. 'pdf'); value: pathdir														IN: defined; OUT: unchanged
	my	$href_value_of_scr_par									=	$_[19] ;	# key1: sealion user available parameter options; value: optional parameter setting																 									IN: defined; OUT: unchanged
	my	$href_hoh_sup_dist_of_clan4_type_of_treepair			=	$_[20] ;	# key1: clan4 code; key2: type of support ('nap', 'napr', or 'napv'); key3: tree-pair code ('first2nd' or 'first3rd'); value: distance Nap											IN: undefined; OUT: defined 
	my	$href_hoh_qrej_of_filtertype_of_scale_of_fn				=	$_[21] ;	# key1: type of supportfilter ('nap_risk1', 'napr_risk1', or 'napv_risk1'); key2: scale value; key3: quartet number; value: filter corresponding score-value						IN: undefined; OUT: defined
	my	$href_bestTree_of_fc_of_cm								=	$_[22] ;	# key1: quartet-number ($fc), key2: calculation method (nap, napv, napr); value: quartet related best supported tree																IN: undefined; OUT: defined 
	my	$href_rNcNa_of_cmfm_of_fc_of_tree						=	$_[23] ;	# key1: quartet-number ($fc), key2: calculation+filter method (nap_risk1, napv_dist2, napr_dist3); Key3: tree; value: rNc/na														IN: undefined; OUT: defined 
	my	$href_algorithm_def_methods								=	$_[24] ;	# list of defined support methods (nap, napv ...)																																	IN: defined; OUT: unchanged
	my	$href_hol_scales_of_filtertype							=	$_[25] ;	# key1: type of supportfilter ('risk', 'qdist'); value: list of ftype corresponding threshold scales																				IN: defined; OUT: unchanged
	my	$href_seen_filter_method								=	$_[26] ;	# key1: filter method incl OTU clipping approaches; value: seen +1																													IN: undefined; OUT: defined
	
	
	######################################################################
	## main evaluation of single pattern frequencies of given quartet   ##
	######################################################################
	#print "::&evaluation_quartet_pattern_main\n";
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# &evaluation_quartet_singletons
	# &evaluation_quartet_inv_var
	# &evaluation_quartet_NaNpNc_signal
	# &evaluation_quartet_ratio_NcNa
	# &evaluation_quartet_corrected_NaNpNc
	# &evaluation_quartet_support_napX
	# &evaluation_quartet_napv_pattern_info
	# &print_tsv_Q1
	# &print_tsv_Q2
	# &print_tsv_Q3
	# &print_tsv_Q4
	######################################################################
	
	######################################################################
	####### evaluate SINGLETON split pattern distribution ################
	######################################################################
	# sort distribution of observed singletons (K, L, M, N) from lowest to highest in @sorted_Ns_obs
	# determine lowest number of observed singletons (sNs_obs)
	# calculate total number of observed singletons (tNs_obs)
	# calculate cf value based on observed shortest and total number of singletons ($cf_obs)
	# calculate proportion of shortest number of singletons given the total number of observed pattern (sPs_obs)
	# calculate proportion of the total number of singletons given the total number of observed pattern (tPs_obs)
	# calculate ratio (rsNs_tNs) of lowest numbered singleton pattern (sNs_obs) to the total (tNs_obs) number of singletons
	my	(
		
		%sing_inf_of_type			,	# key: tree independent info types (.e.g 'sNs'); value_ type value
		%sing_inf_of_type_of_tree	,	# key: tree dependent info types (.e.g 'sNs_exp_from_qt1'); value_ type value
		%cf_of_type_of_qt			,	# key1: info type (e.g. ('cf_obs'); key2: rooted quartet-tree; value: correction factor for site patterns of actual tree	
	);
	
	&evaluation_quartet_singletons(
		
		\@$aref_rooted_qtrees			,	# list of rooted quartet trees 																								IN: defined; OUT: unchanged
		\%$href_obs_n_of_pattern		,	# key: split-pattern code (A to V); value N observed pattern																IN: defined; OUT: unchanged
		\%$href_exp_n_of_tree1_pattern 	,	# key: split-pattern code (A to V); value N expected given quartet-tree 1													IN: defined; OUT: unchanged
		\%$href_exp_n_of_tree2_pattern 	,	# key: split-pattern code (A to V); value N expected given quartet-tree 2													IN: defined; OUT: unchanged
		\%$href_exp_n_of_tree3_pattern 	,	# key: split-pattern code (A to V); value N expected given quartet-tree 3													IN: defined; OUT: unchanged
		\$$sref_Npat_total_obs			,	# total number of (gap excluded) site positions																				IN: defined; OUT: unchanged
		\%sing_inf_of_type				,	# key: info type (.e.g 'sNs'); value_ type value																			IN:undefined; OUT: defined
		\%sing_inf_of_type_of_tree		,	# key1: tree dependent info types (.e.g 'sNs_exp_from_qt1'); key2: rooted tree, value_ type value							IN:undefined; OUT: defined
		\%cf_of_type_of_qt				,	# key1: info type (e.g. ('cf_obs'); key2: rooted quartet-tree; value: correction factor for site patterns of actual tree	IN:undefined; OUT: defined
		\%$href_value_of_scr_par		,	# key1: sealion user available parameter options; value: optional parameter setting											IN: defined; OUT: unchanged
	);
	######################################################################
	######################################################################
	
	
	######################################################################
	####### evaluate Inv and Var split pattern distribution ##############
	######################################################################
	# identification of total number (tN), single overall site proportion (tP)
	# and pairwise proportion (rNi_Nv_obs) of alignment observed complete 
	# invariable (i) and variable (v) site pattern distribution
	my	(
		
		%invvar_inf_of_type_all_trees	,	# key: info type (.e.g 'sNs'); value_ type value
		%invvar_inf_of_type_of_tree		,	# key1: tree dependent info types (.e.g 'sNs_exp_from_qt1'); key2: rooted tree, value_ type value
	) ;
	
	&evaluation_quartet_inv_var(
		
		\@$aref_rooted_qtrees			,	# list of rooted quartet trees 																				IN: defined; OUT: unchanged
		\%$href_obs_n_of_pattern		,	# key: split-pattern code (A to V); value N observed pattern												IN: defined; OUT: unchanged
		\%$href_exp_n_of_tree1_pattern 	,	# key: split-pattern code (A to V); value N expected given quartet-tree 1									IN: defined; OUT: unchanged
		\%$href_exp_n_of_tree2_pattern 	,	# key: split-pattern code (A to V); value N expected given quartet-tree 2									IN: defined; OUT: unchanged
		\%$href_exp_n_of_tree3_pattern 	,	# key: split-pattern code (A to V); value N expected given quartet-tree 3									IN: defined; OUT: unchanged
		\$$sref_Npat_total_obs			,	# total number of (gap excluded) site positions																IN: defined, OUT: unchanged
		\%sing_inf_of_type				,	# key: info type (.e.g 'sNs'); value_ type value															IN: defined, OUT: unchanged
		\%sing_inf_of_type_of_tree		,	# key1: tree dependent info types (.e.g 'sNs_exp_from_qt1'); key2: rooted tree, value_ type value			IN: defined, OUT: unchanged
		\%invvar_inf_of_type_all_trees	,	# key: info type (.e.g 'sNs'); value_ type value															IN: undefined; OUT: defined
		\%invvar_inf_of_type_of_tree	,	# key1: tree dependent info types (.e.g 'sNs_exp_from_qt1'); key2: rooted tree, value_ type value			IN: undefined; OUT: defined
		\%$href_value_of_scr_par		,	# key1: sealion user available parameter options; value: optional parameter setting							IN: defined; OUT: unchanged
	) ;
	######################################################################
	######################################################################
	
	
	######################################################################
	####### evaluation of Na, Np, and Nc #################################
	######################################################################
	# assignment of number of... 
	#	...observed root directive split support (Na_obs)
	#	...observed plesiomorphic supported splits (Np_obs)
	#	...mean expected number of convergent split signal (mean_Nc_expected)
	# for each of the three rooted trees (qt1, qt2, qt3) and their reverse counter part
	my	( 
		
		%NaNpNc_inf_of_type_of_tree		,	# key1: tree dependent info types (.e.g 'Na', 'Np'...); key2: rooted tree, value_ type value
	) ; 
	
	&evaluation_quartet_NaNpNc_signal(
		
		\@$aref_rooted_qtrees			,	# list of rooted quartet trees 																				IN: defined; OUT: unchanged
		\@$aref_ptaxa					,	# list of quartet taxon-name in seq-quartet order															IN: defined; OUT: unchanged
		\%$href_clan_of_seq				,	# key: clan-code; value: clan assigned taxon name															IN: defined; OUT: unchanged
		\%$href_obs_n_of_pattern		,	# key: split pattern code ( A to V ); value: N observed in given sequnce quartet							IN: defined; OUT: unchanged
		\%$href_exp_n_of_tree1_pattern	,	# pattern frequencies expected given tree1 (key: pattern code; value: frequency number						IN: defined, OUT: unchanged
		\%$href_exp_n_of_tree2_pattern	,	# pattern frequencies expected given tree2 (key: pattern code; value: frequency number						IN: defined, OUT: unchanged
		\%$href_exp_n_of_tree3_pattern	,	# pattern frequencies expected given tree3 (key: pattern code; value: frequency number						IN: defined, OUT: unchanged
		\$$sref_Npat_total_obs			,	# total number of observed quartet site positions															IN: defined, OUT: unchanged
		\%NaNpNc_inf_of_type_of_tree	,	# key1: tree dependent info types (.e.g 'Na', 'Np'...); key2: rooted tree, value_ type value				IN: undefined; OUT: defined
		\%$href_value_of_scr_par		,	# key1: sea	lion user available parameter options; value: optional parameter setting						IN: defined; OUT: unchanged
	);
	######################################################################
	######################################################################
	
	
	######################################################################
	####### Risk evaluation based on signal ratio rNc/Na observed ########
	######################################################################
	# compare rNc/Na of each rooted quartet tree
	# - count number of rooted trees within a rNc/Na range of 0.975 to 1.000 within risk range ($counter_rNcNa_risk_range)  
	# - count number of rooted trees with a rNc/Na score > upper threshold value ($counter_rNc_Na_positive)
	# - count number of rooted trees with a rNc/Na score < lower threshold value ($counter_rNc_Na_no_risk)
	if (	( $href_value_of_scr_par->{seq_filterRisk}{risk1} == 1 )	||
			( $href_value_of_scr_par->{seq_filterRisk}{risk2} == 1 )	){
		
		&evaluation_quartet_ratio_NcNa(
			
			\$$sref_clan4									, # actual 4clan code																																			IN: defined; OUT: unchanged 
			\$$sref_fc										, # process number of the current quartet																														IN: defined; OUT: unchanged 
			\%NaNpNc_inf_of_type_of_tree					, # key1: type-code (e.g. 'rNcNaRisk1', 'rNcNaRisk2' ...); key2: rooted tree, value_ type value																	IN: defined, OUT: unchanged
			\%$href_value_of_scr_par						, # key1: sealion user available parameter options; value: optional parameter setting																			IN: defined; OUT: unchanged
			\%$href_hoh_qrej_of_filtertype_of_scale_of_fn	, # key1: type of supportfilter ('nap_risk1', 'napr_risk1', or 'napv_risk1'); key2: scale value; key3: quartet number; value: filter corresponding score-value	IN: undefined; OUT: defined
			\'rNcNa_org'									, # type of Nc/Na ratio (uncorrected->raw)																														IN: defined; OUT: unchanged
			\%$href_rNcNa_of_cmfm_of_fc_of_tree				, # key1: quartet-number ($fc), key2: calculation+filter method (nap_risk1, napv_dist2, napr_dist3); Key3: tree; value: rNc/na									IN: undefined; OUT: defined 
			\%$href_algorithm_def_methods					, # list of defined support methods (nap, napv ...)																												IN: defined; OUT: unchanged
			\%$href_hol_scales_of_filtertype				, # key1: type of supportfilter ('risk', 'qdist'); value: list of ftype corresponding threshold scales															IN: defined; OUT: unchanged
			\%$href_seen_filter_method						, # key1: filter method incl OTU clipping approaches; value: seen +1																							IN: undefined; OUT: defined
		);
	}
	######################################################################
	######################################################################
	
	
	######################################################################
	####### Correction of Na, Np, Nc and calculation rNc/Na red ##########
	######################################################################
	# for each of the three rooted trees (qt1, qt2, qt3)...
	# (1) calculate cf corrected number of reduced pattern support (Nb_obs_red, Na_obs_red, Np_obs_red, Nc_exp_red)
	# ...using observed calculated cf_obs value for Na_obs and Np_obs
	# ...using tree related cf_exp value for $mean_Nc_exp
	# (2) calculate proportion of reduced pattern support (pNa_obs_red, pNp_obs_red, pNc_exp_red)
	# (3) ratio Nc/Na (corrected)
	if ( $href_value_of_scr_par->{seq_algorithm}{nap} == 1 ){
	
		&evaluation_quartet_corrected_NaNpNc(
			
			\@$aref_rooted_qtrees			,	# list of rooted quartet trees (qtree1, qtree2, qtree3)																		IN: defined; OUT: unchanged 
			\%cf_of_type_of_qt				,	# key1: info type (e.g. ('cf_obs'); key2: rooted quartet-tree; value: correction factor for site patterns of actual tree	IN:	defined; OUT: unchanged
			\%NaNpNc_inf_of_type_of_tree	,	# key1: tree dependent info types (.e.g 'Na', 'Np'...); key2: rooted tree, value_ type value								IN: defined; OUT: changed
			\$$sref_Npat_total_obs			,	# total number of observed quartet site positions																			IN: defined; OUT: unchanged
			\%$href_value_of_scr_par		,	# key1: sealion user available parameter options; value: optional parameter setting											IN: defined; OUT: unchanged
		) ;
	}
	######################################################################
	######################################################################
	
	
	######################################################################
	####### evaluation of Nap score ######################################
	######################################################################
	my	(
		
		%finalsupport_of_type_of_qt	,	# key1: support type (e.g. 'nap', 'napv'...); key2: rooted quartet tree; value: score
	);
	
	if ( $href_value_of_scr_par->{seq_algorithm}{nap} == 1 ){
		
		my	@signal_code = ("nap", "Na_nap_red_obs", "Nc_nap_red_exp") ;	# array of hashkeys used for specific support calculation (1st: support code; 2nd: pot apomorph signal code; 3rd expected convergency code)
		
		&evaluation_quartet_support_napX(
			
			\$$sref_clan4										,	# actual 4clan code																																IN: defined; OUT: unchanged 
			\$$sref_fc											,	# process number of the current quartet																											IN: defined; OUT: unchanged 
			\@signal_code										,	# array of hashkeys used for specific support calculation (1st: support code; 2nd: pot apomorph signal code; 3rd expected convergency code)		IN: defined; OUT: unchanged 
			\@$aref_rooted_qtrees								,	# list of rooted quartet trees (qtree1, qtree2, qtree3)																							IN: defined; OUT: unchanged 
			\%NaNpNc_inf_of_type_of_tree						,	# key1: tree dependent info types (.e.g 'Na', 'Np'...); key2: rooted tree, value_ type value													IN: defined; OUT: unchanged
			\%$href_hol_sup_of_clan4_of_calc_of_qtree			,	# key1: 4-clan code; key2: type NAP; key3: rooted quartet-tree; value: list of single NAP quartet support										IN: defined; OUT: changed
			\%$href_hoh_sup_of_clan4_of_calc_of_qtree_of_fc		,	# key1: 4-clan code; key2: type NAP; key3: rooted quartet-tree; key4: quartet number; value:NAP quartet support									IN: defined; OUT: changed
			\%finalsupport_of_type_of_qt						,	# key1: support type (e.g. 'nap', 'napv'...); key2: rooted quartet tree; value: score															IN: undefined; OUT: defined 
			\%$href_hoh_sup_dist_of_clan4_type_of_treepair		,	# key1: clan4 code; key2: type of support ('nap', 'napr', or 'napv'); key3: tree-pair code ('first2nd' or 'first3rd'); value: distance Nap		IN: undefined; OUT: defined 
			\%$href_bestTree_of_fc_of_cm						,	# key1: quartet-number ($fc), key2: calculation method (nap, napv, napr); value: quartet related best supported tree							IN: undefined; OUT: defined 
			\%$href_value_of_scr_par							,	# key1: sealion user available parameter options; value: optional parameter setting																IN: defined; OUT: unchanged
		)
	}
	######################################################################
	######################################################################
	
	
	######################################################################
	####### quartet evaluation due to tree support distances          ####
	######################################################################
	if ( $href_value_of_scr_par->{seq_filterDist}{qdist1} == 1 ){
		
		&evaluation_quartet_support_distances(
			
			\$$sref_clan4											,	# actual 4clan code																																				IN: defined; OUT: unchanged 
			\$$sref_fc												,	# process number of the current quartet																															IN: defined; OUT: unchanged 
			\%$href_value_of_scr_par								,	# key1: sealion user available parameter options; value: optional parameter setting																				IN: defined; OUT: unchanged
			\%$href_hoh_sup_dist_of_clan4_type_of_treepair			,	# key1: clan4 code; key2: type of support ('nap', 'napr', or 'napv'); key3: tree-pair code ('first2nd' or 'first3rd'); value: distance Nap						IN: defined; OUT: unchanged 
			\%$href_hoh_qrej_of_filtertype_of_scale_of_fn			,	# key1: type of supportfilter ('nap_risk1', 'napr_risk1', or 'napv_risk1'); key2: scale value; key3: quartet number; value: filter corresponding score-value	IN: defined; OUT: changed
			\%$href_hol_scales_of_filtertype						,	# key1: type of supportfilter ('risk', 'qdist'); value: list of ftype corresponding threshold scales															IN: defined; OUT: unchanged
			\%$href_seen_filter_method								,	# key1: filter method incl OTU clipping approaches; value: seen +1																								IN: undefined; OUT: defined
		);
	}
	######################################################################
	######################################################################
	
	######################################################################################################
	######################################################################################################
	
	
	######################################################################################################
	############### print tsv table files (pattern, pairwise and best tree table)  #######################
	######################################################################################################
	
	######################################################################
	# print data to table TSV Q3                                       ###
	######################################################################
	# table includes all kinds of pattern related calculations, like Nc, Na, Np, CF, 
	# ratio calculations etc... as well as info about correctness for each of the three rooted quartet-trees
	&print_tsv_Q3(
		
		\%$href_subfolder_q 								,	# key: output-type (e.g. 'pdf'); value: pathdir																				IN: defined; OUT: unchanged
		\%$href_hoh_tsvdata_of_tsvfilename_of_printcode		,	# key1: R grafic unique printcode (e.g. 'metafileQ1'); key2: *.tsv filename; value: *tsv filename data addings 				IN: defined; OUT: unchanged
		\$$sref_clan4										,	# actual 4clan code																											IN: defined; OUT: unchanged 
		\$$sref_fc											,	# process number of the current quartet																						IN: defined; OUT: unchanged 
		\%sing_inf_of_type									,	# key: info type (.e.g 'sNs'); value_ type value																			IN: defined; OUT: unchanged
		\%sing_inf_of_type_of_tree							,	# key1: tree dependent info types (.e.g 'sNs_exp_from_qt1'); key2: rooted tree, value_ type value							IN: defined; OUT: unchanged
		\%invvar_inf_of_type_all_trees						,	# key: info type (.e.g 'sNs'); value_ type value																			IN: defined; OUT: unchanged
		\%invvar_inf_of_type_of_tree						,	# key1: tree dependent info types (.e.g 'sNs_exp_from_qt1'); key2: rooted tree, value_ type value							IN: defined; OUT: unchanged		
		\@$aref_rooted_qtrees								,	# list of rooted quartet trees (qtree1, qtree2, qtree3)																		IN: defined; OUT: unchanged 
		\%NaNpNc_inf_of_type_of_tree						,	# key1: tree dependent info types (.e.g 'Na', 'Np'...); key2: rooted tree, value_ type value								IN: defined; OUT: unchanged
		\%cf_of_type_of_qt									,	# key1: info type (e.g. ('cf_obs'); key2: rooted quartet-tree; value: correction factor for site patterns of actual tree	IN:	defined; OUT: unchanged
		\%finalsupport_of_type_of_qt						,	# key1: support type (e.g. 'nap', 'napv'...); key2: rooted quartet tree; value: score										IN:	defined; OUT: unchanged
		\%$href_value_of_scr_par							,	# key1: sealion user available parameter options; value: optional parameter setting											IN: defined; OUT: unchanged
	);
	######################################################################
	######################################################################
}

sub evaluation_quartet_ratio_NcNa{

	my	$sref_cl4_code							=	$_[0] ;	# actual 4clan code																																									IN: defined; OUT: unchanged 
	my	$sref_qn								=	$_[1] ;	# process number of the current quartet																																				IN: defined; OUT: unchanged 
	my	$href_NaNpNc_inf_of_type_of_tree		=	$_[2] ;	# key1: tree dependent info types (.e.g 'Na', 'Np'...); key2: rooted tree, value_ type value																						IN: defined, OUT: unchanged
	my	$href_value_of_scr_par					=	$_[3] ;	# key1: sealion user available parameter options; value: optional parameter setting																									IN: defined; OUT: unchanged
	my	$href_hoh_qrej_of_ft_of_scale_of_fn		=	$_[4] ;	# key1: 4clan code; key2: type of supportfilter ('nap_risk1', 'napr_risk1', or 'napv_risk1'); key3: scale; key4: quartet number: filter corresponding score-value					IN: undefined; OUT: defined 
	my	$sref_rNcNa_type						=	$_[5] ;	# type of Nc/Na ratio																																								IN: defined; OUT: unchanged
	my	$href_rNcNa_of_cmfm_of_fc_of_tr			=	$_[6] ;	# key1: quartet-number ($fc), key2: calculation+filter method (nap, napv, napr); key3: calculation method + filter method (nap_risk1, napv_dist2, napr_dist3); value: rNc/na		IN: undefined; OUT: defined 
	my	$href_alg_def_methods					=	$_[7] ;	# list of defined support methods (nap, napv ...)																																	IN: defined; OUT: unchanged
	my	$href_hol_scales_of_filtertype			=	$_[8] ;	# key1: type of supportfilter ('risk', 'qdist'); value: list of ftype corresponding threshold scales																				IN: defined; OUT: unchanged
	my	$href_seen_fmethod						=	$_[9] ;	# key1: filter method approaches; value: seen +1																																	IN: undefined; OUT: defined
	
	######################################################################
	# called from subroutine(s)...
	# &evaluation_quartet_pattern_main
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	## evaluate for each quartet-tree if ratio rNc/Na is in risk range ###
	######################################################################
	
	###################################
	# define calculation method depending filter codes
	my	%method_code_risk ;
	if ( $$sref_rNcNa_type eq 'rNcNa_org' ){ 
		
		if ( $href_value_of_scr_par->{seq_filterRisk}{risk1} == 1 ){
			
			for my $cm ( keys %{$href_alg_def_methods->{original}} ){ 
				
				$method_code_risk{$cm."_risk1"} 	= 1	; 
				$href_seen_fmethod->{$cm."_risk1"} 	= 1
			}
		}
		
		if ( $href_value_of_scr_par->{seq_filterRisk}{risk2} == 1 ){
			
			for my $cm ( keys %{$href_alg_def_methods->{original}} ){ 
				
				$method_code_risk{$cm."_risk2"} 	= 1	;
				$href_seen_fmethod->{$cm."_risk2"} 	= 1
			}
		}
	}
	
	else { die "\nBUG-ERROR: Cannot assign 'rNcNa'-flag in subroutine &evaluation_quartet_ratio_NcNa!\nPlease, report BUG to software developer!\n\n" }
	###################################
	
	
	###################################
	# Analyse RISK filtering for different threshold scales
	# necessary for RISK threshold optimization
	my	$treeflag = 0 ;
	my	$upper_limit = $href_value_of_scr_par->{seq_filterRisk}{turisk} ; 
	
	for my $tl_scale ( @{$href_hol_scales_of_filtertype->{risk}} ){
		
		###################################
		# evaluation of each of the three quartet-tree rNc/Na values due to risk range definition
		# - count number of rooted trees within risk range ($counter_rNcNa_risk_range)  
		# - count number of rooted trees with a rNc/Na score > risk range ($counter_rNc_Na_positive)
		# - count number of rooted trees with a rNc/Na score < risk range ($counter_rNc_Na_no_risk)
		my	$counter_rNcNa_risk_range	= 0 ;
		my	$counter_rNc_Na_positive 	= 0 ;
		my	$counter_rNc_Na_no_risk 	= 0 ;
		
		for my $tree ( keys %{$href_NaNpNc_inf_of_type_of_tree->{$$sref_rNcNa_type}} ){
			
			$treeflag++ ;
			my	$flag_evaluation ;
	
			if		( 	( $href_NaNpNc_inf_of_type_of_tree->{$$sref_rNcNa_type}{$tree} <  $upper_limit 	) && 
						( $href_NaNpNc_inf_of_type_of_tree->{$$sref_rNcNa_type}{$tree} >= $tl_scale 	) ){ 
						
							$counter_rNcNa_risk_range++	; $flag_evaluation = "risk range" ;
			}
			elsif	(	$href_NaNpNc_inf_of_type_of_tree->{$$sref_rNcNa_type}{$tree} <   $tl_scale	  ){ $counter_rNc_Na_no_risk++	; $flag_evaluation = "no risk" 			}
			elsif	(	$href_NaNpNc_inf_of_type_of_tree->{$$sref_rNcNa_type}{$tree} >=  $upper_limit ){ $counter_rNc_Na_positive++	; $flag_evaluation = "positive range" 	}
			
			###################################
			# sample rNcNa only once, not each time again, fiving a different scale threshold
			if ( $treeflag <= 3 ){
				
				for my $methodcode ( keys %method_code_risk ){ 
					
					$href_rNcNa_of_cmfm_of_fc_of_tr->{$$sref_qn}{$methodcode}{$tree} = $href_NaNpNc_inf_of_type_of_tree->{$$sref_rNcNa_type}{$tree}
				}
			}
			###################################
			
			###################################
			# test print
			if ( $href_value_of_scr_par->{sys_parameter}{screen_prt} == 1 ){
				
				print	"::&evaluation_quartet_ratio_NcNa of pot. ", $$sref_rNcNa_type, " risk-range(s)\n" ,
						"\n\t", $$sref_rNcNa_type, "\tlt-scale: ", $tl_scale ,"\t", $tree, "\t", $href_NaNpNc_inf_of_type_of_tree->{$$sref_rNcNa_type}{$tree}, "\t",  $flag_evaluation ; #exit;
			}
			###################################
		}
		###################################
		
		
		###################################
		# set quartet to RISK1 if either...
		# ...1) rNc/Na of two trees are in the risk range and the third tree resulting in rNc/Na > $upper_limit
		# ...2) all three trees result in rNc/Na > $lower_limit
		# allows 1 tree in risk range and two trees above risk range without setting risk1=1
		#my	$risktype_r1	= "no_risk" ;
		my	$risktype_r1	= "no_risk" ;
		
		if 		( ( $counter_rNcNa_risk_range >= 2 ) && ( $counter_rNc_Na_no_risk  == 0	) ){ $risktype_r1 = "risk" }
		elsif 	( 	$counter_rNc_Na_positive  == 3 										)  { $risktype_r1 = "risk" }
		###################################
	
		
		###################################
		# set quartet to RISK2 if either...
		# ...1) rNc/Na of all trees is >= $lower_limit
		# disables 1 tree in risk range and two trees above risk range without setting risk2=1
		#my	$risktype_r2	= "no_risk" ;
		my	$risktype_r2	= "no_risk" ;
		
		if 	( $counter_rNc_Na_no_risk == 0 ){ $risktype_r2 = "risk" }
		###################################
		
		
		###################################
		# assign risk value to a deeper going hash to evaluate 4clan risky quartets overall single quartet analyses of the 4clan
		for my $methodcode ( keys %method_code_risk ){ 
			
			if (	( $href_value_of_scr_par->{seq_filterRisk}{risk1} == 1 ) && ( $risktype_r1 eq "risk" ) && ( $methodcode =~ /risk1/ ) ){ 
				
					$href_hoh_qrej_of_ft_of_scale_of_fn->{$methodcode}{$tl_scale}{$$sref_qn} = 1 ; 
			}
			
			if (	( $href_value_of_scr_par->{seq_filterRisk}{risk2} == 1 ) && ( $risktype_r2 eq "risk" ) && ( $methodcode =~ /risk2/ ) ){ 
					
					$href_hoh_qrej_of_ft_of_scale_of_fn->{$methodcode}{$tl_scale}{$$sref_qn} = 1
			}
		}
		###################################
	}
	######################################################################
	
	#exit;
}

sub evaluation_quartet_singletons{

	my	$aref_rooted_qtrees					= $_[0] ; # list of rooted quartet trees 																								IN: defined; OUT: unchanged
	my	$href_Nsites_obs_of_pattern			= $_[1] ; # key: split-pattern code (A to V); value N observed pattern																	IN: defined; OUT: unchanged
	my	$href_Nsites_exp_qt1_of_pattern		= $_[2] ; # key: split-pattern code (A to V); value N expected given quartet-tree 1														IN: defined; OUT: unchanged
	my	$href_Nsites_exp_qt2_of_pattern		= $_[3] ; # key: split-pattern code (A to V); value N expected given quartet-tree 2														IN: defined; OUT: unchanged
	my	$href_Nsites_exp_qt3_of_pattern		= $_[4] ; # key: split-pattern code (A to V); value N expected given quartet-tree 3														IN: defined; OUT: unchanged
	my	$sref_tot_N_pattern_obs				= $_[5] ; # total number of (gap excluded) site positions																				IN: defined; OUT: unchanged
	my	$href_sing_inf_of_type				= $_[6] ; # key: info type (.e.g 'sNs'); value_ type value																				IN:undefined; OUT: defined
	my	$href_sing_inf_of_type_of_tree		= $_[7] ; # key1: tree dependent info types (.e.g 'sNs_exp_from_qt1'); key2: rooted tree, value_ type value								IN:undefined; OUT: defined
	my	$href_cf_of_type_of_qt				= $_[8] ; # key1: info type (e.g. ('cf_obs'); key2: rooted quartet-tree; value: correction factor for site patterns of actual tree		IN:undefined; OUT: defined
	my	$href_val_of_scr_par				= $_[9] ; # key1: sealion user available parameter options; value: optional parameter setting											IN: defined; OUT: unchanged
	
	######################################################################
	# called from...
	# &evaluation_quartet_pattern_main
	#
	############
	# using 
	# &calc_sum
	# &calc_cf
	# &calc_percent
	# &calc_mean
	######################################################################
	
	######################################################################
	## analysis of singleton frequencies (incl. CFexp and CFobs)        ##
	######################################################################
	
	###################################
	# sort alignment observed singleton pattern distribution ('Ns_obs') 
	my	@sorted_Ns_obs				=	sort {$a<=>$b} (
				
				$href_Nsites_obs_of_pattern->{K}, 
				$href_Nsites_obs_of_pattern->{L}, 
				$href_Nsites_obs_of_pattern->{M}, 
				$href_Nsites_obs_of_pattern->{N} 
	) ;
	###################################
	
	###################################
	# sort single tree (qt1, qt2, qt3) expected singleton pattern distribution (Ns_exp) 
	my	@sorted_Ns_exp_from_qt1		=	sort {$a<=>$b} ( 
	
				$href_Nsites_exp_qt1_of_pattern->{K} , 
				$href_Nsites_exp_qt1_of_pattern->{L} , 
				$href_Nsites_exp_qt1_of_pattern->{M} , 
				$href_Nsites_exp_qt1_of_pattern->{N}	
	) ;

	my	@sorted_Ns_exp_from_qt2		=	sort {$a<=>$b} ( 
				
				$href_Nsites_exp_qt2_of_pattern->{K} , 
				$href_Nsites_exp_qt2_of_pattern->{L} , 
				$href_Nsites_exp_qt2_of_pattern->{M} , 
				$href_Nsites_exp_qt2_of_pattern->{N}	
	) ;
	
	my	@sorted_Ns_exp_from_qt3		=	sort {$a<=>$b} ( 
				
				$href_Nsites_exp_qt3_of_pattern->{K} , 
				$href_Nsites_exp_qt3_of_pattern->{L} , 
				$href_Nsites_exp_qt3_of_pattern->{M} , 
				$href_Nsites_exp_qt3_of_pattern->{N}	
	) ;
	###################################
	
	###################################
	# assignment of alignment observed singleton calculations incl cf observed
	# lowest number of observed singleton-pattern 'sNs'
	$href_sing_inf_of_type->{sNs_obs}		=	$sorted_Ns_obs[0] ;
	
	# total number of observed singletons 'tNs'
	$href_sing_inf_of_type->{tNs_obs}		=	&calc_sum		( \@sorted_Ns_obs );
	
	# proportion of lowest singleton number in relation to all singletons
	$href_sing_inf_of_type->{sPs_obs}		=	&calc_percent	( $href_sing_inf_of_type->{sNs_obs}, $href_sing_inf_of_type->{tNs_obs} 	) ;
	
	# proportion of total N singletons in relation to all site positions of the actual quartet 'tPs_obs'
	$href_sing_inf_of_type->{tPs_obs}		=	&calc_percent	( $href_sing_inf_of_type->{tNs_obs}, $$sref_tot_N_pattern_obs 	) ;
	
	# proportion of observed lowest N of singleton pattern to total number of singletons 'rsNs_tNs_obs'
	$href_sing_inf_of_type->{rsNs_tNs_obs}	=	$href_sing_inf_of_type->{sNs_obs} / $href_sing_inf_of_type->{tNs_obs} ;
	###################################
	
	###################################
	# test print
	if ( $href_val_of_scr_par->{sys_parameter}{screen_prt} == 1 ){
		
		print "::&evaluation_quartet_singletons\n";
		
		for my $s (@sorted_Ns_obs){ print "\n\tsingleton distr observed\t", $s} #exit; 
		
		for my $t ( @sorted_Ns_exp_from_qt1 ){ print "\n\tsingleton distr expected from qt1\t", $t }
		
		print	"\n\tN shortest singleton observed (sNs_obs):\t" 				, $href_sing_inf_of_type->{sNs_obs}			,
				"\n\tN total singleton observed (tNs_obs):\t" 					, $href_sing_inf_of_type->{tNs_obs}			,
				"\n\tP shortest singleton observed (sPs_obs):\t" 				, $href_sing_inf_of_type->{sPs_obs}			,
				"\n\tP total singleton observed (tPs_obs):\t" 					, $href_sing_inf_of_type->{tPs_obs}			,
				"\n\tratio 'N_shortest/N_total' observed (rsNs_tNs_obs)\t"		, $href_sing_inf_of_type->{rsNs_tNs_obs}	;# exit;
	}
	###################################
	
	###################################
	# assignment of quartet-tree depending singleton calculations 'sNs_exp' and 'tNs_exp' (needed for following for-loop calculations)
	for my $tr_no ( 0 .. 2 ){
		
		# for given rooted tree assign lowest number of singleton pattern 'sNs_exp'
		if 		( $tr_no == 0 ){ $href_sing_inf_of_type_of_tree->{sNs_exp}{$aref_rooted_qtrees->[$tr_no]}	=	$sorted_Ns_exp_from_qt1[0] }
		elsif 	( $tr_no == 1 ){ $href_sing_inf_of_type_of_tree->{sNs_exp}{$aref_rooted_qtrees->[$tr_no]}	=	$sorted_Ns_exp_from_qt2[0] }
		elsif 	( $tr_no == 2 ){ $href_sing_inf_of_type_of_tree->{sNs_exp}{$aref_rooted_qtrees->[$tr_no]}	=	$sorted_Ns_exp_from_qt3[0] }
		
		# for given rooted tree calculate the total number of expected singletons 'tNs_exp'
		if 		( $tr_no == 0 ){ $href_sing_inf_of_type_of_tree->{tNs_exp}{$aref_rooted_qtrees->[$tr_no]}	=	&calc_sum( \@sorted_Ns_exp_from_qt1 ) }
		elsif 	( $tr_no == 1 ){ $href_sing_inf_of_type_of_tree->{tNs_exp}{$aref_rooted_qtrees->[$tr_no]}	=	&calc_sum( \@sorted_Ns_exp_from_qt2 ) }
		elsif 	( $tr_no == 2 ){ $href_sing_inf_of_type_of_tree->{tNs_exp}{$aref_rooted_qtrees->[$tr_no]}	=	&calc_sum( \@sorted_Ns_exp_from_qt3 ) }
		
		# test print
		if ( $href_val_of_scr_par->{sys_parameter}{screen_prt} == 1 ){
			
			print 	"\n\tN shortest singleton expected (sNs_exp) from ", $aref_rooted_qtrees->[$tr_no] ,"\t", $href_sing_inf_of_type_of_tree->{sNs_exp}{$aref_rooted_qtrees->[$tr_no]} ,
					"\n\tN total singleton expected (tNs_exp) from "   , $aref_rooted_qtrees->[$tr_no] ,"\t", $href_sing_inf_of_type_of_tree->{tNs_exp}{$aref_rooted_qtrees->[$tr_no]} ; #exit;
		}			
	}
	###################################
	
	###################################
	# singleton calculations for each quartet-tree
	for my $tr_no ( 0 .. 2 ){
		
		# calculation of cf correction factor for observed split pattern 'cf_obs'
		$href_cf_of_type_of_qt->{cf_obs}{$aref_rooted_qtrees->[$tr_no]} =	&calc_cf( $href_sing_inf_of_type->{sNs_obs}, $href_sing_inf_of_type->{tNs_obs} 	) ;
		
		# calculate correction factor expected for each tree, including variable sites in singletons 'cf_exp'
		my	@scores ;
		if 		( $tr_no == 0 ){ 
			
			$scores[0] = &calc_cf( $href_sing_inf_of_type_of_tree->{sNs_exp}{$aref_rooted_qtrees->[1]}, $href_sing_inf_of_type_of_tree->{tNs_exp}{$aref_rooted_qtrees->[1]} ) ;
			$scores[1] = &calc_cf( $href_sing_inf_of_type_of_tree->{sNs_exp}{$aref_rooted_qtrees->[2]}, $href_sing_inf_of_type_of_tree->{tNs_exp}{$aref_rooted_qtrees->[2]} ) ;
		}
		elsif 	( $tr_no == 1 ){ 
			
			$scores[0] = &calc_cf( $href_sing_inf_of_type_of_tree->{sNs_exp}{$aref_rooted_qtrees->[0]}, $href_sing_inf_of_type_of_tree->{tNs_exp}{$aref_rooted_qtrees->[0]} ) ;
			$scores[1] = &calc_cf( $href_sing_inf_of_type_of_tree->{sNs_exp}{$aref_rooted_qtrees->[2]}, $href_sing_inf_of_type_of_tree->{tNs_exp}{$aref_rooted_qtrees->[2]} ) ;
		}	
		elsif 	( $tr_no == 2 ){ 
			
			$scores[0] = &calc_cf( $href_sing_inf_of_type_of_tree->{sNs_exp}{$aref_rooted_qtrees->[0]}, $href_sing_inf_of_type_of_tree->{tNs_exp}{$aref_rooted_qtrees->[0]} ) ;
			$scores[1] = &calc_cf( $href_sing_inf_of_type_of_tree->{sNs_exp}{$aref_rooted_qtrees->[1]}, $href_sing_inf_of_type_of_tree->{tNs_exp}{$aref_rooted_qtrees->[1]} ) ;
		}
		$href_cf_of_type_of_qt->{cf_exp}{$aref_rooted_qtrees->[$tr_no]}							=	&calc_mean( \@scores ) ;
		
		# for given rooted tree calculate ratio of expected cf to observed cf 'rcf_exp_obs'
		$href_cf_of_type_of_qt->{rcf_ExpObs}{$aref_rooted_qtrees->[$tr_no]}						=	sprintf ( "%.4f", $href_cf_of_type_of_qt->{cf_exp}{$aref_rooted_qtrees->[$tr_no]} / $href_cf_of_type_of_qt->{cf_obs}{$aref_rooted_qtrees->[$tr_no]} );
		
		# for given rooted tree calculate proportion of lowest estimated N of singletons in relation to the total number of observed sites 'sPs_exp'
		$href_sing_inf_of_type_of_tree->{sPs_exp}{$aref_rooted_qtrees->[$tr_no]}				=	&calc_percent( $href_sing_inf_of_type_of_tree->{sNs_exp}{$aref_rooted_qtrees->[$tr_no]}, $href_sing_inf_of_type_of_tree->{tNs_exp}{$aref_rooted_qtrees->[$tr_no]} ) ;
		
		# for given rooted tree calculate the difference between lowest expected N singletons and lowest number of observed singletons 'dsNs_ExpObs'
		$href_sing_inf_of_type_of_tree->{dsNs_ExpObs}{$aref_rooted_qtrees->[$tr_no]}			=	$href_sing_inf_of_type_of_tree->{sNs_exp}{$aref_rooted_qtrees->[$tr_no]} - $href_sing_inf_of_type->{sNs_obs} ;
	
		# for given rooted tree calculate the difference between lowest expected P singletons and lowest P of observed singletons 'dsNs_ExpObs'
		$href_sing_inf_of_type_of_tree->{dsPs_ExpObs}{$aref_rooted_qtrees->[$tr_no]}			=	$href_sing_inf_of_type_of_tree->{sPs_exp}{$aref_rooted_qtrees->[$tr_no]} - $href_sing_inf_of_type->{sPs_obs} ;
		
		# for given rooted tree calculate the total proportion of expected singletons to the total nmber of observed quartet sites 'tPs_exp'
		$href_sing_inf_of_type_of_tree->{tPs_exp}{$aref_rooted_qtrees->[$tr_no]}				=	&calc_percent( $href_sing_inf_of_type_of_tree->{tNs_exp}{$aref_rooted_qtrees->[$tr_no]}, $$sref_tot_N_pattern_obs ) ;
		
		# for given rooted tree calculate the ratio of lowest N expected singletons to the total number of expected singletons 'sNs_exp'
		$href_sing_inf_of_type_of_tree->{rsNs_tNs_exp}{$aref_rooted_qtrees->[$tr_no]}			=	$href_sing_inf_of_type_of_tree->{sNs_exp}{$aref_rooted_qtrees->[$tr_no]} / $href_sing_inf_of_type_of_tree->{tNs_exp}{$aref_rooted_qtrees->[$tr_no]} ;
		
		# for given rooted tree calculate the distance of the ratio expected lowest N singletons/total number singletons and the total number of observed singletons 'drsNs_tNs_exp'
		$href_sing_inf_of_type_of_tree->{drsNstNs_ExpObs}{$aref_rooted_qtrees->[$tr_no]}		=	$href_sing_inf_of_type_of_tree->{rsNs_tNs_exp}{$aref_rooted_qtrees->[$tr_no]} - $href_sing_inf_of_type->{rsNs_tNs_obs} ;
		
		# for given rooted tree calculate the mean number of lowest convergent expected number of singletons (expected from the other two trees) 'sNs_exp_for'
		@scores = () ;
		if 		( $tr_no == 0 ){ @scores = ( $href_sing_inf_of_type_of_tree->{sNs_exp}{$aref_rooted_qtrees->[1]}, $href_sing_inf_of_type_of_tree->{sNs_exp}{$aref_rooted_qtrees->[2]} ) }
		elsif 	( $tr_no == 1 ){ @scores = ( $href_sing_inf_of_type_of_tree->{sNs_exp}{$aref_rooted_qtrees->[0]}, $href_sing_inf_of_type_of_tree->{sNs_exp}{$aref_rooted_qtrees->[2]} ) }
		elsif 	( $tr_no == 2 ){ @scores = ( $href_sing_inf_of_type_of_tree->{sNs_exp}{$aref_rooted_qtrees->[0]}, $href_sing_inf_of_type_of_tree->{sNs_exp}{$aref_rooted_qtrees->[1]} ) }
		$href_sing_inf_of_type_of_tree->{sNs_exp_for}{$aref_rooted_qtrees->[$tr_no]}			=	&calc_mean	 ( \@scores ) ;
		
		# for given rooted tree calculate the mean number of total convergent expected number of singletons (expected from the other two trees) 'tNs_exp_for'
		@scores = () ;
		if 		( $tr_no == 0 ){ @scores = ( $href_sing_inf_of_type_of_tree->{tNs_exp}{$aref_rooted_qtrees->[1]}, $href_sing_inf_of_type_of_tree->{tNs_exp}{$aref_rooted_qtrees->[2]} ) }
		elsif 	( $tr_no == 1 ){ @scores = ( $href_sing_inf_of_type_of_tree->{tNs_exp}{$aref_rooted_qtrees->[0]}, $href_sing_inf_of_type_of_tree->{tNs_exp}{$aref_rooted_qtrees->[2]} ) }
		elsif 	( $tr_no == 2 ){ @scores = ( $href_sing_inf_of_type_of_tree->{tNs_exp}{$aref_rooted_qtrees->[0]}, $href_sing_inf_of_type_of_tree->{tNs_exp}{$aref_rooted_qtrees->[1]} ) }
		$href_sing_inf_of_type_of_tree->{tNs_exp_for}{$aref_rooted_qtrees->[$tr_no]}			=	&calc_mean	 ( \@scores ) ;
		
		# for given rooted tree calculate the proportion of lowest convergent number of singletons to the total number of quartet observed sites expected from the other two trees 'sPs_exp_for'
		$href_sing_inf_of_type_of_tree->{sPs_exp_for}{$aref_rooted_qtrees->[$tr_no]}			=	&calc_percent( $href_sing_inf_of_type_of_tree->{sNs_exp_for}{$aref_rooted_qtrees->[$tr_no]}, $href_sing_inf_of_type_of_tree->{tNs_exp_for}{$aref_rooted_qtrees->[$tr_no]} ) ;
		
		# for given rooted tree calculate the expected proportion of total convergent number of singletons to the total number of quartet observed sites 'tPs_exp_for'
		$href_sing_inf_of_type_of_tree->{tPs_exp_for}{$aref_rooted_qtrees->[$tr_no]}			=	&calc_percent( $href_sing_inf_of_type_of_tree->{tNs_exp_for}{$aref_rooted_qtrees->[$tr_no]}, $$sref_tot_N_pattern_obs ) ;
		
		# for given rooted tree calculate the expected ratio of convergent expected lowest number of singletons to total number of singletons 'rsNs_tNs_exp_for'
		$href_sing_inf_of_type_of_tree->{rsNs_tNs_exp_for}{$aref_rooted_qtrees->[$tr_no]}		=	$href_sing_inf_of_type_of_tree->{sNs_exp_for}{$aref_rooted_qtrees->[$tr_no]} / $href_sing_inf_of_type_of_tree->{tNs_exp_for}{$aref_rooted_qtrees->[$tr_no]} ;
		
		# for given rooted tree calculate the distance of ratio sNs/tNs expected for given tree and ratio sNs/tNs of observed pattern 'drsNs_tNs_ExpObs_for'
		$href_sing_inf_of_type_of_tree->{drsNstNs_ExpObs_for}{$aref_rooted_qtrees->[$tr_no]}	=	$href_sing_inf_of_type_of_tree->{rsNs_tNs_exp_for}{$aref_rooted_qtrees->[$tr_no]} - $href_sing_inf_of_type->{rsNs_tNs_obs} ;
		
		# for given rooted tree calculate the difference between lowest number of singletons expected for actual tree to lowest number of observed singletons 'dsNs_ExpObs_for'
		$href_sing_inf_of_type_of_tree->{dsNs_ExpObs_for}{$aref_rooted_qtrees->[$tr_no]}		=	$href_sing_inf_of_type_of_tree->{sNs_exp_for}{$aref_rooted_qtrees->[$tr_no]} - $href_sing_inf_of_type->{sNs_obs} ;
		
		# for given rooted tree calculate the difference between lowest P of singletons expected for actual tree to lowest Pof observed singletons 'dsPs_ExpObs_for'
		$href_sing_inf_of_type_of_tree->{dsPs_ExpObs_for}{$aref_rooted_qtrees->[$tr_no]}		=	$href_sing_inf_of_type_of_tree->{sPs_exp_for}{$aref_rooted_qtrees->[$tr_no]} - $href_sing_inf_of_type->{sPs_obs} ;
	}
	###################################
	
	###################################
	# test print
	if ( $href_val_of_scr_par->{sys_parameter}{screen_prt} == 1 ){
		
		print	"\n\tcf observed all trees:\t"	, 	$href_cf_of_type_of_qt->{cf_obs}{$aref_rooted_qtrees->[0]} ,
				"\n\tcf expected for tree1:\t"	,	$href_cf_of_type_of_qt->{cf_exp}{$aref_rooted_qtrees->[0]} ,
				"\n\tcf expected for tree2:\t"	,	$href_cf_of_type_of_qt->{cf_exp}{$aref_rooted_qtrees->[1]} ,
				"\n\tcf expected for tree3:\t"	,	$href_cf_of_type_of_qt->{cf_exp}{$aref_rooted_qtrees->[2]} , "\n"
	}
	######################################################################
}

sub evaluation_quartet_support_distances{
	
	my	$sref_cl4_code											= $_[0] ;	# actual 4clan code																																					IN: defined; OUT: unchanged 
	my	$sref_qn												= $_[1] ;	# process number of the current quartet																																IN: defined; OUT: unchanged 
	my	$href_value_of_scr_par									= $_[2] ;	# key1: sealion user available parameter options; value: optional parameter setting																					IN: defined; OUT: unchanged
	my	$href_Nap_distance_of_clan4_type_of_treepair			= $_[3] ;	# key1: clan4 code; key2: type of support ('nap', 'napr', or 'napv'); key3: tree-pair code ('first2nd' or 'first3rd'); value: distance Nap							IN: defined; OUT: unchanged 
	my	$href_hoh_qrej_of_ft_of_scale_of_fn						= $_[4] ;	# key1: 4clan code; key2: type of supportfilter ('nap_risk1', 'napr_risk1', or 'napv_risk1'); key2: scale; key3: quartet number: filter corresponding score-value	IN: defined; OUT: changed
	my	$href_hol_scales_of_filtertype							= $_[5] ;	# key1: type of supportfilter ('risk', 'qdist'); value: list of ftype corresponding threshold scales																IN: defined; OUT: unchanged
	my	$href_seen_fmethod										= $_[6] ;	# key1: filter method approaches; value: seen +1																													IN: undefined; OUT: defined
	
	#print "::&evaluation_quartet_support_distances\n";
	
	######################################################################
	# called from...
	# &evaluation_quartet_pattern_main
	#
	############
	# using 
	# ...
	######################################################################
	
	######################################################################
	## analysing tree support distances of gicen quartet combination    ##
	######################################################################
	
	###################################
	# for each quartet method (nap, napv ...), compare optionally specified support distances
	# and code each quartet with less/equal supprt distance as the distance-method (qdist1, qdist2, or qdist3) specified threshold to 1
	for my $tmethod ( sort keys %{$href_Nap_distance_of_clan4_type_of_treepair->{$$sref_cl4_code}} ){
		
		for my $tl_scale ( @{$href_hol_scales_of_filtertype->{qdist}} ){
			
			if ( $href_Nap_distance_of_clan4_type_of_treepair->{$$sref_cl4_code}{$tmethod}{$$sref_qn}{first2nd} <= $tl_scale ){
				
				my	$fmethod = $tmethod."_qdist1" ;
				$href_seen_fmethod->{$fmethod} = 1 ;
				
				$href_hoh_qrej_of_ft_of_scale_of_fn->{$fmethod}{$tl_scale}{$$sref_qn} = 1
			}
		}
	} 
	######################################################################
}

sub evaluation_quartet_support_napX{
	
	my	$sref_clan4_combination							= $_[0]  ; # actual 4clan code																																IN: defined; OUT: unchanged 
	my	$sref_quartet_number							= $_[1]  ; # process number of the current quartet																											IN: defined; OUT: unchanged 
	my	$aref_signal_code								= $_[2]  ; # array of hashkeys used for specific support calculation (1st: support code; 2nd: pot apomorph signal code; 3rd expected convergency code)		IN: defined; OUT: unchanged 
	my	$aref_rooted_qtrees								= $_[3]  ; # list of rooted quartet trees (qtree1, qtree2, qtree3)																							IN: defined; OUT: unchanged 
	my	$href_NaNpNc_inf_of_type_of_tree				= $_[4]  ; # key1: tree dependent info types (.e.g 'Na', 'Np'...); key2: rooted tree, value_ type value														IN: defined, OUT: unchanged
	my	$href_hol_Nap_of_clan4_of_calc_of_qtree			= $_[5]  ; # key1: 4-clan code; key2: type NAP; key3: rooted quartet-tree; value: list of single NAP quartet support										IN: defined; OUT: changed
	my	$href_hoh_Nap_of_clan4_of_calc_of_qtree_of_fc	= $_[6]  ; # key1: 4-clan code; key2: type NAP; key3: rooted quartet-tree; key4: quartet number; value:NAP quartet support									IN: defined; OUT: changed
	my	$href_finalsupport_of_type_of_qt				= $_[7]  ; # key1: support type (e.g. 'nap', 'napv'...); key2: rooted quartet tree; value: score															IN: undefined; OUT: defined 
	my	$href_Nap_distance_of_clan4_type_of_treepair	= $_[8]  ; # key1: clan4 code; key2: type of support ('nap', 'napr', or 'napv'); key3: tree-pair code ('first2nd' or 'first3rd'); value: distance Nap		IN: undefined; OUT: defined 
	my	$href_bestTr_of_fc_of_cm						= $_[9]  ; # key1: quartet-number ($fc), key2: calculation method (nap, napv, napr); value: quartet related best supported tree								IN: undefined; OUT: defined 
	my	$href_val_of_scr_par							= $_[10] ; # key1: sealion user available parameter options; value: optional parameter setting																IN: defined; OUT: unchanged
	
	######################################################################
	# called from...
	# &evaluation_quartet_pattern_main
	#
	############
	# using 
	# &norm_values
	######################################################################
	
	######################################################################
	## calculate final quartet-tree support and support distances       ##
	######################################################################
	
	###################################
	# Nap = Na_red - Nc_red (final, raw signal result, can be negative)
	my	$cNs		=	$aref_signal_code->[0]	 ;	# code for final support type (e.g. 'nap', 'nav'...)
	my	$cNa		=	$aref_signal_code->[1]	 ;	# code for potent. apomorph split signal (e.g. 'Na_red_obs', 'Nav_red_obs'...)
	my	$cNc		=	$aref_signal_code->[2]	 ;	# code for potent. convergent split-signal (e.g. 'Nc_red_exp', 'Ncv_red_exp'...)
	my	$code_n0	=	$cNs."_n0" ;
	
	$href_finalsupport_of_type_of_qt->{$code_n0}{$aref_rooted_qtrees->[0]}	=	$href_NaNpNc_inf_of_type_of_tree->{$cNa}{$aref_rooted_qtrees->[0]} - $href_NaNpNc_inf_of_type_of_tree->{$cNc}{$aref_rooted_qtrees->[0]}; 
	$href_finalsupport_of_type_of_qt->{$code_n0}{$aref_rooted_qtrees->[1]}	=	$href_NaNpNc_inf_of_type_of_tree->{$cNa}{$aref_rooted_qtrees->[1]} - $href_NaNpNc_inf_of_type_of_tree->{$cNc}{$aref_rooted_qtrees->[1]}; 
	$href_finalsupport_of_type_of_qt->{$code_n0}{$aref_rooted_qtrees->[2]}	=	$href_NaNpNc_inf_of_type_of_tree->{$cNa}{$aref_rooted_qtrees->[2]} - $href_NaNpNc_inf_of_type_of_tree->{$cNc}{$aref_rooted_qtrees->[2]}; 
	
	my	$cNa_rev	=	$cNa."_rev";
	my	$cNc_rev	=	$cNc."_rev";
	
	my	%final_support_reverse_of_qt ;
		$final_support_reverse_of_qt{$aref_rooted_qtrees->[0]}	=	$href_NaNpNc_inf_of_type_of_tree->{$cNa_rev}{$aref_rooted_qtrees->[0]} - $href_NaNpNc_inf_of_type_of_tree->{$cNc_rev}{$aref_rooted_qtrees->[0]};
		$final_support_reverse_of_qt{$aref_rooted_qtrees->[1]}	=	$href_NaNpNc_inf_of_type_of_tree->{$cNa_rev}{$aref_rooted_qtrees->[1]} - $href_NaNpNc_inf_of_type_of_tree->{$cNc_rev}{$aref_rooted_qtrees->[1]};
		$final_support_reverse_of_qt{$aref_rooted_qtrees->[2]}	=	$href_NaNpNc_inf_of_type_of_tree->{$cNa_rev}{$aref_rooted_qtrees->[2]} - $href_NaNpNc_inf_of_type_of_tree->{$cNc_rev}{$aref_rooted_qtrees->[2]};
	
	#############
	# test print
	if ( $href_val_of_scr_par->{sys_parameter}{screen_prt} == 1 ){
		
		print 	"::&evaluation_quartet_support_napX\n",
				"\n\t", $cNs, " N0\t"			, $aref_rooted_qtrees->[0], "\t", $href_finalsupport_of_type_of_qt->{$code_n0}{$aref_rooted_qtrees->[0]}	, 
				"\n\t", $cNs, " N0\t"			, $aref_rooted_qtrees->[1], "\t", $href_finalsupport_of_type_of_qt->{$code_n0}{$aref_rooted_qtrees->[1]}	, 
				"\n\t", $cNs, " N0\t"			, $aref_rooted_qtrees->[2], "\t", $href_finalsupport_of_type_of_qt->{$code_n0}{$aref_rooted_qtrees->[2]}	,  "\n" ,
				"\n\t", $cNs, " N0 reverse\t"	, $aref_rooted_qtrees->[0], "\t", $final_support_reverse_of_qt{$aref_rooted_qtrees->[0]}					, 
				"\n\t", $cNs, " N0 reverse\t"	, $aref_rooted_qtrees->[1], "\t", $final_support_reverse_of_qt{$aref_rooted_qtrees->[1]}					, 
				"\n\t", $cNs, " N0 reverse\t"	, $aref_rooted_qtrees->[2], "\t", $final_support_reverse_of_qt{$aref_rooted_qtrees->[2]}					,  "\n" ;
	}
	###################################
	
	
	###################################
	# Normaliation step 1: set lowest rooted tree signal to zero (reverse rooted trees incl)
	# identification of lowest rooted tree signal
	# adaption of remaining two rooted tree scores due to lowest score change
	my	$code_n1	=	$cNs."_n1" ;
	
	my	@sorted_Nap	=	sort {$a<=>$b} ( 
		
		$href_finalsupport_of_type_of_qt->{$code_n0}{$aref_rooted_qtrees->[0]}	, 
		$href_finalsupport_of_type_of_qt->{$code_n0}{$aref_rooted_qtrees->[1]}	, 
		$href_finalsupport_of_type_of_qt->{$code_n0}{$aref_rooted_qtrees->[2]}	, 
		$final_support_reverse_of_qt{$aref_rooted_qtrees->[0]}					,
		$final_support_reverse_of_qt{$aref_rooted_qtrees->[1]}					,
		$final_support_reverse_of_qt{$aref_rooted_qtrees->[2]}					,
	) ;
	
	
	if	( $sorted_Nap[0] < 1 ){ 
		
		$href_finalsupport_of_type_of_qt->{$code_n1}{$aref_rooted_qtrees->[0]}	=	sprintf "%.4f", $href_finalsupport_of_type_of_qt->{$code_n0}{$aref_rooted_qtrees->[0]}	+ abs($sorted_Nap[0]) ;
		$href_finalsupport_of_type_of_qt->{$code_n1}{$aref_rooted_qtrees->[1]}	=	sprintf "%.4f", $href_finalsupport_of_type_of_qt->{$code_n0}{$aref_rooted_qtrees->[1]}	+ abs($sorted_Nap[0]) ;
		$href_finalsupport_of_type_of_qt->{$code_n1}{$aref_rooted_qtrees->[2]}	=	sprintf "%.4f", $href_finalsupport_of_type_of_qt->{$code_n0}{$aref_rooted_qtrees->[2]}	+ abs($sorted_Nap[0]) ;
	}
	else{
		
		$href_finalsupport_of_type_of_qt->{$code_n1}{$aref_rooted_qtrees->[0]}	=	sprintf "%.4f", $href_finalsupport_of_type_of_qt->{$code_n0}{$aref_rooted_qtrees->[0]}	- abs($sorted_Nap[0]) ;
		$href_finalsupport_of_type_of_qt->{$code_n1}{$aref_rooted_qtrees->[1]}	=	sprintf "%.4f", $href_finalsupport_of_type_of_qt->{$code_n0}{$aref_rooted_qtrees->[1]}	- abs($sorted_Nap[0]) ;
		$href_finalsupport_of_type_of_qt->{$code_n1}{$aref_rooted_qtrees->[2]}	=	sprintf "%.4f", $href_finalsupport_of_type_of_qt->{$code_n0}{$aref_rooted_qtrees->[2]}	- abs($sorted_Nap[0]) ;
	}
	
	#############
	# test print
	if ( $href_val_of_scr_par->{sys_parameter}{screen_prt} == 1 ){
		
		print 	"\n\t", $cNs, " N1\t", $aref_rooted_qtrees->[0], "\t", $href_finalsupport_of_type_of_qt->{$code_n1}{$aref_rooted_qtrees->[0]} , 
				"\n\t", $cNs, " N1\t", $aref_rooted_qtrees->[1], "\t", $href_finalsupport_of_type_of_qt->{$code_n1}{$aref_rooted_qtrees->[1]} , 
				"\n\t", $cNs, " N1\t", $aref_rooted_qtrees->[2], "\t", $href_finalsupport_of_type_of_qt->{$code_n1}{$aref_rooted_qtrees->[2]} ,  "\n\n" ;
	}
	###################################
	
	
	###################################
	# Normaliation step 2: set norm1 values in relation to each other (norm2 support is between 0 and 1)
	# Norm2 values are stored in list @Nap_norm2_values using the same order as in @Nap_norm1_values
	# Assign norm2 values as value to corresponding tree key in %nap2_of_qt
	# Sort tree keys with corresponding norm2 values with descending order in @ranked_trees
	my	$code_n2	=	$cNs."_n2" ;
	
	my	@Nap_norm1_values	=	( 
		
		$href_finalsupport_of_type_of_qt->{$code_n1}{$aref_rooted_qtrees->[0]} , 
		$href_finalsupport_of_type_of_qt->{$code_n1}{$aref_rooted_qtrees->[1]} ,
		$href_finalsupport_of_type_of_qt->{$code_n1}{$aref_rooted_qtrees->[2]} 
	) ;
	
	
	my	@Nap_norm2_values	=	&norm_values( \@Nap_norm1_values ) ; #for (@Nap_norm2_values){print "\nnorm value\t", $_}
	
	$href_finalsupport_of_type_of_qt->{$code_n2}{$aref_rooted_qtrees->[0]}	=	$Nap_norm2_values[0]	;
	$href_finalsupport_of_type_of_qt->{$code_n2}{$aref_rooted_qtrees->[1]}	=	$Nap_norm2_values[1]	;
	$href_finalsupport_of_type_of_qt->{$code_n2}{$aref_rooted_qtrees->[2]}	=	$Nap_norm2_values[2]	;
		
	push @{$href_hol_Nap_of_clan4_of_calc_of_qtree->{$$sref_clan4_combination}{$cNs}{$aref_rooted_qtrees->[0]}}, $Nap_norm2_values[0]	;
	push @{$href_hol_Nap_of_clan4_of_calc_of_qtree->{$$sref_clan4_combination}{$cNs}{$aref_rooted_qtrees->[1]}}, $Nap_norm2_values[1]	;
	push @{$href_hol_Nap_of_clan4_of_calc_of_qtree->{$$sref_clan4_combination}{$cNs}{$aref_rooted_qtrees->[2]}}, $Nap_norm2_values[2]	;
	
	my %napN2_of_tree = ( $aref_rooted_qtrees->[0] => $Nap_norm2_values[0], $aref_rooted_qtrees->[1] => $Nap_norm2_values[1], $aref_rooted_qtrees->[2] => $Nap_norm2_values[2] ); 
	$href_hoh_Nap_of_clan4_of_calc_of_qtree_of_fc->{$$sref_clan4_combination}{$cNs}{$aref_rooted_qtrees->[0]}{$$sref_quartet_number} = $Nap_norm2_values[0]	;
	$href_hoh_Nap_of_clan4_of_calc_of_qtree_of_fc->{$$sref_clan4_combination}{$cNs}{$aref_rooted_qtrees->[1]}{$$sref_quartet_number} = $Nap_norm2_values[1]	;
	$href_hoh_Nap_of_clan4_of_calc_of_qtree_of_fc->{$$sref_clan4_combination}{$cNs}{$aref_rooted_qtrees->[2]}{$$sref_quartet_number} = $Nap_norm2_values[2]	;
	
	#############
	# test print
	if ( $href_val_of_scr_par->{sys_parameter}{screen_prt} == 1 ){
	
		print	"\n\t", $cNs, " N2\t", $aref_rooted_qtrees->[0], "\t", $Nap_norm2_values[0]	,
				"\n\t", $cNs, " N2\t", $aref_rooted_qtrees->[1], "\t", $Nap_norm2_values[1]	,
				"\n\t", $cNs, " N2\t", $aref_rooted_qtrees->[2], "\t", $Nap_norm2_values[2]	, "\n"	;
	}
	###################################
	
	
	###################################
	# rank rooted quartet-trees in @ranked_trees in order of strongest to lowest support
	my	$code_1stbesttree = $cNs."_1stbestTree" ;
	my	$code_2ndbesttree = $cNs."_2ndbestTree" ;
	my	$code_3rdbesttree = $cNs."_3rdbestTree" ;
	
	my	@ranked_trees	=	sort { $napN2_of_tree{$b} <=> $napN2_of_tree{$a} } keys %napN2_of_tree ;
	
	$href_finalsupport_of_type_of_qt->{$code_1stbesttree} = $ranked_trees[0] ;
	$href_finalsupport_of_type_of_qt->{$code_2ndbesttree} = $ranked_trees[1] ;
	$href_finalsupport_of_type_of_qt->{$code_3rdbesttree} = $ranked_trees[2] ;
	
	#############
	# test print
	if ( $href_val_of_scr_par->{sys_parameter}{screen_prt} == 1 ){
		
		print 	"\n\t", $cNs, " N2 ranked 1st\t", $ranked_trees[0], "\t", $napN2_of_tree{$ranked_trees[0]} , 
				"\n\t", $cNs, " N2 ranked 2nd\t", $ranked_trees[1], "\t", $napN2_of_tree{$ranked_trees[1]} , 
				"\n\t", $cNs, " N2 ranked 3rd\t", $ranked_trees[2], "\t", $napN2_of_tree{$ranked_trees[2]} , "\n" ;
	}
	#############
	
	$href_bestTr_of_fc_of_cm->{$$sref_quartet_number}{$cNs} = $ranked_trees[0].":".$napN2_of_tree{$ranked_trees[0]} ;
	###################################
	
	
	###################################
	# calculate final support distance between best and 2nd best as well as between best and 3rd best tree
	$href_Nap_distance_of_clan4_type_of_treepair->{$$sref_clan4_combination}{$cNs}{$$sref_quartet_number}{first2nd} 	= $napN2_of_tree{$ranked_trees[0]} - $napN2_of_tree{$ranked_trees[1]} ;
	$href_Nap_distance_of_clan4_type_of_treepair->{$$sref_clan4_combination}{$cNs}{$$sref_quartet_number}{first3rd}		= $napN2_of_tree{$ranked_trees[0]} - $napN2_of_tree{$ranked_trees[2]} ;
	$href_Nap_distance_of_clan4_type_of_treepair->{$$sref_clan4_combination}{$cNs}{$$sref_quartet_number}{second3rd} 	= $napN2_of_tree{$ranked_trees[1]} - $napN2_of_tree{$ranked_trees[2]} ;
	
	my	$code_1st2nd	= $cNs."_d1st2nd" ;
	my	$code_1st3rd	= $cNs."_d1st3rd" ;
	my	$code_2nd3rd	= $cNs."_d2nd3rd" ;
	
	$href_finalsupport_of_type_of_qt->{$code_1st2nd} = $href_Nap_distance_of_clan4_type_of_treepair->{$$sref_clan4_combination}{$cNs}{$$sref_quartet_number}{first2nd} 	;
	$href_finalsupport_of_type_of_qt->{$code_1st3rd} = $href_Nap_distance_of_clan4_type_of_treepair->{$$sref_clan4_combination}{$cNs}{$$sref_quartet_number}{first3rd} 	;
	$href_finalsupport_of_type_of_qt->{$code_2nd3rd} = $href_Nap_distance_of_clan4_type_of_treepair->{$$sref_clan4_combination}{$cNs}{$$sref_quartet_number}{second3rd} ;
	
	#############
	# test print
	if ( $href_val_of_scr_par->{sys_parameter}{screen_prt} == 1 ){
		
		print	"\n\t", $code_1st2nd, "\t", $href_Nap_distance_of_clan4_type_of_treepair->{$$sref_clan4_combination}{$cNs}{$$sref_quartet_number}{first2nd} ,
				"\n\t", $code_1st3rd, "\t", $href_Nap_distance_of_clan4_type_of_treepair->{$$sref_clan4_combination}{$cNs}{$$sref_quartet_number}{first3rd} , "\n" ;
	}
	######################################################################
	######################################################################
}

sub extract_4taxon_combination_from_fileline{

	my	$sref_line_in					= $_[0]	; # fileline with 4-taxon info (must be the first line in '*_detailed_split_calc_*'							IN: defined; OUT: unchanged
	my	$aref_4taxa						= $_[1]	; # list of quartet depending taxa in the same order as in the infile										IN: undefined; OUT: defined
	my	$href_taxa4_of_fc				= $_[2]	; # key1	: quartet-number ($fc); value : quartet assigned taxa (alphabetically ordered, joined by ':'	IN: undefined; OUT: defined
	my	$sref_fc						= $_[3]	; # actual quartet-number																					IN: defined; OUT: unchanged
	my	$href_hoh_seen_q_of_t_of_clan4	= $_[4]	; # key1: 4clan combination code; key2: taxon; value: N analysed quartets									IN: undefined; OUT: defined
	my	$sref_clan4						= $_[5]	; # Clan4 of given quartet combination																		IN: defined; OUT: unchanged
	my	$href_clan_of_seq				= $_[6]	; # hash with sequence name (key) assigned clan (value)														IN: defined; OUT: unchanged
	my	$href_value_of_scr_par			= $_[7]	; # key1: sealion user available parameter options; value: optional parameter setting						IN: defined; OUT: unchanged
	my	$sref_p_qcalc_file				= $_[8] ; # actual spd infile																						IN: defined; OUT: unchanged
	
	######################################################################
	## extracting taxa from the 1st line of '*_detailed_split_calc_*'   ##
	######################################################################
	#print 	"::&extract_4taxon_combination_from_fileline:\n\n",
	#		"\t<", $$sref_line_in , ">\n\n" ; 
			
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	#########################################
	# check for correct line structure
	unless ( $$sref_line_in =~ /^Q\d+:\t\w+:\w+:\w+:\w+/ ){ die "\nSPDFILE-ERROR: Cannot identify OTUs from first line of ", $$sref_p_qcalc_file, "!\n\n"  }
	#########################################
	
	#########################################
	# delete quartetnumber at line start (e.g. Q76\t)
	( my $taxa4	= $$sref_line_in ) =~ s/^Q\d+:\t// ;
	#print "taxa4 combination: ", $taxa4, "\n" ;
	#########################################
	
	#########################################
	# assign set of taxa to corresponding quartet-number
	$href_taxa4_of_fc->{$$sref_fc} = $taxa4 ;
	#########################################
	
	#########################################
	# store each taxon of 4taxa in @ptaxa
	@$aref_4taxa  = split ":", $taxa4 ;
	for my $t ( @$aref_4taxa ){ 
		
		$href_hoh_seen_q_of_t_of_clan4->{$$sref_clan4}{$t}++ ;
		
		if ( $href_value_of_scr_par->{sys_parameter}{screen_prt} == 1 ){ print "\tfound taxon: ", $t, " ", $href_clan_of_seq->{$t}, "\n" }
	} 
	#########################################
	#exit;
	######################################################################
}

sub filter_quartets_dist{
	
	my	$sref_clan4												=	$_[0]  ;	# clan4 combination of given quartet combination																																	IN: defined; OUT: unchanged
	my	$href_hol_sup_of_clan4_of_calc_of_qtree					=	$_[1]  ;	# multiple hash, listing for each of the three rooted quartet-trees of a given 4clan combination all single quartet inferred support (NAP or NAPV)									IN: defined; OUT: changed
	my	$href_hoh_sup_of_clan4_of_calc_of_qtree_of_fc			=	$_[2]  ;	# multiple hash, storing the final support (nap or napv) related to each rooted quartet-tree																						IN: defined; OUT: changed
	my	$href_hoh_seen_rejected_of_clan4_of_filtertype_of_fn	=	$_[3]  ;	# key1: 4clan code; key2: type of supportfilter ('nap_risk1', 'napr_risk1', or 'napv_risk1'); key3: actual sealion quartet number; value: filter corresponding score-value			IN: defined; OUT: unchanged 
	my	$href_value_of_scr_par									=	$_[4]  ;	# key1: sealion user available parameter options; value: optional parameter setting																 									IN: defined; OUT: unchanged
	my	$href_specified_support_method_code						=	$_[5]  ;	# list of defined support methods (nap, napv ...)																																	IN: defined; OUT: changed	
	my	$href_hoh_seen_filter_of_type							=	$_[6]  ;	# key1: filter ('dist' or 'risk'); key2: combi-code method+filter ('nap_dist1', 'nap_dist2'...); value: number of filter passing quartets											IN: defined; OUT: changed
	my	$href_qn_of_qfiltered_otu_of_4clan_of_ft				=	$_[7]  ;	# key1: 4seq filtermethod (e.g. nap_risk1); key2: taxon; key3: 4clan|total; value: list of rejected qn's																			IN: undefined; OUT: defined
	my	$href_hoh_qremaining_of_otu_of_4clan_of_ft				=	$_[8]  ;	# key1: 4seq filtermethod (e.g. nap_risk1); key2: taxon; key3: 4clan|total; value: N remaining quartets																				IN: undefined; OUT: defined
	my	$href_taxa4_of_fc										=	$_[9]  ;	# key1: quartet-number ($fc); value : quartet assigned taxa (alphabetically ordered, joined by ':')																					IN: defined; OUT: unchanged
	my	$href_hoh_N_qrej_of_clan4_of_ft							=	$_[10] ;	# key1: 4seq filtermethod (e.g. nap_risk1); key2: 4clan combination code; value: N rejected quartets*3 (counted three times, each time for each quartet tree)						IN: undefined; OUT: defined
	my	$href_seen_remain_of_cmfm_of_fc							=	$_[11] ;	# key1: calculation method + filter method (nap_risk1, napv_dist2, napr_dist3); key2: quartet-number ($fc); value: seen +1															IN: undefined; OUT: defined
	my	$href_seen_filter_method								=	$_[12] ;	# key1: filter method incl OTU clipping approaches; value: seen +1																													IN: undefined; OUT: defined
	
	#print "::&filter_quartets_dist\n";
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# ...
	######################################################################
	
	######################################################################
	# filter quartets for each quartet calculation method                # 
	######################################################################
	for	my $cmethod ( keys %{$href_specified_support_method_code->{original}} ){
		
		############################
		# for each tree of given calculation method...
		for my $tree ( sort keys %{$href_hoh_sup_of_clan4_of_calc_of_qtree_of_fc->{$$sref_clan4}{$cmethod}} ){
			
			############################
			# go through all quartet numbers... 
			for my $qn ( sort {$a<=>$b} keys %{$href_hoh_sup_of_clan4_of_calc_of_qtree_of_fc->{$$sref_clan4}{$cmethod}{$tree}} ){
				
				my	@taxa = split ":", $href_taxa4_of_fc->{$qn} ;
				
				############################
				# filter quartets due to qdist1
				if ( $href_value_of_scr_par->{seq_filterDist}{qdist1} == 1 ) {
					
					my	$fm = $cmethod."_qdist1" ;
					
					unless ( $href_hoh_seen_rejected_of_clan4_of_filtertype_of_fn->{$$sref_clan4}{$fm}{$qn} ){#print "qn risk1 ", $qn, "\t" , $tree, "\n";
						
						push @{$href_hol_sup_of_clan4_of_calc_of_qtree->{$$sref_clan4}{$fm}{$tree}} , $href_hoh_sup_of_clan4_of_calc_of_qtree_of_fc->{$$sref_clan4}{$cmethod}{$tree}{$qn} ;
						
						for my $t ( @taxa ){ $href_hoh_qremaining_of_otu_of_4clan_of_ft->{$fm}{$t}{$$sref_clan4}++ }
						
						$href_seen_remain_of_cmfm_of_fc->{$fm}{$qn}++
					}
					
					else{ 
						
						for my $t ( @taxa ){ push @{$href_qn_of_qfiltered_otu_of_4clan_of_ft->{$fm}{$t}{$$sref_clan4}}, $qn }  
						
						$href_hoh_N_qrej_of_clan4_of_ft->{$fm}{$$sref_clan4}++ ;
					}
					
					$href_hoh_seen_filter_of_type->{dist}{$fm}++ ;
				}
				############################
			}
			############################
		}
		############################
	} #exit;
	######################################################################
}

sub filter_quartets_risk{
	
	my	$sref_clan4												=	$_[0]  ;	# clan4 combination of given quartet combination																																	IN: defined; OUT: unchanged
	my	$href_hol_sup_of_clan4_of_calc_of_qtree					=	$_[1]  ;	# multiple hash, listing for each of the three rooted quartet-trees of a given 4clan combination all single quartet inferred support (NAP or NAPV)									IN: defined; OUT: changed
	my	$href_hoh_sup_of_clan4_of_calc_of_qtree_of_fc			=	$_[2]  ;	# multiple hash, storing the final support (nap or napv) related to each rooted quartet-tree																						IN: defined; OUT: changed
	my	$href_hoh_seen_rejected_of_clan4_of_filtertype_of_fn	=	$_[3]  ;	# key1: 4clan code; key2: type of supportfilter ('nap_risk1', 'napr_risk1', or 'napv_risk1'); key3: actual sealion quartet number; value: filter corresponding score-value			IN: defined; OUT: unchanged 
	my	$href_value_of_scr_par									=	$_[4]  ;	# key1: sealion user available parameter options; value: optional parameter setting																 									IN: defined; OUT: unchanged
	my	$href_specified_support_method_code						=	$_[5]  ;	# list of defined support methods (nap, napv ...)																																	IN: defined; OUT: unchanged	
	my	$href_hoh_seen_filter_of_type							=	$_[6]  ;	# key1: filter ('dist' or 'risk'); key2: combi-code method+filter ('nap_risk1', 'napv_risk2'...); value: number of filter passing quartets											IN: undefined; OUT: defined
	my	$href_qn_of_qfiltered_otu_of_4clan_of_ft				=	$_[7]  ;	# key1: 4seq filtermethod (e.g. nap_risk1); key2: taxon; key3: 4clan|total; value: list of rejected qn's																			IN: undefined; OUT: defined
	my	$href_hoh_qremaining_of_otu_of_4clan_of_ft				=	$_[8]  ;	# key1: 4seq filtermethod (e.g. nap_risk1); key2: taxon; key3: 4clan|total; value: N remaining quartets																				IN: undefined; OUT: defined
	my	$href_taxa4_of_fc										=	$_[9]  ;	# key1: quartet-number ($fc); value : quartet assigned taxa (alphabetically ordered, joined by ':')																					IN: defined; OUT: unchanged
	my	$href_hoh_N_qrej_of_clan4_of_ft							=	$_[10] ;	# key1: 4seq filtermethod (e.g. nap_risk1); key2: 4clan combination code; value: N rejected quartets*3 (counted three times, each time for each quartet tree)						IN: undefined; OUT: defined
	my	$href_seen_remain_of_cmfm_of_fc							=	$_[11] ;	# key1: calculation method + filter method (nap_risk1, napv_dist2, napr_dist3); key2: quartet-number ($fc); value: seen +1															IN: undefined; OUT: defined
	
	#print "::&filter_quartets_risk\n";
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# ...
	######################################################################
	
	######################################################################
	# filter quartets for each quartet calculation method                # 
	######################################################################
	for	my $cmethod ( keys %{$href_specified_support_method_code->{original}} ){
		
		############################
		# for each tree of given calculation method...
		for my $tree ( sort keys %{$href_hoh_sup_of_clan4_of_calc_of_qtree_of_fc->{$$sref_clan4}{$cmethod}} ){
			
			############################
			# go through all quartet numbers... 
			for my $qn ( sort {$a<=>$b} keys %{$href_hoh_sup_of_clan4_of_calc_of_qtree_of_fc->{$$sref_clan4}{$cmethod}{$tree}} ){
				
				my	@taxa = split ":", $href_taxa4_of_fc->{$qn} ;
				
				############################
				# filter quartets due to risk1 nap (uncorrected)
				if (( $href_value_of_scr_par->{seq_filterRisk}{risk1} == 1 ) &&
					( $cmethod eq 'nap' )									 ){
					
					my	$fm = 'nap_risk1' ;
					
					unless ( $href_hoh_seen_rejected_of_clan4_of_filtertype_of_fn->{$$sref_clan4}{$fm}{$qn} ){
						
						push @{$href_hol_sup_of_clan4_of_calc_of_qtree->{$$sref_clan4}{$fm}{$tree}} , $href_hoh_sup_of_clan4_of_calc_of_qtree_of_fc->{$$sref_clan4}{$cmethod}{$tree}{$qn} ;
						
						for my $t ( @taxa ){ $href_hoh_qremaining_of_otu_of_4clan_of_ft->{$fm}{$t}{$$sref_clan4}++ }
						
						$href_seen_remain_of_cmfm_of_fc->{$fm}{$qn}++
					}
					
					else{ 
						
						for my $t ( @taxa ){ push @{$href_qn_of_qfiltered_otu_of_4clan_of_ft->{$fm}{$t}{$$sref_clan4}}, $qn }
						
						$href_hoh_N_qrej_of_clan4_of_ft->{$fm}{$$sref_clan4}++ ;
					}
					
					$href_hoh_seen_filter_of_type->{risk}{$fm}++ ;
				}
				############################
				
				############################
				# filter quartets due to risk2 nap (uncorrected)
				if ( ( $href_value_of_scr_par->{seq_filterRisk}{risk2} == 1 ) &&
					 ( $cmethod eq 'nap' )									  ){
					
					my	$fm = 'nap_risk2' ;	
					
					unless ( $href_hoh_seen_rejected_of_clan4_of_filtertype_of_fn->{$$sref_clan4}{$fm}{$qn} ){
						
						push @{$href_hol_sup_of_clan4_of_calc_of_qtree->{$$sref_clan4}{$fm}{$tree}} , $href_hoh_sup_of_clan4_of_calc_of_qtree_of_fc->{$$sref_clan4}{$cmethod}{$tree}{$qn} ;
						
						for my $t ( @taxa ){ $href_hoh_qremaining_of_otu_of_4clan_of_ft->{$fm}{$t}{$$sref_clan4}++ }
						
						$href_seen_remain_of_cmfm_of_fc->{$fm}{$qn}++
					}
					
					else{ 
						
						for my $t ( @taxa ){ push @{$href_qn_of_qfiltered_otu_of_4clan_of_ft->{$fm}{$t}{$$sref_clan4}}, $qn }
						
						$href_hoh_N_qrej_of_clan4_of_ft->{$fm}{$$sref_clan4}++ ;
					}
					
					$href_hoh_seen_filter_of_type->{risk}{$fm}++ ;
				}
				############################
			}
			############################
		}
		############################
	}
	######################################################################
}

sub filter_quartets_riskdist{
	
	my	$sref_clan4												=	$_[0]  ;	# clan4 combination of given quartet combination																																	IN: defined; OUT: unchanged
	my	$href_hol_sup_of_clan4_of_calc_of_qtree					=	$_[1]  ;	# multiple hash, listing for each of the three rooted quartet-trees of a given 4clan combination all single quartet inferred support (NAP or NAPV)									IN: defined; OUT: changed
	my	$href_hoh_sup_of_clan4_of_calc_of_qtree_of_fc			=	$_[2]  ;	# multiple hash, storing the final support (nap or napv) related to each rooted quartet-tree																						IN: defined; OUT: changed
	my	$href_hoh_seen_rejected_of_clan4_of_filtertype_of_fn	=	$_[3]  ;	# key1: 4clan code; key2: type of supportfilter ('nap_risk1', 'napr_risk1', or 'napv_risk1'); key3: actual sealion quartet number; value: filter corresponding score-value			IN: defined; OUT: unchanged 
	my	$href_value_of_scr_par									=	$_[4]  ;	# key1: sealion user available parameter options; value: optional parameter setting																 									IN: defined; OUT: unchanged
	my	$href_specified_support_method_code						=	$_[5]  ;	# list of defined support methods (nap, napv ...)																																	IN: defined; OUT: changed	
	my	$href_hoh_seen_filter_of_type							=	$_[6]  ;	# key1: filter ('dist' or 'risk'); key2: combi-code method+filter ('nap_dist1', 'nap_risk2'...); value: number of filter passing quartets											IN: defined; OUT: changed
	my	$href_qn_of_qfiltered_otu_of_4clan_of_ft				=	$_[7]  ;	# key1: 4seq filtermethod (e.g. nap_risk1); key2: taxon; key3: 4clan|total; value: list of rejected qn's																			IN: undefined; OUT: defined
	my	$href_hoh_qremaining_of_otu_of_4clan_of_ft				=	$_[8]  ;	# key1: 4seq filtermethod (e.g. nap_risk1); key2: taxon; key3: 4clan|total; value: N remaining quartets																				IN: undefined; OUT: defined
	my	$href_taxa4_of_fc										=	$_[9]  ;	# key1: quartet-number ($fc); value : quartet assigned taxa (alphabetically ordered, joined by ':')																					IN: defined; OUT: unchanged
	my	$href_hoh_N_qrej_of_clan4_of_ft							=	$_[10] ;	# key1: 4seq filtermethod (e.g. nap_risk1); key2: 4clan combination code; value: N rejected quartets*3 (counted three times, each time for each quartet tree)						IN: defined; OUT: changed
	my	$href_seen_remain_of_cmfm_of_fc							=	$_[11] ;	# key1: calculation method + filter method (nap_risk1, napv_dist2, napr_dist3); key2: quartet-number ($fc); value: seen +1															IN: undefined; OUT: defined
	my	$href_seen_filter_method								=	$_[12] ;	# key1: filter method; value: seen +1																																				IN: defined; OUT: changed
	
	#print "::&filter_quartets_riskdist\n";
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# ...
	######################################################################
	
	######################################################################
	# filter quartets for each quartet calculation method                # 
	######################################################################
	
	############################
	# for each risk-filter-method, check for calculation-method corresponding dist-filter-methods
	# and sample commonly non-filtered quartet support as new paired-filter-method in %href_hol_sup_of_clan4_of_calc_of_qtree
	for	my $rfm ( sort keys %{$href_hoh_seen_filter_of_type->{risk}} ){
		
		############################
		# assign risk1 to risk4 codings ($rc) used in %href_hoh_seen_rejected_of_clan4_of_filtertype_of_fn
		# coding of risk1 and risk2 different from $rfm
		unless	( $rfm =~ /risk1|risk2/ ){ die "\nBUG-ERROR: Cannot assign 'risk' code", $rfm ," in subroutine &filter_quartets_riskdist!\nPlease, report BUG to software developer!\n\n" }
		############################
		
		############################
		# extract calc-method ($cm_r: e.g. 'nap') and assigned filter-method ($rm: e.g. 'risk1') from risk-filter-code ($rfm: e.g. 'nap_risk1')
		my	@prts_r	=	split "_", $rfm ;
		my	$cm_r 	=	$prts_r[0] ;
		my	$rm 	=	$prts_r[1] ;	
		############################
		
		############################
		# check each distance-filter method against actual risk-filter methods
		for my	$dfm ( sort keys %{$href_hoh_seen_filter_of_type->{dist}} ){
			
			############################
			# extract calc-method ($cm_d: e.g. 'nap') and assigned filter-method ($dm: e.g. 'dist1') from dist-filter-code ($dfm: e.g. 'nap_dist1')
			my	@prts_d	=	split "_", $dfm ;
			my	$cm_d 	=	$prts_d[0] ;
			my	$dm 	=	$prts_d[1] ;
			############################
			
			############################
			# if calculation-methods of risk and dist filter-method match (e.g. both are based on 'nap')
			if ( $cm_d eq $cm_r ){
				
				############################
				# generate common filter-code for paired risk and dist filtering
				my	$fm = $cm_d."_".$rm."_".$dm ;
				$href_seen_filter_method->{$fm}++ ;
				
				#print	"\n\trisk-filter-method code:\t"	, $rfm	,
				#		"\n\trisk-filter-method code new:\t", $rc	,
				#		"\n\tdist-filter-method code:\t"	, $dfm	,
				#		"\n\tfilter-combination code\t"		, $fm	, 
				#		"\n\t-----------------------\n" 	; #exit;
				############################
				
				############################
				# store each tree support of given calculation method ($cm) which is not filtered by 'risk' or 'dist' in %href_hol_sup_of_clan4_of_calc_of_qtree
				for my $tree ( sort keys %{$href_hoh_sup_of_clan4_of_calc_of_qtree_of_fc->{$$sref_clan4}{$cm_d}} ){
					
					############################
					# go through all quartet numbers of calculation-method... 
					for my $qn ( sort {$a<=>$b} keys %{$href_hoh_sup_of_clan4_of_calc_of_qtree_of_fc->{$$sref_clan4}{$cm_d}{$tree}} ){
						
						$href_hoh_seen_filter_of_type->{riskdist}{$fm}++ ;
						my	@taxa = split ":", $href_taxa4_of_fc->{$qn} ;
						
						############################
						# if quartet is filtered by at least one of the two methods, mark quartet-number as rejected for mixed filter-method
						if (	( $href_hoh_seen_rejected_of_clan4_of_filtertype_of_fn->{$$sref_clan4}{$rfm}{$qn} )	||
								( $href_hoh_seen_rejected_of_clan4_of_filtertype_of_fn->{$$sref_clan4}{$dfm}{$qn} )	){
							
							$href_hoh_seen_rejected_of_clan4_of_filtertype_of_fn->{$$sref_clan4}{$fm}{$qn} = 1;
							
							for my $t ( @taxa ){ push @{$href_qn_of_qfiltered_otu_of_4clan_of_ft->{$fm}{$t}{$$sref_clan4}}, $qn }
							$href_hoh_N_qrej_of_clan4_of_ft->{$fm}{$$sref_clan4}++ ;
						}
						############################
						
						############################
						# sample only support of non actual dist-method or risk-method filtered quartets
						else{ 
							
							push @{$href_hol_sup_of_clan4_of_calc_of_qtree->{$$sref_clan4}{$fm}{$tree}} , $href_hoh_sup_of_clan4_of_calc_of_qtree_of_fc->{$$sref_clan4}{$cm_d}{$tree}{$qn} ;
							
							for my $t ( @taxa ){ $href_hoh_qremaining_of_otu_of_4clan_of_ft->{$fm}{$t}{$$sref_clan4}++ }
							
							$href_seen_remain_of_cmfm_of_fc->{$fm}{$qn}++
						}
						############################
					}
					############################
				}
				############################
			}
			############################
		}
		############################
	}
	############################
}

sub filter_threshold_optimization{
	
	my	$href_printpath							= $_[0]  ; # key1: R grafic unique printcode (e.g. 'metafileQ1'); key2: *.tsv filename; value: *tsv filename data addings 																IN: defined; OUT: unchanged
	my	$href_hoh_data_of_tsvfile_of_prcode		= $_[1]  ; # key: output-type (e.g. 'pdf'); value: pathdir																																IN: defined; OUT: changed
	my	$href_hoh_data_of_texfile_of_prcode		= $_[2]  ; # key1: Latex table unique printcode (e.g. 'latexC'); key2: *.tex filename; value: @list of assigned table lines																IN: defined; OUT: changed
	my	$href_textfile_of_prcode				= $_[3]  ; # key1: textfile unique printcode (e.g. 'textA'); value: filename for print																									IN: defined; OUT: changed
	my	$href_value_of_scr_par					= $_[4]	 ; # key1: sealion user available parameter options; value: optional parameter setting																 							IN: defined, OUT: unchanged
	my	$sref_clan4								= $_[5]	 ; # Clan4 of given quartet combination																																			IN: defined, OUT: unchanged
	my	$href_hol_fc_of_clan4					= $_[6]	 ; # key1: 4clan code; value: list of 4clan analysed quartet numbers																											IN: defined, OUT: unchanged
	my	$href_hoh_rejq_of_ftype_of_scale_of_fn	= $_[7]	 ; # key1: type of supportfilter ('nap_risk1', 'napr_risk1', or 'napv_risk1'); key2: quartet number; ; key3: scale value; value: 1 (reject)										IN: defined, OUT: undefined
	my	$href_hoh_rejq_of_4clan_of_ftype_of_fn	= $_[8]	 ; # key1: 4clan code; key2: type of supportfilter ('nap_risk1', 'napr_risk1', or 'napv_risk1'); key3: actual sealion quartet number; value: filter corresponding score-value	IN: undefined; OUT: defined 
	my	$href_hol_scales_of_filtertype			= $_[9]  ; # key1: type of supportfilter ('risk', 'qdist'); value: list of ftype corresponding threshold scales																			IN: defined; OUT: unchanged
	my	$href_seen_fmethod						= $_[10] ; # key1: filter method approaches; value: seen +1																																IN: defined; OUT: unchanged
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# &print_tsv_MQ6
	# &latex_data_sample_LQ6
	######################################################################
	
	######################################################################
	# Find the best threshold value for quartet-filtering (RISK || DIST) #
	######################################################################
	
	#######################
	# total number of quartets for clan4
	my $Nq_clan4 = @{$href_hol_fc_of_clan4->{$$sref_clan4}} ;
	#######################
	
	LOOP:
	#for my $supportfilter ( sort keys %$href_hoh_rejq_of_ftype_of_scale_of_fn ){ #print "\n\tft\t", $supportfilter, "\n" ; exit;
	for my $supportfilter ( sort keys %$href_seen_fmethod ){ #print "\n\tft\t", $supportfilter, "\n" ; exit;
		
		#######################
		# assign RISK and DIST individual defined alpha for matching filter
		my @scale_values ;
		if    ( $supportfilter =~ /risk._qdist./ ){ next LOOP  }
		elsif ( $supportfilter =~ /risk/ ){ @scale_values = @{$href_hol_scales_of_filtertype->{risk}}  }
		elsif ( $supportfilter =~ /dist/ ){ @scale_values = @{$href_hol_scales_of_filtertype->{qdist}} }
		else  { die "\nBUG-ERROR: Cannot assign 'supportfilter' in subroutine &filter_threshold_optimization!\nPlease, report BUG to software developer!\n\n" }
		#######################
		
		#######################
		# for each scale-value, count number of remaining quartets ($ni) and rejected quartets ($N_qrej)
		my (
				%ni_of_scale 		, # key1: scale value; value: N remaining quartets
				%hol_scale_of_ni	, # key1: N remaining quartets; value: list of scales
				@prtlines_MQ6		, # samples datalines for tsv table MQ6
		) ;
		
		#for my $scale ( sort keys %{$href_hoh_rejq_of_ftype_of_scale_of_fn->{$supportfilter}} ){
		for my $scale ( @scale_values ){
				
			my $N_qrej = 0 ;
			
			if ( keys %{$href_hoh_rejq_of_ftype_of_scale_of_fn->{$supportfilter}{$scale}} ){
				
				$N_qrej = keys %{$href_hoh_rejq_of_ftype_of_scale_of_fn->{$supportfilter}{$scale}} ;
			}
			
			my $ni = $Nq_clan4 - $N_qrej ;
			#$ni_of_scale{$scale} = $N_qrej ;
			
			$ni_of_scale{$scale} = $ni ;
			push @{$hol_scale_of_ni{$ni}}, $scale ;
		}
		#######################
		
		#######################
		# sort scale-values of identical number of remaining quartets ($ni) from low to high for better graphical visualization
		my @sorted_scales ;
		for my $sorted_ni ( sort {$a<=>$b} keys %hol_scale_of_ni ){
			
			for my $sorted_scale ( sort {$a<=>$b} @{$hol_scale_of_ni{$sorted_ni}} ){
				
				#print "\nni\t", $sorted_ni, "\tscale\t", $sorted_scale ; 
				
				push @sorted_scales, $sorted_scale
			}
		}
		#for (@sorted_scales){print "\nscale\t", $_} exit;
		#######################
		
		##############################################
		# start optimization process
		my $Nmatrixcells_start	= @sorted_scales * $Nq_clan4 ;	# N matrix cells given all scale-values (y-axis) and quartets (x-axis)
		my $fp_0				= 0 ;							# old optimum-score
		my $fp_1				= 0 ; 							# new optimum score
		my $fp_bestscale		= 0 ;							# best optimum score
		
		#######################
		# assign RISK and DIST individual defined alpha for matching filter
		my $alpha ;
		if    ( $supportfilter =~ /risk/ ){ $alpha = $href_value_of_scr_par->{seq_filterRisk}{ralpha} }
		elsif ( $supportfilter =~ /dist/ ){ $alpha = $href_value_of_scr_par->{seq_filterDist}{dalpha} }
		else  { die "\nBUG-ERROR: Cannot assign 'supportfilter' in subroutine &filter_threshold_optimization!\nPlease, report BUG to software developer!\n\n" }
		#######################
		
		#######################
		# assign latex print firendly filter method code
		my	$filter_type	=	&recode_filter_printout( \$supportfilter )	;
		my	@cm_fm			=	split "_", $supportfilter ;
		my	$cm				=	$cm_fm[0] ;
		#######################
		
		#######################
		# while new calculated optimum ($fp_1) > old optimum (fp_0), 
		# ...(1) calculate optimum ($fp_1) for next increasingly sorted scale-value $ni
		# ...(2) if ($fp_1) > ($fp_0), define new optimum ($fp_1) as old optimum ($fp_0) 
		# ...(3) continue with next scale-value
		while ( ( $fp_0 == $fp_1 ) && ( @sorted_scales > 0 ) ){
			
			############
			# calculate remaining quartets ($ni_total) given actual set of unchecked scale-values
			my	$ni_total 	=	0 ;
			for my $scale ( @sorted_scales ){ $ni_total += $ni_of_scale{$scale} }
			############
			
			############
			# Number of unchecked scale-values
			my	$Nscales	=	@sorted_scales ;
			############
			
			############
			# calculate P(B), N remaining quartets / total number quartets, sumed over all unchecked scale-values
			my	$pb			=	$ni_total / ( $Nscales * $Nq_clan4 ) ;
			############
			
			############
			# calculate lamda, the ratio of the actual PB related to the PB value of the previous iteration step
			my	$lambda		=	( $Nscales * $Nq_clan4 ) / ( ( $Nscales + 1 ) * $Nq_clan4 ) ;
			############
			
			############
			# calculate new optimum-score fp_1
			$fp_1			=	1 - abs( $lambda - ( $pb ** ( $alpha * ( 1 - $pb ) ) ) );
			############
			
			############
			# define new optimum as old optimum if fp_1 > fp_0
			# assign scale value of actual scale as new best value
			if ( $fp_1 > $fp_0 ){ 
				
				$fp_0 = $fp_1 ;
				$fp_bestscale = $sorted_scales[0] ; 
			}
			############
			
			############
			# join and sample single meta data as printline for tsv table MQ6
			push @prtlines_MQ6, $cm."\t".$filter_type."\t".$$sref_clan4."\t".$sorted_scales[0]."\t".$ni_of_scale{$sorted_scales[0]}."\t".$alpha."\t".$lambda."\t".$fp_0."\t".$fp_1."\t".$fp_bestscale ;
			############
			
			############
			# remove actual scale-value from the list of yet unchecked values
			shift @sorted_scales
			############
		}
		
		if    ( $supportfilter =~ /risk/ ){ $href_value_of_scr_par->{seq_filterRisk}{$$sref_clan4}{tlrisk}{$supportfilter} = $fp_bestscale }
		elsif ( $supportfilter =~ /dist/ ){ $href_value_of_scr_par->{seq_filterDist}{$$sref_clan4}{tldist}{$supportfilter} = $fp_bestscale }
		#######################
		
		#######################
		# transfer rejected (filtered) quartet-numbers of best scale-value to multiple hash %$href_hoh_rejq_of_4clan_of_ftype_of_fn
		for my $qn ( sort {$a<=>$b} keys %{$href_hoh_rejq_of_ftype_of_scale_of_fn->{$supportfilter}{$fp_bestscale}} ){
			
			$href_hoh_rejq_of_4clan_of_ftype_of_fn->{$$sref_clan4}{$supportfilter}{$qn} = $href_hoh_rejq_of_ftype_of_scale_of_fn->{$supportfilter}{$fp_bestscale}{$qn}
		}
		##############################################
		
		##############################################
		# adding best scale value info to each tsv table line MQ6
		# extract sub metadata from each line MQ6 for tex table LQ6 print
		my $prtlines_LQ6	; 
		for my $prtline ( @prtlines_MQ6 ){ 
			
			$prtline .= "\t".$fp_bestscale ;
			my @prts = split "\t", $prtline ;
			
			if ( $prts[3] == $fp_bestscale ){ 
				
				$prtline .= "\tbest" ;
				
				$prtlines_LQ6 .= $prts[0]."\t".$prts[1]."\t".$prts[2]."\t".$prts[3]."\t".$prts[4]."\t".$Nq_clan4."\n"
			} 
			
			else{ $prtline .= "\tminor" }
		}
		##############################################
		
		##############################################
		# table prints
		&print_tsv_MQ6(
					
			\%$href_printpath						, # key1: R grafic unique printcode (e.g. 'metafileQ1'); key2: *.tsv filename; value: *tsv filename data addings	IN: defined; OUT: unchanged
			\%$href_hoh_data_of_tsvfile_of_prcode	, # key: output-type (e.g. 'pdf'); value: pathdir																	IN: defined; OUT: changed
			\@prtlines_MQ6							, # tsv data lines 																									IN: defined; OUT: unchanged
			\$fp_bestscale							, # overall best scale value of actual filtermethod																	IN: defined; OUT: unchanged
			\$supportfilter							, # actual filtermethod																								IN: defined; OUT: unchanged
		) ;
		
		&latex_data_sample_LQ6(
			
			\%$href_printpath						, # key1: R grafic unique printcode (e.g. 'metafileQ1'); key2: *.tsv filename; value: *tsv filename data addings 	IN: defined; OUT: unchanged
			\%$href_hoh_data_of_texfile_of_prcode	, # key1: Latex table unique printcode (e.g. 'latexC'); key2: *.tex filename; value: @list of assigned table lines	IN: defined; OUT: changed
			\%$href_textfile_of_prcode				, # key1: textfile unique printcode (e.g. 'textA'); value: filename for print										IN: defined; OUT: changed
			\$supportfilter							, # actual filtermethod																								IN: defined; OUT: unchanged
			\$prtlines_LQ6							, # tex data lines																									IN: defined; OUT: unchanged
		) ;
		##############################################
	}
	#exit;
	#######################
	# delete info about rejected quartets from analysed sclade values 
	for my $key ( keys %$href_hoh_rejq_of_ftype_of_scale_of_fn ){ delete $href_hoh_rejq_of_ftype_of_scale_of_fn->{$key} } ;
	#######################
	#exit;
}

sub generate_clan_quartets{

	my	$href_hol_taxa_of_clan		= $_[0]	;	# key1: subclade name; value: list of subclade assigned taxa								IN: defined; OUT: unchanged
	my	$aref_clan_quartets			= $_[1]	;	# list of generated 4-clan combinations														IN: undefined; OUT: defined
	my	$href_value_of_scr_par		= $_[2] ;	# key1: sealion user available parameter options; value: optional parameter setting			IN: defined; OUT: unchanged
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# &help
	######################################################################
	
	######################################################################
	# Build all possible 4clan combinations & extract outgroup combis    #
	######################################################################
	
	#######################
	# Generate all possible quartets of defined clans
	my	@clans	=	sort keys %$href_hol_taxa_of_clan ;		#for (@clans){print "\nclan:", $_}	exit;
	
	my	@start_clans	=	@clans	;
	my	@second_clans	=	()		;
	my	@all_quartets	=	()		;

	while (@clans){

		my $c1			= shift @clans ;
		my @clans_copy	= @clans ;

		while (@clans_copy > 2){

			while (@clans_copy > 2){

				my $c2	= shift @clans_copy ;
				my $c3	= shift @clans_copy ;

				for my $c4 ( @clans_copy ){

					my $q	=	$c1.":".$c2.":".$c3.":".$c4	;
					push @all_quartets, $q ;
				}

				push	@second_clans, $c3 ;
				unshift	@clans_copy, $c2
			}

			push @second_clans, $clans_copy[1] ;
			@clans_copy = @second_clans ;		#print "second: ", @second_clans;
			@second_clans = () ;
		}
	}
	#######################
	
	#######################
	# select only 4-clan combinations of the outgroup-clan if outgroup has been defined
	# else take all possible 4-clan combinations
	if	( $href_value_of_scr_par->{sys_parameter}{outgroup}	){

		#### Selection of the 4-clan combinations given the outgroup clan
		my $outgroup_clan = $href_value_of_scr_par->{sys_parameter}{outgroup} ;
		
		for my $clan_q	( @all_quartets	){
			
			if	( $clan_q =~ /$outgroup_clan/ ){ push @$aref_clan_quartets,	$clan_q	 }
		}
		####
		
		#### die script process if defined outgroup clan cannot be found in list of all generated 4-clan combinations
		unless	( @$aref_clan_quartets	>=	1 ){

			print "\nPARAMETER-ERROR: Unknown outgroup: ", $outgroup_clan , " !\n\n" ;
			
			&help( \'o', \%$href_value_of_scr_par )
		}
		####
	}
	else { @$aref_clan_quartets	= @all_quartets	}
	#######################

	#######################
	# count and terminal print of final number of 4-clan combinations
	my	$N_clan_quartets =	@$aref_clan_quartets ;
	#######################
	
	#######################
	## test print
	if ( $href_value_of_scr_par->{sys_parameter}{screen_prt} == 1 ){
		
		print "\n\tNumber 4clan combinations: ", $N_clan_quartets, "\n"	;
		
		my	$counter = 1 ;
		for my $clan_q	( @$aref_clan_quartets ){ print "\n\t4clan ", $counter,": ", $clan_q ; $counter++ }
		print "\n" ;
	}
	#######################

	#######################
	# empty variables which are not needed in further processes
	@all_quartets	=	() ;
	#######################
	
	######################################################################
}

sub generate_taxon_quartets{

	my	$href_hol_taxa_of_clan	=	$_[0]	; # key: clan name; value: list of clan assigned taxa											IN: defined; OUT: unchanged
	my	$sref_clan_combination	=	$_[1]	; # 4-clan string, separated by doublepoints (e.g.: c1:c2:c3:c4)								IN: defined; OUT: unchanged
	my	$aref_taxon_quartets	=	$_[2]	; # list of generated 4seq combinations															IN: undefined; OUT: defined
	my	$href_val_of_scr_par	=	$_[3]	; # key1: sealion user available parameter options; value: optional parameter setting			IN: defined; OUT: unchanged
	my	$href_hol_states_of_t	=	$_[4]	; # key1: taxon name; key2: state position number; value: state 								IN: defined; OUT: unchanged
	my	$href_data_of_fileprop	=	$_[5]	; # key1: property; value : property value														IN: defined; OUT: unchanged
	
	print "\n\t::Building quartets and site-pattern check of 4clan ", $$sref_clan_combination ,"...(may take a while)\n\n";
	
	######################################################################
	# called from subroutine(s)...
	# &quartet_building_and_p4_main
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	# Build all sequence-quartets of 4clan combination                   #
	######################################################################
	my	@clans		=	split ":",  $$sref_clan_combination	;
	my	@clantaxa1	=	@{$href_hol_taxa_of_clan->{$clans[0]}}	;
	my	@clantaxa2	=	@{$href_hol_taxa_of_clan->{$clans[1]}}	;
	my	@clantaxa3	=	@{$href_hol_taxa_of_clan->{$clans[2]}}	;
	my	@clantaxa4	=	@{$href_hol_taxa_of_clan->{$clans[3]}}	;
	
	my	%seen_rejected_quartets ;
	my	$quartet_counter  = 0 ;
	my	$quartet_counter2 = 0 ;
	
	##############################################
	# Define allowed sequence states, all other sequence states are rejected from further split analyses
	my %allowed_states ;
	
	if			( $href_data_of_fileprop->{type} =~ 'nuc'	){ for ( qw/A C G T/									){ $allowed_states{$_}++ } } # all nuc states allowed
	elsif		( $href_data_of_fileprop->{type} =~ 'aa'	){ for ( qw/A C G T N Y R W S K M D V H I E L Q F P/	){ $allowed_states{$_}++ } } # all aa states allowed
	else 		{ die "\nBUG-ERROR: Cannot assign data type in subroutine &perform_quartet_analyses!\nPlease, report BUG to software developer!\n\n" }
	##############################################
	
	for my $taxon1 ( @clantaxa1 ){
		
		for my $taxon2 ( @clantaxa2 ){
			
			for my $taxon3 ( @clantaxa3 ){
				
				for my $taxon4 ( @clantaxa4 ){
					
					$quartet_counter++  ;
					$quartet_counter2++ ;
					
					my	@sampled		=	sort ($taxon1,$taxon2,$taxon3,$taxon4); 
					my	$qstring		=	join ":", @sampled ; 
					my	$lengthseq		=	@{$href_hol_states_of_t->{$taxon1}} ;
					my	$counter_inf	=	0 ;
					
					##############################################
					# analyse quartet states of each sequence positions
					# reject sites of unallowed characters
					for my $seq_pos ( 0 .. $href_data_of_fileprop->{length}-1 ){ 
						
						#print "\n\t", $seq_pos, "\t", $aref_sampled->[0], "\t", $aref_sampled->[1], "\t", $aref_sampled->[2], "\t", $aref_sampled->[3], "\t", $href_hol_seqstates_of_t->{$aref_sampled->[0]}[$seq_pos];
						
						#######################
						# check site characters of each of the 4 taxa
						if	( ( $allowed_states{$href_hol_states_of_t->{$sampled[0]}[$seq_pos]} ) &&
							  ( $allowed_states{$href_hol_states_of_t->{$sampled[1]}[$seq_pos]} ) &&
							  ( $allowed_states{$href_hol_states_of_t->{$sampled[2]}[$seq_pos]} ) &&
							  ( $allowed_states{$href_hol_states_of_t->{$sampled[3]}[$seq_pos]} ) ){
							
							$counter_inf++
						}
						#######################
					}
					#print "\n\tquartet\t", $quartet_counter, "\tN inf\t", $counter_inf, "\tN orig\t", $lengthseq ;
					#######################
					
					if ( $counter_inf >= $href_val_of_scr_par->{qlimit}{qlength} ){
					
						push @$aref_taxon_quartets, $sampled[0].":".$sampled[1].":".$sampled[2].":".$sampled[3] ;
					}
					
					else{ $seen_rejected_quartets{$qstring}++ }
					#######################
					
					if ( $quartet_counter == 100 ){ print "\tcheck quartets (q) of clade-quartet\t", $$sref_clan_combination, "...\tq", $quartet_counter2, "\tchecked\n"; $quartet_counter = 0 }
				}
			}
		}
	}
	
	for my $qstring ( sort keys %seen_rejected_quartets ){ print "\n\t", $qstring, " rejected: N site pattern < minimum allowed sites of ", $href_val_of_scr_par->{qlimit}{qlength} }
	##############################################
	
	##############################################
	# die if no quartet above/equal minimum allowed sequence length
	unless ( @$aref_taxon_quartets ){ 
		
		die "\nQUARTET-ERROR: No quartet left with sequence length > minimum sequence limit (".$href_val_of_scr_par->{qlimit}{qlength}.") for 4clan combination ".$$sref_clan_combination." after excluding Indel/Ambiguity site positions!\n\n" 
	}
	##############################################
	#exit;
	######################################################################
}

sub latex_data_sample_LRC1{
	
	my	$href_hoh_data_of_tsvfile_of_prc	= $_[0] ; # key: output-type (e.g. 'pdf'); value: pathdir									IN: defined; OUT: changed
	my	$sref_sm							= $_[1] ; # support method (nap, napv, napr)												IN: defined; OUT: unchanged
	my	$sref_fm							= $_[2] ; # filter method (risk1, strict_all:hardclipping, risk2_strict_all:hardclipping)	IN: defined; OUT: unchanged
	my	$sref_btr							= $_[3] ; # best tree in nwk																IN: defined; OUT: unchanged
	my	$sref_ts							= $_[4] ; # mean|median																		IN: defined; OUT: unchanged
	my	$sref_tr_supp						= $_[5] ; # final best tree support															IN: defined; OUT: unchanged
	my	$href_textfilename_of_prtcode		= $_[6] ; # key1: textfile unique printcode (e.g. 'textA'); value: filename for print		IN: defined; OUT: unchanged
	my	$href_subfolder_q					= $_[7] ; # key1: output-type (e.g. 'pdf'); value: pathdir									IN: defined; OUT: unchanged
	
	#print "::&latex_data_sample_LRC1\n"; 
	
	######################################################################
	# called from subroutine(s)...
	# &treefishing_result_processing
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	# Open TXT outfile and write header line             		         #
	######################################################################
	my		$f_pinfo = $href_subfolder_q->{txt}."/".$href_textfilename_of_prtcode->{txtLRC1} ;
	open 	OUTtxt,	">>$f_pinfo" or die "\nTEXTFILE-ERROR: Cannot write to ", $f_pinfo, "!\n\n" ;
	
	############
	# print header only once
	unless ( $href_textfilename_of_prtcode->{txtLRC1hcode} ){
		
		print	OUTtxt	"Score", "\t", "Filter","\t", "Best_Rooted-Clade_Tree", "\t", "Average", "\t", "Support", "\n" ;  
		$href_textfilename_of_prtcode->{txtLRC1hcode}++
	}
	######################################################################
	
	######################################################################
	# Generate latex data line, listing best final tree support data    ##
	######################################################################
	
	###############
	# data sampling for latextable LT3 (ternaryplot-legend)
	my	$topscore_rounded	=	sprintf("%.2f", $$sref_tr_supp	) ;
	my	$filter				=	&recode_filter_printout( \$$sref_fm )	;
						
	my	$latex_info_line	=	$$sref_sm."\t".$filter."\t".$$sref_btr."\t".$$sref_ts."\t".$topscore_rounded."\n" ;
	push @{$href_hoh_data_of_tsvfile_of_prc->{latexLRC1}{$$sref_sm}}, $latex_info_line	;
	######################################################################
	
	######################################################################
	# write latex dataline to text file                                  #
	######################################################################
	print	OUTtxt	$latex_info_line ;
	close	OUTtxt ;
	######################################################################
}

sub latex_data_sample_LRC2{
	
	my	$href_hoh_data_of_tsvfile_of_prc	= $_[0] ; # key: output-type (e.g. 'pdf'); value: pathdir									IN: defined; OUT: changed
	my	$sref_sm							= $_[1] ; # support method (nap, napv, napr)												IN: defined; OUT: unchanged
	my	$sref_fm							= $_[2] ; # filter method (risk1, strict_all:hardclipping, risk2_strict_all:hardclipping)	IN: defined; OUT: unchanged
	my	$aref_trees							= $_[3] ; # best three trees in nwk	sorted by support										IN: defined; OUT: unchanged
	my	$sref_ts							= $_[4] ; # mean|median																		IN: defined; OUT: unchanged
	my	$aref_tree_supp						= $_[5] ; # sorted support	of the three best trees											IN: defined; OUT: unchanged
	my	$href_textfilename_of_prtcode		= $_[6] ; # key1: textfile unique printcode (e.g. 'textA'); value: filename for print		IN: defined; OUT: unchanged
	my	$href_subfolder_q					= $_[7] ; # key1: output-type (e.g. 'pdf'); value: pathdir									IN: defined; OUT: unchanged
	my	$sref_dist_best_2nd					= $_[8] ; # distance support best rooted-clade tree to 2nd best								IN: defined; OUT: unchanged
	
	#print "::&latex_data_sample_LRC2\n"; 
	
	######################################################################
	# called from subroutine(s)...
	# &treefishing_result_processing
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	# Open TXT outfile and write header line             		         #
	######################################################################
	my		$f_pinfo = $href_subfolder_q->{txt}."/".$href_textfilename_of_prtcode->{txtLRC2} ;
	open 	OUTtxt,	">>$f_pinfo" or die "\nTEXTFILE-ERROR: Cannot write to ", $f_pinfo, "!\n\n" ;
	
	############
	# print header only once
	unless ( $href_textfilename_of_prtcode->{txtLRC2hcode} ){
		
		print	OUTtxt	"Score", "\t", "Filter", "\t", "Average", "\t", "Best", "\t", "Support_Best" ,"\t", "2nd", "\t", "Support_2nd" ,"\t", "3rd", "\t", "Support_3rd" ,"\t", "Distance_Best_2nd", "\n" ;
		
		$href_textfilename_of_prtcode->{txtLRC2hcode}++
	}
	######################################################################
	
	######################################################################
	# Generate latex data line, listing best final tree support data    ##
	######################################################################
	
	###############
	# data sampling for latextable LT3 (ternaryplot-legend)
	my	$treescore_1_rounded	=	sprintf("%.2f", $aref_tree_supp->[0]	) ;
	my	$treescore_2_rounded	=	sprintf("%.2f", $aref_tree_supp->[1]	) ;
	my	$treescore_3_rounded	=	sprintf("%.2f", $aref_tree_supp->[2]	) ;
	
	my	$filter				=	&recode_filter_printout( \$$sref_fm )	;
						
	my	$latex_info_line	=	$$sref_sm."\t".$filter."\t".$$sref_ts."\t".$aref_trees->[0]."\t".$treescore_1_rounded."\t".$aref_trees->[1]."\t".$treescore_2_rounded."\t".$aref_trees->[2]."\t".$treescore_3_rounded."\t".$$sref_dist_best_2nd."\n" ;
	push @{$href_hoh_data_of_tsvfile_of_prc->{latexLRC2}{$$sref_sm}}, $latex_info_line	;
	######################################################################
	
	######################################################################
	# write latex dataline to text file                                  #
	######################################################################
	print	OUTtxt	$latex_info_line ;
	close	OUTtxt ;
	######################################################################
}

sub latex_data_sample_LP1{
	
	my	$href_val_of_scr_par	= $_[0] ; # key1: options; value: optional parameter setting		IN: defined; OUT: unchanged
	my	$href_ldata_of_flatex	= $_[1] ; # key1: printcode; key2: fname; value: table lines		IN: defined; OUT: changed
	
	######################################################################
	# called from subroutine(s)...
	# &print_script_settings
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	# Generate latex data line, listing specified sealion parameter      #
	######################################################################
	
	###############
	# data sampling for latextable LP1
	@{$href_ldata_of_flatex->{latexLP1}} = (
			
			"MSA infile:\t".$href_val_of_scr_par->{menu_parameter}{msa_infile}."\n" 			,
			"Clade infile:\t".$href_val_of_scr_par->{menu_parameter}{clfile}."\n"				,
			"Outgroup:\t".$href_val_of_scr_par->{menu_parameter}{outgroup}."\n"					,
			"Nmax quartets:\t".$href_val_of_scr_par->{menu_parameter}{maxquartets}."\n"			, 
			"Nmin quartet sites:\t".$href_val_of_scr_par->{menu_parameter}{minlengths}."\n"		, 
			"Support-average:\t".$href_val_of_scr_par->{menu_parameter}{average}."\n"			, 
			"SPD-infolder:\t".$href_val_of_scr_par->{menu_parameter}{spd_folder}."\n"			,
			"\t\n"																				, 
			"QUARTET-FILTER 'RISK'\t \n" 														,
			"RISK1:\t".$href_val_of_scr_par->{menu_parameter}{risk1}."\n" 						,
			"RISK2:\t".$href_val_of_scr_par->{menu_parameter}{risk2}."\n" 						,
			"Upper RISK limit:\t".$href_val_of_scr_par->{menu_parameter}{upperlimit}."\n"		,
			"Lower RISK limit:\t".$href_val_of_scr_par->{menu_parameter}{lowerlimit}."\n"		,
			"Lower RISK scale:\t".$href_val_of_scr_par->{menu_parameter}{rscale}."\n"			,
			"Lower RISK alpha:\t".$href_val_of_scr_par->{menu_parameter}{ralpha}."\n"			,
			"\t\n"																				, 
			"QUARTET-FILTER 'DIST'\t \n" 														,
			"DIST:\t".$href_val_of_scr_par->{menu_parameter}{dist1}."\n" 						,
			"Upper DIST limit:\t".$href_val_of_scr_par->{menu_parameter}{upperlimit_d}."\n"		,
			"Lower DIST limit:\t".$href_val_of_scr_par->{menu_parameter}{lowerlimit_d}."\n"		,
			"Lower DIST scale:\t".$href_val_of_scr_par->{menu_parameter}{dscale}."\n"			,
			"Lower DIST alpha:\t".$href_val_of_scr_par->{menu_parameter}{dalpha}."\n"			,
			"\t\n"																				, 
			"P4-PARAMETER:\t \n" 																,
			"Model (aa):\t".$href_val_of_scr_par->{menu_parameter}{model_aa}."\n"				, 
			"Model (nu):\t".$href_val_of_scr_par->{menu_parameter}{model_nu}."\n"				, 
			"Start alpha:\t".$href_val_of_scr_par->{menu_parameter}{start_alpha}."\n"			, 
			"Start pINV:\t".$href_val_of_scr_par->{menu_parameter}{start_pinv}."\n"				, 
			"\t\n"																				, 
			"OUTPUT:\n\t " 																		,
			"Output-folder:\t".$href_val_of_scr_par->{menu_parameter}{output_folder}."\n"		, 
			"Latex-table Print\t".$href_val_of_scr_par->{menu_parameter}{prt_latex}."\n"		,
			"R-plotting\t".$href_val_of_scr_par->{menu_parameter}{prt_Rplot}."\n"				,
	) ;
	###############
	#exit;
	######################################################################
}

sub latex_data_sample_LP2{
	
	my	$href_ldata_of_flatex				= $_[0] ; # key1: printcode; key2: fname; value: table lines								IN: defined; OUT: changed
	my	$href_code_of_clan					= $_[1]	; # key1: clanname; value assigned code												IN: defined; OUT: unchanged
	my	$href_textfilename_of_printcode		= $_[2] ; # key1: textfile unique printcode (e.g. 'textA'); value: filename for print		IN: defined; OUT: unchanged
	my	$href_subfolder_q					= $_[3] ; # key1: output-type (e.g. 'pdf'); value: pathdir									IN: defined; OUT: unchanged
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	# Add clan names and their codes to TXT outfile                      #
	######################################################################
	my		$f_pinfo = $href_subfolder_q->{txt}."/".$href_textfilename_of_printcode->{txtA} ;
	open 	OUTtxt,	">>$f_pinfo" or die "\nTEXTFILE-ERROR: Cannot write to ", $f_pinfo, "!\n\n" ;
	
	print	OUTtxt	"\n\n"														,
					"\t---------------------------------------------------\n"	,
					"\tClade\tCode\n"											,
					"\t---------------------------------------------------\n"	;
	
	######################################################################
	
	
	######################################################################
	# Generate latex data line, listing clan names and their codes       #
	######################################################################
	for my $clan ( sort keys %$href_code_of_clan ){
		
		my	$latex_info_line = $clan."\t".$href_code_of_clan->{$clan}."\n" ;
		push @{$href_ldata_of_flatex->{latexLP2}}, $latex_info_line	;
		
		print OUTtxt "\t".$latex_info_line ;
	}
	
	print OUTtxt "\t---------------------------------------------------\n" ;
	close OUTtxt ;
	######################################################################
}

sub latex_data_sample_LP3{
	
	my	$href_ldata_of_flatex				= $_[0] ; # key1: printcode; key2: fname; value: table lines								IN: defined; OUT: changed
	my	$href_code_of_4clan					= $_[1]	; # key1: 4clan assigned clan names;  value assigned code							IN: defined; OUT: unchanged
	my	$href_textfilename_of_printcode		= $_[2] ; # key1: textfile unique printcode (e.g. 'textA'); value: filename for print		IN: defined; OUT: unchanged
	my	$href_subfolder_q					= $_[3] ; # key1: output-type (e.g. 'pdf'); value: pathdir									IN: defined; OUT: unchanged
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	# Add 4clan combinations and their codes to TXT outfile              #
	######################################################################
	my		$f_pinfo = $href_subfolder_q->{txt}."/".$href_textfilename_of_printcode->{txtA} ;
	open 	OUTtxt,	">>$f_pinfo" or die "TEXTFILE-ERROR: Cannot write to ", $f_pinfo, "!\n\n" ;
	
	print	OUTtxt	"\n\n"														,
					"\t---------------------------------------------------\n"	,
					"\tClade-Quartet\tCode\n"								,
					"\t---------------------------------------------------\n"	;
	
	######################################################################
	
	######################################################################
	# Generate latex data line, listing 4clan combi and their codes      #
	######################################################################
	for my $clan ( sort keys %$href_code_of_4clan ){
		
		my	$latex_info_line = $clan."\t".$href_code_of_4clan->{$clan}."\n" ;
		push @{$href_ldata_of_flatex->{latexLP3}}, $latex_info_line	;
		
		print OUTtxt "\t".$latex_info_line ;
	}
	
	print OUTtxt "\t---------------------------------------------------\n" ;
	close OUTtxt ;
	######################################################################
}

sub latex_data_sample_LQ1{
	
	my	$href_hoh_data_of_lfile_of_prtcode	= $_[0]	; # key1: Latex table unique printcode (e.g. 'latexC'); key2: *.tex filename; value: @list of assigned table lines				IN: defined; OUT: changed
	my	$sref_cl4							= $_[1]	; # clan4-combination																											IN: defined; OUT: unchanged
	my	$sref_cm							= $_[2]	; # calcultaion method (nap, napv, napr)																						IN: defined; OUT: unchanged
	my	$sref_fm							= $_[3]	; # filter method (risk1, risk2...)																								IN: defined; OUT: unchanged
	my	$sref_qtcode						= $_[4]	; # 4clan tree code (QT1, QT2, QT3)																								IN: defined; OUT: unchanged
	my	$sref_tree							= $_[5]	; # 4clan tree																													IN: defined; OUT: unchanged
	my	$sref_average						= $_[6]	; # mean|median																													IN: defined; OUT: unchanged
	my	$sref_support_rounded				= $_[7]	; # mean|median supported rounded																								IN: defined; OUT: unchanged
	my	$href_subfolder_q					= $_[8] ; # hashpath to output-type specific subdirs of single quartet-pattern analyses; key: output-type (e.g. 'pdf'); value: pathdir	IN: defined; OUT: unchanged
	my	$href_textfilename_of_prtcode		= $_[9] ; # key1: textfile unique printcode (e.g. 'textA'); value: filename for print													IN: defined; OUT: unchanged
	
	#print "::&latex_data_sample_LQ1\n"; 
	
	######################################################################
	# called from subroutine(s)...
	# &calc_4clan_tree_support
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	# Add 4clan combinations and their codes to TXT outfile              #
	######################################################################
	my		$f_pinfo = $href_subfolder_q->{txt}."/".$href_textfilename_of_prtcode->{txtLQ1} ;
	open 	OUTtxt,	">>$f_pinfo" or die "\nTXTFILE-ERROR: Cannot write to ", $f_pinfo, "!\n\n" ;
	
	############
	# print header only once
	unless ( $href_textfilename_of_prtcode->{txtLQ1hcode} ){
		
		print	OUTtxt	"Clade-Quartet", "\t", "Score", "\t", "Filter","\t", "Tree", "\t", "Tree Code", "\t", "Average", "\t", "Support", "\n" ;  
		$href_textfilename_of_prtcode->{txtLQ1hcode}++
	}
	######################################################################
	
	######################################################################
	# Generate latex data line, listing LQ1 of defined triangle corners  #
	######################################################################
	my	$filter			 =	&recode_filter_printout( \$$sref_fm )	;
	my	$latex_info_line =	$$sref_cl4."\t".$$sref_cm."\t".$filter."\t".$$sref_tree."\t".$$sref_qtcode."\t".$$sref_average."\t".$$sref_support_rounded."\n" ;
	push @{$href_hoh_data_of_lfile_of_prtcode->{latexLQ1}{$$sref_cm}}, $latex_info_line	;
	
	print	OUTtxt	$latex_info_line ;
	close	OUTtxt ;
	######################################################################
}

sub latex_data_sample_LQ2{
	
	my	$href_hoh_data_of_lfile_of_prtcode	= $_[0]	 ; # key1: Latex table unique printcode (e.g. 'latexC'); key2: *.tex filename; value: @list of assigned table lines				IN: defined; OUT: changed
	my	$sref_cl4							= $_[1]	 ; # clan4-combination																											IN: defined; OUT: unchanged
	my	$sref_cm							= $_[2]	 ; # calcultaion method (nap, napv, napr)																						IN: defined; OUT: unchanged
	my	$sref_fm							= $_[3]	 ; # filter method (risk1, risk2...)																							IN: defined; OUT: unchanged
	my	$sref_besttree						= $_[4]	 ; # 4clan tree																													IN: defined; OUT: unchanged
	my	$sref_qtcode						= $_[5]	 ; # 4clan tree code (QT1, QT2, QT3)																							IN: defined; OUT: unchanged
	my	$sref_average						= $_[6]	 ; # mean|median																												IN: defined; OUT: unchanged
	my	$sref_support_value					= $_[7]	 ; # best tree support																											IN: defined; OUT: unchanged
	my	$sref_Nquartets_4clan				= $_[8]	 ; # total number of single 4clan quartet analyses (unfiltered)																	IN: defined; OUT: unchanged
	my	$sref_Nquartets_remain				= $_[9]	 ; # number of remaining quartets																								IN: defined; OUT: unchanged
	my	$href_subfolder_q					= $_[10] ; # hashpath to output-type specific subdirs of single quartet-pattern analyses; key: output-type (e.g. 'pdf'); value: pathdir	IN: defined; OUT: unchanged
	my	$href_textfilename_of_prtcode		= $_[11] ; # key1: textfile unique printcode (e.g. 'textA'); value: filename for print													IN: defined; OUT: unchanged
	
	#print "::&latex_data_sample_LQ2\n"; 
	
	######################################################################
	# called from subroutine(s)...
	# &calc_4clan_tree_support
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	# Add 4clan combinations and their codes to TXT outfile              #
	######################################################################
	my		$f_pinfo = $href_subfolder_q->{txt}."/".$href_textfilename_of_prtcode->{txtLQ2} ;
	open 	OUTtxt,	">>$f_pinfo" or die "\nTEXTFILE-ERROR: Cannot write to ", $f_pinfo, "!\n\n" ;
	
	############
	# print header only once
	unless ( $href_textfilename_of_prtcode->{txtLQ2hcode} ){
		
		print	OUTtxt	"Clade-Quartet", "\t", "Score", "\t", "Filter", "\t", "Best_Tree", "\t", "Tree_Code", "\t", "Average", "\t", "Support", "\t", "N_Species-Quartets_Available", "\t", "N_Species-Quartets_Analysed", "\n" ;  
		$href_textfilename_of_prtcode->{txtLQ2hcode}++
	}
	######################################################################
	
	######################################################################
	# Generate latex data line, listing LQ2 best 4clan tree support      #
	######################################################################
	my	$filter			  =	&recode_filter_printout( \$$sref_fm )	;
	my	$latex_info_line  =	$$sref_cl4."\t".$$sref_cm."\t".$filter."\t".$$sref_besttree."\t".$$sref_qtcode."\t".$$sref_average."\t".$$sref_support_value."\t".$$sref_Nquartets_4clan."\t".$$sref_Nquartets_remain."\n" ;
    push @{$href_hoh_data_of_lfile_of_prtcode->{latexLQ2}{$$sref_cm}}, $latex_info_line	; 
    
    print	OUTtxt	$latex_info_line ;
	close	OUTtxt ;
	######################################################################
}

sub latex_data_sample_LQ3{
	
	my	$href_hoh_linedata_of_latex_of_prtcode		= $_[0] ; # key1: Latex table unique printcode (e.g. 'latexC'); key2: *.tex filename; value: @list of assigned table lines						IN: defined; OUT: changed
	my	$href_qn_of_qfiltered_otu_of_4clan_of_ft	= $_[1] ; # key1: 4seq filtermethod (e.g. nap_risk1); key2: taxon; key3: 4clan|total; value: list of rejected qn's								IN: defined; OUT: unchanged
	my	$href_hoh_qremaining_of_otu_of_4clan_of_ft	= $_[2] ; # key1: 4seq filtermethod (e.g. nap_risk1); key2: taxon; key3: 4clan|total; value: N remaining quartets								IN: defined; OUT: unchanged
	my	$href_clan_of_seq							= $_[3] ; # key1: sequence-name; value: clan definition																							IN: defined; OUT: unchanged
	my	$href_seen_clan4							= $_[4] ; # key1: clan4; value: counter-number																									IN: defined; OUT: unchanged
	my	$href_subfolder_q							= $_[5] ; # hashpath to output-type specific subdirs of single quartet-pattern analyses; key: output-type (e.g. 'pdf'); value: pathdir			IN: defined; OUT: unchanged
	my	$href_textfilename_of_prtcode				= $_[6] ; # key1: textfile unique printcode (e.g. 'textA'); value: filename for print															IN: defined; OUT: unchanged
	my	$href_value_of_scr_par						= $_[7] ; # key1: sealion user available parameter options; value: optional parameter setting													IN: defined; OUT: unchanged
	my	$href_seen_filter_method					= $_[8] ; # key1: filter method incl OTU clipping approaches; value: seen +1																	IN: defined; OUT: unchanged
	
	#print "::&latex_data_sample_LQ3\n"; 
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	# Add 4clan combinations and their codes to TXT outfile              #
	######################################################################
	my		$f_pinfo = $href_subfolder_q->{txt}."/".$href_textfilename_of_prtcode->{txtLQ3} ;
	open 	OUTtxt,	">>$f_pinfo" or die "\nTEXTFILE-ERROR: Cannot write to ", $f_pinfo, "!\n\n" ;
	
	############
	# print header only once
	unless ( $href_textfilename_of_prtcode->{txtLQ3hcode} ){
		
		print	OUTtxt	"Species", "\t", "Clade", "\t", "Score", "\t", "Filter", "\t", "Clade-Quartet","\t", "N_Rejecetd_Species-Quartets", "\t", "N_Analysed_Species-Quartets", "\t", "N_Species-Quartets", "\t", "N_Total_Rejected_Species-Quartets", "\t", "N_Total_Remaining_Species-Quartets", "\t", "N_Total_Species-Quartets", "\n" ;  
		$href_textfilename_of_prtcode->{txtLQ3hcode}++
	}
	######################################################################
	
	######################################################################
	# Generate latex data line, listing N filtered taxa in quartets     ##
	######################################################################
	
	####################################
	# sample for each 4seq filter ($ftype) the number of taxon ($t) related quartet-rejections
	for	my $ftype ( sort keys %{$href_seen_filter_method} ){ #print "\n\t", $ftype, "\n"; #exit;
	#for	my $ftype ( sort keys %{$href_qn_of_qfiltered_otu_of_4clan_of_ft} ){ 
		
		my	@prts	=	split "_", $ftype ;
		my	$cm		=	$prts[0] ;
		my	$filter	=	&recode_filter_printout( \$ftype )	;
		
		for my $t ( sort keys %$href_clan_of_seq ){ 
			
			############
			# total number of rejected and remaining quartets
			my	$total_rej = 0 ;
			my	$total_rem = 0 ;
			
			for my $clan4 ( sort keys %{$href_qn_of_qfiltered_otu_of_4clan_of_ft->{$ftype}{$t}}){
			
				$total_rej += @{$href_qn_of_qfiltered_otu_of_4clan_of_ft->{$ftype}{$t}{$clan4}} ;
			}
			
			for my $clan4 ( sort keys %{$href_hoh_qremaining_of_otu_of_4clan_of_ft->{$ftype}{$t}}){
			
				$total_rem += $href_hoh_qremaining_of_otu_of_4clan_of_ft->{$ftype}{$t}{$clan4}
			}
			
			my $total_all = $total_rej + $total_rem ;
			############
			
			############
			# test print
			if ( $href_value_of_scr_par->{sys_parameter}{screen_prt} == 1 ){
				
				print "\n\t", $t, "\t", $total_all
			}
			############
			
			########################
			# print number of rejected and remaining quartets only if taxon is part of the actual 4clan combination
			for my $clan4 ( sort keys %$href_seen_clan4 ){
				
				if ( ( $href_qn_of_qfiltered_otu_of_4clan_of_ft->{$ftype}{$t}{$clan4}  )	||
					 ( $href_hoh_qremaining_of_otu_of_4clan_of_ft->{$ftype}{$t}{$clan4} )	){
					
					############
					# N rejected taxon in 4clan...
					my $n_rej_qn = 0 ;
					
					if ( $href_qn_of_qfiltered_otu_of_4clan_of_ft->{$ftype}{$t}{$clan4} ){
						
						############
						# 4clan specific rejections
						$n_rej_qn = @{$href_qn_of_qfiltered_otu_of_4clan_of_ft->{$ftype}{$t}{$clan4}} / 3 ; # counted for each tree, thus division by 3
						############
					}
					############
					
					############
					# N remaining taxon in 4clan...
					my $n_rem_qn = 0 ;
					
					if ( $href_hoh_qremaining_of_otu_of_4clan_of_ft->{$ftype}{$t}{$clan4} ){
						
						############
						# 4clan specific rejections
						$n_rem_qn = $href_hoh_qremaining_of_otu_of_4clan_of_ft->{$ftype}{$t}{$clan4} / 3 ; # counted for each tree, thus division by 3
						############
					}
					############
					
					############
					# N total quartets of 4clan
					my $n_all_qn = $n_rem_qn + $n_rej_qn ;
					############
					
					############
					# line preparation and storage
					my	$latex_info_line =	$t."\t".$href_clan_of_seq->{$t}."\t".$cm."\t".$filter."\t".$clan4."\t".$n_rej_qn."\t".$n_rem_qn."\t".$n_all_qn."\t".$total_rej."\t".$total_rem."\t".$total_all	;
					push @{$href_hoh_linedata_of_latex_of_prtcode->{latexLQ3}{$ftype}}, $latex_info_line	;
					
					print	OUTtxt	$latex_info_line, "\n";
					#print	"\n\t", $latex_info_line ; exit;
					############
				} 
				############
			}
			########################
		} #exit;
	}
	
	close OUTtxt ;
	######################################################################
}

sub latex_data_sample_LQ4{
	
	my	$href_hoh_linedata_of_latexfile_of_prtcode	= $_[0]	; # key1: Latex table unique printcode (e.g. 'latexC'); key2: *.tex filename; value: @list of assigned table lines					IN: defined; OUT: changed
	my	$href_seen_tax_of_4clan						= $_[1]	; # key1: 4clan code; key2: taxon; value: counter number																			IN: defined; OUT: unchanged
	my	$href_clan_of_seq							= $_[2]	; # key1: taxon; value: assigned clan																								IN: defined; OUT: unchanged 
	my	$href_algorithm_def_methods					= $_[3] ; # key1: methode type; key2: method code; value: integer																			IN: defined; OUT: unchanged
	my	$href_subfolder_q							= $_[4] ; # hashpath to output-type specific subdirs of single quartet-pattern analyses; key: output-type (e.g. 'pdf'); value: pathdir		IN: defined; OUT: unchanged
	my	$href_textfilename_of_prtcode				= $_[5] ; # key1: textfile unique printcode (e.g. 'textA'); value: filename for print														IN: defined; OUT: unchanged
			
	#print "::&latex_data_sample_LQ4\n"; 
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	# Add 4clan combinations and their codes to TXT outfile              #
	######################################################################
	my		$f_pinfo = $href_subfolder_q->{txt}."/".$href_textfilename_of_prtcode->{txtLQ4} ;
	open 	OUTtxt,	">>$f_pinfo" or die "\nTEXTFILE-ERROR: Cannot write to ", $f_pinfo, "!\n\n" ;
	
	############
	# print header only once
	unless ( $href_textfilename_of_prtcode->{txtLQ4hcode} ){
		
		print	OUTtxt	"Species", "\t", "Clade", "\t", "Score", "\t","Clade-Quartet","\t","N_Species-Quartet_Participations", "\t", "N_Total_Species-Quartet_Participations", "\n" ;  
		$href_textfilename_of_prtcode->{txtLQ4hcode}++
	}
	######################################################################
	
	######################################################################
	# Generate latex data line, listing N unfiltert taxa in quartets    ##
	######################################################################
	
	###############
	# data sampling for latextable LQT2, for each quartet-calculation method (nap, napv, napr) separately collected
	for my $cm ( sort keys %{$href_algorithm_def_methods->{original}} ){
		
		for my $t ( sort keys %$href_clan_of_seq ){
			
			my $qn_total = 0 ;
			for my $clan4 ( sort keys %$href_seen_tax_of_4clan ){
				
				if ( $href_seen_tax_of_4clan->{$clan4}{$t} ){ $qn_total += $href_seen_tax_of_4clan->{$clan4}{$t} }
			}
		
			unless ( $qn_total == 0 ){
				
				for my $clan4 ( sort keys %$href_seen_tax_of_4clan ){
					
					if ( $href_seen_tax_of_4clan->{$clan4}{$t} ){ 
						
						my $qn_4clan = $href_seen_tax_of_4clan->{$clan4}{$t} ;
						push @{$href_hoh_linedata_of_latexfile_of_prtcode->{latexLQ4}{$cm}}, $t."\t".$href_clan_of_seq->{$t}."\t".$cm."\t".$clan4."\t".$qn_4clan."\t".$qn_total ;
						
						print	OUTtxt $t, "\t", $href_clan_of_seq->{$t}, "\t", $cm, "\t", $clan4, "\t", $qn_4clan, "\t", $qn_total, "\n" ;
					}
				}
			}
		
			else{ 
				
				push @{$href_hoh_linedata_of_latexfile_of_prtcode->{latexLQ4}{$cm}}, $t."\t".$href_clan_of_seq->{$t}."\t".$cm."\t"."NA"."\t"."NA"."\t"."NA" ;
				
				print	OUTtxt $t, "\t", $href_clan_of_seq->{$t}, "\t", $cm, "\t", "NA", "\t", "NA", "\t", "NA", "\n" ;
			}
		}
	}
	###############
	
	######################################################################
}

sub latex_data_sample_LQ5{
	
	my	$href_hoh_linedata_of_latexfile_of_prtcode	= $_[0]	; # key1: Latex table unique printcode (e.g. 'latexC'); key2: *.tex filename; value: @list of assigned table lines					IN: defined; OUT: changed
	my	$href_seen_tax_of_4clan						= $_[1]	; # key1: 4clan code; key2: taxon; value: counter number																			IN: defined; OUT: unchanged
	my	$href_clan_of_seq							= $_[2]	; # key1: taxon; value: assigned clan																								IN: defined; OUT: unchanged 
	my	$href_algorithm_def_methods					= $_[3] ; # key1: methode type; key2: method code; value: integer																			IN: defined; OUT: unchanged
	my	$href_subfolder_q							= $_[4] ; # hashpath to output-type specific subdirs of single quartet-pattern analyses; key: output-type (e.g. 'pdf'); value: pathdir		IN: defined; OUT: unchanged
	my	$href_textfilename_of_prtcode				= $_[5] ; # key1: textfile unique printcode (e.g. 'textA'); value: filename for print														IN: defined; OUT: unchanged
			
	#print "::&latex_data_sample_LQ5\n"; 
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	# Add 4clan combinations and their codes to TXT outfile              #
	######################################################################
	my		$f_pinfo = $href_subfolder_q->{txt}."/".$href_textfilename_of_prtcode->{txtLQ5} ;
	open 	OUTtxt,	">>$f_pinfo" or die "\nTEXTFILE-ERROR: Cannot write to ", $f_pinfo, "!\n\n" ;
	
	############
	# print header only once
	unless ( $href_textfilename_of_prtcode->{txtLQ5hcode} ){
		
		print	OUTtxt	"Species", "\t", "Clade", "\t", "Score", "\n" ;  
		$href_textfilename_of_prtcode->{txtLQ5hcode}++
	}
	######################################################################
	
	######################################################################
	# Generate latex data line, listing non quartet participating taxa  ##
	######################################################################
	
	###############
	# data sampling for latextable LQ5, for each quartet-calculation method (nap, napv, napr) separately collected
	for my $cm ( sort keys %{$href_algorithm_def_methods->{original}} ){
		
		my $flag = 0 ;
		
		TAXON:
		for my $t ( sort keys %$href_clan_of_seq ){
			
			for my $clan4 ( sort keys %$href_seen_tax_of_4clan ){
				
				if ( $href_seen_tax_of_4clan->{$clan4}{$t} ){ next TAXON }
			}	
			
			push @{$href_hoh_linedata_of_latexfile_of_prtcode->{latexLQ5}{$cm}}, $t."\t".$href_clan_of_seq->{$t}."\t".$cm ;
			
			print	OUTtxt $t, "\t", $href_clan_of_seq->{$t}, "\t", $cm, "\n" ; $flag = 1
		}
		
		if ( $flag == 0 ){
			
			 print	OUTtxt "NA", "\t", "NA", "\t", $cm, "\n" ;
			 push @{$href_hoh_linedata_of_latexfile_of_prtcode->{latexLQ5}{$cm}}, "NA"."\t"."NA"."\t".$cm  ;
		}
	}
	###############
	
	######################################################################
}

sub latex_data_sample_LQ6{
	
	my	$href_prtpath						= $_[0] ; # key1: R grafic unique printcode (e.g. 'metafileQ1'); key2: *.tsv filename; value: *tsv filename data addings 	IN: defined; OUT: unchanged
	my	$href_hoh_data_of_tex_of_prcode		= $_[1] ; # key1: Latex table unique printcode (e.g. 'latexC'); key2: *.tex filename; value: @list of assigned table lines	IN: defined; OUT: changed
	my	$href_text_of_prcode				= $_[2] ; # key1: textfile unique printcode (e.g. 'textA'); value: filename for print										IN: defined; OUT: changed
	my	$sref_supportfilter					= $_[3] ; # actual filtermethod																								IN: defined; OUT: unchanged
	my	$sref_prtlines_LQ6					= $_[4] ; # tex data lines																									IN: defined; OUT: unchanged
	
	######################################################################
	# called from subroutine(s)...
	# &filter_threshold_optimization
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	# Add 4clan combinations and optimized thresholds to TXT outfile     #
	######################################################################
	my		$f_pinfo = $href_prtpath->{txt}."/".$href_text_of_prcode->{txtLQ6} ;
	open 	OUTtxt,	">>$f_pinfo" or die "\nTEXTFILE-ERROR: Cannot write to ", $f_pinfo, "!\n\n" ;
	
	############
	# print header only once
	unless ( $href_text_of_prcode->{txtLQ6hcode} ){
		
		print	OUTtxt	"Score", "\t", "Filter", "\t","Clade-Quartet", "\t", "Optimized_Filter_Threshold_Tl", "\t","N_Remaining_Species-Quartets", "\t", "N_Total_Species-Quartets", "\n" ;  
		$href_text_of_prcode->{txtLQ6hcode}++
	}
	######################################################################
	
	######################################################################
	# Generate latex data line, and print txtfile line                   #
	######################################################################
	print OUTtxt $$sref_prtlines_LQ6 ;
	
	push @{$href_hoh_data_of_tex_of_prcode->{latexLQ6}{$$sref_supportfilter}}, $$sref_prtlines_LQ6 ;
	######################################################################
}

sub latex_data_sample_LQ7{
	
	my	$href_hoh_data_of_tex_of_prcode		= $_[0] ; # key1: Latex table unique printcode (e.g. 'latexC'); key2: *.tex filename; value: @list of assigned table lines	IN: defined; OUT: changed
	my	$href_text_of_prcode				= $_[1] ; # key1: textfile unique printcode (e.g. 'textA'); value: filename for print										IN: defined; OUT: changed
	my	$href_prtpath						= $_[2] ; # key1: R grafic unique printcode (e.g. 'metafileQ1'); key2: *.tsv filename; value: *tsv filename data addings 	IN: defined; OUT: unchanged
	my	$href_hoh_N_qrej_of_clan4_of_ft		= $_[3] ; # key1: 4seq filtermethod (e.g. nap_risk1); key2: 4clan combination code; value:  N rejected quartets				IN: defined; OUT: unchanged
	my	$href_seen_nq_of_clan4				= $_[4] ; # key1: clan4; value: N analyseed quartets*3 (counted three times, each time for each quartet tree)				IN: defined; OUT: unchanged
	my	$href_seen_filter_method			= $_[5] ; # key1: filter method incl OTU clipping approaches; value: seen +1												IN: defined; OUT: unchanged
	
	#print "::&print_tsv_Q6\n";
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	#####################################################################################
	###### print N rejected quartets for each quartet-filter to TEX & TXT LQ7         ###
	#####################################################################################
	for my $cm_fm ( sort keys %$href_seen_filter_method ){
		
		my	$filter	=	&recode_filter_printout( \$cm_fm )	;
		my	@prts 	=	split "_", $cm_fm ;
		my	$cm		=	shift @prts ;
		my	$fm		=	join "_", @prts ;
		
		########
		( my	$f_pinfo = $href_prtpath->{txt}."/".$href_text_of_prcode->{txtLQ7} ) =~ s/.txt$/_${cm}.txt/;
		open 	OUTtxt,	">>$f_pinfo" or die "\nTEXTFILE-ERROR: Cannot write to ", $f_pinfo, "!\n\n" ;
		########
		
		########
		# print header only once to each clanfile corresponding tsv table file
		unless ( $href_text_of_prcode->{txtLQ7hcode} ){
		
			print 	OUTtxt		"Score"						,
								"\tFilter"					,
								"\tClade_Quartet"			,
								"\tStatus"					,
								"\tNumber_Species_Quartets"	,
								"\n"						;
		
			$href_text_of_prcode->{txtLQ7hcode}++
		}
		########
		
		########
		# print table data
		for my $clan4 ( sort keys %$href_seen_nq_of_clan4 ){
			
			my $qtotal = $href_seen_nq_of_clan4->{$clan4} ;
			
			my $qrejec = 0 ;
			if ( $href_hoh_N_qrej_of_clan4_of_ft->{$cm_fm}{$clan4} ){ $qrejec = $href_hoh_N_qrej_of_clan4_of_ft->{$cm_fm}{$clan4} / 3 }
			
			print	OUTtxt		$cm				,	# Score
								"\t",$filter	,	# Filter
								"\t",$clan4		,	# Clan4
								"\t","initially",	# Nq_c
								"\t",$qtotal	,	# Nq_v
								"\n",
								$cm				,	# Score
								"\t",$filter	,	# Filter
								"\t",$clan4		,	# Clan4
								"\t","rejected"	,	# Nq_c
								"\t",$qrejec	,	# Nq_v
								"\n"			;
			
			my $latex_prt_line1 = $cm."\t".$filter."\t".$clan4."\tinitially"."\t".$qtotal ;
			push @{$href_hoh_data_of_tex_of_prcode->{latexLQ7}{$cm}}, $latex_prt_line1 ;
			
			my $latex_prt_line2 = $cm."\t".$filter."\t".$clan4."\trejected"."\t".$qrejec ;
			push @{$href_hoh_data_of_tex_of_prcode->{latexLQ7}{$cm}}, $latex_prt_line2 ;
		}
		########
		
		########
		# close tsv filehandle
		close	OUTtxt ; #exit;
		########
	}
	#####################################################################################
}

sub latex_data_sample_LQT1{
	
	my	$sref_ftype									= $_[0]	; # cm+fm method combination (e.g. nap_risk1)																								IN: defined; OUT: unchanged
	my	$prtinf										= $_[1]	; # original filter code (e.g. 'nap_risk1')																									IN: defined; OUT: unchanged
	my	$href_hoh_linedata_of_latexfile_of_prtcode	= $_[2]	; # key1: Latex table unique printcode (e.g. 'latexC'); key2: *.tex filename; value: @list of assigned table lines							IN: defined; OUT: changed
	my	$sref_clan4									= $_[3]	; # clan4-combination																														IN: defined; OUT: unchanged
	my	$aref_trees_q								= $_[4]	; # list of 4clan-trees in same order as support in @normed_supp																			IN: defined; OUT: unchanged
	my	$href_subf_q								= $_[5]	; # hashpath to output-type specific subdirs of single quartet-pattern analyses; key: output-type (e.g. 'pdf'); value: pathdir				IN: defined; OUT: unchanged
	my	$href_textfile_of_prtcode					= $_[6] ; # key1: textfile unique printcode (e.g. 'textA'); value: filename for print																IN: defined; OUT: unchanged
	
	#print "::&latex_data_sample_LQT1\n"; 
	
	######################################################################
	# called from subroutine(s)...
	# &calc_4clan_tree_support
	# &taxon_clipping
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	# Add 4clan combinations and their codes to TXT outfile              #
	######################################################################
	my		$f_pinfo = $href_subf_q->{txt}."/".$href_textfile_of_prtcode->{txtLQT1} ;
	open 	OUTtxt,	">>$f_pinfo" or die "\nTEXTFILE-ERROR: Cannot write to ", $f_pinfo, "!\n\n" ;
	
	############
	# print header only once
	unless ( $href_textfile_of_prtcode->{txtLQT1hcode} ){
		
		print	OUTtxt	"Score", "\t", "Filter", "\t", "Clade-Quartet", "\t", "Tree 'QT1'", "\t", "Tree 'QT2'", "\t", "Tree 'QT3'", "\n" ;  
		$href_textfile_of_prtcode->{txtLQT1hcode}++
	}
	######################################################################
	
	######################################################################
	# Generate latex data line, listing QT of defined triangle corners  ##
	######################################################################
	
	###############
	# data sampling for latextable LQT1 (ternaryplot-legend)
	my	@prts			 =	split "_", $$sref_ftype ;
	my	$cm				 =	$prts[0] ;
	my	$filter			 =	&recode_filter_printout( \$$sref_ftype )	;
			
	my	$latex_info_line = $cm."\t".$filter."\t".$$sref_clan4."\t".$aref_trees_q->[0]."\t".$aref_trees_q->[1]."\t".$aref_trees_q->[2]."\n" ;
	push 	@{$href_hoh_linedata_of_latexfile_of_prtcode->{latexLQT1}{$$prtinf}}, $latex_info_line	;
	
	print	OUTtxt $latex_info_line ;
	close	OUTtxt ;
	#print "\n\tlinfo\t", $latex_info_line ;
	###############
	
	######################################################################
}

sub latextable_default{

	my $aref_lol_value_line_column	= $_[0]								; # two-dimensional matrix: $mulitarray [lines][columns]
	my $aref_list_headers			= $_[1]								; # have to be sorted before giving to subroutine -> (column1, column2, column3, column4) header of table
	my $sref_filename				= $_[2]								; # new filename for latex table
	my $sref_caption				= $_[3]								; # contain header of the table

	my $N_cols						= @$aref_list_headers				; # saves N columns
	#my $header						= join " & ", @$aref_list_headers	;

	my @temporary_splitted	=	split ('/', $$sref_filename		);
	my $temporary			=	pop			@temporary_splitted	 ;
	my $path				=	join  ("/", @temporary_splitted	);
	#  print "Amount ",$amount, "\n";exit;
	open (my $FH, '>', $temporary) or die "Cannot open $$sref_filename : $!\n";

	print {$FH}		"\\documentclass[10pt,oneside,a4paper]{article}\n"			,
					"\n"														,
					"\\usepackage{longtable}\n"									,
					"\\usepackage[table]{xcolor}\n"								,
					"\\usepackage{helvet}\n"									,
					"\\usepackage{pdflscape}\n"									,
					"\\usepackage{booktabs}\n"									,
					"\n"														,
					"\\renewcommand{\\familydefault}{\\sfdefault}\n"			,
					"\n"														,
					"\n"														,
					"\\definecolor{cell}{RGB}{220,230,240}\n"					,
					"\\definecolor{line}{RGB}{220,20,60}\n"						,
					"\n"														,
					"\n"														,
					"\\begin{document}\\small\n"								,
					"\n"														,
					"\n"														,
					"\\arrayrulecolor{line}\n"									,
					"\n"														,
					"\\begin{longtable}[c]{"	, "l" x $N_cols	, "}\n"			,
					"\n"														,
					"\\caption[VARIABEL]{\\label{tab:parameter}",$$sref_caption,"}\\\\\n"	,
					"\\hline \\toprule\n"											,
					"\n"														;

	for	my $i (	0 .. $N_cols-1 ){
		(my $current_header =$aref_list_headers->[$i]) =~ s/_/\\_/g;
		print {$FH}	'&' if $i > 0	;
		print {$FH}	'\textbf{', $current_header,'}';
	}

	print {$FH}		"\\\\\n", '\cmidrule(lr){1-', $N_cols, '}'									,
					"\n"																		,
					'\\endfirsthead'															,
					"\n"																		,
					"\n"																		,
					'\\multicolumn{', $N_cols, '}{r}{Continues on next page(s).}\\\\'			,
					"\n"																		,
					'\\endfoot'																	,
					"\n"																		,
					'\\hline'																	,
					"\n"																		,
					"\n"																		,
					'\\multicolumn{', $N_cols, '}{r}{End of table \\ref{tab:parameter}}\\\\'	,
					"\n"																		,
					'\\endlastfoot'																,
					"\n"																		,
					"\n"																		;


	my $N_lines = @$aref_lol_value_line_column ; #saves the amount of lines

	for	my $i	( 0 .. $N_lines-1 ){

		print {$FH} '\rowcolor{cell}' if $i % 2 == 0;

		for	my $x (	0 .. $N_cols-1	){
			
			(my $character = $aref_lol_value_line_column->[$i][$x]) =~ s/_/\\_/g;
			
			if ( $x == 0){
				
				print{$FH} '\textcolor{black}{', $character, '}';
			}
			
			else{
				
				print{$FH} ' & \textcolor{black}{', $character, '}';
			}
		}

		print{$FH} '\\\\', "\n"
	}

	print {$FH}		'\bottomrule'		,
					"\n"				,
					'\end{longtable}'	,
					"\n"				,
					'\end{document}'	;
	close $FH	;

	
	####################################
	# must be executed twice, otherwise table format can be incorrect
	my $ftemp = "latex.tmp" ;
	system ("pdflatex -interaction=batchmode $temporary >$ftemp") ;
	system ("pdflatex -interaction=batchmode $temporary >$ftemp") ;
	####################################
	
	(my $current = $temporary) =~ s/\.tex$//;
	for my $new_files (<$current*>){
		
		if ( $new_files =~ /.tex$|.pdf$/ ){
			
			my $new_new_files = $path ."/". $new_files;
			move ($new_files, $new_new_files);
		}
		else{unlink $new_files}
	}
	#
	
}

sub latextable_landscape{

	my $aref_lol_value_line_column	= $_[0]								; # two-dimensional matrix: $mulitarray [lines][columns]
	my $aref_list_headers			= $_[1]								; # have to be sorted before giving to subroutine -> (column1, column2, column3, column4) header of table
	my $sref_filename				= $_[2]								; # new filename for latex table
	my $sref_caption				= $_[3]								; # contain header of the table

	my $N_cols						= @$aref_list_headers				; # saves N columns
	#my $header						= join " & ", @$aref_list_headers	;

	my @temporary_splitted	=	split ('/', $$sref_filename		);
	my $temporary			=	pop			@temporary_splitted	 ;
	my $path				=	join  ("/", @temporary_splitted	);
	#  print "Amount ",$amount, "\n";exit;
	open (my $FH, '>', $temporary) or die "Cannot open $$sref_filename : $!\n";

	print {$FH}		"\\documentclass[10pt,oneside,a4paper]{article}\n"			,
					"\n"														,
					"\\usepackage{longtable}\n"									,
					"\\usepackage[table]{xcolor}\n"								,
					"\\usepackage{helvet}\n"									,
					"\\usepackage{pdflscape}\n"									,
					"\\usepackage{booktabs}\n"									,
					"\n"														,
					"\\renewcommand{\\familydefault}{\\sfdefault}\n"				,
					"\n"														,
					"\n"														,
					"\\definecolor{cell}{RGB}{220,230,240}\n"					,
					"\\definecolor{line}{RGB}{220,20,60}\n"						,
					"\n"														,
					"\n"														,
					"\\begin{document}\n"										,
					"\\begin{landscape}\\small"									,
					"\n"														,
					"\n"														,
					"\\arrayrulecolor{line}\n"									,
					"\n"														,
					"\\begin{longtable}[c]{"	, "l" x $N_cols	, "}\n"				,
					"\n"														,
					"\\caption[VARIABEL]{\\label{tab:parameter}",$$sref_caption,"}\\\\\n"	,
					"\\hline \\toprule\n"											,
					"\n"														;

	for	my $i (	0 .. $N_cols-1 ){
		(my $current_header =$aref_list_headers->[$i]) =~ s/_/\\_/g;
		print {$FH}	'&' if $i > 0	;
		print {$FH}	'\textbf{', $current_header,'}';
	}

	print {$FH}		"\\\\\n", '\cmidrule(lr){1-', $N_cols, '}'									,
					"\n"																		,
					'\\endfirsthead'															,
					"\n"																		,
					"\n"																		,
					'\\multicolumn{', $N_cols, '}{r}{Continues on next page(s).}\\\\'			,
					"\n"																		,
					'\\endfoot'																	,
					"\n"																		,
					'\\hline'																	,
					"\n"																		,
					"\n"																		,
					'\\multicolumn{', $N_cols, '}{r}{End of table \\ref{tab:parameter}}\\\\'	,
					"\n"																		,
					'\\endlastfoot'																,
					"\n"																		,
					"\n"																		;


	my $N_lines = @$aref_lol_value_line_column ; #saves the amount of lines

	for	my $i	( 0 .. $N_lines-1 ){

		print {$FH} '\rowcolor{cell}' if $i % 2 == 0;

		for	my $x (	0 .. $N_cols-1	){
				
			(my $character = $aref_lol_value_line_column->[$i][$x]) =~ s/_/\\_/g;
			
			if ($x == 0){
				
				print{$FH} '\textcolor{black}{', $character, '}';
			}
			
			else{
				
				print{$FH} ' & \textcolor{black}{', $character, '}';
			}
		}

		print{$FH} '\\\\', "\n"
	}

	print {$FH}		'\bottomrule'		,
					"\n"				,
					'\end{longtable}'	,
					'\end{landscape}'	,
					"\n"				,
					'\end{document}'	;
	close $FH	;

	
	####################################
	# must be executed twice, otherwise table format can be incorrect
	my $ftemp = "latex.tmp" ;
	system ("pdflatex -interaction=batchmode $temporary >$ftemp ") ;
	system ("pdflatex -interaction=batchmode $temporary >$ftemp ") ;
	####################################
	
	(my $current = $temporary) =~ s/\.tex$//;
	for my $new_files (<$current*>){
		
		if ( $new_files =~ /.tex$|.pdf$/ ){
			
			my $new_new_files = $path ."/". $new_files;
			move ($new_files, $new_new_files);
		}
		else{unlink $new_files}
	}
	#
	
}

sub latex_quartet_analyses{
	
	my	$href_print_latex									= $_[0]	; # Print control of single latex tables																						IN: defined; OUT: unchanged	
	my	$href_hoh_linedata_of_latexfilename_of_printcode	= $_[1]	; # key1: Latex table unique printcode (e.g. 'latexC'); key2: *.tex filename; value: @list of assigned table lines				IN: defined; OUT: changed
	my	$href_hol_pdf_of_type								= $_[2]	; # key: result type; value: list of associated pdf's in order of construction date												IN: defined; OUT: changed
	my	$aref_pdf_pattern									= $_[3]	; # list variable of single  pdf's																								IN: defined; OUT: changed	
	my	$href_subfolder_q									= $_[4]	; # key1: file type ('svg', 'pdf', 'tex', 'tsv') ; value: printdir																IN: defined; OUT: unchanged	
	my	$href_subsubfolder_q								= $_[5] ; # key1: pdf related subfolder codes ('pdf_single', 'pdf_sumed'...)															IN: defined; OUT: unchanged	
	my	$href_value_of_scr_par								= $_[6] ; # key1: sealion user available parameter options; value: optional parameter setting											IN: defined; OUT: unchanged	
	my	$href_specified_support_method_code					= $_[7] ; # list of defined support methods (nap, napv ...)																				IN: defined; OUT: changed	
	
	######################################################################
	## print latex table info about quartet data                        ##
	######################################################################
	print	"\n\t---------------------------------------------------\n",
			"\twriting LATEX file(s)...\n" ;
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	# definition of graphic-info assigned output-subcategories for pdf summarization
	my	$main	= 'mainTable'		;	# file_suffix for pdf summarization of main overall 4clan result-files
	######################################################################
	
	##################################################################################### LP1 works ## main
	## latex landscape table of sealion defined parameter                              ##
	#####################################################################################
	if	( 	( $href_print_latex->{table_parameter_sealion}						== 1 ) &&
			( $href_value_of_scr_par->{pdf_prtResult}{main} 					== 1 ) ){
		
		my	$latex_filename		=	$href_subfolder_q->{tex}."/LP1_SeaLion_parameter_setup.tex" ;
		print	"\n\t\t...", $latex_filename ;
		
		my	$latex_caption		=	"LP1 -- SeaLion specified parameter." ;
		my	@latex_fileheader	=	( "SeaLion Option\t", "Parameter" ) ;
		
		my	@lol_latex_linedata ;
		my	$ycounter	= 0;
		for my $latexline ( @{$href_hoh_linedata_of_latexfilename_of_printcode->{latexLP1}} ){
				
			my	@latex_data = split "\t", $latexline ;
			
			my	$xcounter	= 0 ;
			for	my $single_d ( @latex_data ){ $lol_latex_linedata[$ycounter][$xcounter] = $single_d ; $xcounter++ }
			$ycounter++
		}
		
		&latextable_default(
		
			\@lol_latex_linedata	,	# two-dimensional matrix: $mulitarray [lines][columns]
			\@latex_fileheader		,	# latex header line (column sorted)
			\$latex_filename		,	# latex output name
			\$latex_caption			,	# latex table caption
		) ;
		
		for my $ft ( sort keys %{$href_specified_support_method_code->{original}} ){
			
			my	$ftype	=	$main."_".$ft ;
			push @{$href_hol_pdf_of_type->{$ftype}}, $href_subfolder_q->{tex}."/LP1_SeaLion_parameter_setup.pdf" ;
		}
	}
	#####################################################################################
	#####################################################################################
	
	##################################################################################### LP2 works ## main
	## latex landscape table of defined clans and their assigned codes                 ##
	#####################################################################################
	if	( 	( $href_print_latex->{table_clancodes_sealion}						== 1 ) &&
			( $href_value_of_scr_par->{pdf_prtResult}{main} 					== 1 ) ){
		
		my	$latex_filename		=	$href_subfolder_q->{tex}."/LP2_SeaLion_clade_codes.tex" ;
		print	"\n\t\t...", $latex_filename ;
		
		my	$latex_caption		=	"LP2 -- Defined clades and assigned codes." ;
		my	@latex_fileheader	=	( "Clade\t", "Code" ) ;
		
		my	@lol_latex_linedata ;
		my	$ycounter	= 0;
		for my $latexline ( @{$href_hoh_linedata_of_latexfilename_of_printcode->{latexLP2}} ){
				
			my	@latex_data = split "\t", $latexline ;
			
			my	$xcounter	= 0 ;
			for	my $single_d ( @latex_data ){ $lol_latex_linedata[$ycounter][$xcounter] = $single_d ; $xcounter++ }
			$ycounter++
		}
		
		&latextable_default(
		
			\@lol_latex_linedata	,	# two-dimensional matrix: $mulitarray [lines][columns]
			\@latex_fileheader		,	# latex header line (column sorted)
			\$latex_filename		,	# latex output name
			\$latex_caption			,	# latex table caption
		) ;
		
		for my $ft ( sort keys %{$href_specified_support_method_code->{original}} ){
			
			my	$ftype	=	$main."_".$ft ;
			push @{$href_hol_pdf_of_type->{$ftype}}, $href_subfolder_q->{tex}."/LP2_SeaLion_clade_codes.pdf" ;
		}
	}
	#####################################################################################
	#####################################################################################
	
	##################################################################################### LP3 works ## main
	## latex landscape table of defined 4clan combi and their assigned codes           ##
	#####################################################################################
	if	( 	( $href_print_latex->{table_4clancodes_sealion}						== 1 ) &&
			( $href_value_of_scr_par->{pdf_prtResult}{main} 					== 1 ) ){
		
		my	$latex_filename		=	$href_subfolder_q->{tex}."/LP3_SeaLion_clade-quartet_codes.tex" ;
		print	"\n\t\t...", $latex_filename ;
		
		my	$latex_caption		=	"LP3 -- Analysed clade-quartets and assigned codes." ;
		my	@latex_fileheader	=	( "Clade-Quartet\t", "Code" ) ;
		
		my	@lol_latex_linedata ;
		my	$ycounter	= 0;
		for my $latexline ( @{$href_hoh_linedata_of_latexfilename_of_printcode->{latexLP3}} ){
				
			my	@latex_data = split "\t", $latexline ;
			
			my	$xcounter	= 0 ;
			for	my $single_d ( @latex_data ){ $lol_latex_linedata[$ycounter][$xcounter] = $single_d ; $xcounter++ }
			$ycounter++
		}
		
		&latextable_default(
		
			\@lol_latex_linedata	,	# two-dimensional matrix: $mulitarray [lines][columns]
			\@latex_fileheader		,	# latex header line (column sorted)
			\$latex_filename		,	# latex output name
			\$latex_caption			,	# latex table caption
		) ;
		
		for my $ft ( sort keys %{$href_specified_support_method_code->{original}} ){
			
			my	$ftype	=	$main."_".$ft ;
			push @{$href_hol_pdf_of_type->{$ftype}}, $href_subfolder_q->{tex}."/LP3_SeaLion_clade-quartet_codes.pdf" ;
		}
	}
	#####################################################################################
	#####################################################################################
	
	##################################################################################### LRC1 works ## main
	## latex landscape table of best mclan tree support (filtered|unfiltered)          ##
	#####################################################################################
	if	( 	( $href_print_latex->{table_besttreesupport_mclan}					== 1 ) &&
			( $href_value_of_scr_par->{pdf_prtResult}{main} 					== 1 ) ){
		
		for my $ft ( sort keys %{$href_hoh_linedata_of_latexfilename_of_printcode->{latexLRC1}} ){
			
			my	$latex_filename		=	$href_subfolder_q->{tex}."/LRC1_SeaLion_rooted-clades_best-tree_support_".$ft.".tex" ;
			print	"\n\t\t...", $latex_filename ;
			
			my	$latex_caption		=	"LRC1 -- Best rooted clade-tree support for specified quartet analyses." ;
			my	@latex_fileheader	=	(
				
				"Score\t", "Filter\t", "Best Rooted-Clade Tree\t", "Average\t", "Support"
			) ;
		
			my	@lol_latex_linedata ;
			my	$ycounter	= 0;
			for my $latexline ( @{$href_hoh_linedata_of_latexfilename_of_printcode->{latexLRC1}{$ft}} ){
				
				my	@latex_data = split "\t", $latexline ;
				
				my	$xcounter	= 0 ;
				for	my $single_d ( @latex_data ){ $lol_latex_linedata[$ycounter][$xcounter] = $single_d ; $xcounter++ }
				$ycounter++
			}
			
			&latextable_landscape(
				
				\@lol_latex_linedata	,	# two-dimensional matrix: $mulitarray [lines][columns]
				\@latex_fileheader		,	# latex header line (column sorted)
				\$latex_filename		,	# latex output name
				\$latex_caption			,	# latex table caption
			) ;
			
			my	$ftype	=	$main."_".$ft ;
			push @{$href_hol_pdf_of_type->{$ftype}}, $href_subfolder_q->{tex}."/LRC1_SeaLion_rooted-clades_best-tree_support_".$ft.".pdf" ;
		}
	}
	#####################################################################################
	#####################################################################################
	
	##################################################################################### LRC2 works ## main
	## latex landscape table of best mclan tree support (filtered|unfiltered)          ##
	#####################################################################################
	if	( 	( $href_print_latex->{table_singletreesupport_mclan}				== 1 ) &&
			( $href_value_of_scr_par->{pdf_prtResult}{main} 					== 1 ) ){
		
		for my $ft ( sort keys %{$href_hoh_linedata_of_latexfilename_of_printcode->{latexLRC2}} ){
			
			my	$latex_filename		=	$href_subfolder_q->{tex}."/LRC2_SeaLion_rooted-clades_single-tree_support_".$ft.".tex" ;
			print	"\n\t\t...", $latex_filename ;
			
			my	$latex_caption		=	"LRC2 -- Best-three rooted clade-tree support for specified quartet analyses." ;
			my	@latex_fileheader	=	(
				
				"Score\t", "Filter\t", "Average\t", "Best\t", "Support Best\t", "2nd Best\t", "Support 2nd\t" ,"3rd Best\t" ,"Support 3rd\t" , "Distance Best 2nd"
			) ;
		
			my	@lol_latex_linedata ;
			my	$ycounter	= 0;
			for my $latexline ( @{$href_hoh_linedata_of_latexfilename_of_printcode->{latexLRC2}{$ft}} ){
				
				my	@latex_data = split "\t", $latexline ;
				
				my	$xcounter	= 0 ;
				for	my $single_d ( @latex_data ){ $lol_latex_linedata[$ycounter][$xcounter] = $single_d ; $xcounter++ }
				$ycounter++
			}
			
			&latextable_landscape(
				
				\@lol_latex_linedata	,	# two-dimensional matrix: $mulitarray [lines][columns]
				\@latex_fileheader		,	# latex header line (column sorted)
				\$latex_filename		,	# latex output name
				\$latex_caption			,	# latex table caption
			) ;
			
			my	$ftype	=	$main."_".$ft ;
			push @{$href_hol_pdf_of_type->{$ftype}}, $href_subfolder_q->{tex}."/LRC2_SeaLion_rooted-clades_single-tree_support_".$ft.".pdf" ;
		}
	}
	#####################################################################################
	#####################################################################################
	
	##################################################################################### LQ7 ## main
	## latex landscape table of rejected species-quartets (filtered|unfiltered)        ##
	#####################################################################################
	if	( 	( $href_print_latex->{table_number_analysed_quartets}				== 1 ) &&
			( $href_value_of_scr_par->{pdf_prtResult}{main} 					== 1 ) ){
		
		for my $ft ( sort keys %{$href_hoh_linedata_of_latexfilename_of_printcode->{latexLQ7}} ){
			
			my	$latex_filename		=	$href_subfolder_q->{tex}."/LQ7_SeaLion_species-quartet_numbers_per_clade-quartet_".$ft.".tex" ;
			print	"\n\t\t...", $latex_filename ;
			
			my	$latex_caption		=	"LQ7 -- Number of rejected species-quartets in single clade-quartets after quartet filtering." ;
			my	@latex_fileheader	=	(
				
				"Score\t", "Filter\t", "Clade_Quartet\t" , "Status\t" , "Number_Species_Quartets"
			) ;
		
			my	@lol_latex_linedata ;
			my	$ycounter	= 0;
			for my $latexline ( @{$href_hoh_linedata_of_latexfilename_of_printcode->{latexLQ7}{$ft}} ){
				
				my	@latex_data = split "\t", $latexline ;
				
				my	$xcounter	= 0 ;
				for	my $single_d ( @latex_data ){ $lol_latex_linedata[$ycounter][$xcounter] = $single_d ; $xcounter++ }
				$ycounter++
			}
			
			&latextable_landscape(
				
				\@lol_latex_linedata	,	# two-dimensional matrix: $mulitarray [lines][columns]
				\@latex_fileheader		,	# latex header line (column sorted)
				\$latex_filename		,	# latex output name
				\$latex_caption			,	# latex table caption
			) ;
			
			my	$ftype	=	$main."_".$ft ;
			push @{$href_hol_pdf_of_type->{$ftype}}, $href_subfolder_q->{tex}."/LQ7_SeaLion_species-quartet_numbers_per_clade-quartet_".$ft.".pdf" ;
		}
	}
	#####################################################################################
	#####################################################################################
	
	##################################################################################### LQ1 works ## main
	## latex landscape table of single (filtered|unfiltered) 4clan tree support        ##
	#####################################################################################
	if	( 	( $href_print_latex->{table_treesupport_4clan}						== 1 ) &&
			( $href_value_of_scr_par->{pdf_prtResult}{main} 					== 1 ) ){
		
		for my $ft ( sort keys %{$href_hoh_linedata_of_latexfilename_of_printcode->{latexLQ1}} ){
			
			my	$latex_filename		=	$href_subfolder_q->{tex}."/LQ1_SeaLion_clade-quartet_single-tree_support_".$ft.".tex" ;
			print	"\n\t\t...", $latex_filename ;
			
			my	$latex_caption		=	"LQ1 -- Support values for single clade-quartet trees in specified quartet analyses." ;
			my	@latex_fileheader	=	(
			
				"Clade-Quartet\t", "Score\t", "Filter\t", "Tree\t", "Tree_Code\t", "Average", "Support"
			) ;
			
			my	@lol_latex_linedata ;
			my	$ycounter	= 0;
			for my $latexline ( @{$href_hoh_linedata_of_latexfilename_of_printcode->{latexLQ1}{$ft}} ){
				
				my	@latex_data = split "\t", $latexline ;
				
				my	$xcounter	= 0 ;
				for	my $single_d ( @latex_data ){ $lol_latex_linedata[$ycounter][$xcounter] = $single_d ; $xcounter++ }
				$ycounter++
			}
			
			&latextable_landscape(
				
				\@lol_latex_linedata	,	# two-dimensional matrix: $mulitarray [lines][columns]
				\@latex_fileheader		,	# latex header line (column sorted)
				\$latex_filename		,	# latex output name
				\$latex_caption			,	# latex table caption
			) ;
			
			my	$ftype	=	$main."_".$ft ;
			push @{$href_hol_pdf_of_type->{$ftype}}, $href_subfolder_q->{tex}."/LQ1_SeaLion_clade-quartet_single-tree_support_".$ft.".pdf" ;
		}
	}
	#####################################################################################
	#####################################################################################
	
	##################################################################################### LQ2 works ## main
	## latex landscape table of best (filtered|unfiltered) 4clan tree support          ##
	#####################################################################################
	if	( 	( $href_print_latex->{table_besttreesupport_4clan}					== 1 ) &&
			( $href_value_of_scr_par->{pdf_prtResult}{main} 					== 1 ) ){
		
		for my $ft ( sort keys %{$href_hoh_linedata_of_latexfilename_of_printcode->{latexLQ2}} ){
			
			my	$latex_filename		=	$href_subfolder_q->{tex}."/LQ2_SeaLion_clade-quartet_best-tree_support_".$ft.".tex" ;
			print	"\n\t\t...", $latex_filename ;
			
			my	$latex_caption		=	"LQ2 -- Best clade-quartet tree support for specified quartet analyses, including 'Nq TOT' (total count of clade-quartet assigned species-quartets) and 'Nq IN' (number of remaining analysed species-quartets after filtering)." ;
			my	@latex_fileheader	=	(
			
				"Clade-Quartet\t", "Score\t", "Filter\t", "Best Tree\t", "Tree Code", "Average\t", "Support\t", "Nq TOT", "Nq IN"
			) ;
		
			my	@lol_latex_linedata ;
			my	$ycounter	= 0;
			for my $latexline ( @{$href_hoh_linedata_of_latexfilename_of_printcode->{latexLQ2}{$ft}} ){
				
				my	@latex_data = split "\t", $latexline ;
				
				my	$xcounter	= 0 ;
				for	my $single_d ( @latex_data ){ $lol_latex_linedata[$ycounter][$xcounter] = $single_d ; $xcounter++ }
				$ycounter++
			}
			
			&latextable_landscape(
				
				\@lol_latex_linedata	,	# two-dimensional matrix: $mulitarray [lines][columns]
				\@latex_fileheader		,	# latex header line (column sorted)
				\$latex_filename		,	# latex output name
				\$latex_caption			,	# latex table caption
			) ;
			
			my	$ftype	=	$main."_".$ft ;
			push @{$href_hol_pdf_of_type->{$ftype}}, $href_subfolder_q->{tex}."/LQ2_SeaLion_clade-quartet_best-tree_support_".$ft.".pdf" ;
		}
	}
	#####################################################################################
	#####################################################################################
	
	##################################################################################### LQ5 works ## main
	## latex landscape table of taxa without 4seq-quartet participations (unfiltered)  ##
	#####################################################################################
	if	( 	( $href_print_latex->{table_4seq_unassociated_taxa_unfiltered}		== 1 ) &&
			( $href_value_of_scr_par->{pdf_prtResult}{main} 					== 1 ) ){
		
		for my $ft ( sort keys %{$href_hoh_linedata_of_latexfilename_of_printcode->{latexLQ5}} ){
			
			( my	$fm_latex		=	$ft )	=~ s/_/\\_/g ;
			my	$latex_filename		=	$href_subfolder_q->{tex}."/LQ5_SeaLion_species_unassignable_".$ft.".tex" ;
			print	"\n\t\t...", $latex_filename ;
			
			my	$latex_caption		=	"LQ5 -- Summary of species defined in the clade-file that cannot be identified in the multiple-sequence alignment." ;
			my	@latex_fileheader	=	(
				
				"Species"			, 
				"Clade"				, 
				"Score"						
			) ;
			
			my	@lol_latex_linedata ;
			my	$ycounter	= 0;
			for my $latexline ( @{$href_hoh_linedata_of_latexfilename_of_printcode->{latexLQ5}{$ft}} ){
				
				my	@latex_data = split "\t", $latexline ;
				
				my	$xcounter	= 0 ;
				for	my $single_d ( @latex_data ){ $lol_latex_linedata[$ycounter][$xcounter] = $single_d ; $xcounter++ }
				$ycounter++
			}
			
			&latextable_default(
				
				\@lol_latex_linedata	,	# two-dimensional matrix: $mulitarray [lines][columns]
				\@latex_fileheader		,	# latex header line (column sorted)
				\$latex_filename		,	# latex output name
				\$latex_caption			,	# latex table caption
			) ;
			
			my	$ftype	=	$main."_".$ft ;
			push @{$href_hol_pdf_of_type->{$ftype}}, $href_subfolder_q->{tex}."/LQ5_SeaLion_species_unassignable_".$ft.".pdf" ;
		}
	}
	#####################################################################################
	#####################################################################################
	
	##################################################################################### LQ4 works ## main
	## latex landscape table of taxon related 4seq-quartet participations (unfiltered) ##
	#####################################################################################
	if	( 	( $href_print_latex->{table_4seq_associated_taxa_unfiltered}		== 1 ) &&
			( $href_value_of_scr_par->{pdf_prtResult}{main} 					== 1 ) ){
		
		for my $ft ( sort keys %{$href_hoh_linedata_of_latexfilename_of_printcode->{latexLQ4}} ){
			
			( my	$fm_latex		=	$ft )	=~ s/_/\\_/g ;
			my	$latex_filename		=	$href_subfolder_q->{tex}."/LQ4_SeaLion_clade-quartet_species_participations_unfiltered_".$ft.".tex" ;
			print	"\n\t\t...", $latex_filename ;
			
			my	$latex_caption		=	"LQ4 -- Overview of species-related quartet participations (unfiltered) in single clade-quartets, including 'N' (Number of species-related clade-quartet participations) and 'N TOT' (Total number of species-related quartet participations across all clade-quartets)." ;
			my	@latex_fileheader	=	(
				
				"Species"			, 
				"Clade"				, 
				"Score"				, 
				"Clade-Quartet"		, 
				"N"					, 
				"N TOT"			
			) ;
			
			my	@lol_latex_linedata ;
			my	$ycounter	= 0;
			for my $latexline ( @{$href_hoh_linedata_of_latexfilename_of_printcode->{latexLQ4}{$ft}} ){
				
				my	@latex_data = split "\t", $latexline ;
				
				my	$xcounter	= 0 ;
				for	my $single_d ( @latex_data ){ $lol_latex_linedata[$ycounter][$xcounter] = $single_d ; $xcounter++ }
				$ycounter++
			}
			
			&latextable_default(
				
				\@lol_latex_linedata	,	# two-dimensional matrix: $mulitarray [lines][columns]
				\@latex_fileheader		,	# latex header line (column sorted)
				\$latex_filename		,	# latex output name
				\$latex_caption			,	# latex table caption
			) ;
			
			my	$ftype	=	$main."_".$ft ;
			push @{$href_hol_pdf_of_type->{$ftype}}, $href_subfolder_q->{tex}."/LQ4_SeaLion_clade-quartet_species_participations_unfiltered_".$ft.".pdf" ;
		}
	}
	#####################################################################################
	#####################################################################################
	
	##################################################################################### LQ3 works ## main
	## latex landscape table of taxon related 4seq-quartet participations (filtered)   ##
	#####################################################################################
	if	( 	( $href_print_latex->{table_4seq_associated_taxa_filtered}			== 1 ) &&
			( $href_value_of_scr_par->{pdf_prtResult}{main} 					== 1 ) ){
		
		for my $ft ( sort keys %{$href_hoh_linedata_of_latexfilename_of_printcode->{latexLQ3}} ){
			
			( my	$fm_latex		=	$ft )	=~ s/_/\\_/g ;
			my	@prts				=	split "_", $ft ;
			my	$cm					=	$prts[0] ;
			my	$filtertype			=	&recode_filter_printout( \$ft )	;
			my	$latex_filename		=	$href_subfolder_q->{tex}."/LQ3_SeaLion_clade-quartet_species_participations_".$filtertype."_filtered_".$cm.".tex" ;
			print	"\n\t\t...", $latex_filename ;
			
			
			( my $ftype_latex = $filtertype )	=~	s/_/\\_/g ;
			my	$latex_caption		=	"LQ3 -- Overview of single species participations in quartets before and after '".$ftype_latex."'-filtering. 
			For each species participating clade-quartet, the number of rejected species-quartets (Nrej), the number of remaining species-quartets (Nrem), and the total number of species related quartet participations (N) are listed together with the total number (overall clade-quartets) of species-quartet participations (Nrej TOT, Nrem TOT, N TOT)" ;
			my	@latex_fileheader	=	(
				
				"Species"		, 
				"Clade"			, 
				"Score"			,
				"Filter"		, 
				"Clade-Quartet"	, 
				"Nrej"			, 
				"Nrem"			,
				"N"				, 
				"Nrej TOT"		, 
				"Nrem TOT"		,
				"N TOT"			
			) ;
			
			my	@lol_latex_linedata ;
			my	$ycounter	= 0;
			for my $latexline ( @{$href_hoh_linedata_of_latexfilename_of_printcode->{latexLQ3}{$ft}} ){
				
				my	@latex_data = split "\t", $latexline ;
				
				my	$xcounter	= 0 ;
				for	my $single_d ( @latex_data ){ $lol_latex_linedata[$ycounter][$xcounter] = $single_d ; $xcounter++ }
				$ycounter++
			}
			
			&latextable_landscape(
				
				\@lol_latex_linedata	,	# two-dimensional matrix: $mulitarray [lines][columns]
				\@latex_fileheader		,	# latex header line (column sorted)
				\$latex_filename		,	# latex output name
				\$latex_caption			,	# latex table caption
			) ;
			
			my	$ftype	=	$main."_".$ft ;
			push @{$href_hol_pdf_of_type->{$ftype}}, $href_subfolder_q->{tex}."/LQ3_SeaLion_clade-quartet_species_participations_".$filtertype."_filtered_".$cm.".pdf" ;
		}
	}
	#####################################################################################
	#####################################################################################
	
	##################################################################################### LQ6 works ## main
	## latex landscape table of 4clan conflicted OTU's                                 ##
	#####################################################################################
	if	(	( $href_print_latex->{table_optimized_4seq_filter_threshold} 		== 1 ) &&
			( $href_value_of_scr_par->{pdf_prtResult}{main} 					== 1 ) ){
		
		for my $filter ( sort keys %{$href_hoh_linedata_of_latexfilename_of_printcode->{latexLQ6}} ){
			
			my	$ftype				=	$main."_".$filter ;
			( 	my	$filter_latex	=	$filter ) =~	s/_/\\_/g ;
			my	@prts				=	split "_", $filter ;
			
			my	$filtertype			 =	&recode_filter_printout( \$filter )	;
			my	$latex_filename		=	$href_subfolder_q->{tex}."/LQ6_SeaLion_clade-quartet_optimized_".$filtertype."_filter_threshold_".$prts[0].".tex" ;
			print	"\n\t\t...", $latex_filename ;
			
			my	$latex_caption		=	"LQ6 -- Number of remaining species-quartets (Nq IN) for clade-quartet individually optimized ".$filtertype."-threshold (Tl), along with the total count of remaining species-quartets (Nq TOT) across all clade-quartets."   ;
			my	@latex_fileheader	=	(
				
				"Score"				,
				"Filter"			,
				"Clade-Quartet"		,
				"Tl"				,
				"Nq IN"				,
				"Nq TOT"
			) ;
			
			my	@data_lines	= @{$hoh_linedata_of_latexfilename_of_printcode{latexLQ6}{$filter}} ;
			
			my	@lol_latex_linedata ;
			my	$ycounter	= 0;
			for my $latexline ( @data_lines ){
				
				my	@latex_data = split "\t", $latexline ;
				
				my	$xcounter	= 0 ;
				for	my $single_d ( @latex_data ){ $lol_latex_linedata[$ycounter][$xcounter] = $single_d ; $xcounter++ }
				$ycounter++
			}
			
			&latextable_default(
				
				\@lol_latex_linedata	,	# two-dimensional matrix: $mulitarray [lines][columns]
				\@latex_fileheader		,	# latex header line (column sorted)
				\$latex_filename		,	# latex output name
				\$latex_caption			,	# latex table caption
			) ;
			
			push @{$href_hol_pdf_of_type->{$ftype}}, $href_subfolder_q->{tex}."/LQ6_SeaLion_clade-quartet_optimized_".$filtertype."_filter_threshold_".$prts[0].".pdf" ;
		}
	}
	#####################################################################################
	#####################################################################################
	
	##################################################################################### LQT1 works ## main
	## latex landscape table of nwk tree-legend QT1, QT2, QT3                          ##
	#####################################################################################
	if	( 	( $href_print_latex->{table_4clan_quartet_trees}					== 1 ) &&
			( $href_value_of_scr_par->{pdf_prtResult}{main} 					== 1 ) ){
		
		for my $ft ( sort keys %{$hoh_linedata_of_latexfilename_of_printcode{latexLQT1}} ){
			
			( my	$fm_latex		=	$ft )	=~ s/_/\_/g ;
			my	@prts				=	split "_", $ft ;
			my	$qfilter_prtype		=	&recode_filter_printout( \$ft ) ;
			
			my	$latex_filename		=	$href_subfolder_q->{tex}."/LQT1_SeaLion_clade-quartet_tree_codes_".$qfilter_prtype."_".$prts[0].".tex" ;
			print	"\n\t\t...", $latex_filename ;
			
			my	$latex_caption		=	"LQT1 -- List of coded clade-quartet trees, as utilized in result plots such as triangle plots and tables, including 'QT1', 'QT2', and 'QT3'." ;
			my	@latex_fileheader	=	(
			
				"Score",
				"Filter",
				"Clade-Quartet",
				"QT1",
				"QT2",
				"QT3"
			) ;
			
			my	@data_lines	= @{$hoh_linedata_of_latexfilename_of_printcode{latexLQT1}{$ft}} ;
			
			my	@lol_latex_linedata ;
			my	$ycounter	= 0;
			for my $latexline ( sort @data_lines ){
				
				$latexline =~ s/_/\_/g ;
				
				my	@latex_data = split "\t", $latexline ;
				
				my	$xcounter	= 0 ;
				for	my $single_d ( @latex_data ){ $lol_latex_linedata[$ycounter][$xcounter] = $single_d ; $xcounter++ }
				$ycounter++
			}
			
			&latextable_landscape(
				
				\@lol_latex_linedata	,	# two-dimensional matrix: $mulitarray [lines][columns]
				\@latex_fileheader		,	# latex header line (column sorted)
				\$latex_filename		,	# latex output name
				\$latex_caption			,	# latex table caption
			) ;
			
			my $ftype = $main."_".$ft ;
			push @{$href_hol_pdf_of_type->{$ftype}}, $href_subfolder_q->{tex}."/LQT1_SeaLion_clade-quartet_tree_codes_".$qfilter_prtype."_".$prts[0].".pdf" ;
		}
	}
	#####################################################################################
	#####################################################################################
}

sub norm_values{
	
	my	$aref_values	=	$_[0] ;	# list of single values		IN: defined; OUT: unchanged
	
	######################################################################
	## normalize a list of values to each other and return new list     ##
	######################################################################
	#print "::&norm_values\n";
	
	######################################################################
	# called from subroutine(s)...
	# &evaluation_quartet_support_napX
	# &filter_quartets
	# &taxon_clipping
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	#########################################
	# calculate total of single values
	my	$sum_values	= &calc_sum( \@$aref_values ) ;
	#########################################
	
	#########################################
	# normalize each value and return normalized values as new list
	my	@normed_values ;
	
	for	( @$aref_values ){ 
		
		unless ( $sum_values == 0 ){
			
			push @normed_values, ( sprintf "%.4f", $_ / $sum_values ) 
		}
		
		else{ push @normed_values, 0 }
	}
	#########################################
	
	#########################################
	# return new list of normalized values
	return @normed_values ;
	#########################################
	
	######################################################################
}

sub plotsR_quartet_analyses{
	
	my	$href_print_Rplot_quartets						= $_[0]	; # Print control of single R plots due to single split-pattern quartet analyses																			IN: defined; OUT: unchanged	
	my	$href_hoh_tsvdata_of_tsvfilename_of_printcode	= $_[1]	; # key1: R grafic unique printcode (e.g. 'metafileT'); key2: *.tsv filename; value: *tsv filename data addings (e.g. clanfilename, type of support...)  	IN: defined; OUT: unchanged	
	my	$href_hol_pdf_of_type							= $_[2]	; # key: result type; value: list of associated pdf's in order of construction date																			IN: defined; OUT: changed
	my	$aref_pdf_pattern								= $_[3]	; # list variable of single quartet inference denerated pdf's																								IN: defined; OUT: changed	
	my	$href_subfolder_q								= $_[4]	; # key1: file type ('svg', 'pdf', 'tex', 'tsv') ; value: printdir																							IN: defined; OUT: unchanged	
	my	$href_subsubfolder_q							= $_[5] ; # key1: pdf related subfolder codes ('pdf_single', 'pdf_sumed'...)																						IN: defined; OUT: unchanged	
	my	$href_value_of_scr_par							= $_[6] ; # key1: sealion user available parameter options; value: optional parameter setting																 		IN: defined; OUT: unchanged	
	my	$href_specified_support_method_code				= $_[7] ; # list of defined support methods (nap, napv ...)																											IN: defined; OUT: changed	
	
	######################################################################
	## print R graphics of extracted quartet data                       ##
	######################################################################
	print	"\n\n\t---------------------------------------------------\n",
			"\twriting R-PLOT file(s)...\n" ;
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# &r_lineplot
	# &r_ternplot
	# &r_boxplot
	# &r_pointplot
	# &r_density
	# &r_barplot
	######################################################################
	
	
	######################################################################
	# Definition of different color-lists used for graphical output    ###
	######################################################################
	
	###################################################
	# 12 colors
	my	@colorTable_13c	=	(	"cornflowerblue", "salmon"	, "mediumseagreen"	, "purple", "gold"	, "gray",	"black"		,
								"tomato"	, "orange"	, "moccasin", "linen"	, "pink",	"violet"	) ;
	###################################################
	
	###################################################
	# 12 colors
	my	@colorTable_12c	=	(	"cornflowerblue", "salmon"	, "mediumseagreen"	, "gold"	, "gray",	"black"		,
								"tomato"	, "orange"	, "moccasin", "linen"	, "pink",	"violet"	) ;
	###################################################
	
	###################################################
	# 6 colors
	#my	@colorTable_6c	=	(	"royalblue4", "red2", "green4", "gold", "gray", "black"	) ;
	my	@colorTable_6c	=	(	"mediumpurple", "cornflowerblue", "lightsteelblue", "chocolate", "sandybrown", "slategray"	) ;
	###################################################
	
	###################################################
	# 4 colors
	my	@colorTable_4c	=	(	"chocolate", "cadetblue", "royalblue1", "black"		) ;
	my	@colorTable_4cb	=	(	"darkgray", "darkorange", "royalblue1", "red"		) ;
	my	@colorTable_4cc	=	(	"royalblue4", "red", "black","deepskyblue"			) ;
	my	@colorTable_4cd	=	(	"royalblue4", "red","gold1", "gold4"				) ;
	my	@colorTable_4ce	=	(	"royalblue4", "darkred", "black","deepskyblue"		) ;
	my	@colorTable_4cf	=	(	"royalblue4", "darkred","gold1", "gold4"			) ;
	my	@colorTable_4cg	=	(	"royalblue1", "tomato1", "springgreen4", "black"	) ;
	###################################################
	
	###################################################
	# 3 colors
	#my	@colorTable_3c	=	(	"royalblue4", "red2", "green4"	) ;
	my	@colorTable_3c	=	(	"cornflowerblue", "red2", "green4"			) ;
	my	@colorTable_3cb	=	(	"royalblue4", "red","orange"				) ;
	my	@colorTable_3cc	=	(	"royalblue1", "tomato1", "springgreen4"		) ;
	my	@colorTable_3cd	=	(	"royalblue4", "red", "darkgreen"			) ;
	my	@colorTable_3ce	=	(	"darkorange", "dimgray", "darkgreen"		) ;
	my	@colorTable_3cf	=	(	"royalblue1", "tomato1", "gold1"			) ;
	my	@colorTable_3cg	=	(	"orange", "darkgray", "springgreen4"		) ;
	my	@colorTable_3ch	=	(	"lightsalmon", "lightcoral", "springgreen4"	) ;
	my	@colorTable_3ci	=	(	"lightgray", "darkgray", "springgreen4"		) ;
	my	@colorTable_3cj	=	(	"coral", "burlywood", "rosybrown"			) ;
	my	@colorTable_3ck	=	(	"cadetblue", "chocolate", "mediumseagreen"	) ;
	my	@colorTable_3cl	=	(	"dimgray", "lightcoral", "springgreen4"		) ;
	my	@colorTable_3cm	=	(	"royalblue4", "rosybrown", "red"			) ;
	###################################################
	
	###################################################
	# 2 colors
	my	@colorTable_2c	=	(	"gold", "gray"						) ;
	my	@colorTable_2cb	=	(	"mediumslateblue", "darkmagenta"	) ;
	my	@colorTable_2cc	=	(	"lightsteelblue", "chocolate"		) ;
	my	@colorTable_2cd	=	(	"royalblue1", "tomato1"				) ;
	my	@colorTable_2ce	=	(	"royalblue4", "red"					) ;
	my	@colorTable_2cf	=	(	"black", "darkorange1"				) ;
	my	@colorTable_2cg	=	(	"royalblue4", "chocolate"			) ;
	my	@colorTable_2ch	=	(	"tomato1", "black"					) ;
	###################################################
	
	###################################################
	# 1 color
	my	@colorTable_1c	=	(	"black"	) ;
	###################################################
	
	
	######################################################################
	# if string length of newick 4clan trees is to long for grafical R prints, use tree corresponding table code ('QT_code') instead
	my $topo = 'Topo' ;
	if ( $href_value_of_scr_par->{sys_parameter}{flag_qtrPrt} == 1 ){ $topo = 'QT_code' }
	######################################################################
	
	######################################################################
	# definition of graphic-info assigned output-subcategories for pdf summarization
	my	$main	= 'main'		;	# file_suffix for pdf summarization of main overall 4clan result-files
	my	$dev1	= 'dev1_PD'		;	# results about single pattern distributions for further script developments
	my	$dev2	= 'dev2_NcNa'	;	# results about Nc/Na pattern analyses
	my	$dev3	= 'dev3_NsCf'	;	# results about singleton and cf analyses
	my	$dev4	= 'dev4_VI'		;	# results about complete variable (V) and invariable (I) site-pattern frequencies
	######################################################################
	
	###########################################################################################################################
	############### Plot MAIN results of single clade-quartet and rooted-clade analyses                                 #######
	###########################################################################################################################
	
	##################################################################################### MQ3 
	## R barchart-plot of best rooted-clade trees                                      ## works ## main userversion (adapted)
	##################################################################################### 
	if	( 	( $href_print_Rplot_quartets->{pbarplot_treesuccess_mclan}			== 1 ) &&
			( $href_value_of_scr_par->{pdf_prtResult}{main} 					== 1 ) ){
		
		for my	$tsv_mclan_metafileMQ3	( keys %{$href_hoh_tsvdata_of_tsvfilename_of_printcode->{metafileMQ3}} ){
			
			my	$clanfile	=	$hoh_tsvdata_of_tsvfilename_of_printcode{metafileMQ3}{$tsv_mclan_metafileMQ3}	;
			
			my	$titl	=	"Plot MQ3 -- Best rooted-clade tree support for unfiltered and filtered analyses of ".$clanfile ;
				
			&r_barplot(	
				
				\$tsv_mclan_metafileMQ3						,	# name of infile table (*.tsv)
				\"Filter"									,	# table header of x-plot categories
				\"Support"									,	# table header of y-plot values
				\"Best_Rooted_Clade_Tree"					,	# table header (factor) for value separation given different fill colors
				\"Best_Rooted_Clade_Tree"					,	# table header (factor) for boxplot frames given different or a single line color(s)
				\"facet"									,	# if defined, facet grid is activated in respect of the next two group parameters
				\"Average"									,	# table header plot-splitting vertically (defined as dot "." if splitting not needed)
				\"Score"									,	# table header plot-splitting horizontally (undefined "" if splitting not needed)
				\"fixed"									,	# default (undef) all the panels have the same scales ('fixed') can be made independent if defined 'free', 'free_x', or 'free_y'
				\""											,	# highlight horizontal line at axis (multiple lines comma separated) interception; neglected if undef
				\"dashed"									,	# horizontal linetype: solid, longdash, dashed, blank, dotdash, dotted, twodash, default if undef
				\1											,	# horizontal linesize; default if undef
				\"gray49"									,	# horizontal line color, default if undef
				\0											,	# 0 -> a stacked barplot is created by default, 1 -> non-stagged barplot with each value of the same x-plot category in a own x-axis position
				\0											,	# 0 -> no labels at barchart, 1 -> label (geom_text) at barchart
				\"Qdata"									,	# table header of barchart label
				\"-0.1"										,	# vertical shift value (integer, float, positive or negative) of barchart label
				\"black"									,	# barchart label color
				\"2.5"										,	# barchart label size 
				\$titl										,	# table title
				\@colorTable_13c							,	# list of specified fill colors, or if undefined using a grey-scaling
				\@colorTable_1c								,	# list of specified line color(s)
				\"Analysis"									,	# label x-axis
				\1											,	# if value == 1 -> rotate x-axis labels via 90°, undef (0) keeps labels horizontally
				\"Final Support"							,	# labe y-axis
				\"Best Rooted-Clade Tree"					,	# box colour label
				\40											,	# width of pdf graphic, keep undef for def width: 40
				\15											,	# height of pdf graphic, keep undef for def height: 15
				\"barchart"									,	# pdf file suffix to avoid different plots of identic names
				\$href_subfolder_q->{svg}					,	# svg print folder
				\$href_subsubfolder_q->{pdf_single}			,	# pdf print folder
				\@$aref_pdf_pattern							,	# pdf sampling for pdfunite
			);
			
			# go through each defined quartet calculation method (nap, napv, napr)
			for my $cm ( sort keys %{$href_specified_support_method_code->{original}} ){
				
				my	$ftype	=	$main."_".$cm	;
				push @{$href_hol_pdf_of_type->{$ftype}}, $aref_pdf_pattern->[-1] ;
			}
		}
	}
	#####################################################################################
	#####################################################################################
	
	##################################################################################### R1 
	## R barchart-plot of best rooted clade newick trees                               ## works ## main userversion (adapted)
	##################################################################################### 
	if	( 	( $href_print_Rplot_quartets->{ggtree_rooted_clade_trees}			== 1 ) &&
			( $href_value_of_scr_par->{pdf_prtResult}{main} 					== 1 ) ){
		
		my	@trees = @{$href_hoh_tsvdata_of_tsvfilename_of_printcode->{metafileR1}{besttree}}	;
		
		my	@fm ; my $cm ;
		for my $av_cm_fm ( @{$href_hoh_tsvdata_of_tsvfilename_of_printcode->{metafileR1}{method}} ){
			
			my	@prts  			=	split "_strict_", $av_cm_fm ;
			my	@av_cm 			=	split "_", $av_cm_fm ;
				$cm				=	$av_cm[1]	;
			my	$qfilter_type	=	&recode_filter_printout( \$prts[0] )	;
			my	$tree_name		;	
			if ( $prts[1] ){
				
				my	$tfilter_type	= &recode_filter_printout( \$av_cm_fm )	;
					$tree_name		= $qfilter_type."_".$tfilter_type."_".$av_cm[0]."_".$cm
			}
			else{	$tree_name		= $qfilter_type."_".$av_cm[0] }
			
			push @fm, $tree_name ;
		}
		
		my	$outfile	= "R1_SeaLion_metaData_best_rooted-clade_trees_".$cm	;
		
		#my	$titl	=	"Plot MQ3 -- Best rooted-clade tree support for unfiltered and filtered analyses of ".$clanfile ;
				
		&r_tree(	
			
			\$outfile									,	# outfile name
			\@trees										,	# best rooted-clade trees of each filetered and unfiltered analysis
			\@fm										,	# best tree corresponding filter method (same list order as @trees)
			\15											,	# height of pdf graphic, keep undef for def height: 15
			\40											,	# width of pdf graphic, keep undef for def width: 40
			\12											,	# label letter-size
			\$href_subfolder_q->{svg}					,	# svg print folder
			\$href_subsubfolder_q->{pdf_single}			,	# pdf print folder
			\@$aref_pdf_pattern							,	# pdf sampling for pdfunite
			\@colorTable_13c							,	# list of specified fill colors, or if undefined using a grey-scaling
		);
		
		my	$ftype	=	$main."_".$cm	;
		push @{$href_hol_pdf_of_type->{$ftype}}, $aref_pdf_pattern->[-1] ;
	}
	#exit;
	#####################################################################################
	#####################################################################################
	
	##################################################################################### MQ2
	# R ternary plot of single, mean, and median quartet support of different 4clan    ##
	# relationships (nap, napr, napv) unfiltered                                       ## works ## main userversion (adapted)
	#####################################################################################
	if	( 	( $href_print_Rplot_quartets->{ptern_seqsupport_clan4}				== 1 ) &&
			( $href_value_of_scr_par->{pdf_prtResult}{main} 					== 1 ) ){
		
		for my	$tsv_mclan_metafileMQ2	( keys %{$href_hoh_tsvdata_of_tsvfilename_of_printcode->{metafileMQ2}} ){
			
			my	@filters_used	=	split "::", $href_hoh_tsvdata_of_tsvfilename_of_printcode->{metafileMQ2}{$tsv_mclan_metafileMQ2} ;
			my	$cm 			=	shift @filters_used ; #print "\n\t", $cm, "\n\n" ;
			my	$ftype			=	$main."_".$cm	;
			my	$qfilter_type	=	&recode_filter_printout( \$cm )	;
			my	$pdf_suffix 	= 	$qfilter_type."_ternaryplot" ;
			my	$ttitl			=	"Plot MQ2 -- ".$qfilter_type." species-quartet support for single clade-quartet trees QT1, QT2, and QT3 (".$cm.")" ;
			
			&r_ternplot(	
				
				\$tsv_mclan_metafileMQ2							,	# name of infile table (*.tsv)
				\"QT1"											,	# x-axis
				\"QT2"											,	# y-axis
				\"QT3"											,	# z-axis
				\"Support_Type"									,	# table header (factor) for value separation given different point colors within each plot
				\"facet"										,	# if defined, facet grid is activated in respect of the next two group parameters
				\"Clade_Quartet"								,	# table header plot-splitting vertically (defined as dot "." if splitting not needed)
				\"Score"										,	# table header plot-splitting horizontally (undefined "" if splitting not needed)
				\$ttitl											,	# plot title
				\4												,	# point size within triangle (if undef = 6)
				\20												,	# height of pdf graphic, keep undef for def height: 15
				\40												,	# width of pdf graphic, keep undef for def width: 40
				\"unfiltered_ternaryplot"						,	# pdf file suffix to avoid different plots of identic names
				\$href_subfolder_q->{svg}						,	# svg print folder
				\$href_subsubfolder_q->{pdf_single}				,	# pdf print folder
				\@$aref_pdf_pattern								,	# pdf sampling for pdfunite
			) ;
			
			push @{$href_hol_pdf_of_type->{$ftype}}, $aref_pdf_pattern->[-1] ;
		}
	}
	#####################################################################################
	#####################################################################################
	
	##################################################################################### MQ6
	# R lineplot of quartet-filter optimization of each method                         ## works ## main userversion (adapted)
	#####################################################################################
	if	( 	( $href_print_Rplot_quartets->{pline_qfilter_optimzation}			== 1 ) &&
			( $href_value_of_scr_par->{pdf_prtResult}{main} 					== 1 ) ){
		
		for my	$tsv_mclan_metafileMQ6	( keys %{$href_hoh_tsvdata_of_tsvfilename_of_printcode->{metafileMQ6}} ){
			
			my	@prts				=	split "::", $href_hoh_tsvdata_of_tsvfilename_of_printcode->{metafileMQ6}{$tsv_mclan_metafileMQ6} ;
			my	$type_support		=	$prts[0] ;
			my	@prts_support		=	split "_", $type_support ; 
			my	$vline				=	$prts[1] ;
			my	$filter				=	&recode_filter_printout( \$type_support )	;
			my	$titl				=	"Plot MQ6 -- Change of f(P) values during the optimization of '".$filter."'-thresholds for species-quartet-filtering related to each clade-quartet (".$prts_support[0].")" ;
				
			&r_lineplot(	
					
					\$tsv_mclan_metafileMQ6						,	# name of infile table (*.tsv)
					\"Scale"									,	# table header of x-plot categories
					\"fPScale"									,	# table header of y-plot values
					\"Rank"										,	# table header (factor) for value separation given different point colors within each plot
					\""											,	# table header (factor) for value separation given different line colors within each plot
					\""											,	# table header (factor) for showing different points in different shapes (neglegted if undef)
					\1											,	# do lines between plot points -> 1; otherwise, no line connections between points
					\4											,	# size lines, default size if undef
					\12											,	# size points, default size if undef
					\"facet"									,	# if defined, facet grid is activated in respect of the next two group parameters
					\"Clade_Quartet"							,	# table header plot-splitting vertically (defined as dot "." if splitting not needed)
					\""											,	# table header plot-splitting horizontally (undefined "" if splitting not needed)
					\"fixed"									,	# default (undef) all the panels have the same scales ('fixed') can be made independent if defined 'free', 'free_x', or 'free_y'
					\"1"										,	# highlight horizontal line at axis interception; neglected if undef
					\"dashed"									,	# horizontal linetype: solid, longdash, dashed, blank, dotdash, dotted, twodash, default if undef
					\1											,	# horizontal linesize; default if undef
					\"gray49"									,	# horizontal line color, default if undef
					\""											,	# highlight vertical line at axis interception; neglected if undef
					\"dashed"									,	# vertical linetype: solid, longdash, dashed, blank, dotdash, dotted, twodash, default if undef
					\1											,	# vertical linesize; default if undef
					\"gray49"									,	# vertical line color, default if undef
					\$titl										,	# table title
					\@colorTable_2ch							,	# list of specified line colors
					\"Filter-Threshold"							,	# label x-axis
					\0											,	# if value == 1 -> rotate x-axis labels via 90°, 0 -> keeps labels horizontally
					\"f(P)"										,	# labe y-axis
					\""											,	# point colour label
					\""											,	# line (outline) colour label
					\""											,	# size shape label 
					\15											,	# height of pdf graphic, keep undef for def height: 15
					\40											,	# width of pdf graphic, keep undef for def width: 40
					\"change_of_optimum_lineplot"				,	# pdf file suffix to avoid different plots of identic names
					\$href_subfolder_q->{svg}					,	# svg print folder
					\$href_subsubfolder_q->{pdf_single}			,	# pdf print folder
					\@$aref_pdf_pattern							,	# pdf sampling for pdfunite
			) ;
			
			my	$ftype	=	$main."_".$type_support ;
			push @{$href_hol_pdf_of_type->{$ftype}}, $aref_pdf_pattern->[-1] ;
		}
	}
	#####################################################################################
	#####################################################################################
	
	##################################################################################### MQ6
	# R lineplot of quartet-filter optimization of each method                         ## works ## main userversion (adapted)
	#####################################################################################
	if	( 	( $href_print_Rplot_quartets->{pline_qnumber_optimzation}			== 1 ) &&
			( $href_value_of_scr_par->{pdf_prtResult}{main} 					== 1 ) ){
		
		for my	$tsv_mclan_metafileMQ6	( keys %{$href_hoh_tsvdata_of_tsvfilename_of_printcode->{metafileMQ6}} ){
			
			my	@prts				=	split "::", $href_hoh_tsvdata_of_tsvfilename_of_printcode->{metafileMQ6}{$tsv_mclan_metafileMQ6} ;
			my	$type_support		=	$prts[0] ;
			my	@prts_support		=	split "_", $type_support ; 
			my	$vline				=	$prts[1] ;
			my	$filter				=	&recode_filter_printout( \$type_support )	;
			my	$titl				=	"Plot MQ6 -- Count of retained species-quartets throughout the '".$filter."'-threshold optimization process for each clade-quartet (".$prts_support[0].")" ;
				
			&r_lineplot(	
					
					\$tsv_mclan_metafileMQ6						,	# name of infile table (*.tsv)
					\"Scale"									,	# table header of x-plot categories
					\"NqRemain"									,	# table header of y-plot values
					\"Rank"										,	# table header (factor) for value separation given different point colors within each plot
					\""											,	# table header (factor) for value separation given different line colors within each plot
					\""											,	# table header (factor) for showing different points in different shapes (neglegted if undef)
					\1											,	# do lines between plot points -> 1; otherwise, no line connections between points
					\4											,	# size lines, default size if undef
					\12											,	# size points, default size if undef
					\"facet"									,	# if defined, facet grid is activated in respect of the next two group parameters
					\"Clade_Quartet"							,	# table header plot-splitting vertically (defined as dot "." if splitting not needed)
					\""											,	# table header plot-splitting horizontally (undefined "" if splitting not needed)
					\"fixed"									,	# default (undef) all the panels have the same scales ('fixed') can be made independent if defined 'free', 'free_x', or 'free_y'
					\""											,	# highlight horizontal line at axis interception; neglected if undef
					\"dashed"									,	# horizontal linetype: solid, longdash, dashed, blank, dotdash, dotted, twodash, default if undef
					\2											,	# horizontal linesize; default if undef
					\"red"										,	# horizontal line color, default if undef
					\""											,	# highlight vertical line at axis interception; neglected if undef
					\"dashed"									,	# vertical linetype: solid, longdash, dashed, blank, dotdash, dotted, twodash, default if undef
					\1											,	# vertical linesize; default if undef
					\"gray49"									,	# vertical line color, default if undef
					\$titl										,	# table title
					\@colorTable_2ch							,	# list of specified line colors
					\"Filter-threshold"							,	# label x-axis
					\0											,	# if value == 1 -> rotate x-axis labels via 90°, 0 -> keeps labels horizontally
					\"N remaining quartets"						,	# labe y-axis
					\""											,	# point colour label
					\""											,	# line (outline) colour label
					\""											,	# size shape label 
					\15											,	# height of pdf graphic, keep undef for def height: 15
					\40											,	# width of pdf graphic, keep undef for def width: 40
					\"remaining_quartets_lineplot"				,	# pdf file suffix to avoid different plots of identic names
					\$href_subfolder_q->{svg}					,	# svg print folder
					\$href_subsubfolder_q->{pdf_single}			,	# pdf print folder
					\@$aref_pdf_pattern							,	# pdf sampling for pdfunite
			) ;
			
			my	$ftype	=	$main."_".$type_support ;
			push @{$href_hol_pdf_of_type->{$ftype}}, $aref_pdf_pattern->[-1] ;
		}
	}
	#####################################################################################
	#####################################################################################
	
	##################################################################################### MQ1
	# R lineplot of Mean & Median 4clan support with a single plot for each method     ## works ## main userversion (adapted)
	#####################################################################################
	if	( 	( $href_print_Rplot_quartets->{pline_clan4_MeanMedian}				== 1 ) &&
			( $href_value_of_scr_par->{pdf_prtResult}{main} 					== 1 ) ){
		
		for my	$tsv_mclan_metafileMQ1	( keys %{$href_hoh_tsvdata_of_tsvfilename_of_printcode->{metafileMQ1}} ){
			
			my	$supportline 	=	$href_hoh_tsvdata_of_tsvfilename_of_printcode->{metafileMQ1}{$tsv_mclan_metafileMQ1} ;
			my	@supporttypes	=	split "::", $supportline ;
			
			for my $type_support ( @supporttypes ){
				
				my	$header_str 	;
				my	$pdf_suffix		;
				
				my	@prts			=	split "_strict_", $type_support ;
				my	@cm_fm			=	split "_", $prts[0] ;
				my	$qfilter_type	=	&recode_filter_printout( \$prts[0] ) ;
					
				if		( $type_support =~ /_strict_/  		){ 
					
					my	$tfilter_type	=	&recode_filter_printout( \$type_support ) ;
						$header_str 	=	"Plot MQ1 -- '".$tfilter_type."'-filtered support values for individual clade-quartet trees using ".$qfilter_type."-filtered species-quartets (".$cm_fm[0].")" ;
						$pdf_suffix		=	"clade-quartet_single-tree_support_".$qfilter_type."_".$tfilter_type."_".$cm_fm[0]."_lineplot" ;
				}
				elsif	( $qfilter_type =~ /unfiltered/  	){ 
						
						$header_str 	=	"Plot MQ1 -- Unfiltered support values for individual clade-quartet trees (".$cm_fm[0].")" ;
						$pdf_suffix		=	"clade-quartet_single-tree_support_".$qfilter_type."_".$cm_fm[0]."_lineplot" ;
				}
				else	{
						$header_str 	=	"Plot MQ1 -- '".$qfilter_type."'-filtered support values for individual clade-quartet trees (".$cm_fm[0].")" ;
						$pdf_suffix		=	"clade-quartet_single-tree_support_".$qfilter_type."_".$cm_fm[0]."_lineplot" ;
				}
				
				&r_lineplot(	
					
					\$tsv_mclan_metafileMQ1						,	# name of infile table (*.tsv)
					\"Tree"										,	# table header of x-plot categories
					\$qfilter_type								,	# table header of y-plot values
					\"Average"									,	# table header (factor) for value separation given different point colors within each plot
					\""											,	# table header (factor) for value separation given different line colors within each plot
					\""											,	# table header (factor) for showing different points in different shapes (neglegted if undef)
					\""											,	# do lines between plot points -> 1; otherwise, no line connections between points
					\4											,	# size lines, default size if undef
					\12											,	# size points, default size if undef
					\"facet"									,	# if defined, facet grid is activated in respect of the next two group parameters
					\"Clade_Quartet"							,	# table header plot-splitting vertically (defined as dot "." if splitting not needed)
					\""											,	# table header plot-splitting horizontally (undefined "" if splitting not needed)
					\"free"										,	# default (undef) all the panels have the same scales ('fixed') can be made independent if defined 'free', 'free_x', or 'free_y'
					\"0.6,0.4"									,	# highlight horizontal line at axis interception; neglected if undef
					\"dashed"									,	# horizontal linetype: solid, longdash, dashed, blank, dotdash, dotted, twodash, default if undef
					\2											,	# horizontal linesize; default if undef
					\"gray49"									,	# horizontal line color, default if undef
					\""											,	# highlight vertical line at axis interception; neglected if undef
					\"dashed"									,	# vertical linetype: solid, longdash, dashed, blank, dotdash, dotted, twodash, default if undef
					\2											,	# vertical linesize; default if undef
					\"gray49"									,	# vertical line color, default if undef
					\$header_str								,	# table title
					\@colorTable_2cd							,	# list of specified line colors
					\"Clade-Quartet Tree"						,	# label x-axis
					\0											,	# if value == 1 -> rotate x-axis labels via 90°, 0 -> keeps labels horizontally
					\"Support"									,	# labe y-axis
					\"Type of Support"							,	# point colour label
					\""											,	# line (outline) colour label
					\""											,	# size shape label 
					\15											,	# height of pdf graphic, keep undef for def height: 15
					\40											,	# width of pdf graphic, keep undef for def width: 40
					\$pdf_suffix								,	# pdf file suffix to avoid different plots of identic names
					\$href_subfolder_q->{svg}					,	# svg print folder
					\$href_subsubfolder_q->{pdf_single}			,	# pdf print folder
					\@$aref_pdf_pattern							,	# pdf sampling for pdfunite
				) ;
				
				if ( $type_support =~ /strict/ ){ $type_support =~ s/_strict.*$// }
				my	$ftype	=	$main."_".$type_support ;
				push @{$href_hol_pdf_of_type->{$ftype}}, $aref_pdf_pattern->[-1] ;
			}
		}
	}
	#####################################################################################
	#####################################################################################
	
	##################################################################################### MQ4
	# R ternary plot of single, mean, and median quartet support of different 4clan    ##
	# relationships (nap, napr, napv), highlighting filtered + unfiltered quartets     ## works ## main userversion (adapted)
	#####################################################################################
	if	( 	( $href_print_Rplot_quartets->{ptern_seqfilter_clan4_fil_unfil}		== 1 ) &&
			( $href_value_of_scr_par->{pdf_prtResult}{main} 					== 1 ) ){
		
		for my	$tsv_mclan_metafileMQ4	( keys %{$href_hoh_tsvdata_of_tsvfilename_of_printcode->{metafileMQ4}} ){
			
			my	$cm_fm			=	$href_hoh_tsvdata_of_tsvfilename_of_printcode->{metafileMQ4}{$tsv_mclan_metafileMQ4} ;
			my	@prts_support	=	split "_", $cm_fm ; 
			my	$ftype			=	$main."_".$cm_fm ;
			my	$filter			=	&recode_filter_printout( \$cm_fm )	;
			my	$ttitl			=	"Plot MQ4 -- Unfiltered (bottom) and '".$filter."'-filtered (top) species-quartet support for single clade-quartet trees QT1, QT2, and QT3 (".$prts_support[0].")" ;
			my	$pdf_suffix 	= 	"ternaryplot" ;
			
			&r_ternplot(	
				
				\$tsv_mclan_metafileMQ4							,	# name of infile table (*.tsv)
				\"QT1"											,	# x-axis
				\"QT2"											,	# y-axis
				\"QT3"											,	# z-axis
				\"Quartet_Status"								,	# table header (factor) for value separation given different point colors within each plot
				\"facet"										,	# if defined, facet grid is activated in respect of the next two group parameters
				\"Clade_Quartet"								,	# table header plot-splitting vertically (defined as dot "." if splitting not needed)
				\"Filter_Status"								,	# table header plot-splitting horizontally (undefined "" if splitting not needed)
				\$ttitl											,	# plot title
				\4												,	# point size within triangle (if undef = 6)
				\20												,	# height of pdf graphic, keep undef for def height: 15
				\40												,	# width of pdf graphic, keep undef for def width: 40
				\$pdf_suffix									,	# pdf file suffix to avoid different plots of identic names
				\$href_subfolder_q->{svg}						,	# svg print folder
				\$href_subsubfolder_q->{pdf_single}				,	# pdf print folder
				\@$aref_pdf_pattern								,	# pdf sampling for pdfunite
			) ;
			
			push @{$href_hol_pdf_of_type->{$ftype}}, $aref_pdf_pattern->[-1] ;
		}
	}
	#####################################################################################
	#####################################################################################
	
	##################################################################################### MQ2
	# R ternary plot of single, mean, and median quartet support of different 4clan    ##
	# relationships (nap, napr, napv), highlighting filtered quartets                  ## works ## main userversion (adapted)
	#####################################################################################
	if	( 	( $href_print_Rplot_quartets->{ptern_seqfilter_clan4}				== 1 ) &&
			( $href_value_of_scr_par->{pdf_prtResult}{main} 					== 1 ) ){
		
		for my	$tsv_mclan_metafileMQ2	( keys %{$href_hoh_tsvdata_of_tsvfilename_of_printcode->{metafileMQ2}} ){
			
			my	@filters_used	=	split "::", $href_hoh_tsvdata_of_tsvfilename_of_printcode->{metafileMQ2}{$tsv_mclan_metafileMQ2} ;
			my	$cm 			=	shift @filters_used ;
			
			for my $f ( @filters_used ){# print "\n\tf\t", $f;
				
				my	$filter			=	&recode_filter_printout( \$f )	;
				my	$ftype			=	$main."_".$cm."_".$f ;
				my	$ttitl			=	"Plot MQ2 -- '".$filter."'-filtered species-quartet support of single clade-quartet trees QT1, QT2, QT3 (".$cm.")" ;
				my	$pdf_suffix 	= 	$filter."_ternaryplot" ;
				
				&r_ternplot(	
					
					\$tsv_mclan_metafileMQ2							,	# name of infile table (*.tsv)
					\"QT1"											,	# x-axis
					\"QT2"											,	# y-axis
					\"QT3"											,	# z-axis
					\$f												,	# table header (factor) for value separation given different point colors within each plot
					\"facet"										,	# if defined, facet grid is activated in respect of the next two group parameters
					\"Clade_Quartet"								,	# table header plot-splitting vertically (defined as dot "." if splitting not needed)
					\"Score"										,	# table header plot-splitting horizontally (undefined "" if splitting not needed)
					\$ttitl											,	# plot title
					\4												,	# point size within triangle (if undef = 6)
					\20												,	# height of pdf graphic, keep undef for def height: 15
					\40												,	# width of pdf graphic, keep undef for def width: 40
					\$pdf_suffix									,	# pdf file suffix to avoid different plots of identic names
					\$href_subfolder_q->{svg}						,	# svg print folder
					\$href_subsubfolder_q->{pdf_single}				,	# pdf print folder
					\@$aref_pdf_pattern								,	# pdf sampling for pdfunite
				) ;
				
				push @{$href_hol_pdf_of_type->{$ftype}}, $aref_pdf_pattern->[-1] ;
			}
			#exit;
		}
	}
	#####################################################################################
	#####################################################################################
	
	##################################################################################### Q3
	## R pointplot ratio CF exp/obs versus ratio Nc/Na orig. (unfiltered)              ## works ## main userversion (adapted)
	#####################################################################################
	if	( 	( $href_print_Rplot_quartets->{ppoint_rNcNaOrg_CFobsCFexp_unfilt} 	== 1 ) &&
			( $href_value_of_scr_par->{pdf_prtResult}{main} 					== 1 ) ){
		
		for my	$tsv_pwout	( keys %{$href_hoh_tsvdata_of_tsvfilename_of_printcode->{metafileQ3}} ){
			
			my	$clan4 		=	$href_hoh_tsvdata_of_tsvfilename_of_printcode->{metafileQ3}{$tsv_pwout} ;
			my	$tl			=	$href_value_of_scr_par->{seq_filterRisk}{tlrisk} ;
			my	$tu			=	$href_value_of_scr_par->{seq_filterRisk}{turisk} ;
			my	$hlines		=	$tl.",".$tu ;
			
			# go through each defined quartet calculation method (nap, napv, napr)
			for my $cm ( sort keys %{$href_specified_support_method_code->{original}} ){
				
				# napr method does not reduce the original split-pattern
				unless ( $cm eq 'napr' ){ 
					
					my	$qfilter_type			=	&recode_filter_printout( \$cm )	;	
					my	$ptitl_rNcNa_CFexpCFobs	=	"Plot Q3 -- Nc/Na ratio for species-quartet best and alternative-supported (minor) trees in clade-quartet ".$clan4." (".$cm.")"  ;
					my	$x_par					=	"rCF_ExpObs_".$cm."_v";
					my	$highlight				=	"tree_rank_".$cm."_c" ;
					my	$pdf_suffix				=	"Nc-to-Na-ratio_clouds_".$qfilter_type."_".$cm."_pointplot" ;
					my	$hlines					=	$href_value_of_scr_par->{seq_filterRisk}{tlrisk}.",".$href_value_of_scr_par->{seq_filterRisk}{turisk} ;
					my	$ftype					=	$main."_".$cm 			;
					
					&r_pointplot(	
						
						\$tsv_pwout										,	# name of infile table (*.tsv)
						\$x_par											,	# table header of x-plot categories
						\"rNcNa_original"								,	# table header of y-plot values
						\$highlight										,	# table header (factor) for value separation given different fill colors within each plot (neglegted if undef)
						\""												,	# table header (factor) for value separation given different point shapes within each plot (neglegted if undef)
						\""												,	# table header (factor) for showing different elements in different sizes (neglegted if undef)
						\""												,	# if defined, smooth regular expression line is activated (method "lm"), if undef expression line is deactivated
						\""												,	# table header (factor) for smooth regression line colour separation (neglegted if undef, smooth must be activated)
						\""												,	# if 1 -> error intervals of smooth line are colored due to smooth factor, undef -> default is used
						\"facet"										,	# if defined, facet grid is activated in respect of the next two group parameters
						\"Clade_Quartet"								,	# table header facet grid plot-splitting vertically (undefined "" if splitting not needed)
						\"Tree"											,	# table header facet grid plot-splitting horizontally (undefined "" if splitting not needed)
						\"free"											,	# default (undef) all the panels have the same scales ('fixed') can be made independent if defined 'free', 'free_x', or 'free_y'
						\$hlines										,	# highlight horizontal line at axis interception; neglected if undef
						\"dashed"										,	# horizontal linetype: solid, longdash, dashed, blank, dotdash, dotted, twodash, default if undef
						\2												,	# horizontal linesize; default if undef
						\"gray49"										,	# horizontal line color, default if undef
						\""												,	# highlight vertical line at axis interception; neglected if undef
						\"dashed"										,	# vertical linetype: solid, longdash, dashed, blank, dotdash, dotted, twodash, default if undef
						\2												,	# vertical linesize; default if undef
						\"gray49"										,	# vertical line color, default if undef
						\$ptitl_rNcNa_CFexpCFobs						,	# table title
						\@colorTable_3cd								,	# list of specified fill colors
						\@colorTable_3cd								,	# list of specified edge colors
						\"Ratio CF exp/obs"								,	# label x-axis
						\"Ratio Nc/Na"									,	# labe y-axis
						\"Best Tree"									,	# label legend title color fill (area), neglected if undef
						\"Best Tree"									,	# label legend title color scale (outlines), neglected if undef
						\""												,	# label legend title shape
						\""												,	# size shape label, neglected if undef
						\45												,	# height of pdf graphic, keep undef for def height: 15
						\""												,	# width of pdf graphic, keep undef for def width: 40
						\$pdf_suffix									,	# pdf file suffix to avoid different plots of identic names
						\$href_subfolder_q->{svg}						,	# svg print folder
						\$href_subsubfolder_q->{pdf_single}				,	# pdf print folder
						\@$aref_pdf_pattern								,	# pdf sampling for pdfunite
					) ;
			
					push @{$href_hol_pdf_of_type->{$ftype}}, $aref_pdf_pattern->[-1] ;
				}
			}
		}
	}
	#####################################################################################
	#####################################################################################
	
	##################################################################################### Q3
	## R pointplot ratio CF exp/obs versus ratio Nc/Na orig. (quartet-filtered)        ## works ## main userversion (adapted)
	#####################################################################################
	if	( 	( $href_print_Rplot_quartets->{ppoint_rNcNa_CFobsCFexp_filt}		== 1 ) &&
			( $href_value_of_scr_par->{pdf_prtResult}{main} 					== 1 ) ){
		
		my	$tl						=	$href_value_of_scr_par->{seq_filterRisk}{tlrisk} ;
		my	$tu						=	$href_value_of_scr_par->{seq_filterRisk}{turisk} ;
		my	$hlines					=	$tl.",".$tu ;					
		
		for my	$tsv_pwout	( keys %{$href_hoh_tsvdata_of_tsvfilename_of_printcode->{metafileQ3}} ){
			
			my	$clan4 =	$href_hoh_tsvdata_of_tsvfilename_of_printcode->{metafileQ3}{$tsv_pwout} ;
			
			for my $cm ( sort keys %{$href_specified_support_method_code->{expanded}} ){
				
				my	$qfilter_type			=	&recode_filter_printout( \$cm )	;	
				my	@cmfm					=	split "_", $cm 			;
				my	$x_par					=	"rCF_ExpObs_".$cmfm[0]."_v";
				my	$highlight				=	"tree_rank_".$cmfm[0]."_c" ;
				my	$pdf_suffix				=	"Nc-to-Na-ratio_clouds_".$qfilter_type."_".$cmfm[0]."_pointplot" ;
				#my	$hlines					=	$href_value_of_scr_par->{seq_filterRisk}{tlrisk}.",".$href_value_of_scr_par->{seq_filterRisk}{turisk} ;
				my	$ftype					=	$main."_".$cm 			;
				my	$ptitl_rNcNa_CFexpCFobs	=	"Plot Q3 -- Nc/Na ratio of 'RISK'-rejected and remaining species-quartets for best and alternative-supported (minor) trees in clade-quartet ".$clan4." (".$cmfm[0].")" ;
				my	$y_par 					=	"rNcNa_original" ;
					
				# napr method does not reduce the original split-pattern
				if (  $cmfm[1] =~ /^risk/  ){ 
					
					#my	$hlines	=	$href_value_of_scr_par->{seq_filterRisk}{$clan4}{tlrisk}{$cm}.",".$href_value_of_scr_par->{seq_filterRisk}{turisk} ;
					my	$hlines	=	$href_value_of_scr_par->{seq_filterRisk}{$clan4}{tlrisk}{$cm} ;
					
					&r_pointplot(	
						
						\$tsv_pwout										,	# name of infile table (*.tsv)
						\$x_par											,	# table header of x-plot categories
						\$y_par											,	# table header of y-plot values
						\$highlight										,	# table header (factor) for value separation given different fill colors within each plot (neglegted if undef)
						\""												,	# table header (factor) for value separation given different point shapes within each plot (neglegted if undef)
						\""												,	# table header (factor) for showing different elements in different sizes (neglegted if undef)
						\""												,	# if defined, smooth regular expression line is activated (method "lm"), if undef expression line is deactivated
						\""												,	# table header (factor) for smooth regression line colour separation (neglegted if undef, smooth must be activated)
						\""												,	# if 1 -> error intervals of smooth line are colored due to smooth factor, undef -> default is used
						\"facet"										,	# if defined, facet grid is activated in respect of the next two group parameters
						\$cm											,	# table header facet grid plot-splitting vertically (undefined "" if splitting not needed)
						\"Tree"											,	# table header facet grid plot-splitting horizontally (undefined "" if splitting not needed)
						\"free"											,	# default (undef) all the panels have the same scales ('fixed') can be made independent if defined 'free', 'free_x', or 'free_y'
						\$hlines										,	# highlight horizontal line at axis interception; neglected if undef
						\"dashed"										,	# horizontal linetype: solid, longdash, dashed, blank, dotdash, dotted, twodash, default if undef
						\2												,	# horizontal linesize; default if undef
						\"gray49"										,	# horizontal line color, default if undef
						\""												,	# highlight vertical line at axis interception; neglected if undef
						\"dashed"										,	# vertical linetype: solid, longdash, dashed, blank, dotdash, dotted, twodash, default if undef
						\2												,	# vertical linesize; default if undef
						\"gray49"										,	# vertical line color, default if undef
						\$ptitl_rNcNa_CFexpCFobs						,	# table title
						\@colorTable_3cd								,	# list of specified fill colors
						\@colorTable_3cd								,	# list of specified edge colors
						\"Ratio CFexp/CFobs"							,	# label x-axis
						\"Ratio Nc/Na"									,	# labe y-axis
						\"Best Tree"									,	# label legend title color fill (area), neglected if undef
						\"Best Tree"									,	# label legend title color scale (outlines), neglected if undef
						\""												,	# label legend title shape
						\""												,	# size shape label, neglected if undef
						\45												,	# height of pdf graphic, keep undef for def height: 15
						\""												,	# width of pdf graphic, keep undef for def width: 40
						\$pdf_suffix									,	# pdf file suffix to avoid different plots of identic names
						\$href_subfolder_q->{svg}						,	# svg print folder
						\$href_subsubfolder_q->{pdf_single}				,	# pdf print folder
						\@$aref_pdf_pattern								,	# pdf sampling for pdfunite
					) ;
			
					push @{$href_hol_pdf_of_type->{$ftype}}, $aref_pdf_pattern->[-1] ;
				}
			}
		}
	}
	#####################################################################################
	#####################################################################################
	
	##################################################################################### Q3
	## R pointplot of quartet individual ratio Nc/Na orig.                             ##
	## for each 4clan tree RISK1 & RISK2 (Risk3 & Risk4 are based on cf corr. Nc/Na)   ## works ## main DEV version (adapted)
	#####################################################################################
	if	( 	( $href_print_Rplot_quartets->{ppoint_rNcNaOrgRisk_FN}				== 1 ) &&
			( $href_value_of_scr_par->{pdf_prtResult}{main} 					== 1 ) ){
		
		for my	$tsv_mclan_metafileQ3	( keys %{$href_hoh_tsvdata_of_tsvfilename_of_printcode->{metafileQ3}} ){
			
			my	$clan4 	=	$href_hoh_tsvdata_of_tsvfilename_of_printcode->{metafileQ3}{$tsv_mclan_metafileQ3} ;
			
			# go through each defined quartet calculation method (nap, napv, napr)
			for my $cm ( sort keys %{$href_specified_support_method_code->{original}} ){
				
				for my $fm ( "risk1", "risk2" ){
					
					if ( $href_value_of_scr_par->{seq_filterRisk}{$fm} == 1 ){
						
						my	@cmfm					=	split "_", $cm 			;
						my	$filter					=	&recode_filter_printout( \$fm )	;
						my	$ptitl_rNcNaOrgRisk_FN	=	"Plot Q3 -- Nc/Na ratio of '".$filter."'-rejected and remaining species-quartets for best and alternative-supported (minor) trees in clade-quartet ".$clan4." (".$cm.")" ;
						#my	$highlighting 			=	"rNcNa_org_".$fm."_c" ;
						my	$highlighting 			=	$cm."_".$fm ;
						my	$highlighting2 			=	"tree_rank_".$cm."_c" ;
						my	$hlines					=	$href_value_of_scr_par->{seq_filterRisk}{$clan4}{tlrisk}{$highlighting}.",".$href_value_of_scr_par->{seq_filterRisk}{turisk} ;
						my	$pdf_suffix				=	"Nc-to-Na-ratio_per_quartet-number_".$filter."_".$cmfm[0]."_pointplot" ;
					
						&r_pointplot(	
							
							\$tsv_mclan_metafileQ3							,	# name of infile table (*.tsv)
							\"Quartet_Number"								,	# table header of x-plot categories
							\"rNcNa_original"								,	# table header of y-plot values
							\$highlighting2									,	# table header (factor) for value separation given different fill colors within each plot (neglegted if undef)
							\$highlighting									,	# table header (factor) for value separation given different point shapes within each plot (neglegted if undef)
							\$highlighting									,	# table header (factor) for showing different elements in different sizes (neglegted if undef)
							\""												,	# if defined, smooth regular expression line is activated (method "lm"), if undef expression line is deactivated
							\$highlighting									,	# table header (factor) for smooth regression line colour separation (neglegted if undef, smooth must be activated)
							\""												,	# if 1 -> error intervals of smooth line are colored due to smooth factor, undef -> default is used
							\"facet"										,	# if defined, facet grid is activated in respect of the next two group parameters
							\"Clade_Quartet"								,	# table header facet grid plot-splitting vertically (undefined "" if splitting not needed)
							\"Tree"											,	# table header facet grid plot-splitting horizontally (undefined "" if splitting not needed)
							\"fixed"										,	# default (undef) all the panels have the same scales ('fixed') can be made independent if defined 'free', 'free_x', or 'free_y'
							\$hlines										,	# highlight horizontal line at axis interception; neglected if undef
							\"dashed"										,	# horizontal linetype: solid, longdash, dashed, blank, dotdash, dotted, twodash, default if undef
							\2												,	# horizontal linesize; default if undef
							\"gray49"										,	# horizontal line color, default if undef
							\""												,	# highlight vertical line at axis interception; neglected if undef
							\"dashed"										,	# vertical linetype: solid, longdash, dashed, blank, dotdash, dotted, twodash, default if undef
							\2												,	# vertical linesize; default if undef
							\"gray49"										,	# vertical line color, default if undef
							\$ptitl_rNcNaOrgRisk_FN							,	# table title
							\@colorTable_4c									,	# list of specified fill colors
							\@colorTable_4c									,	# list of specified edge colors
							\"Species-Quartet Number"						,	# label x-axis
							\"Ratio Nc/Na"									,	# labe y-axis
							\"Tree-Rank"									,	# label legend title color fill (area), neglected if undef
							\"Tree-Rank"									,	# label legend title color scale (outlines), neglected if undef
							\"Ratio Nc/Na"									,	# label legend title shape
							\"Ratio Nc/Na"									,	# size shape label, neglected if undef
							\30												,	# height of pdf graphic, keep undef for def height: 15
							\""												,	# width of pdf graphic, keep undef for def width: 40
							\$pdf_suffix									,	# pdf file suffix to avoid different plots of identic names
							\$href_subfolder_q->{svg}						,	# svg print folder
							\$href_subsubfolder_q->{pdf_single}				,	# pdf print folder
							\@$aref_pdf_pattern								,	# pdf sampling for pdfunite
						) ;
						
						my	$ftype	=	$main."_".$cm."_".$fm	;
						push @{$href_hol_pdf_of_type->{$ftype}}, $aref_pdf_pattern->[-1] ;
					}
				}
			}
		}
	}
	#####################################################################################
	#####################################################################################
	
	##################################################################################### Q5
	# R lineplot of 4-clan individual taxon support (mean or median)                   ## works ## main userversion (adapted)
	#####################################################################################
	if	( 	( $href_print_Rplot_quartets->{pline_taxreject_clan4}				== 1 ) &&
			( $href_value_of_scr_par->{pdf_prtResult}{main} 					== 1 ) ){
		
		for my $tsv_4clan_metafileQ5 ( sort keys %{$href_hoh_tsvdata_of_tsvfilename_of_printcode->{metafileQ5}} ){
			
			my	$cm_fm		=	$href_hoh_tsvdata_of_tsvfilename_of_printcode->{metafileQ5}{$tsv_4clan_metafileQ5} ;
			my	@prts		=	split "_", $cm_fm ;	
			my	$filter		=	&recode_filter_printout( \$cm_fm )	;
			my	$ftype		=	$main."_".$href_hoh_tsvdata_of_tsvfilename_of_printcode->{metafileQ5}{$tsv_4clan_metafileQ5} ;
			my	$titl		=	"Plot Q5 -- Count of '".$filter."'-rejected quartet participations of single species for each clade-quartet (".$prts[0].")";
				
			&r_lineplot(	
				
				\$tsv_4clan_metafileQ5						,	# name of infile table (*.tsv)
				\"Species"									,	# table header of x-plot categories
				\"N_Rejected_Quartets"						,	# table header of y-plot values
				\"Clade"									,	# table header (factor) for value separation given different point colors within each plot
				\""											,	# table header (factor) for value separation given different line colors within each plot
				\""											,	# table header (factor) for showing different points in different shapes (neglegted if undef)
				\1											,	# do lines between plot points -> 1; otherwise, no line connections between points
				\1											,	# size lines, default size if undef
				\1											,	# size points, default size if undef
				\"facet"									,	# if defined, facet grid is activated in respect of the next two group parameters
				\"Clade"									,	# table header plot-splitting vertically (defined as dot "." if splitting not needed)
				\"Clade_Quartet"							,	# table header plot-splitting horizontally (undefined "" if splitting not needed)
				\"free"										,	# default (undef) all the panels have the same scales ('fixed') can be made independent if defined 'free', 'free_x', or 'free_y'
				\""											,	# highlight horizontal line at axis interception; neglected if undef
				\"dashed"									,	# horizontal linetype: solid, longdash, dashed, blank, dotdash, dotted, twodash, default if undef
				\2											,	# horizontal linesize; default if undef
				\"gray49"									,	# horizontal line color, default if undef
				\""											,	# highlight vertical line at axis interception; neglected if undef
				\"dashed"									,	# vertical linetype: solid, longdash, dashed, blank, dotdash, dotted, twodash, default if undef
				\2											,	# vertical linesize; default if undef
				\"gray49"									,	# vertical line color, default if undef
				\$titl										,	# table title
				\@colorTable_12c							,	# list of specified line colors
				\"Species"									,	# label x-axis
				\1											,	# if value == 1 -> rotate x-axis labels via 90°, undef (0) keeps labels horizontally
				\"N Rejected Quartet Particpations"			,	# labe y-axis
				\"Clade"									,	# point colour label
				\"Clade"									,	# line (outline) colour label
				\""											,	# size shape label 
				\15											,	# height of pdf graphic, keep undef for def height: 15
				\40											,	# width of pdf graphic, keep undef for def width: 40
				\"rejected_in_total_lineplot"				,	# pdf file suffix to avoid different plots of identic names
				\$href_subfolder_q->{svg}					,	# svg print folder
				\$href_subsubfolder_q->{pdf_single}			,	# pdf print folder
				\@$aref_pdf_pattern							,	# pdf sampling for pdfunite
			) ; 
			
			push @{$href_hol_pdf_of_type->{$ftype}}, $aref_pdf_pattern->[-1] ;
		}
	}
	#####################################################################################
	#####################################################################################
	
	##################################################################################### Q5
	# R lineplot of 4-clan individual taxon support (mean or median)                   ## works ## main userversion (adapted)
	#####################################################################################
	if	( 	( $href_print_Rplot_quartets->{pline_taxreject_total}				== 1 ) &&
			( $href_value_of_scr_par->{pdf_prtResult}{main} 					== 1 ) ){
		
		for my $tsv_4clan_metafileQ5 ( sort keys %{$href_hoh_tsvdata_of_tsvfilename_of_printcode->{metafileQ5}} ){
			
			my	$fm		=	$href_hoh_tsvdata_of_tsvfilename_of_printcode->{metafileQ5}{$tsv_4clan_metafileQ5} ;
			my	@prts	=	split "_", $fm ;	
			my	$filter	=	&recode_filter_printout( \$fm )	;
			my	$ftype	=	$main."_".$href_hoh_tsvdata_of_tsvfilename_of_printcode->{metafileQ5}{$tsv_4clan_metafileQ5} ;
			my	$titl	=	"Plot Q5: Total count of '".$filter."'-rejected quartet participations of single species across all clade-quartets (".$prts[0].")";
				
			&r_lineplot(	
				
				\$tsv_4clan_metafileQ5						,	# name of infile table (*.tsv)
				\"Species"										,	# table header of x-plot categories
				\"N_Rejected_Quartets_Total"				,	# table header of y-plot values
				\"Clade"									,	# table header (factor) for value separation given different point colors within each plot
				\""											,	# table header (factor) for value separation given different line colors within each plot
				\""											,	# table header (factor) for showing different points in different shapes (neglegted if undef)
				\1											,	# do lines between plot points -> 1; otherwise, no line connections between points
				\1											,	# size lines, default size if undef
				\1											,	# size points, default size if undef
				\"facet"									,	# if defined, facet grid is activated in respect of the next two group parameters
				\"Clade"									,	# table header plot-splitting vertically (defined as dot "." if splitting not needed)
				\""											,	# table header plot-splitting horizontally (undefined "" if splitting not needed)
				\"free"										,	# default (undef) all the panels have the same scales ('fixed') can be made independent if defined 'free', 'free_x', or 'free_y'
				\""											,	# highlight horizontal line at axis interception; neglected if undef
				\"dashed"									,	# horizontal linetype: solid, longdash, dashed, blank, dotdash, dotted, twodash, default if undef
				\2											,	# horizontal linesize; default if undef
				\"gray49"									,	# horizontal line color, default if undef
				\""											,	# highlight vertical line at axis interception; neglected if undef
				\"dashed"									,	# vertical linetype: solid, longdash, dashed, blank, dotdash, dotted, twodash, default if undef
				\2											,	# vertical linesize; default if undef
				\"gray49"									,	# vertical line color, default if undef
				\$titl										,	# table title
				\@colorTable_12c							,	# list of specified line colors
				\"Species"									,	# label x-axis
				\1											,	# if value == 1 -> rotate x-axis labels via 90°, undef (0) keeps labels horizontally
				\"N Total Rejected Quartet Particpations"	,	# labe y-axis
				\"Clade"									,	# point colour label
				\"Clade"									,	# line (outline) colour label
				\""											,	# size shape label 
				\15											,	# height of pdf graphic, keep undef for def height: 15
				\40											,	# width of pdf graphic, keep undef for def width: 40
				\"rejected_lineplot"						,	# pdf file suffix to avoid different plots of identic names
				\$href_subfolder_q->{svg}					,	# svg print folder
				\$href_subsubfolder_q->{pdf_single}			,	# pdf print folder
				\@$aref_pdf_pattern							,	# pdf sampling for pdfunite
			) ; 
			
			push @{$href_hol_pdf_of_type->{$ftype}}, $aref_pdf_pattern->[-1] ;
		}
	}
	#####################################################################################
	#####################################################################################
	
	##################################################################################### Q5
	# R lineplot of 4-clan individual taxon support (mean or median)                   ## works ## main userversion (adapted)
	#####################################################################################
	if	( 	( $href_print_Rplot_quartets->{pline_taxremain_clan4}				== 1 ) &&
			( $href_value_of_scr_par->{pdf_prtResult}{main} 					== 1 ) ){
		
		for my $tsv_4clan_metafileQ5 ( sort keys %{$href_hoh_tsvdata_of_tsvfilename_of_printcode->{metafileQ5}} ){
			
			my	$cm_fm		=	$href_hoh_tsvdata_of_tsvfilename_of_printcode->{metafileQ5}{$tsv_4clan_metafileQ5} ;
			my	@prts		=	split "_", $cm_fm ;	
			
			my	$filter	=	&recode_filter_printout( \$cm_fm )	;
			my	$ftype	=	$main."_".$href_hoh_tsvdata_of_tsvfilename_of_printcode->{metafileQ5}{$tsv_4clan_metafileQ5} ;
			my	$titl	=	"Plot Q5 -- Count of '".$filter."'-remaining quartet participations of single species for each clade-quartet (".$prts[0].")";
				
			&r_lineplot(	
				
				\$tsv_4clan_metafileQ5						,	# name of infile table (*.tsv)
				\"Species"									,	# table header of x-plot categories
				\"N_Retained_Quartets"						,	# table header of y-plot values
				\"Clade"									,	# table header (factor) for value separation given different point colors within each plot
				\""											,	# table header (factor) for value separation given different line colors within each plot
				\""											,	# table header (factor) for showing different points in different shapes (neglegted if undef)
				\1											,	# do lines between plot points -> 1; otherwise, no line connections between points
				\1											,	# size lines, default size if undef
				\1											,	# size points, default size if undef
				\"facet"									,	# if defined, facet grid is activated in respect of the next two group parameters
				\"Clade"									,	# table header plot-splitting vertically (defined as dot "." if splitting not needed)
				\"Clade_Quartet"							,	# table header plot-splitting horizontally (undefined "" if splitting not needed)
				\"free"										,	# default (undef) all the panels have the same scales ('fixed') can be made independent if defined 'free', 'free_x', or 'free_y'
				\""											,	# highlight horizontal line at axis interception; neglected if undef
				\"dashed"									,	# horizontal linetype: solid, longdash, dashed, blank, dotdash, dotted, twodash, default if undef
				\2											,	# horizontal linesize; default if undef
				\"gray49"									,	# horizontal line color, default if undef
				\""											,	# highlight vertical line at axis interception; neglected if undef
				\"dashed"									,	# vertical linetype: solid, longdash, dashed, blank, dotdash, dotted, twodash, default if undef
				\2											,	# vertical linesize; default if undef
				\"gray49"									,	# vertical line color, default if undef
				\$titl										,	# table title
				\@colorTable_12c							,	# list of specified line colors
				\"Species"									,	# label x-axis
				\1											,	# if value == 1 -> rotate x-axis labels via 90°, undef (0) keeps labels horizontally
				\"N Remaining Quartet Participations"		,	# labe y-axis
				\"Clade"									,	# point colour label
				\"Clade"									,	# line (outline) colour label
				\""											,	# size shape label 
				\15											,	# height of pdf graphic, keep undef for def height: 15
				\40											,	# width of pdf graphic, keep undef for def width: 40
				\"remaining_lineplot"						,	# pdf file suffix to avoid different plots of identic names
				\$href_subfolder_q->{svg}					,	# svg print folder
				\$href_subsubfolder_q->{pdf_single}			,	# pdf print folder
				\@$aref_pdf_pattern							,	# pdf sampling for pdfunite
			) ; 
			
			push @{$href_hol_pdf_of_type->{$ftype}}, $aref_pdf_pattern->[-1] ;
		}
	}
	#####################################################################################
	#####################################################################################
	
	##################################################################################### Q5
	# R lineplot of 4-clan individual taxon support (mean or median)                   ## works ## main userversion (adapted)
	#####################################################################################
	if	( 	( $href_print_Rplot_quartets->{pline_taxremain_total}				== 1 ) &&
			( $href_value_of_scr_par->{pdf_prtResult}{main} 					== 1 ) ){
		
		for my $tsv_4clan_metafileQ5 ( sort keys %{$href_hoh_tsvdata_of_tsvfilename_of_printcode->{metafileQ5}} ){#former L
			
			my	$fm				=	$href_hoh_tsvdata_of_tsvfilename_of_printcode->{metafileQ5}{$tsv_4clan_metafileQ5} ;
			my	@prts			=	split "_", $fm ;	
			my	$qfilter_type	=	&recode_filter_printout( \$fm )	;
			my	$ftype			=	$main."_".$href_hoh_tsvdata_of_tsvfilename_of_printcode->{metafileQ5}{$tsv_4clan_metafileQ5} ;
			my	$titl			=	"Plot Q5 -- Total count of '".$qfilter_type."'-remaining quartet participations of single species across all clade-quartets (".$prts[0].")";
				
			&r_lineplot(	
				
				\$tsv_4clan_metafileQ5						,	# name of infile table (*.tsv)
				\"Species"									,	# table header of x-plot categories
				\"N_Retained_Quartets_Total"				,	# table header of y-plot values
				\"Clade"									,	# table header (factor) for value separation given different point colors within each plot
				\""											,	# table header (factor) for value separation given different line colors within each plot
				\""											,	# table header (factor) for showing different points in different shapes (neglegted if undef)
				\1											,	# do lines between plot points -> 1; otherwise, no line connections between points
				\1											,	# size lines, default size if undef
				\1											,	# size points, default size if undef
				\"facet"									,	# if defined, facet grid is activated in respect of the next two group parameters
				\"Clade"									,	# table header plot-splitting vertically (defined as dot "." if splitting not needed)
				\""											,	# table header plot-splitting horizontally (undefined "" if splitting not needed)
				\"free"										,	# default (undef) all the panels have the same scales ('fixed') can be made independent if defined 'free', 'free_x', or 'free_y'
				\""											,	# highlight horizontal line at axis interception; neglected if undef
				\"dashed"									,	# horizontal linetype: solid, longdash, dashed, blank, dotdash, dotted, twodash, default if undef
				\2											,	# horizontal linesize; default if undef
				\"gray49"									,	# horizontal line color, default if undef
				\""											,	# highlight vertical line at axis interception; neglected if undef
				\"dashed"									,	# vertical linetype: solid, longdash, dashed, blank, dotdash, dotted, twodash, default if undef
				\2											,	# vertical linesize; default if undef
				\"gray49"									,	# vertical line color, default if undef
				\$titl										,	# table title
				\@colorTable_12c							,	# list of specified line colors
				\"Species"									,	# label x-axis
				\1											,	# if value == 1 -> rotate x-axis labels via 90°, undef (0) keeps labels horizontally
				\"N Total Remaining Quartet Participations"	,	# labe y-axis
				\"Clade"									,	# point colour label
				\"Clade"									,	# line (outline) colour label
				\""											,	# size shape label 
				\15											,	# height of pdf graphic, keep undef for def height: 15
				\40											,	# width of pdf graphic, keep undef for def width: 40
				\"remaining_in_total_lineplot"				,	# pdf file suffix to avoid different plots of identic names
				\$href_subfolder_q->{svg}					,	# svg print folder
				\$href_subsubfolder_q->{pdf_single}			,	# pdf print folder
				\@$aref_pdf_pattern							,	# pdf sampling for pdfunite
			) ; 
			
			push @{$href_hol_pdf_of_type->{$ftype}}, $aref_pdf_pattern->[-1] ;
		}
	}
	#####################################################################################
	#####################################################################################
	
	##################################################################################### Q6
	## R barchart-plot of Penguin individually filtered N quartets                     ## works ## main userversion (adapted)
	#####################################################################################
	if	(	( $href_print_Rplot_quartets->{pbarplot_filtered_Nquartets}			== 1 ) &&
			( $href_value_of_scr_par->{pdf_prtResult}{main} 					== 1 ) ){
		
		for my $tsv_4clan_metafileQ6 ( sort keys %{$hoh_tsvdata_of_tsvfilename_of_printcode{metafileQ6}} ){
		
			my	$ctype		= $hoh_tsvdata_of_tsvfilename_of_printcode{metafileQ6}{$tsv_4clan_metafileQ6}	;
			my	$titl		= "Plot Q6 -- Count of initial and rejected species-quartets in each clade-quartet (".$ctype.")" ;
			my	$ftype	=	$main."_".$ctype ;
	
			&r_barplot(	
			
				\$tsv_4clan_metafileQ6						,	# name of infile table (*.tsv)
				\"Filter"									,	# table header of x-plot categories
				\"Number"									,	# table header of y-plot values
				\"Species_Quartet"							,	# table header (factor) for value separation given different fill colors
				\"Species_Quartet"							,	# table header (factor) for boxplot frames given different or a single line color(s)
				\"facet"									,	# if defined, facet grid is activated in respect of the next two group parameters
				\"Clade_Quartet"							,	# table header plot-splitting vertically (defined as dot "." if splitting not needed)
				\""											,	# table header plot-splitting horizontally (undefined "" if splitting not needed)
				\"fixed"									,	# default (undef) all the panels have the same scales ('fixed') can be made independent if defined 'free', 'free_x', or 'free_y'
				\""											,	# highlight horizontal line at axis (multiple lines comma separated) interception; neglected if undef
				\"dashed"									,	# horizontal linetype: solid, longdash, dashed, blank, dotdash, dotted, twodash, default if undef
				\1											,	# horizontal linesize; default if undef
				\"gray49"									,	# horizontal line color, default if undef
				\1											,	# 0 -> a stacked barplot is created by default, 1 -> non-stagged barplot with each value of the same x-plot category in a own x-axis position
				\1											,	# 0 -> no labels at barchart, 1 -> label (geom_text) at barchart
				\"Number"									,	# table header of barchart label
				\"-0.1"										,	# vertical shift value (integer, float, positive or negative) of barchart label
				\"black"									,	# barchart label color
				\"5"										,	# barchart label size 
				\$titl										,	# table title
				\@colorTable_2cc							,	# list of specified fill colors, or if undefined using a grey-scaling
				\@colorTable_1c								,	# list of specified line color(s)
				\"Filter Method"							,	# label x-axis
				\0											,	# if value == 1 -> rotate x-axis labels via 90°, undef (0) keeps labels horizontally
				\"N Species-Quartets"						,	# labe y-axis
				\"Species-Quartet"							,	# box colour label
				\40											,	# width of pdf graphic, keep undef for def width: 40
				\15											,	# height of pdf graphic, keep undef for def height: 15
				\"barchart"									,	# pdf file suffix to avoid different plots of identic names
				\$href_subfolder_q->{svg}					,	# svg print folder
				\$href_subsubfolder_q->{pdf_single}			,	# pdf print folder
				\@$aref_pdf_pattern							,	# pdf sampling for pdfunite
			) ;
	
			#push @{$href_hol_pdf_of_type->{$main}}, $aref_pdf_pattern->[-1] ;
			push @{$href_hol_pdf_of_type->{$ftype}}, $aref_pdf_pattern->[-1] ;
		}
	}
	#####################################################################################
	#####################################################################################
		
	
	###########################################################################################################################
	############### Plot MAIN taxon results of 4clan-quartet analyses                                                   #######
	###########################################################################################################################
	
	##################################################################################### T1
	## R ternary plot of 4-clan individual taxon support (mean and median)             ## works ## main
	#####################################################################################
	if	( 	( $href_print_Rplot_quartets->{ptern_taxsupport_clan4}				== 1 ) &&
			( $href_value_of_scr_par->{pdf_prtResult}{main} 					== 1 ) ){
		
		for my $tsv_4clan_metafileT1 ( sort keys %{$hoh_tsvdata_of_tsvfilename_of_printcode{metafileT1}} ){
			
			my	$ftype	=	$main."_".$href_hoh_tsvdata_of_tsvfilename_of_printcode->{metafileT1}{$tsv_4clan_metafileT1}	;
			my	$type	=	$href_hoh_tsvdata_of_tsvfilename_of_printcode->{metafileT1}{$tsv_4clan_metafileT1}				;
			my	@cm_fm	=	split "_", $type  ;
			
			my	$qfilter_type	=	&recode_filter_printout( \$type ) ; my $titl ;
			if	( $qfilter_type =~ /unfiltered/ ){ $titl	=	"Plot T1 -- Unfiltered species support contributions for individual species-quartet-related clade-quartet trees (".$cm_fm[0].")" }
			else								 { $titl	=	"Plot T1 -- '".$qfilter_type."'-filtered species support contributions for individual species-quartet-related clade-quartet trees (".$cm_fm[0].")" }
			
			&r_ternplot(	
			
				\$tsv_4clan_metafileT1							,	# name of infile table (*.tsv)
				\"QT1"											,	# x-axis
				\"QT2"											,	# y-axis
				\"QT3"											,	# z-axis
				\"Clade"										,	# table header (factor) for value separation given different point colors within each plot
				\"facet"										,	# if defined, facet grid is activated in respect of the next two group parameters
				\"Clade_Quartet"								,	# table header plot-splitting vertically (defined as dot "." if splitting not needed)
				\"Average"										,	# table header plot-splitting horizontally (undefined "" if splitting not needed)
				\$titl											,	# plot title
				\4												,	# point size within triangle (if undef = 6)
				\20												,	# height of pdf graphic, keep undef for def height: 15
				\40												,	# width of pdf graphic, keep undef for def height: 40
				\"ternaryplot"									,	# pdf file suffix to avoid different plots of identic names
				\$href_subfolder_q->{svg}						,	# svg print folder
				\$href_subsubfolder_q->{pdf_single}				,	# pdf print folder
				\@$aref_pdf_pattern								,	# pdf sampling for pdfunite
			) ;
			
			#push @{$href_hol_pdf_of_type->{$main}}, $aref_pdf_pattern->[-1] ;
			push @{$href_hol_pdf_of_type->{$ftype}}, $aref_pdf_pattern->[-1] ;
		}
	}
	#####################################################################################
	#####################################################################################

	##################################################################################### T2
	# R lineplot of 4-clan individual taxon support (mean or median)                   ## works ## main
	#####################################################################################
	if	( 	( $href_print_Rplot_quartets->{pline_taxsupport_clan4}				== 1 ) &&
			( $href_value_of_scr_par->{pdf_prtResult}{main} 					== 1 ) ){
		
		for my $tsv_4clan_metafileT2 ( sort keys %{$hoh_tsvdata_of_tsvfilename_of_printcode{metafileT2}} ){
			
			my	@file_info	=	split ":", $href_hoh_tsvdata_of_tsvfilename_of_printcode->{metafileT2}{$tsv_4clan_metafileT2}	;
			my	$ftype		=	$main."_".$file_info[0]	;
			my	$fm			=	$file_info[0]	;
			my	@cm_fm		=	split "_", $fm  ;
			my	$stype		=	$file_info[1]	;
			
			my	$qfilter_type	=	&recode_filter_printout( \$fm )	; my $titl ;
			if	( $qfilter_type =~ /unfiltered/ ){ $titl	=	"Plot T2 -- Unfiltered species (".$stype.") support contributions for individual species-quartet-related clade-quartet trees (".$cm_fm[0].")" }
			else								 { $titl	=	"Plot T2 -- '".$qfilter_type."'-filtered species (".$stype.") support contributions for individual species-quartet-related clade-quartet trees (".$cm_fm[0].")" }
				
			&r_lineplot(	
				
				\$tsv_4clan_metafileT2						,	# name of infile table (*.tsv)
				\"Species"									,	# table header of x-plot categories
				\"Support"									,	# table header of y-plot values
				\"Tree_Code"								,	# table header (factor) for value separation given different point colors within each plot
				\""											,	# table header (factor) for value separation given different line colors within each plot
				\""											,	# table header (factor) for showing different points in different shapes (neglegted if undef)
				\""											,	# do lines between plot points -> 1; otherwise, no line connections between points
				\1											,	# size lines, default size if undef
				\1											,	# size points, default size if undef
				\"facet"									,	# if defined, facet grid is activated in respect of the next two group parameters
				\"Clade"									,	# table header plot-splitting vertically (defined as dot "." if splitting not needed)
				\"Clade_Quartet"							,	# table header plot-splitting horizontally (undefined "" if splitting not needed)
				\"free"										,	# default (undef) all the panels have the same scales ('fixed') can be made independent if defined 'free', 'free_x', or 'free_y'
				\"0.6, 0.4"									,	# highlight horizontal line at axis interception; neglected if undef
				\"dashed"									,	# horizontal linetype: solid, longdash, dashed, blank, dotdash, dotted, twodash, default if undef
				\2											,	# horizontal linesize; default if undef
				\"gray49"									,	# horizontal line color, default if undef
				\""											,	# highlight vertical line at axis interception; neglected if undef
				\"dashed"									,	# vertical linetype: solid, longdash, dashed, blank, dotdash, dotted, twodash, default if undef
				\2											,	# vertical linesize; default if undef
				\"gray49"									,	# vertical line color, default if undef
				\$titl										,	# table title
				\@colorTable_12c							,	# list of specified line colors
				\"Species"									,	# label x-axis
				\1											,	# if value == 1 -> rotate x-axis labels via 90°, undef (0) keeps labels horizontally
				\"Tree Support"								,	# labe y-axis
				\"Clade-Quartet Tree"						,	# point colour label
				\""											,	# line (outline) colour label
				\""											,	# size shape label 
				\15											,	# height of pdf graphic, keep undef for def height: 15
				\40											,	# width of pdf graphic, keep undef for def width: 40
				\"lineplot"									,	# pdf file suffix to avoid different plots of identic names
				\$href_subfolder_q->{svg}					,	# svg print folder
				\$href_subsubfolder_q->{pdf_single}			,	# pdf print folder
				\@$aref_pdf_pattern							,	# pdf sampling for pdfunite
			) ; 
			
			#push @{$href_hol_pdf_of_type->{$main}}, $aref_pdf_pattern->[-1] ;
			push @{$href_hol_pdf_of_type->{$ftype}}, $aref_pdf_pattern->[-1] ;
		}
	}
	#####################################################################################
	#####################################################################################

	##################################################################################### Q7
	## R pointplot best 4clan-tree support related to single quartet-seq lengths       ## works ## main
	#####################################################################################
	if	( 	( $href_print_Rplot_quartets->{ppoint_QseqLength_best_tree}			== 1 ) &&
			( $href_value_of_scr_par->{pdf_prtResult}{main} 					== 1 ) ){
		
		for my	$tsv_mclan_metafileQ7	( keys %{$href_hoh_tsvdata_of_tsvfilename_of_printcode->{metafileQ7}} ){
			
			my	@prts	=	split "::", $href_hoh_tsvdata_of_tsvfilename_of_printcode->{metafileQ7}{$tsv_mclan_metafileQ7} ;
			my	$clan4	=	shift @prts ;
			my	$nap	=	shift @prts ;
				@prts	=	sort  @prts ;
			
			
			for my $method ( 0 .. @prts-1 ){
				
				my	@prts_fm			=	split "\t", $prts[$method] ;
				
				my	$filter				=	$prts_fm[0] ;
				my	$ptitl_dNap_Nap		=	"Plot Q7: Best tree support and the number of site-patterns analyzed for '".$filter."'-remaining species-quartets of clade-quartet ".$clan4." (".$nap.")" ;
				
				my	$x_par				=	"N_informative_sites"			;
				my	$y_par				=	$filter						 	;
				my	$highlight			=	"Best_Tree" 					;
				my	$pdf_suffix			=	$filter."_pointplot"	;
				
				&r_pointplot(	
							
							\$tsv_mclan_metafileQ7								,	# name of infile table (*.tsv)
							\$x_par												,	# table header of x-plot categories
							\$y_par												,	# table header of y-plot values
							\$highlight											,	# table header (factor) for value separation given different fill colors within each plot (neglegted if undef)
							\""													,	# table header (factor) for value separation given different point shapes within each plot (neglegted if undef)
							\""													,	# table header (factor) for showing different elements in different sizes (neglegted if undef)
							\"smooth"											,	# if defined, smooth regular expression line is activated (method "lm"), if undef expression line is deactivated
							\$highlight											,	# table header (factor) for smooth regression line colour separation (neglegted if undef, smooth must be activated)
							\""													,	# if 1 -> error intervals of smooth line are colored due to smooth factor, undef -> default is used
							\"facet"											,	# if defined, facet grid is activated in respect of the next two group parameters
							\""													,	# table header facet grid plot-splitting vertically (undefined "" if splitting not needed)
							\"Clade_Quartet"									,	# table header facet grid plot-splitting horizontally (undefined "" if splitting not needed)
							\"free"												,	# default (undef) all the panels have the same scales ('fixed') can be made independent if defined 'free', 'free_x', or 'free_y'
							\"0.6,0.4"											,	# highlight horizontal line at axis interception; neglected if undef
							\"dashed"											,	# horizontal linetype: solid, longdash, dashed, blank, dotdash, dotted, twodash, default if undef
							\2													,	# horizontal linesize; default if undef
							\"gray49"											,	# horizontal line color, default if undef
							\""													,	# highlight vertical line at axis interception; neglected if undef
							\"dashed"											,	# vertical linetype: solid, longdash, dashed, blank, dotdash, dotted, twodash, default if undef
							\2													,	# vertical linesize; default if undef
							\"gray49"											,	# vertical line color, default if undef
							\$ptitl_dNap_Nap									,	# table title
							\@colorTable_3cb									,	# list of specified fill colors
							\@colorTable_3cb									,	# list of specified edge colors
							\"N Site-Pattern (Analysed)"						,	# label x-axis
							\"Best Tree Support"								,	# labe y-axis
							\"Best Tree"										,	# label legend title color fill (area), neglected if undef
							\""													,	# label legend title color scale (outlines), neglected if undef
							\""													,	# label legend title shape
							\""													,	# size shape label, neglected if undef
							\15													,	# height of pdf graphic, keep undef for def height: 15
							\""													,	# width of pdf graphic, keep undef for def width: 40
							\$pdf_suffix										,	# pdf file suffix to avoid different plots of identic names
							\$href_subfolder_q->{svg}							,	# svg print folder
							\$href_subsubfolder_q->{pdf_single}					,	# pdf print folder
							\@$aref_pdf_pattern									,	# pdf sampling for pdfunite
				) ;
				
				
				my	$ftype	= $main."_".$prts_fm[1] ;
				#print "\nftype: ",$ftype ,"\n";
				push @{$href_hol_pdf_of_type->{$ftype}}, $aref_pdf_pattern->[-1] ;
			}
		}
	}
	#####################################################################################
	#####################################################################################
	
	###########################################################################################################################
	###########################################################################################################################
}

sub print_and_processing_multiclan_support_matrix{
	
	my	$href_hoh_score_rooted_sister_pair	= $_[0]	; # key 1: rooted clan (y-axis e.g. RO), key 2: sister pair of rooted clan (x-axis e.g. AB), value: phyquart identified support score		IN: defined; OUT: unchanged
	my	$href_seen_sister_pair				= $_[1]	; # key : sister_pair combination (e.g. 'AB', value: counter of appearance in multiple clan analyses										IN: defined; OUT: unchanged
	my	$href_subfolder_q					= $_[2]	; # hashpath to output-type specific subdirs of single quartet-pattern analyses; key: output-type (e.g. 'pdf'); value: pathdir				IN: defined; OUT: unchanged
	my	$sref_sf_type						= $_[3]	; # code of quartet-calculation and filter-method	(e.g. 'nap_risk1')																		IN: defined; OUT: unchanged
	my	$sref_clanfile						= $_[4]	; # name of given clanfile																													IN: defined; OUT: unchanged
	my	$sref_outgroup_clan					= $_[5]	; # outgroup clan-code																														IN: defined; OUT: unchanged
	my	$href_val_of_scr_par				= $_[6] ; # key1: sealion user available parameter options; value: optional parameter setting														IN: defined; OUT: unchanged
	my	$sref_outfile_prefix				= $_[7]	; # filename-prefix of tree-fishing resultfiles																								IN: undefined; OUT: defined
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# &treefishing_with_icebreaker
	######################################################################
	
	######################################################################
	## Build multi-clan support matrix                                  ##
	######################################################################
	
	####################################
	# Matrix preparation
	my	@sister_pair	=	sort keys %{$href_seen_sister_pair}		;
	my	$line_1_xaxis	=	join "\t",	@sister_pair				;
	my	$line_allother	=	()										;
			
	for my $root_clan ( sort keys %{$href_hoh_score_rooted_sister_pair} ){ 													
				
		$line_allother	.=	"\n".$root_clan ;
		for my $pairing	( @sister_pair ){																					
					
			if		( $href_hoh_score_rooted_sister_pair->{$root_clan}{$pairing} ){
						
				#print "\ndef $root_clan\t$pairing:\t", $hoh_score_of_rooted_clan_of_sister_pair{$root_clan}{$pairing};
				$line_allother	.=	"\t".$href_hoh_score_rooted_sister_pair->{$root_clan}{$pairing}
			}
			
			else	{	$line_allother	.=	"\t0" }
		}
	}
	####################################
	
	####################################
	# matrix print-out
	## if defined outgroup, print taxon in first line of the matrix outfile
	my	$text = "\t".$line_1_xaxis.$line_allother;
	if ( $href_val_of_scr_par->{sys_parameter}{screen_prt} == 1 ){ 
		
		print	"\n\t::&print_and_processing_multiclan_support_matrix\n" ,
				"\n\t\t", $text 
	}
	##
	
	## print out
	( my	$clandeffile	=	$$sref_clanfile ) =~ s/.txt$// ;
	
			$$sref_outfile_prefix	=	$href_subfolder_q->{mcm}."/SeaLion_mclansmatrix_".$clandeffile."_".$$sref_sf_type ;
	my		$outfile				=	$$sref_outfile_prefix.".txt" ;
	
	open	OUTmat	,	">$outfile" or die "\nICEBREAKER-ERROR: Cannot open multi-clan support matrix ", $outfile, "!\n\n" ;
	print	OUTmat	$text, "\n" ;
	close	OUTmat	;
	####################################
	
	####################################
	# icebreaker processing of actual mclan-matrix and outgroup-clan
	my $icebin = $href_val_of_scr_par->{sys_parameter}{ice_bin} ;
	&treefishing_with_icebreaker( \$icebin, \$outfile, \$$sref_outgroup_clan ) ;
	####################################
}

sub print_endscreen{
	
	my	$href_value_of_scr_par	= $_[0] ;	# key1: method_code; key2: method_option; value: float or integer		IN: defined; OUT: unchanged
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	#######################
	# print end
	print	"\n\n\t------------------------------------------------------------",
			"\n\n\tSeaLion analysis completed !",
			"\n\tAll results are printed to ", $href_value_of_scr_par->{sys_parameter}{o_main} ,
			"\n\t",
			"\n\tHave a nice day." ,
			"\n\t",
			"\n" ,
			"                     -~~~~-						\n" ,
			"                    |O __ O|						\n" ,
			"                 -__|_\\__/_|__-					\n" ,
			"                 ---(__/\\__)---					\n" ,
			"                   /||~--~||\\						\n" ,
			"                  / ||    || \\					\n" ,
			"                 |  `'    `'  |                    \n" ,
			"               _/ /\\__    __/\\ \\_				\n" ,
			"             ~~~~~    ~~~~    ~~~~~				\n" ,
			"\n\t___________________________________________________\n\n",
	;
	#######################
}

sub print_script_settings{
	
	my	$href_value_of_scr_par				= $_[0] ; # key1: sealion user available parameter options; value: optional parameter setting									IN: defined; OUT: unchanged
	my	$href_subfolder_q					= $_[1] ; # key1: subfolder code ('pdf_single', 'pdf_sumed' ); value: path														IN: defined; OUT: unchanged
	my  $href_ldata_of_flatex_of_printcode	= $_[2] ; # key1: Latex table unique printcode (e.g. 'latexC'); key2: *.tex filename; value: @list of assigned table lines		IN: defined; OUT: changed	
	my	$href_textfilename_of_printcode		= $_[3] ; # key1: textfile unique printcode (e.g. 'textA'); value: filename for print											IN: defined; OUT: unchanged
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# &latex_data_sample_LP1
	######################################################################
	
	######################################################################
	# Print defined parameter setting                                    #
	######################################################################
	
	#######################
	# define print out of parameter settings
	my @printinfo = (	"\n",
						"\t---------------------------------------------------\n",
						"\tSealion Setup\n",
						"\t---------------------------------------------------\n",
						"\tMSA infile:\t\t" 		, $href_value_of_scr_par->{menu_parameter}{msa_infile}		, "\n" ,
						"\tClade infile:\t\t" 		, $href_value_of_scr_par->{menu_parameter}{clfile} 			, "\n",
						"\tOutgroup clade:\t\t" 	, $href_value_of_scr_par->{menu_parameter}{outgroup}		, "\n",
						"\tNmax quartets" 			, $href_value_of_scr_par->{menu_parameter}{maxquartets}		, "\n", 
						"\tNmin quartet sites:\t"	, $href_value_of_scr_par->{menu_parameter}{minlengths}		, "\n", 
						"\tSupport-Type:\t\t"		, $href_value_of_scr_par->{menu_parameter}{average}			, "\n", 
						"\tSPD infile DIR:\t\t"		, $href_value_of_scr_par->{menu_parameter}{spd_folder}		, "/\n", 
						"\n" ,
						"\tQUARTET-FILTER 'RISK':\n" ,
						"\t-----------------------\n",
						"\tRISK1:\t\t\t"			, $href_value_of_scr_par->{menu_parameter}{risk1} 			, "\n" ,
						"\tRISK2:\t\t\t"			, $href_value_of_scr_par->{menu_parameter}{risk2} 			, "\n" ,
						"\tRISK threshold:\t\toptimization\n", 
						"\tUpper RISK limit:\t"		, $href_value_of_scr_par->{menu_parameter}{upperlimit}		, "\n",
						"\tLower RISK limit:\t"		, $href_value_of_scr_par->{menu_parameter}{lowerlimit}		, "\n",
						"\tLower RISK scale:\t+"	, $href_value_of_scr_par->{menu_parameter}{dscale}			, "\n",
						"\n" ,
						"\tQUARTET-FILTER 'DIST':\n" ,
						"\t-----------------------\n",
						"\tDIST:\t\t\t"				, $href_value_of_scr_par->{menu_parameter}{qdist_default}	, "\n" ,
						"\tDIST threshold:\t\toptimization\n", 
						"\tUpper DIST limit:\t"		, $href_value_of_scr_par->{menu_parameter}{upperlimit_d}	, "\n",
						"\tLower DIST limit:\t"		, $href_value_of_scr_par->{menu_parameter}{lowerlimit_d}	, "\n",
						"\tUpper DIST scale:\t-"	, $href_value_of_scr_par->{menu_parameter}{dscale}			, "\n",
						"\n" ,
						"\tP4-Parameter:\n" ,
						"\t-----------------------\n",
						"\tModel (aa):\t\t"			, $href_value_of_scr_par->{menu_parameter}{model_aa}		, "\n", 
						"\tModel (nu):\t\t"			, $href_value_of_scr_par->{menu_parameter}{model_nu}		, "\n", 
						"\tStart alpha:\t\t"		, $href_value_of_scr_par->{menu_parameter}{start_alpha} 	, "\n", 
						"\tStart pINV:\t\t"			, $href_value_of_scr_par->{menu_parameter}{start_pinv}		, "\n", 
						"\n",
						"\tOutput:\n" ,
						"\t-----------------------\n",
						"\tOutput-Folder:\t\t" 		, $href_value_of_scr_par->{menu_parameter}{output_folder}	, "/\n", 
						"\tLatex-Table Print\t"		, $href_value_of_scr_par->{menu_parameter}{prt_latex}		, "\n",
						"\tRplotting\t\t"			, $href_value_of_scr_par->{menu_parameter}{prt_Rplot}		, "\n",
						"\t---------------------------------------------------\n",
	) ;
	
	my $printline = join "", @printinfo ;
	#######################
	
	#######################
	# print to terminal
	print $printline, "\n" ;#exit;
	#######################
	
	#######################
	# print parameter setting to parameter_info.txt
	my		$f_pinfo = $href_subfolder_q->{txt}."/".$href_textfilename_of_printcode->{txtA} ;
	open 	OUTtxt,	">$f_pinfo" or die "TEXTFILE-ERROR: Cannot write to ", $f_pinfo, "!\n\n" ;
	print	OUTtxt	$printline ;
	close	OUTtxt;
	#######################
	
	########################
	# send parameter to latex table print
	&latex_data_sample_LP1(
		
		\%$href_value_of_scr_par				, # key1: options; value: optional parameter setting		IN: defined; OUT: unchanged
		\%$href_ldata_of_flatex_of_printcode	, # key1: printcode; key2: fname; value: table lines		IN: defined; OUT: changed
	);
	#######################
	
	#####################################################################
}

sub print_taxon_reduced_clanfile{
	
	my	$href_prtpath							= $_[0] ; # hashpath to output-type specific subdirs of single quartet-pattern analyses; key: output-type (e.g. 'pdf'); value: pathdir	IN: defined, OUT: unchanged
	my	$href_clan_of_otu						= $_[1] ; # key1: taxon; value: assigned clan																							IN: defined; OUT: unchanged 
	my	$href_hoh_seen_confl_taxon_of_ft_of_st	= $_[2] ; # key1: taxon_filter (stric_out, strict_in, strict_all); key2: type of support (mean|median); key3: sequence-name; 			IN: defined; OUT: unchanged  
	my	$sref_new_ftype							= $_[3] ; # combined code of current quartet+filter-method and taxon-filter																IN: defined; OUT: unchanged 
	my	$href_val_of_scr_par					= $_[4] ; # key1: sealion user available parameter options; value: optional parameter setting											IN: defined; OUT: unchanged
	
	#print "::&print_taxon_reduced_clanfile\n";
	
	######################################################################
	# called from subroutine(s)...
	# &evaluation_otu_support_main
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	# print new clanfile for each filtered set of taxa (mean|median)     #
	######################################################################
	for my $stype ( keys %{$href_hoh_seen_confl_taxon_of_ft_of_st->{$$sref_new_ftype}} ){
		
		####################################
		# store each unconflicted taxon as listelement in taxon corresponding clanfile-key
		my	(
				%hol_taxa_of_clan	# key: clanfilename; value: list of clanfile assigned, 'non-conflicted' taxa
		) ;
		
		for	my $tax ( sort keys %$href_clan_of_otu ){
			
			unless ( $href_hoh_seen_confl_taxon_of_ft_of_st->{$$sref_new_ftype}{$stype}{$tax} ){ push @{$hol_taxa_of_clan{$href_clan_of_otu->{$tax}}}, $tax	}
		}
		####################################
		
		####################################
		# Print of conflict-taxon reduced clanfiles
		( 	my $cf_org = $href_val_of_scr_par->{sys_parameter}{clfile} ) =~ s/.txt$// ;
			my $cf_new = $href_prtpath->{clf}."/SeaLion_".$cf_org."_OTUred_".$$sref_new_ftype."_".$stype.".txt" ;
		
		print "\t\t-", $cf_new ,"\n" ;
		open	OUTclan, ">$cf_new" or die "Can't open $cf_new: $!"	;
		########
		
		########
		# for each clan, join corresponding taxa by commata and print clanline to clanfile
		for	my $clan ( sort keys %hol_taxa_of_clan ){
			
			my	$clantaxa = join ",", @{$hol_taxa_of_clan{$clan}} ;
			print	OUTclan	$clan.",".$clantaxa."\n" ;
		}
		########
		
		########
		# close clanfile
		close OUTclan
		##################################
	}	
	######################################################################
	######################################################################
}

sub print_tsv_MQ1{
	
	my	$href_printpath											= $_[0] ; # key1: R grafic unique printcode (e.g. 'metafileQ1'); key2: *.tsv filename; value: *tsv filename data addings 																IN: defined; OUT: unchanged
	my	$href_hoh_data_of_tsvfile_of_prcode						= $_[1] ; # key: output-type (e.g. 'pdf'); value: pathdir																																IN: defined; OUT: changed
	my	$href_fc_of_4clan										= $_[2] ; # key: quartet numbers of current 4clan combination																															IN: defined; OUT: unchanged
	my	$href_hol_sup_of_type_of_mfilter						= $_[3]	; # key1: type of supportfilter ('nap_risk1', 'napr_risk1', or 'napv_risk1' ...); key2: mean|median; value: list of 'tree\tmean|medianScore'									IN: defined; OUT: unchanged
	my	$href_val_of_scr_par									= $_[4] ; # key1: sealion user available parameter options; value: optional parameter setting																							IN: defined; OUT: unchanged
	my	$href_clan4_of_tree										= $_[5] ; # key1: 4clan tree; value: 4clan code																																			IN: defined; OUT: unchanged
	my	$href_algorithm_def_methods								= $_[6] ; # key1: methode type; key2: method code; value: integer																														IN: defined; OUT: unchanged
	
	#print "::&print_tsv_MQ1\n";
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	## print 4clan-tree support of different methods to tsv table MQ1   ##
	######################################################################
	
	####################################
	# open filehandle
	my	 	$tsv_mclan_metafileMQ1 	= $href_printpath->{tsv}."/MQ1_SeaLion_metaData_clade-quartet_average_tree_support.tsv" ;
	open	OUTtsvMQ1, ">$tsv_mclan_metafileMQ1" || die "Cannot write tsv file ", $tsv_mclan_metafileMQ1, ": $!" 	;
	####################################
	
	####################################
	# print header only once to each clanfile corresponding tsv table file
	unless ( $href_hoh_data_of_tsvfile_of_prcode->{metafileMQ1}{$tsv_mclan_metafileMQ1} ){
		
		my @methods ;
		
		########
		# defined header line of tsv outfile
		print OUTtsvMQ1	"Clade_Quartet"		,
						"\tTree"			,
						"\tAverage"			;
		
		for my $ftype ( sort keys %$href_hol_sup_of_type_of_mfilter ){ 
			
			my	$filter_type = &recode_filter_printout( \$ftype ); 
			
			$ftype =~ s/:/_/; print OUTtsvMQ1 "\t", $filter_type; push @methods, $ftype }
		
		print OUTtsvMQ1	"\n";
		########
					
		########
		# count metafile header +1
		my $fmethods = join "::", @methods ;
		
		$href_hoh_data_of_tsvfile_of_prcode->{metafileMQ1}{$tsv_mclan_metafileMQ1} = $fmethods ;
		########
	}
	####################################
	
	####################################
	# print table data for each type of support
	for my $stype ( sort @{$href_algorithm_def_methods->{qsuppo}} ){
		
		####################################
		# print trees for each 4 clan combination separately
		for my $clan4 ( sort keys %$href_fc_of_4clan ){
			
			####################################
			# go through all 4clan trees (including non $clan4 specified trees)
			for my $tree ( sort keys %$href_clan4_of_tree ){
				
				####################################
				# if a tree assigned 4clan combination matches with actual $clan4 combination...
				if ( $href_clan4_of_tree->{$tree} eq $clan4 ){
					
					####################################
					# ...go through each filtered and unfiltered analysis
					my @data_support ;
					for my $ftype ( sort keys %$href_hol_sup_of_type_of_mfilter ){ 
						
						####################################
						# if list of 'tree\tsupport' strings of method and support combination exists...
						if ( $href_hol_sup_of_type_of_mfilter->{$ftype}{$stype} ){
							
							my $flag_tree = 0 ;
							
							# go through all 'tree\tsupport' strings (including strings which do not match with actual $clan4 code
							for my $data ( @{$href_hol_sup_of_type_of_mfilter->{$ftype}{$stype}} ){
								
								# split string 'tree\tsupport'
								my @tre_sup = split "\t", $data ;
								
								# if tree of 'tree\tsupport' matches with $tree, sample corresponding support and set flag to defined
								if ( $tre_sup[0] eq $tree ){ push @data_support, $tre_sup[1]; $flag_tree = 1 }
							}
							
							# if flag is undefined, meaning that no support exists for $tree, sample 'NA'
							if ( $flag_tree == 0 ){ push @data_support, 'NA' }
						}
						####################################
						
						####################################
						# if no tree match is found, sample 'NA'
						else { push @data_support, 'NA' }
						####################################
					}
					####################################
					
					####################################
					# join data line and print it with additional data in front
					my $data_support_line = join "\t", @data_support ;
					
					print OUTtsvMQ1 $clan4, "\t", $tree, "\t", $stype, "\t", $data_support_line, "\n"
					####################################
				}
				####################################
			}
			####################################
		}
		####################################
	}
	####################################
	
	####################################
	# close filehandle
	close OUTtsvMQ1 ;
	######################################################################
}

sub print_tsv_MQ2{
	
	my	$href_printpath											=	$_[0] ;	# key1: R grafic unique printcode (e.g. 'metafileQ1'); key2: *.tsv filename; value: *tsv filename data addings 																IN: defined; OUT: unchanged
	my	$href_hoh_data_of_tsvfile_of_prcode						=	$_[1] ;	# key: output-type (e.g. 'pdf'); value: pathdir																																IN: defined; OUT: changed
	my	$sref_clan4												=	$_[2] ;	# clan4 combination of given quartet combination																															IN: defined; OUT: unchanged
	my	$href_fc_of_4clan										=	$_[3] ;	# key: quartet numbers of current 4clan combination																															IN: defined; OUT: unchanged
	my	$href_hoh_Nap_of_clan4_of_calc_of_qtree_of_fc			=	$_[4] ;	# key1: 4-clan code; key2: type NAP; key3: rooted quartet-tree; key4: quartet number; value:NAP quartet support																IN: defined; OUT: unchanged
	my	$href_value_of_scr_par									=	$_[5] ;	# key1: sealion user available parameter options; value: optional parameter setting																 							IN: defined; OUT: unchanged
	my	$href_hoh_seen_rej_of_4clan_of_ftype_of_fn				=	$_[6] ;	# key1: 4clan code; key2: type of supportfilter ('nap_risk1', 'napr_risk1', or 'napv_risk1'); key3: actual sealion quartet number; value: filter corresponding score-value	IN: defined; OUT: unchanged 
	my	$href_support_method_code								=	$_[7] ;	# list of defined support methods (nap, napv ...)																															IN: defined; OUT: unchanged
	my	$href_hoh_seen_filter_of_type							=	$_[8] ;	# key1: filter ('dist' or 'risk'); key2: combi-code method+filter ('nap_dist1', 'nap_risk2'...); value: number of filter passing quartets									IN: defined; OUT: unchanged
	
	######################################################################
	## print header and data to tsv table MQ2                           ##
	######################################################################
	#print "::&print_tsv_MQ2\n";
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	## print 4clan-tree support to tsv table MQ2                        ##
	######################################################################
	# tsv table file, listing quartet individual and overall mean and median support
	# (nap and napv) for each rooted quartet-tree, lengths bl2, and seed number
	for my $type_support ( sort keys %{$href_support_method_code->{original}} ){
		
		my	%seen_filter = () ;
		my	@inflines 	 = () ;
		
		####################################
		# check each filter-type of the different filter classes ('dist', 'risk', 'dist+risk')
		# if filter-type assigned calculation-method (nap, napr,napv) match with the actual $type_support
		# store matching filters without $type_support prefix in %seen_filter (for tsv data print) and @infiles (for R plotting)
		for my $fclass ( sort keys %{$href_hoh_seen_filter_of_type} ){
			
			for my $ftype ( sort keys %{$href_hoh_seen_filter_of_type->{$fclass}} ){
				
				#print "\nclass:", $fclass, "\ntype: ", $ftype, "\n";
				
				my	@code_prts	= split "_", $ftype 	;
				my	$cm			= shift @code_prts  	;
				my	$fm			= join "_", @code_prts  ;
				
				if  ( $cm eq  $type_support ){
				
					push 	@inflines, $fm 	;	# store codes for R plotting
					$seen_filter{$fm}++ 	;	# store codes for tsv mean|median support print (see below)
				}
			}
		}#exit; 
		####################################
		
		####################################
		# open filehandle
		my	 	$tsv_mclan_metafileMQ2 	= $href_printpath->{tsv}."/MQ2_SeaLion_metaData_clade-quartet_species-quartet-specific_tree-support_".$type_support.".tsv" ;
		open	OUTtsvMQ2, ">>$tsv_mclan_metafileMQ2" || die "Cannot write tsv file ", $tsv_mclan_metafileMQ2, ": $!" 	;
		####################################
		
		####################################
		# print header only once to each clanfile corresponding tsv table file
		unless ( $href_hoh_data_of_tsvfile_of_prcode->{metafileMQ2}{$tsv_mclan_metafileMQ2} ){
			
			########
			# defined header line of tsv outfile
			print OUTtsvMQ2	"Clade_Quartet"		,
							"\tQN"				,
							"\tScore"			,	
							"\tSupport_Type"	;
			
			for my $ftype ( @inflines ){ print	OUTtsvMQ2 "\t", $ftype }
										
			print OUTtsvMQ2	"\tQT1"				,
							"\tQT2"				,
							"\tQT3"				,
							"\n" 				;
			########
			
			########
			# count metafile header +1
			unshift @inflines, $type_support ;
			$href_hoh_data_of_tsvfile_of_prcode->{metafileMQ2}{$tsv_mclan_metafileMQ2} = join "::", @inflines ;
			########

			#for ( @inflines ){ print "\n\tinf\t", $_ } #exit;
		}
		####################################
		
		####################################
		# print table data for support and risk of single quartets (extra line for each quartet)
		for my $fn ( sort {$a<=>$b} keys %$href_fc_of_4clan ){	#print "\nfilenumber:\t", $fn;
			
			####################################
			# print general table info
			print OUTtsvMQ2		$$sref_clan4			,	# Clan4
								"\t",$fn				,	# FN
								"\t",$type_support		,	# CF_c
								"\t","species-quartet"	;	# Score_c
			####################################
			
			####################################
			# print filter result for each filter-method of defined filter-classes 
			# (single and combined) which matches the actual quartet-calculation-method
			for my $fclass ( sort keys %{$href_hoh_seen_filter_of_type} ){
				
				for my $ftype ( sort keys %{$href_hoh_seen_filter_of_type->{$fclass}} ){ #print "\n\t\ftype\t", $ftype;
					
					####################################
					# extract quartet-calculation-method (e.g. nap) and filter-method (e.g. combined 'risk3_qdist2') from $fytpe
					# e.g. $fytpe = 'nap_risk3_qdist2'
					my	@code_prts	= split "_", $ftype	    ;
					my	$m			= shift 	 @code_prts ; # e.g. 'nap'
					my	$f			= join  "_", @code_prts ; # e.g. 'risk3_qdist2'
					####################################
					
					####################################
					# print table info only if filtertype is based on $type_support coded quartet-calculation
					if ( $type_support eq $m ){ 
						
							if	( $href_hoh_seen_rej_of_4clan_of_ftype_of_fn->{$$sref_clan4}{$ftype}{$fn} ){ print OUTtsvMQ2 "\t", "reject" } else 	{ print OUTtsvMQ2 "\t", "remain" } 
					}
					####################################
				}
			}#exit;
			####################################
			
			####################################
			# print quartet related tree-support in respect of actual $type_support coded quartet-calculation-method (nap, napr, napv)
			for my $tree ( sort keys %{$href_hoh_Nap_of_clan4_of_calc_of_qtree_of_fc->{$$sref_clan4}{$type_support}} ){
					
				print OUTtsvMQ2	"\t", $href_hoh_Nap_of_clan4_of_calc_of_qtree_of_fc->{$$sref_clan4}{$type_support}{$tree}{$fn} ;	# qtree1, qtree2, qtree3
			}
			####################################
			
			####################################
			# exit lineprint with newline
			print OUTtsvMQ2	"\n";
			####################################
		}
		########
		
		########
		# print table data for mean 4clan support (extra line)
		print OUTtsvMQ2		$$sref_clan4		,	# Clan4
							"\t","NA"			,	# FN
							"\t",$type_support	,	# CF_c
							"\t","mean"			;	# Score_c
							
		for my $f_combi ( keys %seen_filter ){ print OUTtsvMQ2 "\tmean"}# print "\n\tcombi:\t", $f_combi }
		
		for my $tree ( sort keys %{$href_hoh_Nap_of_clan4_of_calc_of_qtree_of_fc->{$$sref_clan4}{$type_support}} ){
					
			print OUTtsvMQ2	"\t", $href_hoh_Nap_of_clan4_of_calc_of_qtree_of_fc->{$$sref_clan4}{$type_support}{$tree}{mean}	# qtree1, qtree2, qtree3
		}
		########
		
		########
		# print table data for median 4clan support (extra line)
		print OUTtsvMQ2		"\n"				,
							$$sref_clan4		,	# Clan4
							"\t","NA"			,	# FN
							"\t",$type_support	,	# CF_c
							"\t","median"		;	# Score_c
		
		for my $f_combi ( keys %seen_filter ){ print OUTtsvMQ2 "\tmedian"}# print "\n\tcombi:\t", $f_combi }
		
		for my $tree ( sort keys %{$href_hoh_Nap_of_clan4_of_calc_of_qtree_of_fc->{$$sref_clan4}{$type_support}} ){
					
			print OUTtsvMQ2	"\t", $href_hoh_Nap_of_clan4_of_calc_of_qtree_of_fc->{$$sref_clan4}{$type_support}{$tree}{median}	# qtree1, qtree2, qtree3
		}
		
		print OUTtsvMQ2		"\n" ;
		####################################
		
		####################################
		# close filehandle
		close OUTtsvMQ2 ;
		####################################
	} #exit;
	######################################################################
	######################################################################
}

sub print_tsv_MQ3{
	
	my	$href_prtpath						= $_[0]	; # key1: R grafic unique printcode (e.g. 'metafileQ1'); key2: *.tsv filename; value: *tsv filename data addings 							IN: defined; OUT: unchanged
	my	$href_hoh_data_of_tsvfile_of_prc	= $_[1] ; # key: output-type (e.g. 'pdf'); value: pathdir																							IN: defined; OUT: changed
	my	$sref_clfile						= $_[2] ; # name of given clanfile																													IN: defined; OUT: unchanged
	my	$sref_sm							= $_[3] ; # support method (nap, napv, napr)																										IN: defined; OUT: unchanged
	my	$sref_fm							= $_[4] ; # filter method (risk1, strict_all:hardclipping, risk2_strict_all:hardclipping)															IN: defined; OUT: unchanged
	my	$sref_btr							= $_[5] ; # best tree in nwk																														IN: defined; OUT: unchanged
	my	$sref_ts							= $_[6] ; # mean|median																																IN: defined; OUT: unchanged
	my	$sref_tr_supp						= $_[7] ; # final best tree support																													IN: defined; OUT: unchanged
	
	######################################################################
	## print header and data to tsv table MQ3                           ##
	######################################################################
	#print "::&print_tsv_MQ3\n";
	
	######################################################################
	# called from subroutine(s)...
	# &treefishing_result_processing
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	#####################################################################################
	###### print final tree and support of specified calculation and filter methods   ###
	#####################################################################################
	( my $clafile = $$sref_clfile ) =~  s/.txt$// ;
	my		$tsv_mclan_metafileMQ3	=	$href_prtpath->{tsv}."/MQ3_SeaLion_metaData_best_rooted-clade_tree_support_".$clafile.".tsv" ;
	open	OUTtsvMQ3, ">>$tsv_mclan_metafileMQ3" || die "Cannot write tsv file ", $tsv_mclan_metafileMQ3, ": $!"	;
	
	###############
	# (a) unless best tree support of non-taxon filtered clanfile has not been line printed, print defined header of tsv outfile as first line
	# (b) define header line as printed to avoid additional printing
	unless ( $href_hoh_data_of_tsvfile_of_prc->{metafileMQ3}{$tsv_mclan_metafileMQ3} ){
		
		####################################
		# (a) print defined header of tsv outfile as first line
		print	OUTtsvMQ3	"Cladefile\t", "Score\t", "Filter\t", "Best_Rooted_Clade_Tree\t", "Average\t", "Support" ;
		####################################
		
		####################################
		# (b) define header line as printed to avoid additional printing
		$href_hoh_data_of_tsvfile_of_prc->{metafileMQ3}{$tsv_mclan_metafileMQ3} = $$sref_clfile ;
		####################################
	}
	###############
	
	###############
	# recode filter codes for print out
	my	@prts 			=	split "_strict_", $$sref_fm ;
	my	$filter_type 	=	&recode_filter_printout( \$prts[0]  )	;
	
	if ( $$sref_fm =~ /_strict_/ ){
		
		my	$tfiltertype  =	&recode_filter_printout( \$$sref_fm )	;
		
			$filter_type  .= "_".$tfiltertype
	}
	###############
	
	###############
	# print data line
	print OUTtsvMQ3	"\n", $$sref_clfile, "\t", $$sref_sm, "\t", $filter_type, "\t", $$sref_btr, "\t", $$sref_ts, "\t", $$sref_tr_supp ;
	###############
	
	###############
	# close filehandle
	close OUTtsvMQ3 ;
	#####################################################################################
}

sub print_tsv_MQ4{
	
	my	$href_printpath											=	$_[0] ;	# key1: R grafic unique printcode (e.g. 'metafileQ1'); key2: *.tsv filename; value: *tsv filename data addings 																IN: defined; OUT: unchanged
	my	$href_hoh_data_of_tsvfile_of_prcode						=	$_[1] ;	# key: output-type (e.g. 'pdf'); value: pathdir																																IN: defined; OUT: changed
	my	$sref_clan4												=	$_[2] ;	# clan4 combination of given quartet combination																															IN: defined; OUT: unchanged
	my	$href_fc_of_4clan										=	$_[3] ;	# key: quartet numbers of current 4clan combination																															IN: defined; OUT: unchanged
	my	$href_hoh_Nap_of_clan4_of_calc_of_qtree_of_fc			=	$_[4] ;	# key1: 4-clan code; key2: type NAP; key3: rooted quartet-tree; key4: quartet number; value:NAP quartet support																IN: defined; OUT: unchanged
	my	$href_value_of_scr_par									=	$_[5] ;	# key1: sealion user available parameter options; value: optional parameter setting																 							IN: defined; OUT: unchanged
	my	$href_hoh_seen_rej_of_4clan_of_ftype_of_fn				=	$_[6] ;	# key1: 4clan code; key2: type of supportfilter ('nap_risk1', 'napr_risk1', or 'napv_risk1'); key3: actual sealion quartet number; value: filter corresponding score-value	IN: defined; OUT: unchanged 
	my	$href_support_method_code								=	$_[7] ;	# list of defined support methods (nap, napv ...)																															IN: defined; OUT: unchanged
	my	$href_hoh_seen_filter_of_type							=	$_[8] ;	# key1: filter ('dist' or 'risk'); key2: combi-code method+filter ('nap_dist1', 'nap_risk2'...); value: number of filter passing quartets									IN: defined; OUT: unchanged
	
	######################################################################
	## print header and data to tsv table MQ4                           ##
	######################################################################
	#print "::&print_tsv_MQ4\n";
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	## print 4clan-tree support to tsv table MQ4                        ##
	######################################################################
	# tsv table file, listing quartet individual and overall mean and median support
	# (nap and napv) for each rooted quartet-tree before and fter quartet-filtering
	for my $fclass ( sort keys %{$href_hoh_seen_filter_of_type} ){
		
		for my $ftype ( sort keys %{$href_hoh_seen_filter_of_type->{$fclass}} ){
			
			my	$qfilter_type	=	&recode_filter_printout( \$ftype )	;
			my	@code_prts		=	split "_", $ftype 	;
			my	$cm				=	shift @code_prts  	;
			my	$fm				=	join "_", @code_prts  ;
			
			#print "\nclass:", $fclass, "\ntype: ", $ftype, "\n", $cm, "\n";
			
			####################################
			# open filehandle
			my	 	$tsv_mclan_metafileMQ4 	= $href_printpath->{tsv}."/MQ4_SeaLion_metaData_clade-quartet_species-quartet-specific_tree-support_unfiltered-vs-filtered_".$qfilter_type."_".$cm.".tsv" ;
			open	OUTtsvMQ4, ">>$tsv_mclan_metafileMQ4" || die "Cannot write tsv file ", $tsv_mclan_metafileMQ4, ": $!" 	;
			####################################
			
			####################################
			# print header only once to each clanfile corresponding tsv table file
			unless ( $href_hoh_data_of_tsvfile_of_prcode->{metafileMQ4}{$tsv_mclan_metafileMQ4} ){
				
				########
				# defined header line of tsv outfile
				print OUTtsvMQ4	"Clade_Quartet"		,
								"\tQuartet_Number"	,
								"\tScore"			,
								"\tFilter"			,
								"\tSupport_Type"	,
								"\tFilter_Status"	,
								"\tQuartet_Status"	,
								"\tQT1_topology"	,
								"\tQT1"				,
								"\tQT2_topology"	,
								"\tQT2"				,
								"\tQT3_topology"	,
								"\tQT3"				,
								"\n" 				;
				########
				
				########
				# count metafile header +1
				$href_hoh_data_of_tsvfile_of_prcode->{metafileMQ4}{$tsv_mclan_metafileMQ4} = $ftype ;
				########
			}
			####################################
			
			####################################
			# print table data for support and risk of single quartets (extra line for each quartet)
			my ( $line1, $line2 ) ;
			
			for my $fn ( sort {$a<=>$b} keys %$href_fc_of_4clan ){	#print "\nfilenumber:\t", $fn;
				
				####################################
				# print general table info
				$line1	= $$sref_clan4."\t".$fn."\t".$cm."\t".$qfilter_type."\tsingle"."\t"."before" ;
				$line2	= $$sref_clan4."\t".$fn."\t".$cm."\t".$qfilter_type."\tsingle"."\t"."after"  ;
				####################################
				
				####################################
				# print risk2, risk3, and qdist1 to qdist3 filter results (uncombined and combined)
				if	( $href_hoh_seen_rej_of_4clan_of_ftype_of_fn->{$$sref_clan4}{$ftype}{$fn} ){ 
					
					$line1 .= "\t"."reject"	;
					$line2  = ()			;
				} 
				
				else { 
					
					$line1	.= "\t"."remain"; 
					$line2	.= "\t"."remain"; 
				} 
				####################################
				
				####################################
				# print quartet related tree-support in respect of actual $type_support coded quartet-calculation-method (nap, napr, napv)
				my $treeline ;
				for my $tree ( sort keys %{$href_hoh_Nap_of_clan4_of_calc_of_qtree_of_fc->{$$sref_clan4}{$cm}} ){
			
					$treeline	.= "\t".$tree."\t".$href_hoh_Nap_of_clan4_of_calc_of_qtree_of_fc->{$$sref_clan4}{$cm}{$tree}{$fn} ;	# qtree1, qtree2, qtree3
				}
				####################################
				
				####################################
				# print single quartet data line(s)
								print OUTtsvMQ4	$line1, $treeline, "\n" ;
				if ( $line2 ){  print OUTtsvMQ4	$line2, $treeline, "\n" }
				####################################
			}
			####################################
			
			########
			# print table data for mean 4clan support (extra line)
			$line1	= $$sref_clan4."\t"."NA"."\t".$cm."\t".$qfilter_type."\tmean"."\t"."before"."\t"."mean" ;
			$line2	= $$sref_clan4."\t"."NA"."\t".$cm."\t".$qfilter_type."\tmean"."\t"."after"."\t"."mean"  ;
								
			for my $tree ( sort keys %{$href_hoh_Nap_of_clan4_of_calc_of_qtree_of_fc->{$$sref_clan4}{$cm}} ){
				
				if ( $href_hoh_Nap_of_clan4_of_calc_of_qtree_of_fc->{$$sref_clan4}{$cm}{$tree}{mean} ){
							
					$line1		.= "\t".$tree."\t".$href_hoh_Nap_of_clan4_of_calc_of_qtree_of_fc->{$$sref_clan4}{$cm}{$tree}{mean}		;	# qtree1, qtree2, qtree3
				}
				else{ $line1	.= "\t"."NA"."\t"."NA" }
				
				if ( $href_hoh_Nap_of_clan4_of_calc_of_qtree_of_fc->{$$sref_clan4}{$ftype}{$tree}{mean} ){
					
					$line2		.= "\t".$tree."\t".$href_hoh_Nap_of_clan4_of_calc_of_qtree_of_fc->{$$sref_clan4}{$ftype}{$tree}{mean}	;	# qtree1, qtree2, qtree3
				}
				else{ $line2	.= "\t"."NA"."\t"."NA" }
			}
			
			print OUTtsvMQ4	$line1, "\n", $line2, "\n" ; 
			########
			
			########
			# print table data for mean 4clan support (extra line)
			$line1	= $$sref_clan4."\t"."NA"."\t".$cm."\t".$qfilter_type."\tmedian"."\t"."before"."\t"."median" ;
			$line2	= $$sref_clan4."\t"."NA"."\t".$cm."\t".$qfilter_type."\tmedian"."\t"."after"."\t"."median"  ;
								
			for my $tree ( sort keys %{$href_hoh_Nap_of_clan4_of_calc_of_qtree_of_fc->{$$sref_clan4}{$cm}} ){
						
				if ( $href_hoh_Nap_of_clan4_of_calc_of_qtree_of_fc->{$$sref_clan4}{$cm}{$tree}{median} ){
							
					$line1		.= "\t".$tree."\t".$href_hoh_Nap_of_clan4_of_calc_of_qtree_of_fc->{$$sref_clan4}{$cm}{$tree}{median}		;	# qtree1, qtree2, qtree3
				}
				else{ $line1	.= "\t"."NA"."\t"."NA" }
				
				if ( $href_hoh_Nap_of_clan4_of_calc_of_qtree_of_fc->{$$sref_clan4}{$ftype}{$tree}{median} ){
					
					$line2		.= "\t".$tree."\t".$href_hoh_Nap_of_clan4_of_calc_of_qtree_of_fc->{$$sref_clan4}{$ftype}{$tree}{median}	;	# qtree1, qtree2, qtree3
				}
				else{ $line2	.= "\t"."NA"."\t"."NA" }
			}
			
			print OUTtsvMQ4	$line1, "\n", $line2, "\n" ; 
			########
			
			########
			# close filehandle
			close OUTtsvMQ4 ;
			########
		}
	} #exit;
	######################################################################
	######################################################################
}

sub print_tsv_MQ6{
	
	my	$href_prtpath						= $_[0]  ;	# key1: R grafic unique printcode (e.g. 'metafileQ1'); key2: *.tsv filename; value: *tsv filename data addings	IN: defined; OUT: unchanged
	my	$href_hoh_data_of_tsvf_of_pcode		= $_[1]  ;	# key: output-type (e.g. 'pdf'); value: pathdir																	IN: defined; OUT: changed
	my	$aref_prtlines_MQ6					= $_[2]  ;	# tsv data lines 																								IN: defined; OUT: changed
	my	$sref_fp_bestscale					= $_[3]  ;	# overall best scale value of actual filtermethod																IN: defined; OUT: unchanged
	my 	$sref_supportfilter					= $_[4]  ;	# actual filtermethod																							IN: defined; OUT: unchanged
	
	######################################################################
	# called from subroutine(s)...
	# &filter_threshold_optimization
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	#####################################################################################
	###### print values used for matrix optimization given differen filter-scales   #####
	#####################################################################################
	
	########
	my		$filter				 	=	&recode_filter_printout( \$$sref_supportfilter )	;
	my		@cm_fm					=	split "_", $$sref_supportfilter ;
	
	my	 	$tsv_4clan_metafileMQ6	= $href_prtpath->{tsv}."/MQ6_SeaLion_metaData_clade-quartet_".$filter."_threshold-optimization_".$cm_fm[0].".tsv" ;
	open 	OUTtsvMQ6, ">>$tsv_4clan_metafileMQ6" || die "Cannot write tsv file ", $tsv_4clan_metafileMQ6, ": $!" ;
	########
	
	########
	# print header only once to each clanfile corresponding tsv table file
	unless ( $href_hoh_data_of_tsvf_of_pcode->{metafileMQ6}{$tsv_4clan_metafileMQ6} ){
		
		########
		# defined header line of tsv outfile
		print 	OUTtsvMQ6	"Score"					,
							"\tFilter"				,
							"\tClade_Quartet"		,
							"\tScale"				,
							"\tNqRemain"			,
							"\tAlpha" 				,
							"\tLambda" 				,
							"\tfPbest" 				,
							"\tfPScale" 			,
							"\tBestScaleTmp"		,
							"\tBestScaleTot"		,
							"\tRank"				,
							"\n" 					;
		########
				
		########
		# count metafile header +1
		$href_hoh_data_of_tsvf_of_pcode->{metafileMQ6}{$tsv_4clan_metafileMQ6} = $$sref_supportfilter."::".$$sref_fp_bestscale ;
		########
	}
	########
	
	########
	# print table data
	my $dataline = join "\n", @$aref_prtlines_MQ6 ; 
	print OUTtsvMQ6	$dataline, "\n" ;
	########
	
	########
	# close filehandle
	close OUTtsvMQ6 ;
	########
	
	######################################################################
}

sub print_tsv_Q3{
	
	my	$href_printpath							=	$_[0] ;		# key1: R grafic unique printcode (e.g. 'metafileQ1'); key2: *.tsv filename; value: *tsv filename data addings 				IN: defined; OUT: unchanged
	my	$href_hoh_data_of_tsvfile_of_prcode		=	$_[1] ;		# key: output-type (e.g. 'pdf'); value: pathdir																				IN: defined; OUT: unchanged
	my	$sref_cl4_code							=	$_[2] ;		# actual 4clan code																											IN: defined; OUT: unchanged 
	my	$sref_qn								=	$_[3] ;		# process number of the current quartet																						IN: defined; OUT: unchanged 
	my	$href_sing_inf_of_type					=	$_[4] ;		# key: info type (.e.g 'sNs'); value_ type value																			IN: defined; OUT: unchanged
	my	$href_sing_inf_of_type_of_tree			=	$_[5] ;		# key1: tree dependent info types (.e.g 'sNs_exp_from_qt1'); key2: rooted tree, value_ type value							IN: defined; OUT: unchanged
	my	$href_invvar_inf_of_type				=	$_[6] ;		# key: info type (.e.g 'sNs'); value_ type value																			IN: defined; OUT: unchanged
	my	$href_invvar_inf_of_type_of_tree		=	$_[7] ;		# key1: tree dependent info types (.e.g 'sNs_exp_from_qt1'); key2: rooted tree, value_ type value							IN: defined; OUT: unchanged
	my	$aref_rooted_qtrees						=	$_[8] ;		# list of rooted quartet trees 																								IN: defined; OUT: unchanged
	my	$href_NaNpNc_inf_of_type_of_tree		=	$_[9] ;		# key1: tree dependent info types (.e.g 'Na', 'Np'...); key2: rooted tree, value_ type value								IN: defined; OUT: unchanged
	my	$href_cf_of_type_of_qt					=	$_[10] ;	# key1: info type (e.g. ('cf_obs'); key2: rooted quartet-tree; value: correction factor for site patterns of actual tree	IN: defined; OUT: unchanged
	my	$href_finalsupport_of_type_of_qt		=	$_[11] ;	# key1: support type (e.g. 'nap', 'napv'...); key2: rooted quartet tree; value: score										IN: defined; OUT: unchanged 
	my	$href_value_of_scr_par					=	$_[12] ;	# key1: sealion user available parameter options; value: optional parameter setting											IN: defined; OUT: unchanged
	
	#print "::&print_tsv_Q3\n";
	
	######################################################################
	# called from subroutine(s)...
	# &evaluation_quartet_pattern_main
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	## print header and data to tsv table Q3                            ##
	######################################################################
	# table includes pattern info about each best supported rooted tree of each quartet
	# used for plot:
	#	pbox_sNsExpObs						ppoint_rNcNaRed_CFobsCFexp	
	#	pbox_tNsExpObs						ppoint_rNcvNavOrgRed_CFvobsCFvexp
	#	pbox_sPsExpObsFor					ppoint_Nap_rCFvExpObs
	#	pbox_rsNstNsExpObs					ppoint_Nap_rCFexpObs
	#	pbox_dsNsExpsNsObs					ppoint_NapN2_rCFobsExp
	#	ppoint_rsNstNs_nNaNc				ppoint_NapN2_Seed
	#	ppoint_sNs_tNs						ppoint_NaprN2_Seed
	#	ppoint_CF_NaNc						ppoint_NapvN2_Seed
	#	ppoint_CFexp_CFobs					ppoint_rNcNa_CFobsCFexp
	#	ppoint_rNcNaOrgRed_CFobsCFexp		ppoint_rNcNaOrgRisk1_Seed
	#	pbox_tPsExpObs						ppoint_rNcNaOrgRisk2_Seed
	#	pbox_sPsExpObs						ppoint_rNcrNarOrgRisk1_Seed
	#	pbox_drsNstNsExpObs					ppoint_rNcrNarOrgRisk2_Seed
	#	pbox_dsPsExpsPsObs					ppoint_rNcvNavOrgRisk1_Seed
	#	pbox_tPvExpObs						ppoint_rNcvNavOrgRisk2_Seed
	#	pbox_dsPsvExpsPsObs
	my	 	$tsv_4clan_metafileQ3 = $href_printpath->{tsv}."/Q3_SeaLion_metaData_species-quartet_split-analysis_results_clade-quartet_".$$sref_cl4_code.".tsv" ;
	open 	OUTtsvQ3, ">>$tsv_4clan_metafileQ3" || die "Can't write tsv file ", $tsv_4clan_metafileQ3, ": $!" ;
	
	########
	# print header only once to each clanfile corresponding tsv table file
	unless ( $href_hoh_data_of_tsvfile_of_prcode->{metafileQ3}{$tsv_4clan_metafileQ3} ){
		
		#print "::printing table header to ", $tsv_4clan_metafileQ3 ,"\n" ;
		
		########
		# defined header line of tsv outfile
		print 	OUTtsvQ3		"Clade_Quartet"											, # Clade4
								"\tQuartet_Number"										, # FN
								"\tTree"												, # Topo
								"\trNcNa_original"										; # rNcNa_org_v
		
		if ( $href_value_of_scr_par->{seq_algorithm}{nap} == 1 ){
			
			print	OUTtsvQ3	"\ttree_rank_nap_c"										,
								"\trCF_ExpObs_nap_v"									,
								#
								#"\trCF_ExpObs_nap_c"									;
		}
		
		print 	OUTtsvQ3	"\n";
		########
		
		########
		# count metafile header +1
		$href_hoh_data_of_tsvfile_of_prcode->{metafileQ3}{$tsv_4clan_metafileQ3} = $$sref_cl4_code ;
		########
	}
	########
	
	########
	# print table data
	#print "::printing table data to ", $tsv_4clan_metafileQ3 ,"\n" ;
	
	for my $tree ( @$aref_rooted_qtrees ){
		
		##### LINE1+LINE2 (same content)
		print 	OUTtsvQ3	$$sref_cl4_code																,	# Clan4
							"\t",$$sref_qn																,	# FN
							"\t",$tree																	,	# Topo
							"\t",$href_NaNpNc_inf_of_type_of_tree->{rNcNa_org}{$tree}					;	# rNcNa_org_v
		#####
		
		######################################################################
		# NAP quartet calculations
		if ( $href_value_of_scr_par->{seq_algorithm}{nap} == 1 ){
			
			##### assign method nap identified quartet-tree ranking (1st, 2nd, or 3rd)
			my $tree_ranking = &assign_tree_order( \%$href_finalsupport_of_type_of_qt, \"nap", $tree )	;
			unless ( $tree_ranking eq 'best' ) { $tree_ranking = 'minor' } 
			#####
			
			##### LINE1+LINE2 (same content)
			print 	OUTtsvQ3	"\t",$tree_ranking														,	# tree_rank_nap_c
								"\t",$href_cf_of_type_of_qt->{rcf_ExpObs}{$tree}						;	# rCF_ExpObs_nap_v
			#####
		}
		######################################################################
					
		print 	OUTtsvQ3		"\n";
		
	}			
	########
	
	########
	# close tsv filehandle
	close	OUTtsvQ3 ; 
	#####################################################################################
	#####################################################################################
}

sub print_tsv_Q3b{
	
	my	$href_printpath							=	$_[0] ;	# key1: R grafic unique printcode (e.g. 'metafileQ1'); key2: *.tsv filename; value: *tsv filename data addings 																IN: defined; OUT: unchanged
	my	$href_hoh_data_of_tsvfile_of_prcode		=	$_[1] ;	# key: output-type (e.g. 'pdf'); value: pathdir																																IN: defined; OUT: unchanged
	my	$sref_clan4								=	$_[2] ;	# Clan4 of given quartet combination																																		IN: defined; OUT: unchanged
	my	$href_fc_of_4clan						=	$_[3] ;	# key: quartet numbers of current 4clan combination																															IN: defined; OUT: unchanged
	my	$href_hoh_qrej_of_4clan_of_ftype_of_fn	=	$_[4] ;	# key1: 4clan code; key2: type of supportfilter ('nap_risk1', 'napr_risk1', or 'napv_risk1'); key3: actual sealion quartet number; value: filter corresponding score-value	IN: defined; OUT: unchanged
	my	$href_seen_filter_method				=	$_[5] ;	# key1: filter method approaches; value: seen +1																															IN: defined; OUT: unchanged
	
	#print "::&print_tsv_Q3b\n";
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	my	$tsv_4clan_metafileQ3 = $href_printpath->{tsv}."/Q3_SeaLion_metaData_species-quartet_split-analysis_results_clade-quartet_".$$sref_clan4.".tsv" ;
	
	my	@tablelines		 = &readIN_file( \$tsv_4clan_metafileQ3 ) ;
	my	@ftypes			 = sort keys %{$href_seen_filter_method} ;
	my	$line_header	 = join "\t", @ftypes ;
		$tablelines[0]	.= "\t".$line_header ;
	
	for my $linenumber ( 1 .. @tablelines-1 ){
		
		my @prts = split "\t", $tablelines[$linenumber] ;
		my $fn	 = $prts[1] ;
		
		for my $ftype ( @ftypes ){
			
			if 	 ( $href_hoh_qrej_of_4clan_of_ftype_of_fn->{$$sref_clan4}{$ftype}{$fn} ){ $tablelines[$linenumber] .= "\trejected" }
			else { $tablelines[$linenumber] .= "\tremain" }
		}
	}
	
	my	$table_complete = join "\n", @tablelines ;
	
	open  OUTtsvQ3, ">$tsv_4clan_metafileQ3" || die "Can't write tsv file ", $tsv_4clan_metafileQ3, ": $!" ;
	print OUTtsvQ3	$table_complete."\n" ;
	close OUTtsvQ3;
}

sub print_tsv_Q5{
	
	my	$href_printpath								= $_[0] ;	# key1: R grafic unique printcode (e.g. 'metafileQ1'); key2: *.tsv filename; value: *tsv filename data addings 				IN: defined; OUT: unchanged
	my	$href_hoh_data_of_tsvfile_of_prcode			= $_[1] ;	# key: output-type (e.g. 'pdf'); value: pathdir																				IN: defined; OUT: unchanged
	my	$href_qn_of_qfiltered_otu_of_4clan_of_ft	= $_[2] ;	# key1: 4seq filtermethod (e.g. nap_risk1); key2: taxon; key3: 4clan|total; value: list of rejected qn's					IN: defined; OUT: unchanged
	my	$href_hoh_qremaining_of_otu_of_4clan_of_ft	= $_[3] ;	# key1: 4seq filtermethod (e.g. nap_risk1); key2: taxon; key3: 4clan|total; value: N remaining quartets						IN: defined; OUT: unchanged
	my	$href_clan_of_seq							= $_[4] ;	# key1: sequence-name; value: clan definition																				IN: defined; OUT: unchanged
	my	$href_seen_clan4							= $_[5] ;	# key1: clan4; value: counter-number																						IN: defined; OUT: unchanged
	my	$href_seen_filt_method						= $_[6] ;	# key1: filter method incl OTU clipping approaches; value: seen +1															IN: defined; OUT: unchanged
	#print "::&print_tsv_Q5\n";
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	#####################################################################################
	###### print support in respect of all three rooted trees to table TSV Q5       #####
	#####################################################################################
	#for	my $ftype ( sort keys %{$href_qn_of_qfiltered_otu_of_4clan_of_ft} ){ 
	for	my $ftype ( sort keys %$href_seen_filt_method ){ 
		
		########
		my		$filter				 	=	&recode_filter_printout( \$ftype )	;
		my		@cm_fm					=	split "_", $ftype ;
		my	 	$tsv_4clan_metafileQ5	=	$href_printpath->{tsv}."/Q5_SeaLion_metaData_clade-quartet_species-specific_quartet-particpations_for_".$filter."_".$cm_fm[0].".tsv" ;
		open 	OUTtsvQ5, ">>$tsv_4clan_metafileQ5" || die "Cannot write tsv file ", $tsv_4clan_metafileQ5, ": $!" ;
		########
		
		########
		#print "::printing table header to ", $tsv_4clan_metafileQ5 ,"\n" ;
			
		########
		# defined header line of tsv outfile
		print 	OUTtsvQ5	"Species"						,
							"\tClade"						,
							"\tScore"						,
							"\tFilter"						,
							"\tClade_Quartet"				,
							"\tN_Rejected_Quartets"			,
							"\tN_Rejected_Quartets_Total"	,
							"\tN_Retained_Quartets"			,
							"\tN_Retained_Quartets_Total"	,
							"\n"				;
		########
		
		########
		# count metafile header +1
		$href_hoh_data_of_tsvfile_of_prcode->{metafileQ5}{$tsv_4clan_metafileQ5} = $ftype ;
		########
		
		########
		# print table data
		for my $t ( sort keys %$href_clan_of_seq ){
				
			############
			# total number of rejected quartets
			my	$total_rej = 0 ;
			my	$total_rem = 0 ;
			
			for my $clan4 ( sort keys %$href_seen_clan4 ){
			#for my $clan4 ( sort keys %{$href_qn_of_qfiltered_otu_of_4clan_of_ft->{$ftype}{$t}}){
				
				if ( $href_qn_of_qfiltered_otu_of_4clan_of_ft->{$ftype}{$t}{$clan4}   ){ $total_rej += @{$href_qn_of_qfiltered_otu_of_4clan_of_ft->{$ftype}{$t}{$clan4}} }
				if ( $href_hoh_qremaining_of_otu_of_4clan_of_ft->{$ftype}{$t}{$clan4} ){ $total_rem += $href_hoh_qremaining_of_otu_of_4clan_of_ft->{$ftype}{$t}{$clan4}  }
			}
			$total_rej = $total_rej / 3 ; # counted 3times per quartet (for each tree once), thus division by three
			$total_rem = $total_rem / 3 ; # counted 3times per quartet (for each tree once), thus division by three
			############$href_seen_clan4
				
			############
			for my $clan4 ( sort keys %$href_seen_clan4 ){
				
				my $n_rej_qn = 0 ;
				my $n_rem_qn = 0 ;
				
				if ( $href_qn_of_qfiltered_otu_of_4clan_of_ft->{$ftype}{$t}{$clan4} ){
					
					############
					# 4clan specific rejections
					$n_rej_qn = @{$href_qn_of_qfiltered_otu_of_4clan_of_ft->{$ftype}{$t}{$clan4}} / 3 ; # counted 3times per quartet (for each tree once), thus division by three
					############
				}
				
				if ( $href_hoh_qremaining_of_otu_of_4clan_of_ft->{$ftype}{$t}{$clan4} ){
					
					############
					# 4clan specific remains
					$n_rem_qn = $href_hoh_qremaining_of_otu_of_4clan_of_ft->{$ftype}{$t}{$clan4} / 3 ; # counted 3times per quartet (for each tree once), thus division by three
					############
				}
				
				print	OUTtsvQ5	$t								,	# OTU
									"\t",$href_clan_of_seq->{$t}	,	# Clan
									"\t",$cm_fm[0]					,	# Score
									"\t",$filter					,	# Filter
									"\t",$clan4						,	# Clan4
									"\t",$n_rej_qn					,	# rej4clan
									"\t",$total_rej					,	# rejTotal
									"\t",$n_rem_qn					,	# remain4clan
									"\t",$total_rem					,	# remainTotal
									"\n"							;
			}
			########
		}
		########
		
		########
		# close tsv filehandle
		close	OUTtsvQ5 ; #exit;
		########
	}
	#####################################################################################
	#####################################################################################
}

sub print_tsv_Q6{
	
	my	$href_printpath									= $_[0] ;	# key1: R grafic unique printcode (e.g. 'metafileQ1'); key2: *.tsv filename; value: *tsv filename data addings 																		IN: defined; OUT: unchanged
	my	$href_hoh_data_of_tsvfile_of_prcode				= $_[1] ;	# key: output-type (e.g. 'pdf'); value: pathdir																																		IN: defined; OUT: changed
	my	$href_hoh_N_qrej_of_clan4_of_ft					= $_[2] ;	# key1: 4seq filtermethod (e.g. nap_risk1); key2: 4clan combination code; value:  N rejected quartets																				IN: defined; OUT: unchanged
	my	$href_seen_nq_of_clan4							= $_[3] ;	# key1: clan4; value: N analyseed quartets*3 (counted three times, each time for each quartet tree)																					IN: defined; OUT: unchanged
	my	$href_seen_filter_method						= $_[4] ;	# key1: filter method incl OTU clipping approaches; value: seen +1																													IN: defined; OUT: unchanged
	
	#print "::&print_tsv_Q6\n";
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	#####################################################################################
	###### print N reject and remain quartets for each quartet-filter to table TSV Q6 ###
	#####################################################################################
	for my $cm_fm ( sort keys %$href_seen_filter_method ){
		
		my	$filter	=	&recode_filter_printout( \$cm_fm )	;
		my	@prts 	=	split "_", $cm_fm ;
		my	$cm		=	shift @prts ;
		my	$fm		=	join "_", @prts ;
		
		########
		my	 	$tsv_4clan_metafileQ6 = $href_printpath->{tsv}."/Q6_SeaLion_metaData_clade-quartet_species-quartet_number-rejections_".$cm.".tsv" ;
		open 	OUTtsvQ6, ">>$tsv_4clan_metafileQ6" || die "Cannot write tsv file ", $tsv_4clan_metafileQ6, ": $!" ;
		########
		
		########
		# print header only once to each clanfile corresponding tsv table file
		unless ( $href_hoh_data_of_tsvfile_of_prcode->{metafileQ6}{$tsv_4clan_metafileQ6} ){
			
			#print "::printing table header to ", $tsv_4clan_metafileQ6 ,"\n" ;
			
			########
			# defined header line of tsv outfile
			print 	OUTtsvQ6	"Score"					,
								"\tFilter"				,
								"\tClade_Quartet"		,
								"\tSpecies_Quartet"		,
								"\tNumber"				,
								"\n"					;
			########
			
			########
			# count metafile header +1
			$href_hoh_data_of_tsvfile_of_prcode->{metafileQ6}{$tsv_4clan_metafileQ6} = $cm ;
			########
		}
		########
		
		########
		# print table data
		for my $clan4 ( sort keys %$href_seen_nq_of_clan4 ){
			
			my $qtotal = $href_seen_nq_of_clan4->{$clan4} ;
			
			
			
			my $qrejec = 0 ;
			if ( $href_hoh_N_qrej_of_clan4_of_ft->{$cm_fm}{$clan4} ){ $qrejec = $href_hoh_N_qrej_of_clan4_of_ft->{$cm_fm}{$clan4} / 3 }
			
			print	OUTtsvQ6	$cm				,	# Score
								"\t",$filter	,	# Filter
								"\t",$clan4		,	# Clan4
								"\t","initially",	# Nq_c
								"\t",$qtotal	,	# Nq_v
								"\n",
								$cm				,	# Score
								"\t",$filter	,	# Filter
								"\t",$clan4		,	# Clan4
								"\t","rejected"	,	# Nq_c
								"\t",$qrejec	,	# Nq_v
								"\n"			;
		}
		########
		
		########
		# close tsv filehandle
		close	OUTtsvQ6 ; #exit;
		########
	}
	#####################################################################################
}

sub print_tsv_Q7{
	
	my	$href_printpath									= $_[0] ;	# key1: R grafic unique printcode (e.g. 'metafileQ1'); key2: *.tsv filename; value: *tsv filename data addings 																		IN: defined; OUT: unchanged
	my	$href_hoh_data_of_tsvfile_of_prcode				= $_[1] ;	# key: output-type (e.g. 'pdf'); value: pathdir																																		IN: defined; OUT: changed
	my	$href_seqL_of_fc								= $_[2] ;	# key1: quartet-number ($fc); value: quartet related sequence length (corrected)																									IN: defined; OUT: unchanged
	my	$href_bestTree_of_fc_of_cm						= $_[3] ;	# key1: quartet-number ($fc), key2: calculation method (nap, napv, napr); value: quartet related best supported tree																IN: defined; OUT: unchanged
	my	$href_seen_remain_of_cmfm_of_fc					= $_[4] ;	# key1: calculation method + filter method (nap_risk1, napv_dist2, napr_dist3); key2: quartet-number ($fc); value: seen +1															IN: defined; OUT: unchanged
	my	$href_taxa4_of_fc								= $_[5] ;	# key1: quartet-number ($fc); value : quartet assigned taxa (alphabetically ordered, joined by ':'																					IN: defined; OUT: unchanged
	my	$href_hol_fc_of_clan4							= $_[6] ;	# key1: 4clan combination code; value: list of assigned quartet-numbers ($fc)																										IN: defined; OUT: unchanged
	my	$href_algorithm_def_methods						= $_[7] ;	# key1: methode type; key2: method code; value: integer																																IN: defined; OUT: unchanged
	my	$href_seen_filter_method						= $_[8] ;	# key1: filter method incl OTU clipping approaches; value: seen +1																													IN: defined; OUT: unchanged
	
	#print "::&print_tsv_Q7\n";
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	## print header and data to tsv table Q7                            ##
	######################################################################
	for my $cm_org ( sort keys %{$href_algorithm_def_methods->{original}} ){
		
		for my $clan4 ( sort keys %$href_hol_fc_of_clan4 ){
			
			########
			my	 	$tsv_4clan_metafileQ7 = $href_printpath->{tsv}."/Q7_SeaLion_metaData_species-quartet_N-sites-analysed_best-tree-support_".$clan4."_".$cm_org.".tsv" ;
			open 	OUTtsvQ7, ">>$tsv_4clan_metafileQ7" || die "Cannot write tsv file ", $tsv_4clan_metafileQ7, ": $!" ;
			########
			
			########
			# print header only once to each clanfile corresponding tsv table file
			unless ( $href_hoh_data_of_tsvfile_of_prcode->{metafileQ7}{$tsv_4clan_metafileQ7} ){
				
				#print "::printing table header to ", $tsv_4clan_metafileQ7 ,"\n" ;
				
				########
				# defined header line of tsv outfile
				my	$unfiltered	=	&recode_filter_printout( \$cm_org )	;
				print 	OUTtsvQ7	"Quartet_Number"	  	,
									"\tSpecies_Quartet"	  	,
									"\tClade_Quartet"	  	,
									"\tN_informative_sites" ,
									"\tBest_Tree" 		  	,
									"\t", $unfiltered		;
				
				my $methods = $clan4."::".$cm_org."::".$unfiltered."\t".$cm_org ;
				
				for my $fm ( sort keys %$href_seen_filter_method ){	#print "\nfm:", $fm, "\n";
					
					if ( $fm =~ /^${cm_org}_/ ){ 
						
						my	@filter_prts 	=	split "_strict_", $fm ;
						my	$qfilter_type	=	&recode_filter_printout( \$filter_prts[0] )	;
						
						if ( $fm =~ /strict/ ){
							
							my	$tfilter_type	=	&recode_filter_printout( \$fm )	;
								$qfilter_type 	=	$qfilter_type."_".$tfilter_type ; 
						}
						
						print 	OUTtsvQ7 "\t", $qfilter_type ;
						
						$methods .= "::".$qfilter_type."\t".$fm ;
						#print "\nQ7:\t", $methods, "\n" ;
					}
				}
				
				print 	OUTtsvQ7 "\n" ;
				########
				
				########
				# count metafile header +1
				$href_hoh_data_of_tsvfile_of_prcode->{metafileQ7}{$tsv_4clan_metafileQ7} = $methods ;
				########
			}
			########
			
			########
			# print table data
			my $aref_qnumbers = \@{$href_hol_fc_of_clan4->{$clan4}} ;
			for my $qn ( @$aref_qnumbers ){
				
				my @prts	 = split ":", $href_bestTree_of_fc_of_cm->{$qn}{$cm_org} ;
				my $besttree = $prts[0] ;
				my $support  = $prts[1] ;
				
				print 	OUTtsvQ7		  $qn 						,
									"\t", $href_taxa4_of_fc->{$qn}	,
									"\t", $clan4					,
									"\t", $href_seqL_of_fc->{$qn}	,
									"\t", $besttree					,
									"\t", $support					;
				
				for my $fm ( sort keys %$href_seen_filter_method ){
					
					if ( $fm =~ /^${cm_org}_/ ){
						
						if  ( $href_seen_remain_of_cmfm_of_fc->{$fm}{$qn}	){ print OUTtsvQ7 "\t"  , $support	}
						else												 { print OUTtsvQ7 "\tNA"			}
					}
				}
				
				print 	OUTtsvQ7 "\n"
			}
			########
			
			########
			# close tsv filehandle
			close	OUTtsvQ7 ; #exit;
			########
		}
	}
	#####################################################################################
}

sub print_tsv_T1{
	
	my	$href_prtpath									= $_[0] ;	# key1: R grafic unique printcode (e.g. 'metafileQ1'); key2: *.tsv filename; value: *tsv filename data addings 							IN: defined; OUT: unchanged
	my	$href_hoh_data_of_tsvfile_of_prtcode			= $_[1] ;	# key: output-type (e.g. 'pdf'); value: pathdir																							IN: defined; OUT: unchanged
	my	$href_supp_of_4clan_of_tax_of_type_of_tree		= $_[2]	;	# key1: 4clan code; key2: 4taxon-string; key3: mean|median; key4: quartet-tree; value: support value 									IN: defined; OUT: unchanged
	my	$sref_cmfm										= $_[3] ;	# type of calculation (nap, napv...) and or filter (rsik_1, risk1_dist1...) method														IN: defined; OUT: unchanged 
	my	$sref_clan4										= $_[4] ;	# clan4 code																															IN: defined; OUT: unchanged 
	my	$href_clan_of_sequence							= $_[5]	;	# key1: taxon; value: assigned clan																										IN: defined; OUT: unchanged 
	my	$href_alg_def_methods							= $_[6] ;	# key1: 'ftaxon'; value: list of defined taxon filter																					IN: defined; OUT: unchanged	
	my	$href_seen_alltax_of_4clan						= $_[7] ;	# key1: 4clan code; key2: each 4clan assigned taxon; value: counter number																IN: defined; OUT: unchanged 
	
	#print "::&print_tsv_T1\n";
	
	######################################################################
	# called from subroutine(s)...
	# &evaluation_otu_support_main
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	# print data to table TSV T1                                       ###
	######################################################################
	
	########
	# open filehandle
	#my	 	$tsv_4clan_metafileT1 = $href_prtpath->{tsv}."/T1_Penguin_taxon_metaDataTable_singleTreeSupport_".$$sref_cmfm."_".$stype.".tsv" ;
	my		@cm_fm 			=	split "_", $$sref_cmfm ;
	my		$qfilter_type	=	&recode_filter_printout( \$$sref_cmfm )	;
	
	my	 	$tsv_4clan_metafileT1 = $href_prtpath->{tsv}."/T1_SeaLion_metaData_clade-quartet_species-specific_tree-support_contributions_".$qfilter_type."_".$cm_fm[0].".tsv" ;
	open 	OUTtsvT1, ">>$tsv_4clan_metafileT1" || die "Can't write tsv file ", $tsv_4clan_metafileT1, ": $!" ;
	########
	
	########
	# print header only once to each clanfile corresponding tsv table file
	unless ( $href_hoh_data_of_tsvfile_of_prtcode->{metafileT1}{$tsv_4clan_metafileT1} ){
		
		#print "::printing table header to ", $tsv_4clan_metafileT1 ,"\n" ;
		
		########
		# defined header line of tsv outfile
		print	OUTtsvT1	"Score\t", "Filter\t", "Species\t", "Clade\t", "Clade_Quartet\t", "Average\t", "QT1_topology\t", "QT1\t", "QT2_topology\t", "QT2\t", "QT3_topology\t", "QT3" ;
		########
		
		########
		# count metafile header +1
		$href_hoh_data_of_tsvfile_of_prtcode->{metafileT1}{$tsv_4clan_metafileT1} = $$sref_cmfm ;
		########
	}
	########
	
	########
	# print data line to 'pbarplot_Lanalysed_quartets'
	for my $taxon ( sort keys %{$href_seen_alltax_of_4clan->{$$sref_clan4}} ){
		
		for my $stype ( @{$href_alg_def_methods->{qsuppo}} ){
			
			########
			# print general data line info
			print	OUTtsvT1	"\n", $cm_fm[0]							,	# Score
								"\t", $qfilter_type						,	# Filter
								"\t", $taxon							,	# OTU
								"\t", $href_clan_of_sequence->{$taxon}	, 	# Clan
								"\t", $$sref_clan4						,	# Clan4
								"\t", $stype							;	# SupportType
			
			########
			# print data line if taxon support exists
			if ( $href_supp_of_4clan_of_tax_of_type_of_tree->{$$sref_clan4}{$taxon} ){	
				
				for my $tree ( sort keys %{$href_supp_of_4clan_of_tax_of_type_of_tree->{$$sref_clan4}{$taxon}{$stype}} ){
					
					print	OUTtsvT1	"\t", $tree																					,	
										"\t", $href_supp_of_4clan_of_tax_of_type_of_tree->{$$sref_clan4}{$taxon}{$stype}{$tree} 	;	
				}
			}
			########
			
			########
			# print data line if taxon support does not exist
			else{ print	OUTtsvT1	"\t", "NA", "\t", "NA", "\t", "NA", "\t", "NA", "\t", "NA", "\t", "NA" }
			########
		}
	}
	########
	
	########
	# close tsv filehandle
	close	OUTtsvT1 ;
	######################################################################
	######################################################################
}

sub print_tsv_T2{
	
	my	$href_prtpath									= $_[0] ;	# key1: R grafic unique printcode (e.g. 'metafileQ1'); key2: *.tsv filename; value: *tsv filename data addings 							IN: defined; OUT: unchanged
	my	$href_hoh_data_of_tsvfile_of_prtcode			= $_[1] ;	# key: output-type (e.g. 'pdf'); value: pathdir																							IN: defined; OUT: unchanged
	my	$href_supp_of_4clan_of_tax_of_type_of_tree		= $_[2]	;	# key1: 4clan code; key2: 4taxon-string; key3: mean|median; key4: quartet-tree; value: support value 									IN: defined; OUT: unchanged
	my	$sref_cmfm										= $_[3] ;	# type of calculation (nap, napv...) and or filter (rsik_1, risk1_dist1...) method														IN: defined; OUT: unchanged 
	my	$sref_clan4										= $_[4] ;	# clan4 code																															IN: defined; OUT: unchanged 
	my	$href_clan_of_sequence							= $_[5]	;	# key1: taxon; value: assigned clan																										IN: defined; OUT: unchanged 
	my	$href_alg_def_methods							= $_[6] ;	# key1: 'ftaxon'; value: list of defined taxon filter																					IN: defined; OUT: unchanged	
	my	$href_seen_alltax_of_4clan						= $_[7] ;	# key1: 4clan code; key2: each 4clan assigned taxon; value: counter number																IN: defined; OUT: unchanged 
	
	#print "::&print_tsv_T2\n";
	
	######################################################################
	# called from subroutine(s)...
	# &evaluation_otu_support_main
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	# print data to table TSV T2                                       ###
	######################################################################
	for my $taxon ( sort keys %{$href_seen_alltax_of_4clan->{$$sref_clan4}} ){
		
		for my $stype ( @{$href_alg_def_methods->{qsuppo}} ){
			
			########
			# open filehandle	
			my		@cm_fm 			=	split "_", $$sref_cmfm ;
			my		$qfilter_type	=	&recode_filter_printout( \$$sref_cmfm )	;
	
			my	 	$tsv_4clan_metafileT2 = $href_prtpath->{tsv}."/T2_SeaLion_metaData_clade-quartet_species-specific_support_contributions_".$qfilter_type."_".$stype."_".$cm_fm[0].".tsv" ;
			open 	OUTtsvT2, ">>$tsv_4clan_metafileT2" || die "Can't write tsv file ", $tsv_4clan_metafileT2, ": $!" ;
			########
			
			########
			# print header only once to each clanfile corresponding tsv table file
			unless ( $href_hoh_data_of_tsvfile_of_prtcode->{metafileT2}{$tsv_4clan_metafileT2} ){
				
				#print "::printing table header to ", $tsv_4clan_metafileT2 ,"\n" ;
				
				########
				# defined header line of tsv outfile
				#print	OUTtsvT2	"Method\t", "OTU\t", "Number\t", "Clan\t", "Clan4\t", "SupportType\t", "Qtree_code\t", "QT\t", "QT_support";
				print	OUTtsvT2	"Score\t", "Filter\t", "Species\t", "Clade\t", "Clade_Quartet\t", "Average\t", "Tree_Code\t", "Tree\t", "Support";
				########
				
				########
				# count metafile header +1
				$href_hoh_data_of_tsvfile_of_prtcode->{metafileT2}{$tsv_4clan_metafileT2} = $$sref_cmfm.":".$stype ;
				########
			}
			########
			
			########
			# print data line if taxon support exists
			if ( $href_supp_of_4clan_of_tax_of_type_of_tree->{$$sref_clan4}{$taxon} ){ #print "\ntaxon ", $taxon, "\t4clan ", $$sref_clan4, "\n";
				
				my $qt_counter = 0 ;
				for my $tree ( sort keys %{$href_supp_of_4clan_of_tax_of_type_of_tree->{$$sref_clan4}{$taxon}{$stype}} ){
					
					$qt_counter++ ;
					
					print	OUTtsvT2	"\n", $cm_fm[0]							,	# Score
										"\t", $qfilter_type						,	# Filter
										"\t", $taxon							,	# Species
										"\t", $href_clan_of_sequence->{$taxon}	, 	# Clade
										"\t", $$sref_clan4						,	# Clade_Quartet
										"\t", $stype							,	# Average
										"\t", "QT", $qt_counter					,	# Tree_Code
										"\t", $tree								,	# Tree
										"\t", $href_supp_of_4clan_of_tax_of_type_of_tree->{$$sref_clan4}{$taxon}{$stype}{$tree} ;	# QTx_support (x = 1,2,3) 
				}
			}
			
			########
			# print data line if taxon support does not exist
			else{ 
		
				print	OUTtsvT2	"\n", $cm_fm[0], "\t", $qfilter_type, "\t", $taxon, "\t", $href_clan_of_sequence->{$taxon}, "\t", $$sref_clan4, "\t", $stype , "\t", "QT1", "\t", "NA", "\t", "NA",
									"\n", $cm_fm[0], "\t", $qfilter_type, "\t", $taxon, "\t", $href_clan_of_sequence->{$taxon}, "\t", $$sref_clan4, "\t", $stype , "\t", "QT2", "\t", "NA", "\t", "NA",
									"\n", $cm_fm[0], "\t", $qfilter_type, "\t", $taxon, "\t", $href_clan_of_sequence->{$taxon}, "\t", $$sref_clan4, "\t", $stype , "\t", "QT3", "\t", "NA", "\t", "NA";
			}
			########
			
			########
			# close tsv filehandle
			close	OUTtsvT2 ;
			########
		}
	}#exit;
	######################################################################
	######################################################################
}

sub quartet_building_and_p4_main{
	
	my	$href_value_of_scr_par			= $_[0] ; # key1: sealion user available parameter options; value: optional parameter setting				IN: defined; OUT: unchanged
	my	$href_subfolder_q				= $_[1]	; # key: output-type (e.g. 'pdf'); value: pathdir													IN: defined; OUT: unchanged
	my	$href_hol_taxa_of_clan			= $_[2]	; # key1: subclade name; value: list of subclade assigned taxa									 	IN: defined; OUT: unchanged
	my	$aref_clanquartets				= $_[3]	; # list of generated 4-clan combinations															IN: defined; OUT: unchanged
	my	$href_unred_seq_of_taxa			= $_[4]	; # key1: taxonname; value sequence																	IN: defined; OUT: unchanged
	my	$href_hol_seq_states_of_t		= $_[5] ; # key1: taxon name; key2: state position number; value: state 									IN: defined; OUT: undefined
	my	$href_data_of_infile_property	= $_[6] ; # key1: property; value : property value															IN: defined; OUT: unchanged
	my	$href_clan_of_seq				= $_[7] ; # key1: sequence-name; value: clan definition														IN: defined; OUT: unchanged
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# &generate_taxon_quartets
	# &check_taxon_quartet_number
	# &site_exclusion
	# &build_quartet_alignment
	# &write_p4_script
	# &start_p4
	# &read_in_p4_results
	# &write_p4_result_file
	######################################################################
	
	######################################################################
	# Build and analyse quartets of each 4clan combination               #
	######################################################################
	for my $clan4 ( sort @$aref_clanquartets ){
		
		##############################################
		# Generate all possible quartets between defined subgroup taxa
		my	@taxon_quartets	= ()	;	# list of 4tax quartets above limit sequence length after gap+amb exclusions
		&generate_taxon_quartets(
			
			\%hol_taxa_of_clan					, # key: clan name; value: list of clan assigned taxa											IN: defined; OUT: unchanged
			\$clan4								, # 4-clan string, separated by doublepoints (e.g.: c1:c2:c3:c4)								IN: defined; OUT: unchanged
			\@taxon_quartets					, # list of generated 4-taxon combinations														IN: undefined; OUT: defined
			\%$href_value_of_scr_par			, # key1: sealion user available parameter options; value: optional parameter setting			IN: defined; OUT: unchanged
			\%$href_hol_seq_states_of_t			, # key1: taxon name; key2: state position number; value: state 								IN: defined; OUT: unchanged
			\%$href_data_of_infile_property		, # key1: property; value : property value														IN: defined; OUT: unchanged
		) ;
		##############################################
		
		
		##############################################
		# if N possible taxon quartets > max N of allowed taxon quartets start user query for random selection of quartets
		&check_taxon_quartet_number(
			
			\@taxon_quartets					, # list of generated 4-taxon combinations														IN: defined; OUT: changed
			\%$href_value_of_scr_par			, # key1: sealion user available parameter options; value: optional parameter setting			IN: defined; OUT: unchanged
			\$clan4								, # 4-clan string, separated by doublepoints (e.g.: c1:c2:c3:c4)								IN: defined; OUT: unchanged
		) ;
		##############################################
		
		
		######################################################################
		# MSA and P4 analysis of seq-quartets of each 4clan combination      #
		######################################################################
		my	$qcounter = 0 ;
		for my $tquartet ( @taxon_quartets ){
			
			##############################################
			# site exclusion (gap+amb)
			my (
					%red_seq_of_taxa, # key1: taxon; value: gap&amb reduced sequence
					$seq_length		, # N characters of reduced quartet alignment
			) ;
			
			my @qtaxa = split ":", $tquartet ;
			
			&site_exclusion(
						
				\%$href_hol_seq_states_of_t			, # key1: taxon name; key2: state position number; value: state 							IN: defined; OUT: unchanged
				\%$href_data_of_infile_property		, # key1: property; value : property value													IN: defined; OUT: unchanged
				\@qtaxa								, # list of the 4 taxon names in alphab. order												IN: defined; OUT: unchanged
				\%red_seq_of_taxa					, # key1: taxon; value: gap&amb reduced sequence											IN: undefined; OUT: defined
				\$seq_length						, # common length of reduced quartet-sequences												IN: undefined; OUT: defined
			) ;
			##############################################
			
			
			##############################################
			# build P4 site-excluded alignment infiles
			my (
					$msa_p4, # P4 input alignment
			) ;
			
			&build_quartet_alignment(
				
				\%red_seq_of_taxa					, # key1: taxon; value: gap&amb reduced sequence											IN: undefined; OUT: defined
				\@qtaxa								, # list of actual 4-taxon combinations														IN: defined; OUT: unchanged
				\$msa_p4							, # P4 input alignment																		IN: undefined; OUT: defined
			) ;
			##############################################
			
			
			##############################################
			# write P4 executable scripts 
			my (
					$script_p4, # P4 input commands
					@qtrees_p4, # P4 trees in the same order as defined in P4 input skript
			) ;
			
			&write_p4_script(
				
				\%$href_value_of_scr_par			, # key1: sealion user available parameter options; value: optional parameter setting		IN: defined; OUT: unchanged
				\%$href_data_of_infile_property		, # key1: property; value : property value													IN: defined; OUT: unchanged
				\@qtaxa								, # list of actual 4-taxon combination														IN: defined; OUT: unchanged
				\$msa_p4							, # P4 input alignment																		IN: defined; OUT: unchanged
				\$script_p4							, # P4 input commands																		IN: undefined; OUT: defined
				\@qtrees_p4							, # P4 trees in the same order as defined in P4 input skript								IN: undefined; OUT: defined
			) ;
			##############################################
			
			
			##############################################
			# process P4
			my (
					$outfile_p4, # P4 resultfile
			) ;
			
			$qcounter++ ;
			
			&start_p4(
				
				\$msa_p4							, # P4 input alignment																		IN: defined; OUT: undefined
				\$script_p4							, # P4 input commands																		IN: defined; OUT: undefined
				\$clan4								, # actual 4clan combination																IN: defined; OUT: unchanged
				\$qcounter							, # actual quartet number																	IN: defined; OUT: unchanged
				\$outfile_p4						, # P4 resultfile																			IN: undefined; OUT: defined
				\$tquartet							, # 4 taxon-quartet names in alphab. order													IN: defined; OUT: unchanged
			) ;
			##############################################
			
			
			##############################################
			# read p4 resultfile
			my (
					%hoh_found_N_of_topo_of_pattern    , # key1: quartet-tree; key2: recoded pattern code; value: N observed
					%hoh_expected_N_of_topo_of_pattern , # key1: quartet-tree; key2: recoded pattern code; value: N expected
			) ;
			
			&read_in_p4_results(
				
				\%$href_value_of_scr_par			, # key1: sealion user available parameter options; value: optional parameter setting		IN: defined; OUT: unchanged
				\%$href_data_of_infile_property		, # key1: property; value : property value													IN: defined; OUT: unchanged
				\$outfile_p4						, # P4 resultfile																			IN: defined; OUT: unchanged
				\%hoh_found_N_of_topo_of_pattern    , # key1: quartet-tree; key2: recoded pattern code; value: N observed						IN: undefined; OUT: defined
				\%hoh_expected_N_of_topo_of_pattern , # key1: quartet-tree; key2: recoded pattern code; value: N expected						IN: undefined; OUT: defined
			) ;
			##############################################
			
			
			##############################################
			# print p4 resultfile in shorter format
			&write_spd_file(
				
				\%$href_value_of_scr_par			, # key1: sealion user available parameter options; value: optional parameter setting		IN: defined; OUT: changed
				\%$href_subfolder_q					, # key: output-type (e.g. 'pdf'); value: pathdir											IN: defined; OUT: unchanged
				\%hoh_found_N_of_topo_of_pattern    , # key1: quartet-tree; key2: recoded pattern code; value: N observed						IN: defined; OUT: unchanged
				\%hoh_expected_N_of_topo_of_pattern , # key1: quartet-tree; key2: recoded pattern code; value: N expected						IN: defined; OUT: unchanged
				\$outfile_p4						, # P4 resultfile																			IN: defined; OUT: unchanged
				\@qtrees_p4							, # P4 trees in the same order as defined in P4 input skript								IN: defined; OUT: unchanged
				\$tquartet							, # 4 taxon-quartet names in alphab. order													IN: defined; OUT: unchanged
				\$clan4								, # 4-clan string, separated by doublepoints (e.g.: c1:c2:c3:c4)							IN: defined; OUT: unchanged
				\$qcounter							, # quartet internal number																	IN: defined; OUT: unchanged
				\$seq_length						, # common length of reduced quartet-sequences												IN: defined; OUT: unchanged
				\%$href_clan_of_seq					, # key1: sequence-name; value: clan definition												IN: defined; OUT: unchanged
			) ;
			##############################################
		}
	}
}

sub r_barplot{
	
	my	$sref_intable		=	$_[0]	;	# name of infile table (*.tsv)
	my	$sref_x_theader		=	$_[1]	;	# table header of x-plot categories
	my	$sref_y_theader		=	$_[2]	;	# table header of y-plot values
	my	$sref_group_fill	=	$_[3]	;	# table header (factor) for value separation given different point shapes within each plot
	my	$sref_group_line	=	$_[4]	;	# table header (factor) for value separation given different line colors within each plot
	my	$sref_flag_facet	=	$_[5]	;	# if defined, facet grid is activated in respect of the next two group parameters
	my	$sref_plot_split_v	=	$_[6]	;	# table header plot-splitting vertically (defined as dot "." if splitting not needed)
	my	$sref_plot_split_h	=	$_[7]	;	# table header plot-splitting horizontally (undefined "" if splitting not needed)
	my	$sref_scale_facet	=	$_[8]	;	# default (undef) all the panels have the same scales ('fixed') can be made independent if defined 'free', 'free_x', or 'free_y'
	my	$sref_hline_value	=	$_[9]	;	# highlight horizontal line at axis interception; neglected if undef, for multiple hline use "vhl1,vhl2..."
	my	$sref_hline_type	=	$_[10]	;	# horizontal linetype: solid, longdash, dashed, blank, dotdash, dotted, twodash
	my	$sref_hline_size	=	$_[11]	;	# horizontal linesize; default if undef
	my	$sref_hline_color	=	$_[12]	;	# horizontal line color
	my	$sref_stagged_val	=	$_[13]	;	# 0 -> a stacked barplot is created by default, 1 -> non-stagged barplot with each value of the same x-plot category in a own x-axis position
	my	$sref_do_bar_label	=	$_[14]	;	# 0 -> no labels at barchart, 1 -> label (geom_text) at barchart
	my	$sref_bar_label_c	=	$_[15]	;	# table header of barchart label
	my	$sref_bar_label_vp	=	$_[16]	;	# vertical shift value (integer, float, positive or negative) of barchart label
	my	$sref_bar_label_co	=	$_[17]	;	# barchart label color
	my	$sref_bar_label_si	=	$_[18]	;	# barchart label size 
	my	$sref_ttitle		=	$_[19]	;	# table title
	my	$aref_cbox			=	$_[20]	;	# list of specified box colors, e.g. (E69F00, 56B4E9) or (coral3, coral4, ...), or if undefined using a grey-scaling
	my	$aref_clin			=	$_[21]	;	# list of specified boxline colors, e.g. (E69F00, 56B4E9) or (coral3, coral4, ...)
	my	$sref_xlab			=	$_[22]	;	# label x-axis
	my	$sref_rotate_xlabel	=	$_[23]	;	# if value == 1 -> rotate x-axis labels via 90°, otherwise keep labels horizontally
	my	$sref_ylab			=	$_[24]	;	# labe y-axis
	my	$sref_lab_fill		=	$_[25]	;	# bar fill label
	my	$sref_pdf_width		=	$_[26]	;	# width of pdf graphic, keep undef for def width: 40
	my	$sref_pdf_height	=	$_[27]	;	# height of pdf graphic, keep undef for def height: 15
	my	$sref_filesuffix	=	$_[28]	;	# pdf file suffix to avoid different plots of identic names
	my	$sref_svg_folder	=	$_[29]	;	# subfolder of svg print
	my	$sref_pdf_folder	=	$_[30]	;	# subfolder of pdf print
	my	$aref_pdf_sampled	=	$_[31]	;	# list of pdf prints, subsequently collected via pdfunite
	
	######################################################################
	## print barplot using R                                            ##
	######################################################################
	
	######################################################################
	# called from subroutine(s)...
	# &plotsR_quartet_analyses
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	my	$home_path			=	getcwd() 											;	# define home path
	my	$R					=	Statistics::R->new()								;	# load R packge Statistics
	
	#########################################
	# define individual plot-filename with subfolder path
	my	$filename 			=	$$sref_intable	;
	if	( $$sref_filesuffix ){	$filename  =~ s/.tsv$/_${$sref_filesuffix}/ }
	else					 {	$filename  =~ s/.tsv$// 					}
	
	my	(	$svg_out, $pdf_out 	); 
	my		@path_parts			=	split "/", $filename ;
	if	( $$sref_svg_folder		){	$svg_out	=	$$sref_svg_folder."/".$path_parts[-1].".svg"	} else{	$svg_out = $filename.".svg"	}
	if	( $$sref_pdf_folder		){	$pdf_out	=	$$sref_pdf_folder."/".$path_parts[-1]			} else{	$pdf_out = $filename		}
	
	# store filename with pdf suffix in @pdf_pattern for subsequent pdf unite process
	push						@$aref_pdf_sampled, $pdf_out.".pdf" 			;
	#########################################
	
	#########################################
	# complete color line and fill for each plot following
	# subroutine assigned custom color lists
	my ( $scale_cbox, $scale_clin );
	if ( @$aref_cbox ){	$scale_cbox = join "\", \"", @$aref_cbox }
	if ( @$aref_clin ){	$scale_clin	= join "\", \"", @$aref_clin }
	#########################################
	
	#########################################
	# storing individual command lines of geom_point, geom_smooth, facet_grid
	my	@command_lines ;
	#########################################
	
	######## GGPLOT COMMAND LINE #####################################################
	##################################################################################
	# make ggplot 'aes' command line
	my	$table			=	"pb"	;
	my	$aes_parameter	=	"x=".$$sref_x_theader.", y=".$$sref_y_theader.", fill=".$$sref_group_fill ;
	push @command_lines	, 	"ggplot(data=".$table.",aes(".$aes_parameter."))" ;
	##################################################################################
	##################################################################################
	
	######## SCALE COLOR COMMAND LINE ################################################
	##################################################################################
	# if fill color defined, use fill color. otherwise use grey scale
	if ( @$aref_cbox )	{ push @command_lines, "scale_fill_manual(values=c(\"".$scale_cbox."\"))" }
	else				{ push @command_lines, "scale_fill_grey()" }
	##################################################################################
	##################################################################################
	
	######## GGPLOT COMMAND LINE #####################################################
	##################################################################################
	# make ggplot geom_bar definitions, incl staggered bar definition
	my	$geom_bar_parameter =	"stat=\"identity\"" ;
	if		( @$aref_clin 				){ $geom_bar_parameter	.=	", color=\"".$scale_clin."\"" }
	if		( $$sref_stagged_val == 1	){ $geom_bar_parameter	.=	", position=position_dodge()" }
	
	push @command_lines,	"geom_bar(".$geom_bar_parameter.")" ;
	##################################################################################
	##################################################################################
	
	######## FACET_GRID COMMAND LINE #################################################
	##################################################################################
	# define facet_grid parameter (if faced_grid is defined to be activated)
	if	( $$sref_flag_facet ){
		
		my	@parameter ;
		
		#########################################
		# if defined, define horizontal (h) and vertical (v) facet_grid grouping
		# following definition of $$sref_plot_split_h and $$sref_plot_split_v respectively
		# store respective group separation commands @parameter
		if		( ( $$sref_plot_split_h	) && ( $$sref_plot_split_v) ){ push @parameter, $$sref_plot_split_h." ~ ".$$sref_plot_split_v	}
		elsif	( 	$$sref_plot_split_h	)							 { push @parameter, $$sref_plot_split_h." ~ ."						}
		elsif	( 	$$sref_plot_split_v	)							 { push @parameter, 					". ~ ".$$sref_plot_split_v	}
		#########################################
		
		#########################################
		# change x,y scale command of each panel from fixed (default) to defined parameter ($$sref_scale_facet)
		if	( $$sref_scale_facet =~ /free|free_x|free_y/ ){	push @parameter, "scales=\"".$$sref_scale_facet."\"" }
		#########################################
		
		#########################################
		# complete facet_grid parameter line
		my	$facet_grid_command	=	join "," ,	@parameter;	
		#########################################
		
		#########################################
		# store complete geom_smooth command in @command_lines
		push	@command_lines, "facet_grid(".$facet_grid_command.")" ;
		#########################################
	}
	##################################################################################
	##################################################################################
	
	######## GEOM_HLINE COMMAND LINE #################################################
	##################################################################################
	# highlight a y-axis line at a certain y-value ($$sref_hline_value, float or integer)
	if	( $$sref_hline_value =~ /\d+(\.\d+)?/ ){
		
		my							$geom_hline	 = "yintercept=c(".$$sref_hline_value.")"	;
		if	( $$sref_hline_color ){	$geom_hline	.= ", color=\"".$$sref_hline_color."\""		}
		if	( $$sref_hline_type  ){	$geom_hline	.= ", linetype=\"".$$sref_hline_type."\""	}
		if	( $$sref_hline_size	 ){ $geom_hline	.= ", size=".$$sref_hline_size				}
									$geom_hline	 =	"geom_hline(".$geom_hline.")"			;
		
		push @command_lines, $geom_hline 
	}
	##################################################################################
	##################################################################################
	
	####### GEOM_TEXT COMMAND LINE ###################################################
	##################################################################################
	# defines if specific values are plotted at each bar
	if	( $$sref_do_bar_label == 1 ){
		
		my	$geom_text_parameter =	"geom_text(aes(label=".$$sref_bar_label_c."), vjust=".$$sref_bar_label_vp.", color=\"".$$sref_bar_label_co."\"" ;
		
		if	( $$sref_stagged_val == 1 )	{ $geom_text_parameter .= ", position = position_dodge(0.9), size=".$$sref_bar_label_si.")" }
		else							{ $geom_text_parameter .= ", size=".$$sref_bar_label_si.")" }
		
		push @command_lines, $geom_text_parameter
	}
	##################################################################################
	##################################################################################
	
	######## X,Y LAB COMMAND LINE ####################################################
	##################################################################################
	push @command_lines, "xlab(\"".$$sref_xlab."\") + ylab(\"".$$sref_ylab."\")" ;
	##################################################################################
	##################################################################################
	
	######## LAB TEXT COMMAND LINE ###################################################
	##################################################################################
	# Change legend label of fill groupings if both label and grouping (fill) are defined
	if ( ( $$sref_group_fill ) && ($$sref_lab_fill) ){ push @command_lines, "labs(fill=\"".$$sref_lab_fill."\")"	}
	##################################################################################
	##################################################################################
	
	######## GG TITLE COMMAND LINE ###################################################
	##################################################################################
	push @command_lines, "ggtitle(\"".$$sref_ttitle."\")" ;
	##################################################################################
	##################################################################################
	
	my	@theme_lines ;
	######## Theme Plot Title #### ###################################################
	##################################################################################
	push @theme_lines, "plot.title = element_text(color=\"#993333\", size=28, face=\"bold\")" ;
	##################################################################################
	##################################################################################
	
	######## Theme Legend Title #### #################################################
	##################################################################################
	push @theme_lines, "legend.title = element_text(color=\"black\", size=28, face=\"bold\")" ;
	##################################################################################
	##################################################################################
	
	######## Theme Legend Background #################################################
	##################################################################################
	push @theme_lines, "legend.box.background = element_blank()" ;
	#push @theme_lines, "legend.background = element_rect(fill=\"ivory2\", size=2, linetype=\"solid\")" ;
	##################################################################################
	##################################################################################
	
	######## Theme Legend Text #######################################################
	##################################################################################
	push @theme_lines, "legend.text = element_text(color=\"black\", size=24, face=\"bold\")" ;
	##################################################################################
	##################################################################################
	
	######## Theme Legend Position #######################################################
	##################################################################################
	push @theme_lines, "legend.position = \"bottom\"" ;
	push @theme_lines, "legend.box = \"vertical\"" ;
	##################################################################################
	##################################################################################
	
	######## Theme Axis Title ########################################################
	##################################################################################
	push @theme_lines, "axis.title.x = element_text(color=\"black\", size=24, face=\"bold\")" ;
	push @theme_lines, "axis.title.y = element_text(color=\"black\", size=24, face=\"bold\")" ;
	##################################################################################
	##################################################################################
	
	######## Theme Axis Text #########################################################
	##################################################################################
	if ( $$sref_rotate_xlabel == 1 ){ push @theme_lines, "axis.text.x = element_text(color=\"black\", size=12, angle=90, hjust = 1)"	}
	else							{ push @theme_lines, "axis.text.x = element_text(color=\"black\", size=12)" 			}
	push @theme_lines, "axis.text.y = element_text(color=\"black\", size=12)" ;
	##################################################################################
	##################################################################################
	
	######## Theme Strip Text ########################################################
	##################################################################################
	push @theme_lines, "strip.text.x = element_text(color=\"black\", size=22, face=\"bold\")" ;
	push @theme_lines, "strip.text.y = element_text(color=\"black\", size=22, face=\"bold\")" ;
	##################################################################################
	##################################################################################
	
	#########################################
	# join geom_point, geom_smooth, facet_grid parameter lines by " + "
	my	$ggplot_theme	=	join ", ", 	@theme_lines		;
	#########################################
	
	#########################################
	# join geom_point, geom_smooth, facet_grid parameter lines by " + "
	my	$ggplot_parameter	=	join " + ", 	@command_lines		;
	#########################################
	
	#########################################
	# change pdf output height (default 15)
	my	$pdf_height	=	15 ;
	if ( $$sref_pdf_height =~ /\d+/ ){ $pdf_height = $$sref_pdf_height }
	#########################################
	
	#########################################
	# change pdf/svg output width (default 40)
	my	$pdf_width	=	40 ;
	if ( $$sref_pdf_width =~ /\d+/ ){ $pdf_width = $$sref_pdf_width }
	#########################################

	#########################################
	# Prepare R commands
	my	@R_commands =(
		'library(ggplot2)',
		'library(svglite)',
		'library(reshape)',
		"setwd(\"$home_path\")",
		"$table <- read.table(\"$$sref_intable\", header = TRUE, sep=\"\\t\")",
		"$ggplot_parameter+theme($ggplot_theme)",
		"plotTitle = \"$svg_out\"",
		"ggsave(plotTitle, width=".$pdf_width.", height =".$pdf_height.", limitsize= FALSE)",
		"ggsave(paste(\"$pdf_out\", \".pdf\",sep = \"\"), width = ".$pdf_width.", height = ".$pdf_height.")"
	);
	#########################################
	
	#########################################
	# terminal print of generated command line
	my	$command_line_print	=	join "\n", @R_commands ;
	#print "\n\nR commandline:\n", $command_line_print, "\n\n",
	#########################################
	
	#########################################
	# R execution
	#print"\n", @R_commands,"\n"; exit;
	#print "\nR pattern plot of ", $filename, " (.svg & .pdf)" ;
	print "\n\t\t...", $filename, " (.svg & .pdf)" ;
	$R->run(@R_commands);
	$R->stop();
	##################################################################################
}

sub r_boxplot{
	
	my	$sref_intable		=	$_[0]	;	# name of infile table (*.tsv)
	my	$sref_x_theader		=	$_[1]	;	# table header of x-plot categories
	my	$sref_y_theader		=	$_[2]	;	# table header of y-plot values
	my	$sref_group_fill	=	$_[3]	;	# table header (factor) for value separation using different fill colors within each plot
	my	$sref_group_edge	=	$_[4]	;	# table header (factor) for value separation using different scatter colors within each plot
	my	$sref_group_shape	=	$_[5]	;	# table header (factor) for value separation using different scatter shapes within each plot
	my	$sref_flag_facet	=	$_[6]	;	# if defined, facet grid is activated in respect of the next two group parameters
	my	$sref_plot_split_v	=	$_[7]	;	# table header plot-splitting vertically (defined as dot "." if splitting not needed)
	my	$sref_plot_split_h	=	$_[8]	;	# table header plot-splitting horizontally (undefined "" if splitting not needed)
	my	$sref_scale_facet	=	$_[9]	;	# default (undef) all the panels have the same scales ('fixed') can be made independent if defined 'free', 'free_x', or 'free_y'
	my	$sref_flag_jitter	=	$_[10]	;	# if defined, jitter plot is activated
	my	$sref_hline_value	=	$_[11]	;	# highlight horizontal line at axis interception; neglected if undef, for multiple hline use "vhl1,vhl2..."
	my	$sref_hline_type	=	$_[12]	;	# horizontal linetype: solid, longdash, dashed, blank, dotdash, dotted, twodash
	my	$sref_hline_size	=	$_[13]	;	# horizontal linesize; default if undef
	my	$sref_hline_color	=	$_[14]	;	# horizontal line color
	my	$sref_ttitle		=	$_[15]	;	# table title
	my	$aref_cfill			=	$_[16]	;	# list of specified fill colors, e.g. (E69F00, 56B4E9) or (coral3, coral4, ...)
	my	$aref_ccolor		=	$_[17]	;	# list of specified edge colors
	my	$sref_xlab			=	$_[18]	;	# label x-axis
	my	$sref_ylab			=	$_[19]	;	# label y-axis
	my	$sref_lab_fill		=	$_[20]	;	# fill colour label
	my	$sref_lab_edge		=	$_[21]	;	# edge (outline) colour label
	my	$sref_lab_shape		=	$_[22]	;	# scatter shape label 	
	my	$sref_pdf_height	=	$_[23]	;	# height of pdf graphic, keep undef for def height: 15
	my	$sref_pdf_width		=	$_[24]	;	# width of pdf graphic, keep undef for def width: 40
	my	$sref_filesuffix	=	$_[25]	;	# pdf file suffix to avoid different plots of identic names
	my	$sref_svg_folder	=	$_[26]	;	# subfolder of svg print
	my	$sref_pdf_folder	=	$_[27]	;	# subfolder of pdf print
	my	$aref_pdf_sampled	=	$_[28]	;	# list of pdf prints, subsequently collected via pdfunite
	
	######################################################################
	## print box plot using R                                       ##
	######################################################################
	
	######################################################################
	# called from subroutine(s)...
	# &plotsR_quartet_analyses
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	my	$home_path			=	getcwd() 											;	# define home path
	my	$R					=	Statistics::R->new()								;	# load R packge Statistics
	
	#########################################
	# define individual plot-filename with subfolder path
	my	$filename 			=	$$sref_intable	;
	if	( $$sref_filesuffix ){	$filename  =~ s/.tsv$/_${$sref_filesuffix}/ }
	else					 {	$filename  =~ s/.tsv$// 					}
	
	my	@path_parts			=	split "/", $filename ;
	my	$svg_out			=	$$sref_svg_folder."/".$path_parts[-1].".svg"	;
	my	$pdf_out			=	$$sref_pdf_folder."/".$path_parts[-1]			;
	push						@$aref_pdf_sampled, $pdf_out.".pdf" 						;	# store filename with pdf suffix in @pdf_pattern for subsequent pdf unite process
	#########################################
	
	#########################################
	# complete color line for each plot following
	# subroutine assigned custom color lists
	my	$scale_fill			=	join "\", \"", 	@$aref_cfill 		;
	my	$scale_color		=	join "\", \"", 	@$aref_ccolor 		;
	#########################################
	
	#########################################
	# storing individual command lines of geom_point, geom_smooth, facet_grid
	my	@command_lines ;
	#########################################
	
	######## GGPLOT COMMAND LINE #####################################################
	##################################################################################
	# make ggplot 'aes' command line
	my	$table			=	"bp"	;
	my	$aes_parameter	=	"x=".$$sref_x_theader.", y=".$$sref_y_theader ;
	if	( $$sref_group_fill  ){	$aes_parameter	.= ", fill=".$$sref_group_fill		}
	if	( $$sref_group_edge  ){	$aes_parameter	.= ", colour=".$$sref_group_edge	}
	if	( $$sref_group_shape ){	$aes_parameter	.= ", shape=".$$sref_group_shape	}
	push @command_lines, "ggplot(".$table.",aes(".$aes_parameter."))" ;
	##################################################################################
	##################################################################################
	
	######## GEOM_HLINE COMMAND LINE #################################################
	##################################################################################
	# highlight a y-axis line at a certain y-value ($$sref_hline_value, float or integer)
	if	( $$sref_hline_value =~ /\d+(\.\d+)?/ ){
		
		my							$geom_hline	 = "yintercept=c(".$$sref_hline_value.")"	;
		if	( $$sref_hline_color ){	$geom_hline	.= ", color=\"".$$sref_hline_color."\""		}
		if	( $$sref_hline_type  ){	$geom_hline	.= ", linetype=\"".$$sref_hline_type."\""	}
		if	( $$sref_hline_size	 ){ $geom_hline	.= ", size=".$$sref_hline_size				}
									$geom_hline	 =	"geom_hline(".$geom_hline.")"			;
		
		push @command_lines, $geom_hline 
	}
	##################################################################################
	##################################################################################
		
	######## GEOM_JITTER COMMAND LINE #################################################
	##################################################################################
	# define facet_grid parameter (if faced_grid is defined to be activated)
	if	( $$sref_flag_jitter ){ push	@command_lines, "geom_jitter(position=position_jitter(0.2))" }
	##################################################################################
	##################################################################################
	
	######## GEOM_BOXPLOT COMMAND LINE ###############################################
	##################################################################################
	push @command_lines, "geom_boxplot(color=\"black\", outlier.colour=\"black\", outlier.shape=8,outlier.size=2, position=position_dodge(1))" ;
	##################################################################################
	##################################################################################
	
	######## STAT_SUMMARY COMMAND LINE ###############################################
	##################################################################################
	push @command_lines, "stat_summary(fun.y=median, geom=\"point\", shape=18, color=\"black\", cex=4.0)" ;
	##################################################################################
	##################################################################################
	
	######## SCALE COLOR COMMAND LINE ################################################
	##################################################################################
	push @command_lines, "scale_fill_manual(values=c(\"".$scale_fill."\"))" ;
	push @command_lines, "scale_color_manual(values=c(\"".$scale_color."\"))" ;
	##################################################################################
	##################################################################################
	
	######## GG TITLE COMMAND LINE ###################################################
	##################################################################################
	push @command_lines, "ggtitle(\"".$$sref_ttitle."\")" ;
	##################################################################################
	##################################################################################
	
	######## X,Y LAB COMMAND LINE ####################################################
	##################################################################################
	push @command_lines, "xlab(\"".$$sref_xlab."\") + ylab(\"".$$sref_ylab."\")" ;
	##################################################################################
	##################################################################################
	
	######## FACET_GRID COMMAND LINE #################################################
	##################################################################################
	# define facet_grid parameter (if faced_grid is defined to be activated)
	if	( $$sref_flag_facet ){
		
		my	@parameter ;
		
		#########################################
		# if defined, define horizontal (h) and vertical (v) facet_grid grouping
		# following definition of $$sref_plot_split_h and $$sref_plot_split_v respectively
		# store respective group separation commands @parameter
		if		( ( $$sref_plot_split_h	) && ( $$sref_plot_split_v) ){ push @parameter, $$sref_plot_split_h." ~ ".$$sref_plot_split_v	}
		elsif	( 	$$sref_plot_split_h	)							 { push @parameter, $$sref_plot_split_h." ~ ."						}
		elsif	( 	$$sref_plot_split_v	)							 { push @parameter, 					". ~ ".$$sref_plot_split_v	}
		#########################################
		
		#########################################
		# change x,y scale command of each panel from fixed (default) to defined parameter ($$sref_scale_facet)
		if	( $$sref_scale_facet =~ /free|free_x|free_y/ ){	push @parameter, "scales=\"".$$sref_scale_facet."\"" }
		#########################################
		
		#########################################
		# complete facet_grid parameter line
		my	$facet_grid_command	=	join "," ,	@parameter;	
		#########################################
		
		#########################################
		# store complete geom_smooth command in @command_lines
		push	@command_lines, "facet_grid(".$facet_grid_command.")" ;
		#########################################
	}
	##################################################################################
	##################################################################################
	
	######## LAB LEGEND COMMAND LINE #################################################
	##################################################################################
	# Change legend label of different groupings if both label and grouping (fill, edge, size) are defined
	if ( ( $$sref_group_fill  ) && ($$sref_lab_fill ) ){ push @command_lines, "labs(fill=\"".$$sref_lab_fill."\")"		}
	if ( ( $$sref_group_edge  ) && ($$sref_lab_edge ) ){ push @command_lines, "labs(colour=\"".$$sref_lab_edge."\")"	}
	if ( ( $$sref_group_shape ) && ($$sref_lab_shape) ){ push @command_lines, "labs(shape=\"".$$sref_lab_shape."\")"	}
	##################################################################################
	##################################################################################
	
	#########################################
	# join geom_point, geom_smooth, facet_grid parameter lines by " + "
	my	$ggplot_parameter	=	join " + ", 	@command_lines		;
	#########################################
	
	my	@theme_lines ;
	######## Theme Plot Title #### ###################################################
	##################################################################################
	push @theme_lines, "plot.title = element_text(color=\"#993333\", size=28, face=\"bold\")" ;
	##################################################################################
	##################################################################################
	
	######## Theme Legend Title #### #################################################
	##################################################################################
	push @theme_lines, "legend.title = element_text(color=\"black\", size=28, face=\"bold\")" ;
	##################################################################################
	##################################################################################
	
	######## Theme Legend Background #################################################
	##################################################################################
	push @theme_lines, "legend.box.background = element_blank()" ;
	#push @theme_lines, "legend.background = element_rect(fill=\"ivory2\", size=2, linetype=\"solid\")" ;
	##################################################################################
	##################################################################################
	
	######## Theme Legend Text #######################################################
	##################################################################################
	push @theme_lines, "legend.text = element_text(color=\"black\", size=24, face=\"bold\")" ;
	##################################################################################
	##################################################################################
	
	######## Theme Legend Position #######################################################
	##################################################################################
	push @theme_lines, "legend.position = \"bottom\"" ;
	push @theme_lines, "legend.box = \"vertical\"" ;
	##################################################################################
	##################################################################################
	
	######## Theme Axis Title ########################################################
	##################################################################################
	push @theme_lines, "axis.title.x = element_text(color=\"black\", size=24, face=\"bold\")" ;
	push @theme_lines, "axis.title.y = element_text(color=\"black\", size=24, face=\"bold\")" ;
	##################################################################################
	##################################################################################
	
	######## Theme Axis Text #########################################################
	##################################################################################
	push @theme_lines, "axis.text.x = element_text(color=\"black\", size=22)" ;
	push @theme_lines, "axis.text.y = element_text(color=\"black\", size=22)" ;
	##################################################################################
	##################################################################################
	
	######## Theme Strip Text ########################################################
	##################################################################################
	push @theme_lines, "strip.text.x = element_text(color=\"black\", size=22, face=\"bold\")" ;
	push @theme_lines, "strip.text.y = element_text(color=\"black\", size=22, face=\"bold\")" ;
	##################################################################################
	##################################################################################
	
	#########################################
	# join geom_point, geom_smooth, facet_grid parameter lines by " + "
	my	$ggplot_theme	=	join ", ", 	@theme_lines		;
	#########################################
	
	#########################################
	# change pdf output height (default 15)
	my	$pdf_height	=	15 ;
	if ( $$sref_pdf_height =~ /\d+/ ){ $pdf_height = $$sref_pdf_height }
	#########################################

	#########################################
	# change pdf/svg output width (default 40)
	my	$pdf_width	=	40 ;
	if ( $$sref_pdf_width =~ /\d+/ ){ $pdf_width = $$sref_pdf_width }
	#########################################
	
	#########################################
	# Prepare R commands
	my	@R_commands =(
		# load R needed R libraries
		'library(ggplot2)',
		'library(svglite)',
		'library(reshape)',
		"setwd(\"$home_path\")",
		"$table <- read.table(\"$$sref_intable\", header = TRUE, sep=\"\\t\")",
		"$ggplot_parameter+theme($ggplot_theme)",
		"plotTitle = \"$svg_out\"",
		"ggsave(plotTitle, width=$pdf_width, height =$pdf_height, limitsize= FALSE)",
		"ggsave(paste(\"$pdf_out\", \".pdf\",sep = \"\"), width = $pdf_width, height = $pdf_height)"
	);
	#########################################
	
	#########################################
	# R execution
	#print"\n", @R_commands,"\n"; #exit;
	print "\n\t\t...", $filename, " (.svg & .pdf)" ;
	#print "\nR pattern plot of ", $filename, " (.svg & .pdf)" ;
	$R->run(@R_commands);
	$R->stop();
	##################################################################################
}

sub r_density{
	
	my	$sref_intable		=	$_[0]	;	# name of infile table (*.tsv)
	my	$sref_x_theader		=	$_[1]	;	# table header of x-plot categories
	my	$sref_y_theader		=	$_[2]	;	# table header of y-plot values (keep empty for density plot)
	my	$sref_group_fill	=	$_[3]	;	# table header (factor) for value separation given different fill colors within each plot
	my	$sref_group_edge	=	$_[4]	;	# table header (factor) for value separation given different edge colors within each plot
	my	$sref_flag_facet	=	$_[5]	;	# if defined, facet grid is activated in respect of the next two group parameters
	my	$sref_plot_split_v	=	$_[6]	;	# table header plot-splitting vertically (defined as dot "." if splitting not needed)
	my	$sref_plot_split_h	=	$_[7]	;	# table header plot-splitting horizontally (undefined "" if splitting not needed)
	my	$sref_scale_facet	=	$_[8]	;	# default (undef) all the panels have the same scales ('fixed') can be made independent if defined 'free', 'free_x', or 'free_y'
	my	$sref_ttitle		=	$_[9]	;	# table title
	my	$aref_cfill			=	$_[10]	;	# list of specified fill colors, e.g. (E69F00, 56B4E9) or (coral3, coral4, ...)
	my	$aref_ccolor		=	$_[11]	;	# list of specified edge colors
	my	$sref_xlab			=	$_[12]	;	# label x-axis
	my	$sref_ylab			=	$_[13]	;	# label y-axis
	my	$sref_lab_fill		=	$_[14]	;	# fill colour label
	my	$sref_lab_edge		=	$_[15]	;	# edge (outline) colour label 	
	my	$sref_pdf_height	=	$_[16]	;	# height of pdf graphic, keep undef for def height: 15
	my	$sref_pdf_width		=	$_[17]	;	# width of pdf graphic, keep undef for def width: 40
	my	$sref_filesuffix	=	$_[18]	;	# pdf file suffix to avoid different plots of identic names
	my	$sref_svg_folder	=	$_[19]	;	# subfolder of svg print
	my	$sref_pdf_folder	=	$_[20]	;	# subfolder of pdf print
	my	$aref_pdf_sampled	=	$_[21]	;	# list of pdf prints, subsequently collected via pdfunite
	
	######################################################################
	## print density plot using R                                       ##
	######################################################################
	
	######################################################################
	# called from subroutine(s)...
	# &plotsR_quartet_analyses
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	my	$home_path			=	getcwd() 											;	# define home path
	my	$R					=	Statistics::R->new()								;	# load R packge Statistics
	
	#########################################
	# define individual plot-filename with subfolder path
	my	$filename 			=	$$sref_intable	;
	if	( $$sref_filesuffix ){	$filename  =~ s/.tsv$/_${$sref_filesuffix}/ }
	else					 {	$filename  =~ s/.tsv$// 					}
	
	my	@path_parts			=	split "/", $filename ;
	my	$svg_out			=	$$sref_svg_folder."/".$path_parts[-1].".svg"	;
	my	$pdf_out			=	$$sref_pdf_folder."/".$path_parts[-1]			;
	push						@$aref_pdf_sampled, $pdf_out.".pdf" 						;	# store filename with pdf suffix in @pdf_pattern for subsequent pdf unite process
	#########################################
	
	#########################################
	# complete color line for each plot following
	# subroutine assigned custom color lists
	my	$scale_fill			=	join "\", \"", 	@$aref_cfill 		;
	my	$scale_color		=	join "\", \"", 	@$aref_ccolor 		;
	#########################################
	
	#########################################
	# storing individual command lines of geom_point, geom_smooth, facet_grid
	my	@command_lines ;
	#########################################
	
	######## GGPLOT COMMAND LINE #####################################################
	##################################################################################
	# make ggplot 'aes' command line
	my	$table			=	"dp"	;
	my	$aes_parameter	=	"x=".$$sref_x_theader ;
	if	( $$sref_y_theader 	){	$aes_parameter	.= ", y=".$$sref_y_theader			}
	if	( $$sref_group_fill ){	$aes_parameter	.= ", fill=".$$sref_group_fill		}
	if	( $$sref_group_edge ){	$aes_parameter	.= ", colour=".$$sref_group_edge	}
	push @command_lines, "ggplot(".$table.",aes(".$aes_parameter."))" ;
	##################################################################################
	##################################################################################
	
	######## GEOM_DENSITY COMMAND LINE ################################################
	##################################################################################
	push @command_lines, "geom_density(alpha=.4)" ;
	##################################################################################
	##################################################################################
	
	######## SCALE COLOR COMMAND LINE ################################################
	##################################################################################
	push @command_lines, "scale_fill_manual(values=c(\"".$scale_fill."\"))" ;
	push @command_lines, "scale_color_manual(values=c(\"".$scale_color."\"))" ;
	##################################################################################
	##################################################################################
	
	######## GG TITLE COMMAND LINE ###################################################
	##################################################################################
	push @command_lines, "ggtitle(\"".$$sref_ttitle."\")" ;
	##################################################################################
	##################################################################################
	
	######## X,Y LAB COMMAND LINE ####################################################
	##################################################################################
	push @command_lines, "xlab(\"".$$sref_xlab."\") + ylab(\"".$$sref_ylab."\")" ;
	##################################################################################
	##################################################################################
	
	######## FACET_GRID COMMAND LINE #################################################
	##################################################################################
	# define facet_grid parameter (if faced_grid is defined to be activated)
	if	( $$sref_flag_facet ){
		
		my	@parameter ;
		
		#########################################
		# if defined, define horizontal (h) and vertical (v) facet_grid grouping
		# following definition of $$sref_plot_split_h and $$sref_plot_split_v respectively
		# store respective group separation commands @parameter
		if		( ( $$sref_plot_split_h	) && ( $$sref_plot_split_v) ){ push @parameter, $$sref_plot_split_h." ~ ".$$sref_plot_split_v	}
		elsif	( 	$$sref_plot_split_h	)							 { push @parameter, $$sref_plot_split_h." ~ ."						}
		elsif	( 	$$sref_plot_split_v	)							 { push @parameter, 					". ~ ".$$sref_plot_split_v	}
		#########################################
		
		#########################################
		# change x,y scale command of each panel from fixed (default) to defined parameter ($$sref_scale_facet)
		if	( $$sref_scale_facet =~ /free|free_x|free_y/ ){	push @parameter, "scales=\"".$$sref_scale_facet."\"" }
		#########################################
		
		#########################################
		# complete facet_grid parameter line
		my	$facet_grid_command	=	join "," ,	@parameter;	
		#########################################
		
		#########################################
		# store complete geom_smooth command in @command_lines
		push	@command_lines, "facet_grid(".$facet_grid_command.")" ;
		#########################################
	}
	##################################################################################
	##################################################################################
	
	######## LAB LEGEND COMMAND LINE #################################################
	##################################################################################
	# Change legend label of different groupings if both label and grouping (fill, edge, size) are defined
	if ( ( $$sref_group_fill ) && ($$sref_lab_fill) ){ push @command_lines, "labs(fill=\"".$$sref_lab_fill."\")"	}
	if ( ( $$sref_group_edge ) && ($$sref_lab_edge) ){ push @command_lines, "labs(colour=\"".$$sref_lab_edge."\")"	}
	##################################################################################
	##################################################################################
	
	#########################################
	# join geom_point, geom_smooth, facet_grid parameter lines by " + "
	my	$ggplot_parameter	=	join " + ", 	@command_lines		;
	#########################################
	
	my	@theme_lines ;
	######## Theme Plot Title #### ###################################################
	##################################################################################
	push @theme_lines, "plot.title = element_text(color=\"#993333\", size=28, face=\"bold\")" ;
	##################################################################################
	##################################################################################
	
	######## Theme Legend Title #### #################################################
	##################################################################################
	push @theme_lines, "legend.title = element_text(color=\"black\", size=28, face=\"bold\")" ;
	##################################################################################
	##################################################################################
	
	######## Theme Legend Background #################################################
	##################################################################################
	push @theme_lines, "legend.box.background = element_blank()" ;
	#push @theme_lines, "legend.background = element_rect(fill=\"ivory2\", size=2, linetype=\"solid\")" ;
	##################################################################################
	##################################################################################
	
	######## Theme Legend Text #######################################################
	##################################################################################
	push @theme_lines, "legend.text = element_text(color=\"black\", size=24, face=\"bold\")" ;
	##################################################################################
	##################################################################################
	
	######## Theme Legend Position #######################################################
	##################################################################################
	push @theme_lines, "legend.position = \"bottom\"" ;
	push @theme_lines, "legend.box = \"vertical\"" ;
	##################################################################################
	##################################################################################
	
	######## Theme Axis Title ########################################################
	##################################################################################
	push @theme_lines, "axis.title.x = element_text(color=\"black\", size=24, face=\"bold\")" ;
	push @theme_lines, "axis.title.y = element_text(color=\"black\", size=24, face=\"bold\")" ;
	##################################################################################
	##################################################################################
	
	######## Theme Axis Text #########################################################
	##################################################################################
	push @theme_lines, "axis.text.x = element_text(color=\"black\", size=22)" ;
	push @theme_lines, "axis.text.y = element_text(color=\"black\", size=22)" ;
	##################################################################################
	##################################################################################
	
	######## Theme Strip Text ########################################################
	##################################################################################
	push @theme_lines, "strip.text.x = element_text(color=\"black\", size=22, face=\"bold\")" ;
	push @theme_lines, "strip.text.y = element_text(color=\"black\", size=22, face=\"bold\")" ;
	##################################################################################
	##################################################################################
	
	#########################################
	# join geom_point, geom_smooth, facet_grid parameter lines by " + "
	my	$ggplot_theme	=	join ", ", 	@theme_lines		;
	#########################################
	
	#########################################
	# change pdf output height (default 15)
	my	$pdf_height	=	15 ;
	if ( $$sref_pdf_height =~ /\d+/ ){ $pdf_height = $$sref_pdf_height }
	#########################################
	
	#########################################
	# change pdf/svg output width (default 40)
	my	$pdf_width	=	40 ;
	if ( $$sref_pdf_width =~ /\d+/ ){ $pdf_width = $$sref_pdf_width }
	#########################################
	
	#########################################
	# Prepare R commands
	my	@R_commands =(
		'library(ggplot2)',
		'library(svglite)',
		'library(reshape)',
		"setwd(\"$home_path\")",
		"$table <- read.table(\"$$sref_intable\", header = TRUE, sep=\"\\t\")",
		"$ggplot_parameter+theme($ggplot_theme)",
		"plotTitle = \"$svg_out\"",
		"ggsave(plotTitle, width=$pdf_width, height =$pdf_height, limitsize= FALSE)",
		"ggsave(paste(\"$pdf_out\", \".pdf\",sep = \"\"), width = $pdf_width, height = ".$pdf_height.")"
	);
	#########################################
	
	#########################################
	# R execution
	print "\n\t\t...", $filename, " (.svg & .pdf)" ;
	#print "\nR pattern plot of ", $filename, " (.svg & .pdf)" ;
	$R->run(@R_commands);
	$R->stop();
	##################################################################################
}

sub r_tree{
	
	my	$sref_outfilename	=	$_[0]	;	# name of the R plot outfile
	my	$aref_trees			=	$_[1]	;	# list of trees in newick string
	my	$aref_treelabels	=	$_[2]	;	# list of trees assigned labels in newick string (labels in same list order as trees)
	my	$sref_pdf_height	=	$_[3]	;	# height of pdf graphic, keep undef for def height: 15
	my	$sref_pdf_width		=	$_[4]	;	# width of pdf graphic, keep undef for def width: 40
	my	$sref_label_size	=	$_[5]	;	# letter size of tree labels
	my	$sref_svg_folder	=	$_[6]	;	# subfolder of svg print
	my	$sref_pdf_folder	=	$_[7]	;	# subfolder of pdf print
	my	$aref_pdf_sampled	=	$_[8]	;	# list of pdf prints, subsequently collected via pdfunite
	my	$aref_tcolor		=	$_[9]	;	# list of specified tree colors
	######################################################################
	## print line plot using R                                          ##
	######################################################################
	
	######################################################################
	# called from subroutine(s)...
	# &plotsR_quartet_analyses
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	my	$home_path			=	getcwd() 				;	# define home path
	my	$R					=	Statistics::R->new()	;	# load R packge Statistics
	
	my	$labels				=	join "', '"			,	@$aref_treelabels	;
	my	$string_treecodes	=	join ", "			,	@$aref_trees		;
	my	$ntrees			 	= 	@$aref_trees								;
	
	#########################################
	# define individual plot-filename with subfolder path
	my	$svg_out			=	$$sref_svg_folder."/".$$sref_outfilename.".svg"	; 
	my	$pdf_out			=	$$sref_pdf_folder."/".$$sref_outfilename		;
	push						@$aref_pdf_sampled, $pdf_out.".pdf" 			;
	#########################################
	
	#########################################
	# Install and load required packages
	#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
	#BiocManager::install("ggtree")
	#########################################
	
	#########################################
	# list of Newick trees
	#my	@besttrees_cp ;
	#for my $besttree ( 0 .. @$aref_trees-1 ){ $besttrees_cp[$besttree] = "\"".$aref_trees->[$besttree]."\"" }
	#my	$treeline  = join ",", @besttrees_cp ;
	#print "\n\ttreeline:\t",  $treeline;
	#my	$cmdline_nwk_trees  = "newick_trees <- c(".$treeline.")" ; exit;
	#########################################

	#########################################
	# Convert the Newick strings to phylo objects
	#my	$cmdline_nwk_to_phylo  = "phylo_trees <- lapply(newick_trees, function(tree) read.tree(text = tree))" ;
	#########################################

	#########################################
	# Create a ggtree plot for each tree
	#my	$cmdline_tree_to_ggtree  = "tree_plots <- lapply(phylo_trees, function(tree) ggtree(tree) + geom_tree() + geom_tiplab(size=15, colour = \"blue\") )" ;
	#########################################
	
	#########################################
	# Arrange and display the plots
	#my	$cmdline_arrange_plot  = "grid.arrange(grobs = tree_plots, ncol = ".$$sref_ncols.")" ;
	#my	$cmdline_arrange_plot  = "grid.arrange(grobs = tree_plots, ncol = 4, nrow = 4)" ;
	#########################################
	
	my 	$counter_tree = 0 ;
	my	@tree_data_line ;
	my	@read_tree_line ;
	my	@tree_code_line ;
	my	@labels = @$aref_treelabels ;
	my	$geom_tippoint_color = "black" ;
	my	@tree_colors = (
		
		"red"	,
		"blue"	,
		"darkgreen"	,	
		"orange",
	);
	
	my %seen_tree ;
	for my $tree ( @$aref_trees ){
	#for my $tree ( "(A,((B,C),(D,E)));", "(A,(C,(B,(D,E))));", "(G,(C,(B,(D,E))));", "(H,(C,(B,(D,E))));", "(A,((B,C),(D,E)));", "(A,(C,(B,(D,E))));", "(H,(C,(B,(D,E))));", "(H,(C,(B,(D,E))));" ){
		
		my	$tree_color_branches ;
		unless ( $seen_tree{$tree} ){ 
			
			if ( @$aref_tcolor ){ $tree_color_branches = shift @$aref_tcolor } else { $tree_color_branches = "black" }
			$seen_tree{$tree} =  $tree_color_branches  
		}
		else{ $tree_color_branches = $seen_tree{$tree} }
		
		my $treelabel = shift @labels ;
		
		$counter_tree++ ;
		
		push @tree_data_line , "tree_data".$counter_tree." <- \"".$tree."\"" ; 
		push @read_tree_line , "tree".$counter_tree." <- ggtree(read.tree(text = tree_data".$counter_tree."), color = \"".$tree_color_branches."\") + theme_tree('grey89') + geom_rootpoint(color=\"black\", size=5) + geom_tiplab(size = ".$$sref_label_size.") + labs(title = \"BestTree ".$treelabel."\") + theme(plot.title = element_text(size = 20) )" ; #tree1 <- read.tree(text = tree_data1)
		push @tree_code_line , "tree".$counter_tree ;
	}
	
	my $cmdline_trees  			= join "\n", @tree_data_line ;
	my $cmdline_read_trees		= join "\n", @read_tree_line ;
	my $cmdline_arrange_plot	= join "," , @tree_code_line ;
	
	my $ncols = 4;
	my $nrows = sprintf("%d", (@$aref_trees/4)); 
	
	# Arrange and display the trees
	$cmdline_arrange_plot = "arranged_plots <- grid.arrange(".$cmdline_arrange_plot.", ncol = $ncols, nrow = ".$nrows.")" ;
	
		
	#########################################
	# change pdf output height (default 15)
	my	$pdf_height	=	15 ;
	if ( $$sref_pdf_height =~ /\d+/ ){ $pdf_height = $$sref_pdf_height }
	#########################################
	
	#########################################
	# change pdf/svg output width (default 40)
	my	$pdf_width	=	40 ;
	if ( $$sref_pdf_width =~ /\d+/ ){ $pdf_width = $$sref_pdf_width }
	#########################################
	
	#########################################
	# Prepare R commands
	my	@R_commands =(
		"library(ggtree)\n",
		"library(reshape)\n",
		"library(gridExtra)\n",
		#"library(grid)\n",
		"library(ggplot2)\n",
		"setwd(\"$home_path\")\n",
		$cmdline_trees."\n",
		$cmdline_read_trees."\n",
		$cmdline_arrange_plot."\n",
		"plotTitle1 <- \"$svg_out\"\n",
		"plotTitle2 <- \"$pdf_out.pdf\"\n",
		"ggsave(plotTitle1, arranged_plots, width = $pdf_width, height = $pdf_height)\n",
		"ggsave(plotTitle2, arranged_plots, width = $pdf_width, height = $pdf_height)\n",
		#"ggsave(paste(\"$pdf_out\", \".pdf\",sep = \"\"), width = 20, height = 10, units = \"cm\")\n"
	);
	#########################################
	
	#########################################
	# terminal print of generated command line
	my	$command_line_print	=	join "\n", @R_commands ;
	#print "\n\nR commandline:\n", $command_line_print, "\n\n",
	#########################################
	
	#########################################
	# R execution
	print "\n\t\t...", $$sref_outfilename, " (.svg & .pdf)" ;
	#print "\n\t...R pattern plot of ", $filename, " (.svg & .pdf)" ;
	#print"\n", @R_commands,"\n"; #exit;
	$R->run(@R_commands);
	$R->stop(); 

	#exit;
	##################################################################################
}

sub r_lineplot{
	
	my	$sref_intable		=	$_[0]	;	# name of infile table (*.tsv)
	my	$sref_x_theader		=	$_[1]	;	# table header of x-plot categories
	my	$sref_y_theader		=	$_[2]	;	# table header of y-plot values
	my	$sref_group_point	=	$_[3]	;	# table header (factor) for value separation given different point shapes within each plot
	my	$sref_group_line	=	$_[4]	;	# table header (factor) for value separation given different line colors within each plot
	my	$sref_group_shape	=	$_[5]	;	# table header (factor) for showing different elements in different sizes (neglegted if undef)
	my	$sref_do_line_con	=	$_[6]	;	# do lines between plot points -> 1; otherwise, no line connections between points
	my	$sref_size_line		=	$_[7]	;	# size lineplot, default size if undef
	my	$sref_size_point	=	$_[8]	;	# size lineplot, default size if undef
	my	$sref_flag_facet	=	$_[9]	;	# if defined, facet grid is activated in respect of the next two group parameters
	my	$sref_plot_split_v	=	$_[10]	;	# table header plot-splitting vertically (defined as dot "." if splitting not needed)
	my	$sref_plot_split_h	=	$_[11]	;	# table header plot-splitting horizontally (undefined "" if splitting not needed)
	my	$sref_scale_facet	=	$_[12]	;	# default (undef) all the panels have the same scales ('fixed') can be made independent if defined 'free', 'free_x', or 'free_y'
	my	$sref_hline_value	=	$_[13]	;	# highlight horizontal line at axis interception; neglected if undef, for multiple hline use "vhl1,vhl2..."
	my	$sref_hline_type	=	$_[14]	;	# horizontal linetype: solid, longdash, dashed, blank, dotdash, dotted, twodash
	my	$sref_hline_size	=	$_[15]	;	# horizontal linesize; default if undef
	my	$sref_hline_color	=	$_[16]	;	# horizontal line color
	my	$sref_vline_value	=	$_[17]	;	# highlight vertical line at axis interception; neglected if undef, for multiple vline use "vvl1,vvl2..."
	my	$sref_vline_type	=	$_[18]	;	# vertical linetype: solid, longdash, dashed, blank, dotdash, dotted, twodash
	my	$sref_vline_size	=	$_[19]	;	# vertical linesize; default if undef
	my	$sref_vline_color	=	$_[20]	;	# vertical line color
	my	$sref_ttitle		=	$_[21]	;	# table title
	my	$aref_cline			=	$_[22]	;	# list of specified line colors, e.g. (E69F00, 56B4E9) or (coral3, coral4, ...)
	my	$sref_xlab			=	$_[23]	;	# label x-axis
	my	$sref_rotate_xlabel	=	$_[24]	;	# if value == 1 -> rotate x-axis labels via 90°, otherwise keep labels horizontally
	my	$sref_ylab			=	$_[25]	;	# labe y-axis
	my	$sref_lab_point		=	$_[26]	;	# point group label
	my	$sref_lab_line		=	$_[27]	;	# line group label 	
	my	$sref_lab_size		=	$_[28]	;	# point size group label 	
	my	$sref_pdf_height	=	$_[29]	;	# height of pdf graphic, keep undef for def height: 15
	my	$sref_pdf_width		=	$_[30]	;	# width of pdf graphic, keep undef for def width: 40
	my	$sref_filesuffix	=	$_[31]	;	# pdf file suffix to avoid different plots of identic names
	my	$sref_svg_folder	=	$_[32]	;	# subfolder of svg print
	my	$sref_pdf_folder	=	$_[33]	;	# subfolder of pdf print
	my	$aref_pdf_sampled	=	$_[34]	;	# list of pdf prints, subsequently collected via pdfunite
	
	######################################################################
	## print line plot using R                                          ##
	######################################################################
	
	######################################################################
	# called from subroutine(s)...
	# &plotsR_quartet_analyses
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	my	$home_path			=	getcwd() 				;	# define home path
	my	$R					=	Statistics::R->new()	;	# load R packge Statistics
	
	#########################################
	# define individual plot-filename with subfolder path
	my	$filename 			=	$$sref_intable	;
	if	( $$sref_filesuffix ){	$filename  =~ s/.tsv$/_${$sref_filesuffix}/ }
	else					 {	$filename  =~ s/.tsv$// 					}
	
	my	@path_parts			=	split "/", $filename ;
	my	$svg_out			=	$$sref_svg_folder."/".$path_parts[-1].".svg"	; 
	my	$pdf_out			=	$$sref_pdf_folder."/".$path_parts[-1]			;
	push						@$aref_pdf_sampled, $pdf_out.".pdf" 			;	# store filename with pdf suffix in @pdf_pattern for subsequent pdf unite process
	#########################################
	
	#########################################
	# complete color line and point for each plot following
	# subroutine assigned custom color lists
	my	$scale_cline		=	join "\", \"", 	@$aref_cline ;
	#########################################
	
	#########################################
	# storing individual command lines of geom_point, geom_smooth, facet_grid
	my	@command_lines ;
	#########################################
	
	######## GGPLOT COMMAND LINE #####################################################
	##################################################################################
	# make ggplot 'aes' command line
	my	$table			=	"pl"	;
	my	$aes_parameter	=	"x=".$$sref_x_theader.", y=".$$sref_y_theader ;
	
	# changed during filter grafic optimization 29.06.23
	if	( $$sref_do_line_con )
			
			{ push @command_lines	, 	"ggplot(".$table.",aes(".$aes_parameter.", group=1))"
	}
	else	{ push @command_lines	, 	"ggplot(".$table.",aes(".$aes_parameter."))" }
	##################################################################################
	##################################################################################
	
	######## SCALE COLOR COMMAND LINE ################################################
	##################################################################################
	push @command_lines, "scale_color_manual(values=c(\"".$scale_cline."\"))" ;
	##################################################################################
	##################################################################################
	
	######## GG TITLE COMMAND LINE ###################################################
	##################################################################################
	push @command_lines, "ggtitle(\"".$$sref_ttitle."\")" ;
	##################################################################################
	##################################################################################
	
	######## X,Y LAB COMMAND LINE ####################################################
	##################################################################################
	push @command_lines, "xlab(\"".$$sref_xlab."\") + ylab(\"".$$sref_ylab."\")" ;
	##################################################################################
	##################################################################################
	
	######## GEOM_LINE COMMAND LINE #################################################
	##################################################################################
	# define geom_point aesthetics (aes) in respect of defined/undef fill, shape, and size grouping 
	my	$geom_line_command	;
	if	( ( $$sref_group_line ) || ( $$sref_size_line ) ){
		
			$geom_line_command .= "aes(" ;
			
		# if defined store respective group separation commands in @parameter
		my	@parameter			;
		
		if	( $$sref_group_line ){ push @parameter, "color=".$$sref_group_line } # if line grouping is defined
		if	( $$sref_size_line  ){ push @parameter, "size=".$$sref_size_line  } # if line size is defined
		
		# complete aes parameter line
		my	$lstring			 =	join "," ,	@parameter;	
			$geom_line_command	.=	$lstring.")" ;
	}
	
	# store complete geom_point command in @command_lines
	push	@command_lines, "geom_line(".$geom_line_command.")" ;
	##################################################################################
	##################################################################################
	
	######## GEOM_POINT COMMAND LINE #################################################
	##################################################################################
	# define geom_point aesthetics (aes) in respect of defined/undef fill, shape, and size grouping 
	my	$geom_point_command	=	"" ;
	if		( ( $$sref_group_point ) || ( $$sref_group_shape ) || ( $$sref_size_point ) ){
		
		# open aes parameter line
			$geom_point_command	.=	"aes(" ;
		
		# if defined store respective group separation commands in @parameter
		my	@parameter			;
		
		if	( $$sref_group_point ){ push @parameter, "color=".$$sref_group_point } # if point shape factor separation is defined
		if	( $$sref_group_shape ){ push @parameter, "shape=".$$sref_group_shape } # if point size factor separation is defined
		if	( $$sref_size_point	 ){ push @parameter, "size=".$$sref_size_point	 }
		
		# complete aes parameter line
		my	$pstring			 =	join "," ,	@parameter;	
			$geom_point_command	.=	$pstring.")" ;
	}
	
	# store complete geom_point command in @command_lines
	push	@command_lines, "geom_point(".$geom_point_command.")" ;
	##################################################################################
	##################################################################################
	
	######## FACET_GRID COMMAND LINE #################################################
	##################################################################################
	# define facet_grid parameter (if faced_grid is defined to be activated)
	if	( $$sref_flag_facet ){
		
		my	@parameter ;
		
		#########################################
		# if defined, define horizontal (h) and vertical (v) facet_grid grouping
		# following definition of $$sref_plot_split_h and $$sref_plot_split_v respectively
		# store respective group separation commands @parameter
		if		( ( $$sref_plot_split_h	) && ( $$sref_plot_split_v) ){ push @parameter, $$sref_plot_split_h." ~ ".$$sref_plot_split_v	}
		elsif	( 	$$sref_plot_split_h	)							 { push @parameter, $$sref_plot_split_h." ~ ."						}
		elsif	( 	$$sref_plot_split_v	)							 { push @parameter, 					". ~ ".$$sref_plot_split_v	}
		#########################################
		
		#########################################
		# change x,y scale command of each panel from fixed (default) to defined parameter ($$sref_scale_facet)
		if	( $$sref_scale_facet =~ /free|free_x|free_y/ ){	push @parameter, "scales=\"".$$sref_scale_facet."\"" }
		#########################################
		
		#########################################
		# complete facet_grid parameter line
		my	$facet_grid_command	=	join "," ,	@parameter;	
		#########################################
		
		#########################################
		# store complete geom_smooth command in @command_lines
		push	@command_lines, "facet_grid(".$facet_grid_command.")" ;
		#########################################
	}
	##################################################################################
	##################################################################################
	
	######## GEOM_HLINE COMMAND LINE #################################################
	##################################################################################
	# highlight a y-axis line at a certain y-value ($$sref_hline_value, float or integer)
	if	( $$sref_hline_value =~ /\d+(\.\d+)?/ ){
		
		my							$geom_hline	 = "yintercept=c(".$$sref_hline_value.")"	;
		if	( $$sref_hline_color ){	$geom_hline	.= ", color=\"".$$sref_hline_color."\""		}
		if	( $$sref_hline_type  ){	$geom_hline	.= ", linetype=\"".$$sref_hline_type."\""	}
		if	( $$sref_hline_size	 ){ $geom_hline	.= ", size=".$$sref_hline_size				}
									$geom_hline	 =	"geom_hline(".$geom_hline.")"			;
		
		push @command_lines, $geom_hline 
	}
	##################################################################################
	##################################################################################
	
	######## GEOM_VLINE COMMAND LINE #################################################
	##################################################################################
	# highlight a y-axis line at a certain y-value ($$sref_hline_value, float or integer)
	if	( $$sref_vline_value =~ /\d+(\.\d+)?/ ){
		
		my							$geom_vline	 = "xintercept=c(".$$sref_vline_value.")"	;
		if	( $$sref_vline_color ){	$geom_vline	.= ", color=\"".$$sref_vline_color."\""		}
		if	( $$sref_vline_type  ){	$geom_vline	.= ", linetype=\"".$$sref_vline_type."\""	}
		if	( $$sref_vline_size	 ){ $geom_vline	.= ", size=".$$sref_vline_size				}
									$geom_vline	 = "geom_vline(".$geom_vline.")"			;
		
		push @command_lines, $geom_vline 
	}
	##################################################################################
	##################################################################################
	
	######## LAB TEXT COMMAND LINE ###################################################
	##################################################################################
	# Change legend label of different groupings if both label and grouping (fill, edge, size) are defined
	if ( ( $$sref_group_point ) && ($$sref_lab_point) ){ push @command_lines, "labs(shape=\"".$$sref_lab_point."\")"	}
	if ( ( $$sref_group_line  ) && ($$sref_lab_line)  ){ push @command_lines, "labs(colour=\"".$$sref_lab_line."\")"	}
	if ( ( $$sref_group_shape ) && ($$sref_lab_size)  ){ push @command_lines, "labs(size=\"".$$sref_lab_size."\")" 	}
	##################################################################################
	##################################################################################
	
	#########################################
	# join geom_point, geom_smooth, facet_grid parameter lines by " + "
	my	$ggplot_parameter	=	join " + ", 	@command_lines		;
	#########################################
	
	my	@theme_lines ;
	######## Theme Plot Title #### ###################################################
	##################################################################################
	push @theme_lines, "plot.title = element_text(color=\"#993333\", size=28, face=\"bold\")" ;
	##################################################################################
	##################################################################################
	
	######## Theme Legend Title #### #################################################
	##################################################################################
	push @theme_lines, "legend.title = element_text(color=\"black\", size=28, face=\"bold\")" ;
	##################################################################################
	##################################################################################
	
	######## Theme Legend Background #################################################
	##################################################################################
	push @theme_lines, "legend.box.background = element_blank()" ;
	#push @theme_lines, "legend.background = element_rect(fill=\"ivory2\", size=2, linetype=\"solid\")" ;
	##################################################################################
	##################################################################################
	
	######## Theme Legend Text #######################################################
	##################################################################################
	push @theme_lines, "legend.text = element_text(color=\"black\", size=24, face=\"bold\")" ;
	##################################################################################
	##################################################################################
	
	######## Theme Legend Position #######################################################
	##################################################################################
	push @theme_lines, "legend.position = \"bottom\"" ;
	push @theme_lines, "legend.box = \"vertical\"" ;
	##################################################################################
	##################################################################################
	
	######## Theme Axis Title ########################################################
	##################################################################################
	push @theme_lines, "axis.title.x = element_text(color=\"black\", size=24, face=\"bold\")" ;
	push @theme_lines, "axis.title.y = element_text(color=\"black\", size=24, face=\"bold\")" ;
	##################################################################################
	##################################################################################
	
	######## Theme Axis Text #########################################################
	##################################################################################
	if ( $$sref_rotate_xlabel == 1 ){ push @theme_lines, "axis.text.x = element_text(color=\"black\", size=12, angle=90, hjust = 1)"	}
	else							{ push @theme_lines, "axis.text.x = element_text(color=\"black\", size=12)" 			}
	push @theme_lines, "axis.text.y = element_text(color=\"black\", size=12)" ;
	##################################################################################
	##################################################################################
	
	######## Theme Strip Text ########################################################
	##################################################################################
	push @theme_lines, "strip.text.x = element_text(color=\"black\", size=22, face=\"bold\")" ;
	push @theme_lines, "strip.text.y = element_text(color=\"black\", size=22, face=\"bold\")" ;
	##################################################################################
	##################################################################################
	
	#########################################
	# join geom_point, geom_smooth, facet_grid parameter lines by " + "
	my	$ggplot_theme	=	join ", ", 	@theme_lines		;
	#########################################
	
	#########################################
	# change pdf output height (default 15)
	my	$pdf_height	=	15 ;
	if ( $$sref_pdf_height =~ /\d+/ ){ $pdf_height = $$sref_pdf_height }
	#########################################
	
	#########################################
	# change pdf/svg output width (default 40)
	my	$pdf_width	=	40 ;
	if ( $$sref_pdf_width =~ /\d+/ ){ $pdf_width = $$sref_pdf_width }
	#########################################
	
	#########################################
	# Prepare R commands
	my	@R_commands =(
		'library(ggplot2)',
		'library(svglite)',
		'library(reshape)',
		"setwd(\"$home_path\")",
		"$table <- read.table(\"$$sref_intable\", header = TRUE, sep=\"\\t\")",
		"$ggplot_parameter+theme($ggplot_theme)",
		"plotTitle = \"$svg_out\"",
		"ggsave(plotTitle, width=$pdf_width, height =$pdf_height, limitsize= FALSE)",
		"ggsave(paste(\"$pdf_out\", \".pdf\",sep = \"\"), width = $pdf_width, height = $pdf_height)"
	);
	#########################################
	
	#########################################
	# terminal print of generated command line
	my	$command_line_print	=	join "\n", @R_commands ;
	#print "\n\nR commandline:\n", $command_line_print, "\n\n",
	#########################################
	
	#########################################
	# R execution
	print "\n\t\t...", $filename, " (.svg & .pdf)" ;
	#print "\n\t...R pattern plot of ", $filename, " (.svg & .pdf)" ;
	#print"\n", @R_commands,"\n"; #exit;
	$R->run(@R_commands);
	$R->stop();
	##################################################################################
}

sub r_pointplot{
	
	my	$sref_intable		=	$_[0]	;	# name of infile table (*.tsv)
	my	$sref_x_theader		=	$_[1]	;	# table header of x-plot categories
	my	$sref_y_theader		=	$_[2]	;	# table header of y-plot values (keep empty for density plot)
	my	$sref_group_fill	=	$_[3]	;	# table header (factor) for value separation given different fill colors within each plot (neglegted if undef)
	my	$sref_group_edge	=	$_[4]	;	# table header (factor) for value separation given different edge colors within each plot (neglegted if undef)
	my	$sref_group_size	=	$_[5]	;	# table header (factor) for showing different elements in different sizes (neglegted if undef)
	my	$sref_flag_smooth	=	$_[6]	;	# smooth regular expression line activated if defined, otherwise neglegted
	my	$sref_group_smooth	=	$_[7]	;	# table header (factor) for smooth regression line colour separation (neglegted if undef, smooth must be activated)
	my	$sref_smooth_fill	=	$_[8]	;	# if 1 -> error intervals of smooth line are colored due to smooth factor, otherwise default is used
	my	$sref_flag_facet	=	$_[9]	;	# if defined, facet grid is activated in respect of the next two group parameters
	my	$sref_plot_split_v	=	$_[10]	;	# table header plot-splitting vertically (defined as dot "." if splitting not needed)
	my	$sref_plot_split_h	=	$_[11]	;	# table header plot-splitting horizontally (undefined "" if splitting not needed)
	my	$sref_scale_facet	=	$_[12]	;	# default (undef) all the panels have the same scales ('fixed') can be made independent if defined 'free', 'free_x', or 'free_y'
	my	$sref_hline_value	=	$_[13]	;	# highlight horizontal line at axis interception; neglected if undef, for multiple hline use "vhl1,vhl2..."
	my	$sref_hline_type	=	$_[14]	;	# horizontal linetype: solid, longdash, dashed, blank, dotdash, dotted, twodash
	my	$sref_hline_size	=	$_[15]	;	# horizontal linesize; default if undef
	my	$sref_hline_color	=	$_[16]	;	# horizontal line color
	my	$sref_vline_value	=	$_[17]	;	# highlight vertical line at axis interception; neglected if undef, for multiple vline use "vvl1,vvl2..."
	my	$sref_vline_type	=	$_[18]	;	# vertical linetype: solid, longdash, dashed, blank, dotdash, dotted, twodash
	my	$sref_vline_size	=	$_[19]	;	# vertical linesize; default if undef
	my	$sref_vline_color	=	$_[20]	;	# vertical line color
	my	$sref_ttitle		=	$_[21]	;	# table title
	my	$aref_cfill			=	$_[22]	;	# list of specified fill colors, e.g. (E69F00, 56B4E9) or (coral3, coral4, ...)
	my	$aref_ccolor		=	$_[23]	;	# list of specified edge colors
	my	$sref_xlab			=	$_[24]	;	# label x-axis
	my	$sref_ylab			=	$_[25]	;	# label y-axis
	my	$sref_lab_fill		=	$_[26]	;	# fill colour label
	my	$sref_lab_edge		=	$_[27]	;	# edge (outline) colour label 	
	my	$sref_lab_shape		=	$_[28]	;	# shape label 
	my	$sref_lab_size		=	$_[29]	;	# size label 	
	my	$sref_pdf_height	=	$_[30]	;	# height of pdf graphic, keep undef for def height: 15
	my	$sref_pdf_width		=	$_[31]	;	# width of pdf graphic, keep undef for def width: 40
	my	$sref_filesuffix	=	$_[32]	;	# pdf file suffix to avoid different plots of identic names
	my	$sref_svg_folder	=	$_[33]	;	# subfolder of svg print
	my	$sref_pdf_folder	=	$_[34]	;	# subfolder of pdf print
	my	$aref_pdf_sampled	=	$_[35]	;	# list of pdf prints, subsequently collected via pdfunite
	
	######################################################################
	## print point plot using R                                         ##
	######################################################################
	
	######################################################################
	# called from subroutine(s)...
	# &plotsR_quartet_analyses
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	my	$home_path			=	getcwd() 											;	# define home path
	my	$R					=	Statistics::R->new()								;	# load R packge Statistics
	
	#########################################
	# define individual plot-filename with subfolder path
	my	$filename 			=	$$sref_intable	;
	if	( $$sref_filesuffix ){	$filename  =~ s/.tsv$/_${$sref_filesuffix}/ }
	else					 {	$filename  =~ s/.tsv$// 					}
	
	my	@path_parts			=	split "/", $filename ;
	my	$svg_out			=	$$sref_svg_folder."/".$path_parts[-1].".svg"	;
	my	$pdf_out			=	$$sref_pdf_folder."/".$path_parts[-1]			;
	push						@$aref_pdf_sampled, $pdf_out.".pdf" 			;	# store filename with pdf suffix in @pdf_pattern for subsequent pdf unite process
	#########################################
	
	#########################################
	# complete color line for each plot following
	# subroutine assigned custom color lists
	my	$scale_fill			=	join "\", \"", 	@$aref_cfill 		;
	my	$scale_color		=	join "\", \"", 	@$aref_ccolor 		;
	#########################################
	
	#########################################
	# storing individual command lines of geom_point, geom_smooth, facet_grid
	my	@command_lines ;
	#########################################
	
	######## GGPLOT COMMAND LINE #####################################################
	##################################################################################
	# make ggplot 'aes' command line
	my	$table			=	"pp"	;
	my	$aes_parameter	=	"x=".$$sref_x_theader.", y=".$$sref_y_theader ;
	#if	( $$sref_group_fill ){	$aes_parameter	.= ", fill=".$$sref_group_edge		}
	#if	( $$sref_group_edge ){	$aes_parameter	.= ", shape=".$$sref_group_edge.", size=".$$sref_group_edge	}
	push @command_lines, "ggplot(".$table.",aes(".$aes_parameter."))" ;
	##################################################################################
	##################################################################################
	
	######## SCALE COLOR COMMAND LINE ################################################
	##################################################################################
	push @command_lines, "scale_fill_manual(values=c(\"".$scale_fill."\"))" ;
	push @command_lines, "scale_color_manual(values=c(\"".$scale_color."\"))" ;
	##################################################################################
	##################################################################################
	
	######## GG TITLE COMMAND LINE ###################################################
	##################################################################################
	push @command_lines, "ggtitle(\"".$$sref_ttitle."\")" ;
	##################################################################################
	##################################################################################
	
	######## X,Y LAB COMMAND LINE ####################################################
	##################################################################################
	push @command_lines, "xlab(\"".$$sref_xlab."\") + ylab(\"".$$sref_ylab."\")" ;
	##################################################################################
	##################################################################################
	
	######## GEOM_POINT COMMAND LINE #################################################
	##################################################################################
	# define geom_point aesthetics (aes) in respect of defined/undef fill, shape, and size grouping 
	my	$geom_point_command	=	"" ;
	if		( ( $$sref_group_fill )	||	( $$sref_group_edge ) || ( $$sref_group_size ) ){
		
		# open aes parameter line
			$geom_point_command	.=	"aes(" ;
		
		# if defined store respective group separation commands @parameter
		my	@parameter			;
		
		if	( $$sref_group_fill ){ push @parameter,	"colour=".$$sref_group_fill 	} # if point fill factor separation is defined
		if	( $$sref_group_edge ){ push @parameter,	"shape=".$$sref_group_edge 		} # if point shape factor separation is defined
		if	( $$sref_group_size ){ push @parameter,	"size=".$$sref_group_size 		} # if point size factor separation is defined
		
		# complete aes parameter line
		my	$pstring			 =	join "," ,	@parameter;	
			$geom_point_command	.=	$pstring.")" ;
	}
	
	# store complete geom_point command in @command_lines
	push	@command_lines, "geom_point(".$geom_point_command.")" ;
	##################################################################################
	##################################################################################
	
	######## GEOM_HLINE COMMAND LINE #################################################
	##################################################################################
	# highlight a y-axis line at a certain y-value ($$sref_hline_value, float or integer)
	if	( $$sref_hline_value =~ /\d+(\.\d+)?/ ){
		
		my							$geom_hline	 = "yintercept=c(".$$sref_hline_value.")"	;
		if	( $$sref_hline_color ){	$geom_hline	.= ", color=\"".$$sref_hline_color."\""		}
		if	( $$sref_hline_type  ){	$geom_hline	.= ", linetype=\"".$$sref_hline_type."\""	}
		if	( $$sref_hline_size	 ){ $geom_hline	.= ", size=".$$sref_hline_size				}
									$geom_hline	 =	"geom_hline(".$geom_hline.")"			;
		
		push @command_lines, $geom_hline 
	}
	##################################################################################
	##################################################################################
	
	######## GEOM_VLINE COMMAND LINE #################################################
	##################################################################################
	# highlight a y-axis line at a certain y-value ($$sref_hline_value, float or integer)
	if	( $$sref_vline_value =~ /\d+(\.\d+)?/ ){
		
		my							$geom_vline	 = "xintercept=c(".$$sref_vline_value.")"	;
		if	( $$sref_vline_color ){	$geom_vline	.= ", color=\"".$$sref_vline_color."\""		}
		if	( $$sref_vline_type  ){	$geom_vline	.= ", linetype=\"".$$sref_vline_type."\""	}
		if	( $$sref_vline_size	 ){ $geom_vline	.= ", size=".$$sref_vline_size				}
									$geom_vline	 = "geom_vline(".$geom_vline.")"			;
		
		push @command_lines, $geom_vline 
	}
	##################################################################################
	##################################################################################
	
	######## GEOM_SMOOTH COMMAND LINE ################################################
	##################################################################################
	# define geom_smooth parameter (if smooth is defined to be activated)
	if	( $$sref_flag_smooth ){
		
		# use method=lm as standard parameter
		my	$geom_smooth_command	=	"method=lm" ;
		
		# if defined, extend smooth parameter due to colour separation of defined groups ($$sref_group_smooth)
		if	( $$sref_group_smooth	){	$geom_smooth_command	.=	",aes(colour=".$$sref_group_smooth	}
		if	( $$sref_smooth_fill	){	$geom_smooth_command	.=	",fill=".$$sref_group_smooth		}
		if	( ( $$sref_group_smooth	) || ( $$sref_smooth_fill	) ){ $geom_smooth_command	.=	")" }
		
		# store complete geom_smooth command in @command_lines
		push	@command_lines, "geom_smooth(".$geom_smooth_command.")" ;
	}
	##################################################################################
	##################################################################################
	
	######## FACET_GRID COMMAND LINE #################################################
	##################################################################################
	# define facet_grid parameter (if faced_grid is defined to be activated)
	if	( $$sref_flag_facet ){
		
		my	@parameter ;
		
		#########################################
		# if defined, define horizontal (h) and vertical (v) facet_grid grouping
		# following definition of $$sref_plot_split_h and $$sref_plot_split_v respectively
		# store respective group separation commands @parameter
		if		( ( $$sref_plot_split_h	) && ( $$sref_plot_split_v) ){ push @parameter, $$sref_plot_split_h." ~ ".$$sref_plot_split_v	}
		elsif	( 	$$sref_plot_split_h	)							 { push @parameter, $$sref_plot_split_h." ~ ."						}
		elsif	( 	$$sref_plot_split_v	)							 { push @parameter, 					". ~ ".$$sref_plot_split_v	}
		#########################################
		
		#########################################
		# change x,y scale command of each panel from fixed (default) to defined parameter ($$sref_scale_facet)
		if	( $$sref_scale_facet =~ /free|free_x|free_y/ ){	push @parameter, "scales=\"".$$sref_scale_facet."\"" }
		#########################################
		
		#########################################
		# complete facet_grid parameter line
		my	$facet_grid_command	=	join "," ,	@parameter;	
		#########################################
		
		#########################################
		# store complete geom_smooth command in @command_lines
		push	@command_lines, "facet_grid(".$facet_grid_command.")" ;
		#########################################
	}
	##################################################################################
	##################################################################################
	
	######## LAB TEXT COMMAND LINE ###################################################
	##################################################################################
	# Change legend label of different groupings if both label and grouping (fill, edge, size) are defined
	if ( ( $$sref_group_fill ) && ($$sref_lab_fill) ){ push @command_lines, "labs(colour=\"".$$sref_lab_fill."\")"	}
	if (   $$sref_smooth_fill 						){ push @command_lines, "labs(fill=\"".$$sref_lab_edge."\")"	}
	if ( ( $$sref_group_edge ) && ($$sref_lab_edge) ){ push @command_lines, "labs(shape=\"".$$sref_lab_shape."\")"	}
	if ( ( $$sref_group_size ) && ($$sref_lab_size) ){ push @command_lines, "labs(size=\"".$$sref_lab_size."\")" 	}
													   #push @command_lines, "labs(fill=element_blank())" ;
	##################################################################################
	##################################################################################
	
	#########################################
	# join geom_point, geom_smooth, facet_grid parameter lines by " + "
	my	$ggplot_parameter	=	join " + ", 	@command_lines		;
	#########################################
	
	my	@theme_lines ;
	######## Theme Plot Title #### ###################################################
	##################################################################################
	push @theme_lines, "plot.title = element_text(color=\"#993333\", size=28, face=\"bold\")" ;
	##################################################################################
	##################################################################################
	
	######## Theme Legend Title #### #################################################
	##################################################################################
	push @theme_lines, "legend.title = element_text(color=\"black\", size=28, face=\"bold\")" ;
	##################################################################################
	##################################################################################
	
	######## Theme Legend Background #################################################
	##################################################################################
	push @theme_lines, "legend.box.background = element_blank()" ;
	#push @theme_lines, "legend.background = element_rect(fill=\"ivory2\", size=2, linetype=\"solid\")" ;
	##################################################################################
	##################################################################################
	
	######## Theme Legend Text #######################################################
	##################################################################################
	push @theme_lines, "legend.text = element_text(color=\"black\", size=24, face=\"bold\")" ;
	##################################################################################
	##################################################################################
	
	######## Theme Legend Position #######################################################
	##################################################################################
	push @theme_lines, "legend.position = \"bottom\"" ;
	push @theme_lines, "legend.box = \"vertical\"" ;
	##################################################################################
	##################################################################################
	
	######## Theme Axis Title ########################################################
	##################################################################################
	push @theme_lines, "axis.title.x = element_text(color=\"black\", size=24, face=\"bold\")" ;
	push @theme_lines, "axis.title.y = element_text(color=\"black\", size=24, face=\"bold\")" ;
	##################################################################################
	##################################################################################
	
	######## Theme Axis Text #########################################################
	##################################################################################
	push @theme_lines, "axis.text.x = element_text(color=\"black\", size=22)" ;
	push @theme_lines, "axis.text.y = element_text(color=\"black\", size=22)" ;
	##################################################################################
	##################################################################################
	
	######## Theme Strip Text ########################################################
	##################################################################################
	push @theme_lines, "strip.text.x = element_text(color=\"black\", size=22, face=\"bold\")" ;
	push @theme_lines, "strip.text.y = element_text(color=\"black\", size=22, face=\"bold\")" ;
	##################################################################################
	##################################################################################
	
	#########################################
	# join geom_point, geom_smooth, facet_grid parameter lines by " + "
	my	$ggplot_theme	=	join ", ", 	@theme_lines		;
	#########################################
	
	#########################################
	# change pdf output height (default 15)
	my	$pdf_height	=	15 ;
	if ( $$sref_pdf_height =~ /\d+/ ){ $pdf_height = $$sref_pdf_height }
	#########################################

	#########################################
	# change pdf/svg output width (default 40)
	my	$pdf_width	=	40 ;
	if ( $$sref_pdf_width =~ /\d+/ ){ $pdf_width = $$sref_pdf_width }
	#########################################
	
	#########################################
	# Prepare R commands
	my	@R_commands =(
		'library(ggplot2)',
		'library(svglite)',
		'library(reshape)',
		"setwd(\"$home_path\")",
		"$table <- read.table(\"$$sref_intable\", header = TRUE, sep=\"\\t\")",
		"$ggplot_parameter+theme($ggplot_theme)",
		"plotTitle = \"$svg_out\"",
		"ggsave(plotTitle, width=$pdf_width, height =$pdf_height, limitsize= FALSE)",
		"ggsave(paste(\"$pdf_out\", \".pdf\",sep = \"\"), width = $pdf_width, height = $pdf_height)"
	);
	#########################################
	
	#########################################
	# terminal print of generated command line
	#my	$command_line_print	=	join "\n", @R_commands ;
	#print "\n\nR commandline:\n", $command_line_print, "\n\n",
	#########################################
	
	#########################################
	# R execution
	#print"\n", @R_commands,"\n"; #exit;
	print "\n\t\t...", $filename, " (.svg & .pdf)" ;
	#print "\nR pattern plot of ", $filename, " (.svg & .pdf)" ;
	$R->run(@R_commands);
	$R->stop();
	##################################################################################
}

sub r_ternplot{
	
	my	$sref_intable		=	$_[0]	;	# name of infile table (*.tsv)
	my	$sref_x_theader		=	$_[1]	;	# table header of x-plot categories
	my	$sref_y_theader		=	$_[2]	;	# table header of y-plot values
	my	$sref_z_theader		=	$_[3]	;	# table header of z-plot values
	my	$sref_group_fill	=	$_[4]	;	# table header (factor) for value separation given different fill colors within each plot
	my	$sref_flag_facet	=	$_[5]	;	# if defined, facet grid is activated in respect of the next two group parameters
	my	$sref_plot_split_v	=	$_[6]	;	# table header plot-splitting vertically (defined as dot "." if splitting not needed)
	my	$sref_plot_split_h	=	$_[7]	;	# table header plot-splitting horizontally (undefined "" if splitting not needed)
	my	$sref_ttitle		=	$_[8]	;	# table title
	my	$sref_size			=	$_[9]	;	# point size within triangle
	my	$sref_pdf_height	=	$_[10]	;	# height of pdf graphic, keep undef for def height: 15
	my	$sref_pdf_width		=	$_[11]	;	# width of pdf graphic, keep undef for def width: 40
	my	$sref_filesuffix	=	$_[12]	;	# pdf file suffix to avoid different plots of identic names
	my	$sref_svg_folder	=	$_[13]	;	# subfolder of svg print
	my	$sref_pdf_folder	=	$_[14]	;	# subfolder of pdf print
	my	$aref_pdf_sampled	=	$_[15]	;	# list of pdf prints, subsequently collected via pdfunite
	
	######################################################################
	## print ternary plot using R                                       ##
	######################################################################
	
	######################################################################
	# called from subroutine(s)...
	# &plotsR_quartet_analyses
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	my	$table				=	"tp"	;
	my	$home_path			=	getcwd() 											;	# define home path
	my	$R					=	Statistics::R->new()								;	# load R packge Statistics
	
	#########################################
	# define individual plot-filename with subfolder path
	my	$filename 			=	$$sref_intable	;
	if	( $$sref_filesuffix ){	$filename  =~ s/.tsv$/_${$sref_filesuffix}/ }
	else					 {	$filename  =~ s/.tsv$// 					}
	
	my	(	$svg_out, $pdf_out 	); 
	my		@path_parts			=	split "/", $filename ;
	if	( $$sref_svg_folder		){	$svg_out	=	$$sref_svg_folder."/".$path_parts[-1].".svg"	} else{	$svg_out = $filename.".svg"	}
	if	( $$sref_pdf_folder		){	$pdf_out	=	$$sref_pdf_folder."/".$path_parts[-1]			} else{	$pdf_out = $filename		}
	
	# store filename with pdf suffix in @pdf_pattern for subsequent pdf unite process
	push						@$aref_pdf_sampled, $pdf_out.".pdf" 			;
	#########################################
	
	#########################################
	# change pdf output height (default 15)
	my	$pdf_height	=	15 ;
	if ( $$sref_pdf_height =~ /\d+/ ){ $pdf_height = $$sref_pdf_height }
	#########################################
	
	#########################################
	# change pdf output width (default 40)
	my	$pdf_width	=	40 ;
	if ( $$sref_pdf_width =~ /\d+/ ){ $pdf_width = $$sref_pdf_width }
	#########################################
	
	#########################################
	# storing individual command lines of geom_point, geom_smooth, facet_grid
	my	@command_lines ;
	#########################################
	
	######## GEOM_LINE COMMAND LINE ##################################################
	##################################################################################
	# define ggtern aesthetics (aes) in respect of defined fill grouping 
	my	$ggtern_command = "ggtern(data=".$table.", aes(x=".$$sref_x_theader.",y=".$$sref_y_theader.",z=".$$sref_z_theader."))";
	push	@command_lines, $ggtern_command ;
	##################################################################################
	##################################################################################
	
	######## GEOM_Point COMMAND LINE #################################################
	##################################################################################
	# define ggtern aesthetics (aes) in respect of defined fill grouping 
	my	$tsize	= 6 ;
	if	( $$sref_size =~ /\d+/ ){ $tsize = $$sref_size } 
	my	$geom_point_command = "geom_point(aes(fill=".$$sref_group_fill."),size=".$tsize.",shape=21,color=\"black\")";
	push	@command_lines, $geom_point_command ;
	##################################################################################
	##################################################################################
	
	######## GG TITLE COMMAND LINE ###################################################
	##################################################################################
	push @command_lines, "ggtitle(\"".$$sref_ttitle."\")" ;
	##################################################################################
	##################################################################################
	
	######## FACET_GRID COMMAND LINE #################################################
	##################################################################################
	# define facet_grid parameter (if faced_grid is defined to be activated)
	if	( $$sref_flag_facet ){
		
		my	@parameter ;
		
		#########################################
		# if defined, define horizontal (h) and vertical (v) facet_grid grouping
		# following definition of $$sref_plot_split_h and $$sref_plot_split_v respectively
		# store respective group separation commands @parameter
		if		( ( $$sref_plot_split_h	) && ( $$sref_plot_split_v) ){ push @parameter, $$sref_plot_split_h." ~ ".$$sref_plot_split_v	}
		elsif	( 	$$sref_plot_split_h	)							 { push @parameter, $$sref_plot_split_h." ~ ."						}
		elsif	( 	$$sref_plot_split_v	)							 { push @parameter, 					". ~ ".$$sref_plot_split_v	}
		#########################################
		
		#########################################
		# complete facet_grid parameter line
		my	$facet_grid_command	=	join "," ,	@parameter;	
		#########################################
		
		#########################################
		# store complete geom_smooth command in @command_lines
		push	@command_lines, "facet_grid(".$facet_grid_command.")" ;
		#########################################
	}
	##################################################################################
	##################################################################################
	
	######## X,Y LAB COMMAND LINE ####################################################
	##################################################################################
	push @command_lines, "labs(fill=\"".$$sref_group_fill."\")" ;
	##################################################################################
	##################################################################################
	
	######## THEME ###################################################################
	##################################################################################
	push @command_lines, "theme_rgbg()+theme(legend.position=c(0,1),legend.justification=c(0,1))" ;
	##################################################################################
	##################################################################################
	
	#########################################
	# join ggtern, geom_point, facet_grid parameter lines by " + "
	my	$ggplot_parameter	=	join " + ", 	@command_lines		;
	#########################################
	
	#########################################
	# Prepare R commands
	my	@R_commands =(
		'library(ggplot2)',
		'library(svglite)',
		'library(reshape)',
		'library(ggtern)',
		"setwd(\"$home_path\")",
		"$table <- read.table(\"$$sref_intable\", header = TRUE, sep=\"\\t\")",
		"$ggplot_parameter",
		"plotTitle = \"$svg_out\"",
		"ggsave(plotTitle, width=".$pdf_width.", height =".$pdf_height.", limitsize= FALSE)",
		"ggsave(paste(\"$pdf_out\", \".pdf\",sep = \"\"), width = ".$pdf_width.", height = ".$pdf_height.")"
	);
	#########################################
	
	#########################################
	# terminal print of generated command line
	print "\n\t\t...", $filename, " (.svg & .pdf)" ;
	#print "\nR pattern plot of ", $filename, " (.svg & .pdf)" ;
	my	$command_line_print	=	join "\n", @R_commands ;
	#print "\n\nR commandline:\n", $command_line_print, "\n\n",
	#########################################
	
	#########################################
	# R execution
	#print "\nR pattern plot of ", $filename, " (.svg & .pdf)" ;
	$R->run(@R_commands);
	$R->stop();
	##################################################################################
}

sub readin_fasta{

	my	$href_set_opt		= $_[0]	;	#	# key1: sealion user available parameter options; value: optional parameter setting			IN: defined; OUT: unchanged
	
	######################################################################
	# called from subroutine(s)...
	# &check_msa_infile
	#
	############
	# using subroutine(s)...
	# &help
	######################################################################
	
	######################################################################
	# Read, check, and assign .fas formated alignment data               #
	######################################################################
	
	####################################
	# READ IN fasta interleaved and non-interleaved formated files
	# store single taxa and their associated sequence in has %seq_of_tax
	# key: taxon name ; value: taxon associated sequence
	my $fas_file	= $href_set_opt->{sys_parameter}{msafile} ;
	
	my 	(
			%seq_of_tax 		, 	# key1: taxon; value: sequence
			$taxon 				, 	# sequence name
			%seen_fasta_taxon		# key1: taxon; value: counter
	) ;
	
	open INfas, $fas_file or die "\nMSAFILE-ERROR: Cannot read ", $fas_file, "!\n\n";
	
	while ( my $line = <INfas> ){

		chomp	$line ;
		
		$line =~ s/\s+//g;

		#######################
		# identification of taxon names
		# delete spaces within taxon names
		# check if taxon names appear multiple times
		# check if taxon names consist of only alphanumeric signs
		if ( $line =~ /^\>/ )	{

			( $taxon = $line ) =~ s/^\>|\s+// ;

			if ( ( $taxon =~ /\w+/ ) && ( $seen_fasta_taxon{$taxon} ) ){

				print "\nMSAFILE-ERROR: OTU ", $taxon, " appears multiple times!\n\n"; 
				
				&help( \'i', \%$href_set_opt )
			}

			elsif( $taxon =~ /\w+/ ){ $seen_fasta_taxon{$taxon}++ }

			else { 
				
				print "\nMSAFILE-ERROR: OTU ", $taxon, " includes non-alphanumeric signs!\n\n"; 
					
				&help( \'i', \%$href_set_opt ) 
			}
		}
		#######################

		#######################
		# store single sequence lines of identified taxon as hash value of taxon
		elsif ( $seq_of_tax{$taxon} ){ $seq_of_tax{$taxon} .= $line }
		else { $seq_of_tax{$taxon} = $line }
		#######################
	}

	close INfas ;
	####################################

	#######################
	# return taxa and their associated sequence in has %seq_of_tax
	return %seq_of_tax ;
	#######################
	
	######################################################################
}

sub readIN_file{
	
	my	$sref_file = $_[0] ;
	
	#print "::&readIN_file ", $$sref_file, "\n" ;
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	## Read file and assign each line as list-element                   ##
	######################################################################
	
	####################################
	# Read IN file and chomp lines
	open IN, "<$$sref_file" or die "\nFILE-ERROR: Cannot read file ", $$sref_file, "!\n\n" ;
	chomp ( my @lines = <IN> ); close IN ;
	####################################
	
	####################################
	# return list of file lines
	return @lines ;
	####################################
	
	###################################################
}

sub read_in_p4_results{
	
	my $href_value_of_scr_par					= $_[0] ; # key1: sealion user available parameter options; value: optional parameter setting				IN: defined; OUT: unchanged
	my $href_data_of_fileprop					= $_[1]	; # key1: property; value : property value															IN: defined; OUT: unchanged
	my $sref_p4_outfile_name					= $_[2] ; # name of the p4 outfile																			IN: defined; OUT: unchanged
	my $href_hoh_found_N_of_topo_of_pattern		= $_[3] ; # key1: topology; key2: recoded site pattern; value: N site pattern in original data				IN: undefined; OUT: defined
	my $href_hoh_expected_N_of_topo_of_pattern	= $_[4] ; # key1: topology; key2: recoded site pattern; value: expected ML N site pattern					IN: undefined; OUT: defined
	
	######################################################################
	# called from subroutine(s)...
	# &quartet_building_and_p4_main
	#
	############
	# using subroutine(s)...
	# &recode_pattern
	######################################################################
	
	######################################################################
	# Read and Recode P4 evaluated site pattern distributions            #
	######################################################################
	
	##########################################################
	open IN, "<$$sref_p4_outfile_name" or die "\nINFILE-ERROR:Cannot read P4 result file ", $$sref_p4_outfile_name , "!\n\n" ;
	print "\n\tRead p4 result file\n";

	#############################
	# three topologies in P4 outfile
	my $topology ;

	#############################
	# read in each single line of the p4 outputfile
	while ( my $line = <IN> ){

		#############################
		# reject newline sign
		chomp $line ;

		#############################
		# Extraction of each topology with associated difference between expected and observed pattern frequencies
		# something went wrong if P4 prints 'bad', stop read in
		if	( $line =~ /Bad/ ){ die "\nP4-ERROR:\n", $line, "!\n\n" }

		#############################
		# if line lists a topology without any other info, e.g. '(A, B, (C, D));'
		# assign topology to string $topology
		elsif 	(	( $line =~ /^\(\w/ ) && ( $line !~ /\d+\);$/ )	){ #print "\nline:\t", $line, "\n";

			##############################################
			# Change unrooted quartet tree to rooted tree
			# (A, B, (C, D)); -> ((A,B),(C,D))
			(	$topology = $line )	=~ s/\s+|\n|;//g	;	#print "\tTopology:\t", $topology, "\n" ;
				#$topology			=~ s/^\(/\(\(/g		;	print "\tTopology:\t", $topology, "\n" ;
				#$topology			=~ s/,\(/\),\(/g	;	print "\tTopology:\t", $topology, "\n" ; exit;
		}

		#############################
		# if line consists of a quartet pattern followed by tab or multiple space signs,
		# followed by oberved number of pattern positions, followed again by tab or multiple space signs,
		# followed by number of P4 expected pattern number, e.g. 'aaaa        125       135.96',
		elsif	(	(	$line =~ /^\w\w\w\w\s+\d+\s+\d+/	) ||
					(	$line =~ /^\(.\w\w\w\w.,/			) ){ 

			#############################
			# split line into three blocks by substitution of multiple space signs or tab to a single space sign,
			# and splitting by single space signs
			# @line_blocks = (pattern (XXXX), Nobserved, Nexpected)
				$line			=~	s/\(|\)|,|\'//g 	; #print "\nline:\t", $line, "\n";
				$line			=~	s/\s+/ /g 			; #print "\nline:\t", $line, "\n";
			my	@line_blocks	=	split " ", $line 	; #for (@line_blocks){print "\n", $_, "\n"} #exit;

			#############################
			# recode quartet pattern by splitting pattern into single states
			# via subroutine &recode_pattern, recoded pattern
				$line_blocks[0]		=~	tr/acgt/ACGT/ ;
			my	@site_pattern		=	split "", $line_blocks[0] ;
			my	$recoded_pattern 	=	&recode_pattern(

													\@site_pattern					, # list of pattern states								IN: defined; OUT: unchanged
													\%$href_value_of_scr_par		, # key1: sealion parameter options; value: parametr	IN: defined; OUT: unchanged
													\$href_data_of_fileprop->{type}	, # sequenec state ('nuc' or 'aa' )						IN: defined; OUT: unchanged
										) ;
			#print "\norig.\t", $line_blocks[0], "\nrecoded\t",  $recoded_pattern ; 
			#############################

			############################# (1)
			# Add observed pattern number to %hoh_found_N_of_topo_of_pattern
			$href_hoh_found_N_of_topo_of_pattern->{$topology}{$recoded_pattern}		+= $line_blocks[1] ;
			#############################

			############################# (2)
			# Add expected pattern number to %hoh_expected_N_of_topo_of_pattern
			$href_hoh_expected_N_of_topo_of_pattern->{$topology}{$recoded_pattern}	+= $line_blocks[2] ;
			#############################
		}
		#############################
	}
	close IN ;
	#exit;
	#############################
	# delete p4 result file
	unlink $$sref_p4_outfile_name ;
	##########################################################
	
	######################################################################
}

sub readin_phylip{

	my	$href_set_opt		= $_[0]	;	#	# key1: sealion user available parameter options; value: optional parameter setting			IN: defined; OUT: unchanged
	
	######################################################################
	# called from subroutine(s)...
	# &check_msa_infile
	#
	############
	# using subroutine(s)...
	# &help
	######################################################################
	
	######################################################################
	# Read, check, and assign .phy formated alignment data               #
	######################################################################

	####################################
	# READ IN phylip interleaved and non-interleaved formated files
	my	$phy_file = $href_set_opt->{sys_parameter}{msafile} ;
	open INphy , $phy_file or die "\nMSAFILE-ERROR: Cannot read ", $phy_file, "!\n\n" ;
	
	chomp (my @all_lines_phy = <INphy>) and close INphy ;
	####################################
	
	####################################
	# exclude empty lines (indelible simulated files often include empty lines)
	my @cleaned_lines_phy	;
	
	my $linenumber 		= 1 ;
	for my $line (@all_lines_phy){ if ( $line =~ /\w+|\?|\-/ ){ push @cleaned_lines_phy, $line; $linenumber++} }
	
	@all_lines_phy		= @cleaned_lines_phy ;
	@cleaned_lines_phy 	= () ;
	####################################

	####################################
	# extract the first line (infoline) and determine the number of taxa ($info_line[0])
	( my $infoline	= shift @all_lines_phy ) =~ s/\s+/ / ;
	my	@infos_line	=  split " ", $infoline ;
	
	unless ( ( $infos_line[0] =~ /\d+/ ) && ( $infos_line[1] =~ /\d+/ ) ){

		print "\nMSAFILE-ERROR: No phylip format! Missing first line with OTU and site number!\n\n";
		
		&help( \'i', \%$href_set_opt )
	}
	####################################

	####################################
	# phylip files can be in interleaved and non-interleaved format
	# interleaved:	tax1 ACGT...	# first part of all taxa
	# 				tax2 ACGT...
	#								#''space line''
	#				tax1 CCCC...	# Second part of all taxa
	#				tax2 GGGG...
	# to concatenate sequence parts correctly store single lines in separate hashkeys
	# until number of taxa is reached ($c equal $infos_line[0]). afterwards remove the following spaceline and concatenate
	# next lines to their corresponding taxon sequences inferred from the first round and so on...
	# If phylip file is in non-interleaved format, the while lopp stops automatically after the first foreach loop
	my 	%seq_phy = () ;
	
	while ( @all_lines_phy ){

		my $N_lines = @all_lines_phy ; #print "\nremain: ", $N_lines, "\n";
		for ( my $c=1; $c<=$infos_line[0]; $c++ ){ my $seq_line_phy = shift @all_lines_phy ; push ( @{$seq_phy{$c}} , $seq_line_phy ) }
	}
	####################################
	
	####################################
	# join single sequence parts of each taxon (if interleaved there are multiple key values)
	# taxonnames are in the same line as sequenceinformation (separated by one or multiple whitespaces), therefore
	# substitute multiple whitespaces into one whitespace, join all sequence parts to one string (only important if interleaved format)
	# split taxonnames from sequenceinformation and store both in the hashreference %href_seq_of_tax (key: taxon; value: complete sequence)
	my ( %seq_of_tax, %seen_phylip_taxon ) ;
	
	for my $line_c ( sort {$a<=>$b} keys %seq_phy ){ #print "\n", $line_c;

		my	@seq_single_parts				=	exists($seq_phy{$line_c}) ? @{$seq_phy{$line_c}} :( ) ;
			$seq_single_parts[0]			=~	s/\s+/:::/ ;
		my	$seq_complete					=	join	"", @seq_single_parts ;
			$seq_complete					=~	s/\s+//g ;
			@seq_single_parts				=	split	":::", $seq_complete ;
		my	$taxon							=	shift @seq_single_parts ; #print "\ntaxon: ", $taxon, "\t", length $seq_single_parts[0], "\n" ;

		if ( ( $taxon =~ /\w+/ ) && ( $seen_phylip_taxon{$taxon} ) ){

			print "\nMSAFILE-ERROR: OTU ", $taxon, " appears multiple times in multiple-sequence alignment!\n\n";
			
			&help( \'i', \%$href_set_opt )
		}

		elsif ( $taxon =~ /\w+/ ){ $seen_phylip_taxon{$taxon}++ }

		else { 
			
			print "\nMSAFILE-ERROR: OTU ", $taxon, " includes non-alphanumeric signs in multiple-sequence alignment!\n\n"; 
			
			&help( \'i', \%$href_set_opt ) 
		}

		$seq_of_tax{$taxon}	=	$seq_single_parts[0] ;
		@seq_single_parts	= () ;
	}
	
	%seq_phy = () ; #exit;
	####################################

	####################################
	# return taxa and their associated sequence in has %seq_of_tax
	return %seq_of_tax ;
	####################################
	
	######################################################################
}

sub read_write_spd_files{
	
	my	$aref_l_qcalcfile				= $_[0]  ; # list of infile-lines of '*_detailed_split_calc_*'												IN: defined; OUT: unchanged
	my	$href_obs_n_of_pattern 			= $_[1]  ; # key: pattern code (e.g. 'A'); value number of observed pattern frequency						IN: undefined; OUT: defined
	my	$href_exp_n_of_tree1_pattern 	= $_[2]  ; # key: pattern code (e.g. 'A'); value number of expected pattern frequency given tree 1			IN: undefined; OUT: defined
	my	$href_exp_n_of_tree2_pattern 	= $_[3]  ; # key: pattern code (e.g. 'A'); value number of expected pattern frequency given tree 2			IN: undefined; OUT: defined
	my	$href_exp_n_of_tree3_pattern 	= $_[4]  ; # key: pattern code (e.g. 'A'); value number of expected pattern frequency given tree 3			IN: undefined; OUT: defined
	my	$href_obs_p_of_pattern 			= $_[5]  ; # key: pattern code (e.g. 'A'); value percentage of observed pattern frequency					IN: undefined; OUT: defined
	my	$href_exp_p_of_tree1_pattern 	= $_[6]  ; # key: pattern code (e.g. 'A'); value percentage of expected pattern frequency given tree 1		IN: undefined; OUT: defined
	my	$href_exp_p_of_tree2_pattern 	= $_[7]  ; # key: pattern code (e.g. 'A'); value percentage of expected pattern frequency given tree 2		IN: undefined; OUT: defined
	my	$href_exp_p_of_tree3_pattern 	= $_[8]  ; # key: pattern code (e.g. 'A'); value percentage of expected pattern frequency given tree 3		IN: undefined; OUT: defined
	my	$sref_Npat_total_obs			= $_[9]  ; # Total number of different split pattern														IN: undefined; OUT: defined
	my	$href_value_of_scr_par			= $_[10] ; # key1: sealion user available parameter options; value: optional parameter setting				IN: defined; OUT: unchanged
	my	$href_subfolder_q				= $_[11] ; # key: output-type (e.g. 'pdf'); value: pathdir													IN: defined; OUT: unchanged
	my	$aref_qtaxa						= $_[12] ; # list of quartet depending taxa in the same order as in the infile								IN: defined; OUT: unchanged
	my	$href_clan_of_seq				= $_[13] ; # key1: sequence name; value: assigned clan 														IN: defined; OUT: unchanged
	my	$sref_fc						= $_[14] ; # actual quartet-number																			IN: defined; OUT: unchanged
	my	$href_seqL_of_fc				= $_[15] ; # key1: quartet-number ($fc); value: quartet related sequence length (corrected)					IN: undefined; OUT: defined
	my	$sref_p_qcalc_file				= $_[16] ; # actual spd infile																				IN: defined; OUT: unchanged
	
	#print "::&read_write_spd_files\n";
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	## assign obs and exp pattern frequencies '*_detailed_split_calc_*' ##
	######################################################################
	my $pattern_tree_line ;
	for my $line ( @$aref_l_qcalcfile ){
		
		#########################################
		# extract tree line of following split-pattern frequencies
		if 	  ( $line =~ /^Pattern\s+Observed/ ){ $pattern_tree_line = $line }
		#########################################
		
		#########################################
		# extract data from split-pattern frequency-lines
		elsif ( $line =~ /^Pattern\s+\w\t/ ){
			
			#########################################
			# substitute possible string 'Pattern ' at the beginning of the line
			# older calculation files have additional string "Pattern " at line beginning
			# e.g. Pattern	A	XXYY	:	1077		923.61			571.88			569.47
			$line =~ s/^Pattern\s+//g	;
			#print "\ns1: ", $line ;
			#########################################
			
			#########################################
			# delete pattern explanation (e.g. 'XXYY\t:')
			$line =~ s/^(\w\t)\w\w\w\w\t(:)?/$1/g	; 
			$line =~ s/://g	; 
			#print "\ns2: ", $line ;
			#########################################
			
			#########################################
			# substitute multiple blanks in a row to a tabstop sign
			$line =~ s/\s+/\t/g ;
			#print "\ns3: ", $line ;
			#########################################
			
			#########################################
			# split line at tabstop signs and sample single elements in @lparts
			my	@lparts	= split "\t", $line ;
			#for my $e ( @lparts ){ print "\ne: ", $e }
			#########################################
			
			#########################################
			# add number of observed splitpattern of actual pattern code to total number of pbserved split pattern
			$$sref_Npat_total_obs	+=	$lparts[1] ;
			#########################################
			
			#########################################
			# assign pattern distribution of actual pattern code to correspondinh hash-variables
			$href_obs_n_of_pattern		->{$lparts[0]} = $lparts[1] ;
			$href_exp_n_of_tree1_pattern->{$lparts[0]} = $lparts[2] ;
			$href_exp_n_of_tree2_pattern->{$lparts[0]} = $lparts[3] ;
			$href_exp_n_of_tree3_pattern->{$lparts[0]} = $lparts[4] ;
			
			if ( $href_value_of_scr_par->{sys_parameter}{screen_prt} == 1 ){
				
				print 	"\n\t"				, $lparts[0]									, 
						"\tN_pat_obs:\t"	, $href_obs_n_of_pattern->{$lparts[0]}			, 
						"\tN_pat_exp1:\t"	, $href_exp_n_of_tree1_pattern->{$lparts[0]}	,	 
						"\tN_pat_exp2:\t"	, $href_exp_n_of_tree2_pattern->{$lparts[0]}	,
						"\tN_pat_exp3:\t"	, $href_exp_n_of_tree3_pattern->{$lparts[0]} 	
			}
			#########################################
		}
		#########################################
		
		#########################################
		# check correct clan association of each of the four taxa in spd file
		elsif ( $line =~ /^Clan assignment\t\w+\t\w+/ ){
			
			my @prts = split "\t", $line ;
			if ( $href_clan_of_seq->{$prts[1]} ne $prts[2] ){
				
				die "\nSPDFILE-ERROR: Clanfile uncongruent defined clan '", $prts[2] ,"' of OTU '", $prts[1] , "' in ", $$sref_p_qcalc_file, "!\n\n" ;
			}
		}
		#########################################
	}
	
	$href_seqL_of_fc->{$$sref_fc} = $$sref_Npat_total_obs ;
	#########################################
	
	#########################################
	# for each N pattern frequency, calculate the percentage of pattern frequencies related to the total number of site-patterns of given quartet
	for my $pattern ( qw/A B C D E F G H J K L M N V I/){
		
		$href_obs_p_of_pattern->{$pattern}			=	&calc_percent( $href_obs_n_of_pattern->{$pattern}		, $$sref_Npat_total_obs )	;
		$href_exp_p_of_tree1_pattern->{$pattern}	=	&calc_percent( $href_exp_n_of_tree1_pattern->{$pattern}	, $$sref_Npat_total_obs )	;
		$href_exp_p_of_tree2_pattern->{$pattern}	=	&calc_percent( $href_exp_n_of_tree2_pattern->{$pattern}	, $$sref_Npat_total_obs )	;
		$href_exp_p_of_tree3_pattern->{$pattern}	=	&calc_percent( $href_exp_n_of_tree3_pattern->{$pattern}	, $$sref_Npat_total_obs )	;
		
		if ( $href_value_of_scr_par->{sys_parameter}{screen_prt} == 1 ){
			
			print 	"\n\t"				, $pattern									, 
					"\tP_pat_obs:\t"	, $href_obs_p_of_pattern->{$pattern}		, 
					"\tP_pat_exp1:\t"	, $href_exp_p_of_tree1_pattern->{$pattern}	,	 
					"\tP_pat_exp2:\t"	, $href_exp_p_of_tree2_pattern->{$pattern}	,
					"\tP_pat_exp3:\t"	, $href_exp_p_of_tree3_pattern->{$pattern} 	
		}
	}
	#########################################
	
	#########################################
	# print new SPD file if optionally defined
	if ( $href_value_of_scr_par->{sys_parameter}{write_spd} == 1 ){
		
		#######################
		# check firts line (taxonline) and change Qnumber of first line to actual quartet number
		my $line1_taxa = () ;
		if ( $aref_l_qcalcfile->[0] =~ /^Q\d+:\t\w+:\w+:\w+:\w+/ ){
		
			( $line1_taxa = $aref_l_qcalcfile->[0] ) =~ s/^Q\d+/Q${$sref_fc}/ 
		}
		
		else{ die "\nSPDFILE-ERROR: Cannot identify OTUs in first line of spdfile ", $$sref_p_qcalc_file, "!\n\n"  }
		#######################
		
		#######################
		# define and make new outfile pathdir
		my		@clan4_uncod	= sort ( $href_clan_of_seq->{$aref_qtaxa->[0]}, $href_clan_of_seq->{$aref_qtaxa->[1]}, $href_clan_of_seq->{$aref_qtaxa->[2]}, $href_clan_of_seq->{$aref_qtaxa->[3]} ) ;
		my		$clan4_uncoded	= join "_", @clan4_uncod ;
		my		$spdpath		= $href_subfolder_q->{spd}."/".$clan4_uncoded ;
		mkdir 	$spdpath 		;
		#######################
		
		#######################
		# open filehandle
		my		$outfile = "SeaLion_detailed_split_calc_q".$$sref_fc.".txt" ;
		my	 	$spdfile = $spdpath."/".$outfile ;
		
		print "\t\twriting SPD file\t", $spdfile, "\n" ;

		open 	OUTspd, ">$spdfile" || die "\nSPDFILE-ERROR: Cannot write to spdfile ", $spdfile, "!\n\n" ;
		#######################
		
		#######################
		# print header lines
		print OUTspd	$line1_taxa,
						"\n---------------------------------------------------------------"	,
						"\nAlignment observed and P4 expected pattern frequencies" ,
						"\n\n", $pattern_tree_line ;
		#######################
		
		#######################
		# print observed frequencies
		for my $pattern ( qw/A B C D E F G H J K L M N V I/){
			
			print OUTspd	"\nPattern\t"	, $pattern, 
							"\t:\t"			, $href_obs_n_of_pattern	 ->{$pattern}	,
							"\t\t"			, $href_exp_n_of_tree1_pattern->{$pattern}	,
							"\t\t\t"		, $href_exp_n_of_tree2_pattern->{$pattern}	,
							"\t\t\t"		, $href_exp_n_of_tree3_pattern->{$pattern}	;
		}
		
		print	OUTspd	"\n---------------------------------------------------------------\n" 							,
						"\nSequence length of split informative site pattern (bp):\t", $href_seqL_of_fc->{$$sref_fc}	,
						"\nClan assignment\t", $aref_qtaxa->[0],"\t", $href_clan_of_seq->{$aref_qtaxa->[0]}				,
						"\nClan assignment\t", $aref_qtaxa->[1],"\t", $href_clan_of_seq->{$aref_qtaxa->[1]}				,
						"\nClan assignment\t", $aref_qtaxa->[2],"\t", $href_clan_of_seq->{$aref_qtaxa->[2]}				,
						"\nClan assignment\t", $aref_qtaxa->[3],"\t", $href_clan_of_seq->{$aref_qtaxa->[3]}				;
		#######################
		
		#######################
		# close filehandle
		close OUTspd;
		#######################
	}
	#########################################
	#exit;
	######################################################################
}

sub recode_filter_printout{
	
	my	$sref_filter	=	$_[0] ;	# string with script-coded filter methods (e.g. risk1, nap_risk1, napv_risk1_qdist1...) IN: defined; OUT: unchanged
	
	######################################################################
	# called from subroutine(s)...
	# &plotsR_quartet_analyses
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	# translate script-internal coded filtermethods of input 			 #
	# string to output firendly lettering                                #
	######################################################################
	
	my	$filtertype	; my $flag = 0 ;
	if  	( $$sref_filter =~ /risk1_qdist/  					){ $flag++; $filtertype = "RISK_DIST"  		}
	elsif  	( $$sref_filter =~ /risk2_qdist/  					){ $flag++; $filtertype = "RISK2_DIST"  	}
	elsif  	( $$sref_filter =~ /risk3_qdist/  					){ $flag++; $filtertype = "RISK3_DIST"  	}
	elsif  	( $$sref_filter =~ /risk4_qdist/  					){ $flag++; $filtertype = "RISK4_DIST"  	}
	elsif  	( $$sref_filter =~ /risk1/  						){ $flag++; $filtertype = "RISK" 			}
	elsif  	( $$sref_filter =~ /risk2/  						){ $flag++; $filtertype = "RISK2" 			}
	elsif  	( $$sref_filter =~ /risk3/  						){ $flag++; $filtertype = "RISK3" 			}
	elsif  	( $$sref_filter =~ /risk4/  						){ $flag++; $filtertype = "RISK4" 			}
	elsif   ( $$sref_filter =~ /qdist1/ 						){ $flag++; $filtertype = "DIST"  			}
	elsif	( $$sref_filter =~ /nap|napv|napr|unfiltered|none/ 	){ $flag++; $filtertype = "unfiltered" 		}
	
	if   	( $$sref_filter =~ /strict_out.softclipping/ 	){ $flag++; $filtertype = "OTU-STRICT-OUT-Soft"	}
	elsif   ( $$sref_filter =~ /strict_in.softclipping/ 	){ $flag++; $filtertype = "OTU-STRICT-IN-Soft"	}
	elsif   ( $$sref_filter =~ /strict_all.softclipping/ 	){ $flag++; $filtertype = "OTU-STRICT-ALL-Soft"	}
	elsif   ( $$sref_filter =~ /strict_out.hardclipping/ 	){ $flag++; $filtertype = "OTU-STRICT-OUT-Hard"	}
	elsif   ( $$sref_filter =~ /strict_in.hardclipping/ 	){ $flag++; $filtertype = "OTU-STRICT-IN-Hard"	}
	elsif   ( $$sref_filter =~ /strict_all.hardclipping/ 	){ $flag++; $filtertype = "OTU-STRICT-ALL-Hard"	}
	elsif   ( $$sref_filter =~ /strict_out/ 				){ $flag++; $filtertype = "OTU-STRICT-OUT"	}
	elsif   ( $$sref_filter =~ /strict_in/ 					){ $flag++; $filtertype = "OTU-STRICT-IN"	}
	elsif   ( $$sref_filter =~ /strict_all/ 				){ $flag++; $filtertype = "OTU-STRICT-ALL"	}
	
	unless	( $flag ){ die "\nBUG-ERROR: Cannot assign '".$$sref_filter."' in subroutine &recode_filter_printout!\nPlease, report BUG to software developer!\n\n" }
	
	return $filtertype ;
	######################################################################
}

sub recode_pattern{

	my $aref_pattern_states		= $_[0] ; # list of pattern states								IN: defined; OUT: unchanged
	my $href_val_of_scr_par		= $_[1] ; # key1: sealion parameter options; value: parametr	IN: defined; OUT: unchanged
	my $sref_sequence_state		= $_[2] ; # sequenec state ('nuc' or 'aa' )						IN: defined; OUT: unchanged
	
	#print "::&recode_pattern\n";
	
	######################################################################
	# called from subroutine(s)...
	# &read_in_p4_results
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	# Recode actual site pattern                                         #
	######################################################################
	
	#######################
	# If msa pattern has to be recoded to state classes (like purine or pyrimidine if nucleotide data)
	# substitute single character states of list @$aref_pattern_states to corresponding class states
	if ( $href_val_of_scr_par->{sys_parameter}{rycoding} == 1 ){

		if ( $$sref_sequence_state eq 'nuc' ){

			map { s/A|G/R/; $_ } @$aref_pattern_states;	# substitute A and G to P (Purine)
			map { s/C|T/Y/; $_ } @$aref_pattern_states	# substitute C and T to Y (Pyrimidine)
		}

		elsif ( $$sref_sequence_state eq 'aa' ){

			map { s/A|W|M|I|L|F|P/R/; $_			} @$aref_pattern_states;	# substitute A|W|M|I|L|F|P to H (Hydrophobic)
			map { s/C|G|T|N|Y|R|S|K|D|V|E|Q/Y/; $_	} @$aref_pattern_states;	# substitute C|G|T|N|Y|R|S|K|D|V|E|Q to Y (Hydrophilic)
		}

		else{ die "\nBUG-ERROR: Cannot assign data type in subroutine &recode_pattern!\n\tPlease, report BUG to software developer!\n\n" }
	}
	#######################

	##############################################
	# START quartet pattern recoding

	#######################
	# definition of single quartet pattern codes
	# key state pattern; value: quartet pattern code
	my %code_of_pattern = (
							'XXYY' => 'A',
							'XYXY' => 'B',
							'XYYX' => 'C',
							'XXYW' => 'D',
							'XYWW' => 'E',
							'XYWY' => 'F',
							'XYXW' => 'G',
							'XYYW' => 'H',
							'XYWX' => 'J',
							'XYYY' => 'K',
							'XYXX' => 'L',
							'XXYX' => 'M',
							'XXXY' => 'N',
							'XYWZ' => 'V',
							'XXXX' => 'I'
						) ;
	#######################

	#######################
	# START quartet pattern recoding via if queries
	if ( $aref_pattern_states->[0] eq $aref_pattern_states->[1] ){

		if ( $aref_pattern_states->[1] eq $aref_pattern_states->[2] ){

			if	( $aref_pattern_states->[2] eq $aref_pattern_states->[3] )	{ return $code_of_pattern{XXXX} }	# I
			else 															{ return $code_of_pattern{XXXY} }	# N
		}

		elsif	( $aref_pattern_states->[1] eq $aref_pattern_states->[3] )	{ return $code_of_pattern{XXYX} }	# M
		elsif	( $aref_pattern_states->[2] eq $aref_pattern_states->[3] )	{ return $code_of_pattern{XXYY} }	# A
		else																{ return $code_of_pattern{XXYW} }	# D
	}

	elsif( $aref_pattern_states->[0] eq $aref_pattern_states->[2] ){

		if		( $aref_pattern_states->[2] eq $aref_pattern_states->[3] )	{ return $code_of_pattern{XYXX} }	# L
		elsif	( $aref_pattern_states->[1] eq $aref_pattern_states->[3] )	{ return $code_of_pattern{XYXY} }	# B
		else																{ return $code_of_pattern{XYXW} }	# G
	}

	elsif( $aref_pattern_states->[0] eq $aref_pattern_states->[3] ){

		if		( $aref_pattern_states->[1] eq $aref_pattern_states->[2] )	{ return $code_of_pattern{XYYX} }	# C
		else																{ return $code_of_pattern{XYWX} }	# J
	}

	elsif( $aref_pattern_states->[1] eq $aref_pattern_states->[2] ){

		if		( $aref_pattern_states->[2] eq $aref_pattern_states->[3] )	{ return $code_of_pattern{XYYY} }	# K
		else																{ return $code_of_pattern{XYYW} }	# H
	}

	elsif		( $aref_pattern_states->[2] eq $aref_pattern_states->[3] )	{ return $code_of_pattern{XYWW} }	# E
	elsif		( $aref_pattern_states->[1] eq $aref_pattern_states->[3] )	{ return $code_of_pattern{XYWY} }	# F

	else 																	{ return $code_of_pattern{XYWZ} }	# V
	# END quartet pattern recoding via if queries
	#######################
}

sub sampling_single_4clan_taxon_support{
	
	my	$href_value_of_scr_param								= $_[0]	 ;	# key1: sealion user available parameter options; value: optional parameter setting																							IN: defined; OUT: unchanged
	my	$href_hol_nwk_and_support_of_clan4_of_taxa4				= $_[1]  ;	# key1: 4clan code; key2: 4taxon-string; value: list of 'tree::supp:' of the three 4clan relationships																		IN: undefined; OUT: defined
	my	$sref_cm												= $_[2]  ;	# type of calculation (nap, napv...) and or filter (rsik_1, risk1_dist1...) method																							IN: defined; OUT: unchanged 
	my	$sref_clan4												= $_[3]	 ;	# actual 4clan combi code																																					IN: defined; OUT: unchanged
	my	$href_hoh_sup_of_clan4_of_calc_of_qtree_of_qn			= $_[4]	 ;	# multiple hash, storing the final support (nap or napv) related to each rooted quartet-tree																				IN: defined; OUT: unchanged
	my	$href_hoh_seen_rejected_of_4clan_of_filtertype_of_qn	= $_[5]	 ;	# key1: 4clan code; key2: type of supportfilter ('nap_risk1', 'napr_risk1', or 'napv_risk1'); key3: actual sealion quartet number; value: filter corresponding score-value	IN: defined; OUT: unchanged 
	my	$href_taxa4_of_qn										= $_[6]	 ;	# key1: quartet-number ($fc); value : quartet assigned taxa (alphabetically ordered, joined by ':'																			IN: defined; OUT: unchanged 
	my	$sref_filttype											= $_[7]	 ;	# type of calculation (nap, napv...) and or filter (rsik_1, risk1_dist1...) method																							IN: defined; OUT: unchanged 
	my	$href_hol_sup_of_4clan_of_tax_of_tree					= $_[8]	 ;	# key1: 4clan code; key2: 4taxon-string; key3: quartet-tree; value: list of taxon																							IN: undefined; OUT: defined
	my	$href_seen_qn_of_4clan									= $_[9]	 ;	# key1: 4clan code; key2: quartet-number; value; list of 4clan assigned quartet number																						IN: undefined; OUT: defined
	my	$href_seen_tax_of_4clan									= $_[10] ;	# key1: 4clan code; key2: remaining taxon after OTU-filtering; value: counter number																						IN: undefined; OUT: defined
	my	$href_seen_alltax_of_4clan								= $_[11] ;	# key1: 4clan code; key2: each 4clan assigned taxon; value: counter number																									IN: undefined; OUT: defined
	
	#print	"::&sampling_single_4clan_taxon_support\n";
	
	######################################################################
	# called from subroutine(s)...
	# &evaluation_otu_support_main
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	# sampling of actual filter-method remaining 4clan taxon-support     #
	######################################################################
	for my $tree ( sort keys %{$href_hoh_sup_of_clan4_of_calc_of_qtree_of_qn->{$$sref_clan4}{$$sref_cm}} ){
		
		############################
		# filtering unneeded hashkeys of 
		# %{$href_hoh_sup_of_clan4_of_calc_of_qtree_of_fc->{$clan4}{$cm}{$tree}}
		# which are NOT integers ('mean' or 'median')
		my	@qnumbers ;
		for my $code ( sort keys %{$href_hoh_sup_of_clan4_of_calc_of_qtree_of_qn->{$$sref_clan4}{$$sref_cm}{$tree}} ){
			
			if ( $code =~ /\d+/ ){ push @qnumbers, $code }
		}
		############################
		
		############################
		# go through all quartet numbers and sample for each quartet corresponding taxon (yet) non-filtered support
		my %seen_t4 ;
		for my $qn ( sort {$a<=>$b} @qnumbers ){
		
			############################
			# split qn assigned 4taxon string ('tax1:tax2:tax3:tax4')
			my	@taxa = split ":", $href_taxa4_of_qn->{$qn} ;
			for my $t ( @taxa ){ $href_seen_alltax_of_4clan->{$$sref_clan4}{$t}++ } 
			############################
			
			############################
			# sample tree and assigned support (joined by '::') in 
			# @{$href_hol_nwk_and_support_of_clan4_of_taxa4->{$$sref_clan4}{$href_taxa4_of_qn->{$qn}}}
			if		( $$sref_filttype =~ /^nap$|^napr$|^napv$/ ){
					
					############################
					# without any previously conducted sequence_filtering, sample tree::support of all 4clan assigned quartets
					for my $t ( @taxa ){ 
					
						push @{$href_hol_sup_of_4clan_of_tax_of_tree->{$$sref_clan4}{$t}{$tree}} , $href_hoh_sup_of_clan4_of_calc_of_qtree_of_qn->{$$sref_clan4}{$$sref_cm}{$tree}{$qn}	; 
						$href_seen_tax_of_4clan->{$$sref_clan4}{$t}++ ; 
					}
				
					$href_seen_qn_of_4clan->{$$sref_clan4}{$qn}++ ;
					push @{$href_hol_nwk_and_support_of_clan4_of_taxa4->{$$sref_clan4}{$href_taxa4_of_qn->{$qn}}}, $tree."::".$href_hoh_sup_of_clan4_of_calc_of_qtree_of_qn->{$$sref_clan4}{$$sref_cm}{$tree}{$qn};
					############################
			}
			############################

			############################
			# sample tree::support of risk and qdist filter remaining quartets
			else{
				
				unless ( $href_hoh_seen_rejected_of_4clan_of_filtertype_of_qn->{$$sref_clan4}{$$sref_filttype}{$qn} ){
					
					for my $t ( @taxa ){ 
						
						push @{$href_hol_sup_of_4clan_of_tax_of_tree->{$$sref_clan4}{$t}{$tree}} , $href_hoh_sup_of_clan4_of_calc_of_qtree_of_qn->{$$sref_clan4}{$$sref_cm}{$tree}{$qn} ;
						$href_seen_tax_of_4clan->{$$sref_clan4}{$t}++ ;
					}
					
					$href_seen_qn_of_4clan->{$$sref_clan4}{$qn}++ ;
					push @{$href_hol_nwk_and_support_of_clan4_of_taxa4->{$$sref_clan4}{$href_taxa4_of_qn->{$qn}}}, $tree."::".$href_hoh_sup_of_clan4_of_calc_of_qtree_of_qn->{$$sref_clan4}{$$sref_cm}{$tree}{$qn} ;
				}
			}
			############################
		}#exit;
		############################
		
		############################
		# test print
		if ( $href_value_of_scr_param->{sys_parameter}{screen_prt} == 1 ){
			
			for my $qn ( sort {$a <=> $b} @qnumbers ){ 
			
				for my $tree_supp ( @{$href_hol_nwk_and_support_of_clan4_of_taxa4->{$$sref_clan4}{$href_taxa4_of_qn->{$qn}}}){
					
					#print "\n\tqn:", $qn, "\t4taxa:", $href_taxa4_of_qn->{$qn}, "\t", $tree_supp
				}
			}
		} #exit ;
		############################
	}
	
	
	######################################################################
}

sub sampling_subdirectories{
	
	my	$sref_searchdir		= $_[0] ;	# main directory for filepath search 		IN: defined; OUT: unchanged
	my	$sref_searchfile	= $_[1] ;	# file to sample from subdir(s)				IN: defined; OUT: unchanged
	
	######################################################################
	## calculate and return mean from a list of values                  ##
	######################################################################
	#print "::&sampling_subdirectories\n";
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	#########################################
	# start file search...
	#print "::searching DIR ", $$sref_searchdir, " for ", $$sref_searchfile ,"\n" ;
	
	# store each penguin calculation-file '*_detailed_split_calc_*' file with path as list in @found
	my @found =	File::Find::Rule->file->name($$sref_searchfile)->in($$sref_searchdir) ;
	#my @files	=	File::Find::Rule->file->name('*_detailed_split_calc_*.txt')->in($presult_folder);
	#########################################
	
	#########################################
	# return list of sampled files
	if ( @found ){ return @found } else { die "\nSPDFILE-ERROR: Cannot find SPD file(s) ", $$sref_searchfile, " in pathDIR ", $$sref_searchdir, "!\n\n" }
	#########################################
	
	######################################################################
}

sub start_p4{
	
	my $sref_alignment 		= $_[0]	; # generated quartet alignment				IN: defined; OUT: undefined
	my $sref_script_p4		= $_[1]	; # P4 input commands						IN: defined; OUT: undefined
	my $sref_clan4			= $_[2]	; # actual 4clan combination				IN: defined; OUT: unchanged
	my $sref_qcounter		= $_[3]	; # actual quartet number					IN: defined; OUT: unchanged
	my $sref_outfile_p4		= $_[4]	; # P4 resultfile							IN: undefined; OUT: defined
	my $sref_tquartet		= $_[5]	; # 4 taxon-quartet names in alphab. order	IN: defined; OUT: unchanged
	
	######################################################################
	# called from subroutine(s)...
	# &quartet_building_and_p4_main
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	# Start P4 quartet analysis for given alignment and p4-commands      #
	######################################################################
	print	"\n\tStart P4...\n\tQ", $$sref_qcounter, ": ", $$sref_clan4," (", $$sref_tquartet, ")\n" ;
	
	my		$program			=	"p4" ;
			$$sref_outfile_p4	=	"p4_result_".$$sref_clan4."_q".$$sref_qcounter.".txt" ;
	
	##############################################
	# system start
	system	( "$program $$sref_script_p4 >$$sref_outfile_p4" )	;
	##############################################
	
	##############################################
	# deletion of p4 infiles
	unlink	( $$sref_script_p4, $$sref_alignment );
	######################################################################
}

sub site_exclusion{

	my $href_hol_seqstates_of_t		= $_[0] ; # key1: taxon name; key2: state position number; value: sequence state at that position	IN: defined; OUT: undefined
	my $href_data_of_fileprop		= $_[1]	; # key1: property; value : property value													IN: defined; OUT: unchanged
	my $aref_sampled				= $_[2]	; # list of the 4 taxon names in alphab. order												IN: defined; OUT: unchanged
	my $href_red_seq_of_t			= $_[3]	; # key1: taxon; value: gap&amb reduced sequence											IN: undefined; OUT: defined
	my $sref_seq_length				= $_[4]	; # common length of reduced quartet-sequences												IN: undefined; OUT: defined
	
	######################################################################
	# called from subroutine(s)...
	# &quartet_building_and_p4_main
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	# Build all sequence-quartets of 4clan combination                   #
	######################################################################
	
	##############################################
	# Define allowed sequence states, all other sequence states are rejected from further split analyses
	my %allowed_states ;
	
	if			( $href_data_of_fileprop->{type} =~ 'nuc'	){ for ( qw/A C G T/									){ $allowed_states{$_}++ } } # all nuc states allowed
	elsif		( $href_data_of_fileprop->{type} =~ 'aa'	){ for ( qw/A C G T N Y R W S K M D V H I E L Q F P/	){ $allowed_states{$_}++ } } # all aa states allowed
	else 		{ die "\nBUG-ERROR: Cannot assign data type in subroutine &perform_quartet_analyses!\nPlease, report BUG to software developer!\n\n" }
	##############################################
	
	##############################################
	# analyse quartet states of each sequence positions
	# reject sites of unallowed characters
	for my $seq_pos ( 0 .. $href_data_of_fileprop->{length}-1 ){ 
	
		#print "\n\t", $seq_pos, "\t", $aref_sampled->[0], "\t", $aref_sampled->[1], "\t", $aref_sampled->[2], "\t", $aref_sampled->[3], "\t", $href_hol_seqstates_of_t->{$aref_sampled->[0]}[$seq_pos];
		
		#######################
		# check site characters of each of the 4 taxa
		if	( ( $allowed_states{$href_hol_seqstates_of_t->{$aref_sampled->[0]}[$seq_pos]} ) &&
			  ( $allowed_states{$href_hol_seqstates_of_t->{$aref_sampled->[1]}[$seq_pos]} ) &&
			  ( $allowed_states{$href_hol_seqstates_of_t->{$aref_sampled->[2]}[$seq_pos]} ) &&
			  ( $allowed_states{$href_hol_seqstates_of_t->{$aref_sampled->[3]}[$seq_pos]} ) ){
			
			$href_red_seq_of_t->{$aref_sampled->[0]} .= $href_hol_seqstates_of_t->{$aref_sampled->[0]}[$seq_pos] ;
			$href_red_seq_of_t->{$aref_sampled->[1]} .= $href_hol_seqstates_of_t->{$aref_sampled->[1]}[$seq_pos] ;
			$href_red_seq_of_t->{$aref_sampled->[2]} .= $href_hol_seqstates_of_t->{$aref_sampled->[2]}[$seq_pos] ;
			$href_red_seq_of_t->{$aref_sampled->[3]} .= $href_hol_seqstates_of_t->{$aref_sampled->[3]}[$seq_pos] 
		}
		#######################
	}
	#######################
	
	#######################
	# seq length
	$$sref_seq_length = length $href_red_seq_of_t->{$aref_sampled->[0]}
	##############################################
	
	#####################################################################
}

sub treefishing_with_icebreaker{
	
	my	$sref_ice_bin		= $_[0] ;	# icebreaker executable				IN: defined; OUT:unchanged
	my	$sref_ice_inmatrix	= $_[1] ;	# matrix infile-name (with path)	IN: defined; OUT:unchanged
	my	$sref_outgr			= $_[2] ;	# defined outgroup-clanname			IN: defined; OUT:unchanged
	
	#print "::&treefishing_with_icebreaker\n";
	
	######################################################################
	# called from subroutine(s)...
	# &print_and_processing_multiclan_support_matrix
	#
	############
	# using subroutine(s)...
	# ...
	######################################################################
	
	######################################################################
	## Process icebreaker on matrix-infile and specified outgroup-clan  ##
	######################################################################
	
	####################################
	# Start tree-fishing amalysis using Penguin_icebreaker
	my	$sysc	=	"./".$$sref_ice_bin." ".$$sref_ice_inmatrix." ".$$sref_outgr ; #print $sysc ;
	
	system ( $sysc );	# bash, if not chmod u+x
	######################################################################
}

sub treefishing_result_processing{
	
	my	$href_value_of_scr_par					= $_[0]	; # key1: sealion user available parameter options; value: optional parameter setting										IN: defined; OUT: unchanged
	my	$href_printpath							= $_[1]	; # key1: R grafic unique printcode (e.g. 'metafileQ1'); key2: *.tsv filename; value: *tsv filename data addings 			IN: defined; OUT: unchanged
	my	$href_hoh_data_of_tsvfile_of_prcode		= $_[2]	; # key: output-type (e.g. 'pdf'); value: pathdir																			IN: defined; OUT: changed
	my	$sref_sf_type							= $_[3]	; # code of quartet-calculation and filter-method	(e.g. 'nap_risk1')														IN: defined; OUT: unchanged
	my	$sref_clanfile							= $_[4]	; # name of given clanfile																									IN: defined; OUT: unchanged
	my	$sref_outfile_prefix					= $_[5]	; # filename-prefix of tree-fishing resultfiles																				IN: defined; OUT: unchanged
	my	$href_hoh_latexdata_of_file_of_prcode	= $_[6]	; # key1: Latex table unique printcode (e.g. 'latexC'); key2: *.tex filename; value: @list of assigned table lines			IN: defined; OUT: changed
	my	$href_textfilename_of_printcode			= $_[7] ; # key1: textfile unique printcode (e.g. 'textA'); value: filename for print												IN: defined; OUT: unchanged
	my	$href_algorithm_def_methods				= $_[8] ; # key1: methode type; key2: method code; value: integer																	IN: defined; OUT: unchanged
	
	########################################################################
	# Processing of Icebreaker result file "*top_topology_scores.txt"    ###
	########################################################################
	
	######################################################################
	# called from subroutine(s)...
	# &main
	#
	############
	# using subroutine(s)...
	# &norm_values
	# &calc_mclan_besttree_support
	# &print_tsv_MQ3
	# &latex_data_sample_LMQ1
	######################################################################
	
	####################################
	# read icebreaker resultfile with the three best supported topologies of the current analysis ( "*top_topology_scores.txt")
	# - extraction of the three best multi-clan trees
	# - support normalization of the three best multi-clan tree support values
	# - final support calculation of overall best multi-clan tree
	# print general icebreaker result file showing the three best multi-clan trees with normalized support
	# print best mean|median icebreaker tree result file, just showing the newick tree
	# print final best multi-clan support of given clanfile and filter approach to .tsv table-file $tsv_4clan_metafileMQ3
	my	$icefolder	=	$href_printpath->{tre} ;
	
	####################################
	# assign penguin-setup parameter 
	# (ts: mean|median; sm: support-method (e.g. nap); fm: filter-method (if used))
	my	@ts_sm_fm	=	split "_", $$sref_sf_type ;
	my	$ts			=	shift @ts_sm_fm ;
	my	$sm			=	shift @ts_sm_fm ;
	my	$fm			;	
	if		(  @ts_sm_fm >= 2 )	{ $fm =  join "_", @ts_sm_fm }
	elsif	(  @ts_sm_fm == 1 )	{ $fm =  $ts_sm_fm[0] }
	else						{ $fm =	'unfiltered' }
	
	############
	# test print
	if ( $href_value_of_scr_par->{sys_parameter}{screen_prt} == 1 ){
		
		print	"\n\ttype_support:\t"	, $ts, 
				"\n\tmethod_support:\t"	, $sm,
				"\n\tmethod_filter:\t"	, $fm, "\n"
	}
	####################################
	
	####################################
	#
	#my	( @best_trees_fm ) ;
	for my $iceresult (<${$sref_outfile_prefix}_*top_topology_scores.txt>){ 
		
		#print "\n\tReading...\t", $iceresult, "\n\n";
		
		####################################
		# Read IN "*top_topology_scores.txt" and assign inlines to @result_lines
		# "*top_topology_scores.txt" consists of three lines, each line highlighting one of the three best multi-clan trees (in descending support order) and corresponding (un-normalized) support (tab delimited)
		open INice,	"<$iceresult" or die "\nICEBREAKER-ERROR: Cannot open ", $iceresult, "!\n\n" ;
		chomp ( my @result_lines = <INice> ) ;
		close INice ;
		####################################
		
		####################################
		# extraction of the three best multi-clan trees and sampling in @support_v
		my @results_l1 = split "\t", $result_lines[0] ;
		my @results_l2 = split "\t", $result_lines[1] ;
		my @results_l3 = split "\t", $result_lines[2] ;	
		
		my @support_v		= ( $results_l1[1], $results_l2[1], $results_l3[1] ) ;
		my @trees			= ( $results_l1[0], $results_l2[0], $results_l3[0] ) ;
		my $best_top_score  = sprintf ( "%.2f", ( $results_l1[1] ) );
		my $dist_best_2nd	= sprintf ( "%.2f", ( $results_l1[1] - $results_l2[1] ));
		
		#push @best_trees_fm, $results_l1[0]."::".$$sref_sf_type ;
		push @{$href_hoh_data_of_tsvfile_of_prcode->{metafileR1}{besttree}}, $results_l1[0]  ;
		push @{$href_hoh_data_of_tsvfile_of_prcode->{metafileR1}{method}}  , $$sref_sf_type  ;
		####################################
		
		####################################
		# support normalization of the three best multi-clan tree support values and sampling in @support_n
		#my @support_n = &norm_values( \@support_v ) ;
		####################################
		
		####################################
		# final support calculation of overall best multi-clan tree
		#my $best_top_score = &calc_mclan_besttree_support ( \@support_n ) ;
		####################################
		
		####################################
		# test print
		if ( $href_value_of_scr_par->{sys_parameter}{screen_prt} == 1 ){
			
			for my $line ( @result_lines ){ print "\traw support:\t", $line	, "\n" }
			print "\tbest tree support:\t ", $best_top_score, "\n\n"
		}
		####################################
		
		####################################
		# print general icebreaker result file showing the three best multi-clan trees with normalized support
		( my $clfile_cp = $$sref_clanfile ) =~ s/.txt$// ;
		
		my	$ice_result_info = $icefolder."/SeaLion_besttree_support_info_".$clfile_cp.".txt" ;
		open OUTtxtA,	">>$ice_result_info" or die "\nICEBREAKER-ERROR: Cannot open ", $ice_result_info, "!\n\n" ;
		
		print OUTtxtA	$$sref_sf_type	, "\tBest Rooted-Clade Tree:\t", $results_l1[0]	, "\tSupport\t", $results_l1[1]	, "\tSupport Distance Best to 2nd-Best Tree\t", $dist_best_2nd, "\n" ;
		
		close OUTtxtA ;
		####################################
		
		####################################
		# print best mean|median icebreaker tree result file, just showing the newick tree
		my	$ice_best_tree = $icefolder."/SeaLion_besttree_".$clfile_cp."_".$$sref_sf_type.".tre" ;
		open OUTtxtB,	">$ice_best_tree" or die "\nICEBREAKER-ERROR: Cannot open ", $ice_best_tree, "!\n\n" ;
		
		print OUTtxtB	$results_l1[0] ;
		close OUTtxtB	;
		####################################
		
		####################################
		# print best tree dataline to tsv table MQ3
		&print_tsv_MQ3(
			
			\%$href_printpath						, # key1: R grafic unique printcode (e.g. 'metafileQ1'); key2: *.tsv filename; value: *tsv filename data addings 						IN: defined; OUT: unchanged
			\%$href_hoh_data_of_tsvfile_of_prcode	, # key: output-type (e.g. 'pdf'); value: pathdir																						IN: defined; OUT: changed
			\$$sref_clanfile						, # name of given clanfile																												IN: defined; OUT: unchanged
			\$sm									, # support method (nap, napv, napr)																									IN: defined; OUT: unchanged
			\$fm									, # filter method (risk1, strict_all:hardclipping, risk2_strict_all:hardclipping)														IN: defined; OUT: unchanged
			\$results_l1[0]							, # best tree in nwk																													IN: defined; OUT: unchanged
			\$ts									, # mean|median																															IN: defined; OUT: unchanged
			\$best_top_score						, # final best tree support																												IN: defined; OUT: unchanged
		) ;
		####################################
		
		####################################
		# sample best tree dataline to latex table LRC1
		&latex_data_sample_LRC1(
			
			\%$href_hoh_latexdata_of_file_of_prcode	, # key1: Latex table unique printcode (e.g. 'latexC'); key2: *.tex filename; value: @list of assigned table lines						IN: defined; OUT: changed
			\$sm									, # support method (nap, napv, napr)																									IN: defined; OUT: unchanged
			\$fm									, # filter method (risk1, strict_all:hardclipping, risk2_strict_all:hardclipping)														IN: defined; OUT: unchanged
			\$results_l1[0]							, # best tree in nwk																													IN: defined; OUT: unchanged
			\$ts									, # mean|median																															IN: defined; OUT: unchanged
			\$best_top_score						, # final best tree support																												IN: defined; OUT: unchanged
			\%$href_textfilename_of_printcode		, # key1: textfile unique printcode (e.g. 'textA'); value: filename for print															IN: defined; OUT: unchanged
			\%$href_printpath						, # key1: R grafic unique printcode (e.g. 'metafileQ1'); key2: *.tsv filename; value: *tsv filename data addings 						IN: defined; OUT: unchanged
		) ;
		####################################
		
		####################################
		# sample best tree dataline to latex table LRC1
		&latex_data_sample_LRC2(
			
			\%$href_hoh_latexdata_of_file_of_prcode	, # key1: Latex table unique printcode (e.g. 'latexC'); key2: *.tex filename; value: @list of assigned table lines						IN: defined; OUT: changed
			\$sm									, # support method (nap, napv, napr)																									IN: defined; OUT: unchanged
			\$fm									, # filter method (risk1, strict_all:hardclipping, risk2_strict_all:hardclipping)														IN: defined; OUT: unchanged
			\@trees									, # best three trees in nwk	sorted by support																							IN: defined; OUT: unchanged
			\$ts									, # mean|median																															IN: defined; OUT: unchanged
			\@support_v								, # sorted best three tree's support																									IN: defined; OUT: unchanged
			\%$href_textfilename_of_printcode		, # key1: textfile unique printcode (e.g. 'textA'); value: filename for print															IN: defined; OUT: unchanged
			\%$href_printpath						, # key1: R grafic unique printcode (e.g. 'metafileQ1'); key2: *.tsv filename; value: *tsv filename data addings 						IN: defined; OUT: unchanged
			\$dist_best_2nd							, # distance support best rooted-clade tree to 2nd best																					IN: defined; OUT: unchanged
		) ;
		####################################
	}
	####################################
}

sub write_p4_script{

	my	$href_val_of_scr_p		= $_[0]	; # key1: sealion user available parameter options; value: optional parameter setting			IN: defined; OUT: unchanged
	my	$href_data_of_prop		= $_[1]	; # key1: property; value : property value														IN: defined; OUT: unchanged
	my	$aref_qtaxa				= $_[2]	; # list of generated 4-taxa combination														IN: defined; OUT: unchanged
	my	$sref_alignment 		= $_[3]	; # generated quartet alignment																	IN: defined; OUT: unchanged
	my	$sref_script_p4			= $_[4]	; # P4 input commands																			IN: undefined; OUT: defined
	my	$aref_qtrees_p4			= $_[5]	; # P4 trees in the same order as defined in P4 input skript									IN: undefined; OUT: defined
	
	######################################################################
	# called from subroutine(s)...
	# &quartet_building_and_p4_main
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	# Write P$ command-script for P4 analysis                            #
	######################################################################
		
	##############################################
	# build unrooted quartet tree
	# (A,B,(C,D)) should make P4 analysis a little bit faster
	my	$p4_tree1	=	"(".$aref_qtaxa->[0].",".$aref_qtaxa->[1].",(".$aref_qtaxa->[2].",".$aref_qtaxa->[3]."))" ;
	my	$p4_tree2	=	"(".$aref_qtaxa->[0].",".$aref_qtaxa->[2].",(".$aref_qtaxa->[1].",".$aref_qtaxa->[3]."))" ;
	my	$p4_tree3	=	"(".$aref_qtaxa->[0].",".$aref_qtaxa->[3].",(".$aref_qtaxa->[1].",".$aref_qtaxa->[2]."))" ;
	
	@$aref_qtrees_p4 = ( $p4_tree1, $p4_tree2, $p4_tree3 ) ;
	##############################################
	
	##############################################
	# ML model parameters
	my	$ncat		= 4		;
	my	$alpha		= $href_val_of_scr_p->{p4_parameter}{alpha} ;
	my	$pinv		= $href_val_of_scr_p->{p4_parameter}{pinv}  ;
	my	$seqtype	= $href_data_of_prop->{type} ;
	
	my	$model	;
	if	( $seqtype eq 'nuc'  ){ $model = $href_val_of_scr_p->{p4_parameter}{model_nu} }
	else					  { $model = $href_val_of_scr_p->{p4_parameter}{model_aa} }
	
	my	%p4_parameters = (

			'gamma_cat'	=>	$ncat										,
			'pinv'		=>	$pinv										,
			'pinvfree'	=>	1											,
			'alpha'		=>	$alpha										,
			'asrv'		=>	"    t.newGdasrv(free=1, val=".$alpha.")"	,
			'model'		=>	$model										,
	) ;
	
	if ( $alpha == 100 ){ ######################################## stehen geblieben

		$p4_parameters{asrv}		=	()		;
		$p4_parameters{pinv}		=	'0.0'	;
		$p4_parameters{pinvfree}	=	'0'		;
		$p4_parameters{gamma_cat}	=	1		;
	}
	##############################################



	##############################################
	# Definement python code ML model parameters
	# JC needs another input code as e.g. GTR
	my	(
			$newcomp	,
			$rmatrix	,
	) ;

	if		( $seqtype	eq	'nuc' ){

		if		( $p4_parameters{model}	eq	'JC' ){

			$newcomp	=	"    t.newComp(free=0, spec=\'equal\')" ;
			$rmatrix	=	"    t.newRMatrix(free=0, spec=\'ones\')" ;
		}

		elsif	( $p4_parameters{model}	eq	'F81' ){

			$newcomp	=	"    t.newComp(free=1, spec=\'empirical\')" ;
			$rmatrix	=	"    t.newRMatrix(free=0, spec=\'ones\')" ;

		}

		elsif	( $p4_parameters{model}	eq	'K2P' ){

			$newcomp	=	"    t.newComp(free=0, spec=\'equal\')" ;
			$rmatrix	=	"    t.newRMatrix(free=1, spec=\'2p\')" ;
		}

		elsif	( $p4_parameters{model}	eq 'HKY' ){

			$newcomp	=	"    t.newComp(free=1, spec=\'empirical\')" ;
			$rmatrix	=	"    t.newRMatrix(free=1, spec=\'2p\')" ;
		}

		elsif	( $p4_parameters{model}	eq	'GTR' ){

			$newcomp	=	"    t.newComp(free=1, spec=\'empirical\')" ;
			$rmatrix	=	"    t.newRMatrix(free=1, spec=\'specified\', val=[2., 3., 4., 5., 6., 7.])" ;
		}

		else{ die "\nBUG-ERROR: Cannot assign defined ", $seqtype ," ML model ", $p4_parameters{model}, " in subroutine &write_p4_script!\nPlease, report BUG to software developer!\n\n" }

		$p4_parameters{seqtype}		=	'dna'	;
		$p4_parameters{npattern}	=	'256'	;
	}

	elsif	( $seqtype	eq	'aa' ){

			$newcomp	=	"    t.newComp(free=1, spec=\'empirical\')" ;
			$rmatrix	=	"    t.newRMatrix(free=0, spec=\'".$model."\')" ;

			$p4_parameters{seqtype}		=	'protein'	;
			$p4_parameters{npattern}	=	'160000'	;
	}

	else	{ die "\nBUG-ERROR: Cannot assign sequence type ", $seqtype, " in subroutine &write_p4_script!\nPlease, report BUG to software developer!\n\n" }
	##############################################



	##############################################
	# open p4 input script
	# print script code
	$$sref_script_p4	= "p4_temp_script.py" ;
	open	OUTp4,	">$$sref_script_p4" or die "\nP4-ERROR: Cannot open p4 inputfile ", $$sref_script_p4,"!\n\n" ;

	print	OUTp4	"var.warnReadNoFile = False\n",
					"read(\'", $$sref_alignment  ,"\')\n",
					"a = var.alignments[0]\n" ,
					"d = Data([a])\n" ,
					"aLength = len( var.alignments[0] )\n" ,
					"print (aLength)\n" ,
					"\n" ,
					"var.PIVEC_MIN = 1.e-5\n",
					"var.RATE_MIN = 1.e-5\n",
					"var.BRLEN_MIN = 1.e-5\n",
					"var.GAMMA_SHAPE_MIN = 0.3\n",
					"\n",
					"myTreeString = \"\"\"\n" ,
					$p4_tree1, ";\n" ,
					$p4_tree2, ";\n" ,
					$p4_tree3, ";\n" ,
					"\"\"\"\n" ,
					"read(myTreeString)\n" ,
					"\n" ,
					"for t in var.trees:\n" ,
					"    \n" ,
					"    t.write()\n",
					"    print()\n",
					"\n" ,
					"    t.taxNames = a.taxNames\n" ,
					"    t.data = d\n" ,
					"\n" ,
					$newcomp, "\n" ,
					$rmatrix, "\n" ,
					"    t.setNGammaCat(nGammaCat=",								$p4_parameters{gamma_cat}	, ")\n" ,
					$p4_parameters{asrv}, "\n",
					"    t.setPInvar(free=", $p4_parameters{pinvfree},", val=",		$p4_parameters{pinv}		, ")\n" ,
					"\n" ,
					"    t.optLogLike()\n" ,
					"    a2 = func.newEmptyAlignment(dataType=\'", $p4_parameters{seqtype} ,"\', taxNames=a.taxNames, length=", $p4_parameters{npattern}, ")\n" ,
					"\n" ,
					"    for s in a2.sequences:\n" ,
					"        s.sequence = list(s.sequence)\n" ,
					"\n" ,
					"    posn = 0\n" ,
					"    nn = list(a2.symbols)\n" ,
					"    for p0 in nn:\n" ,
					"        for p1 in nn:\n" ,
					"            for p2 in nn:\n" ,
					"                for p3 in nn:\n" ,
					"                    a2.sequences[0].sequence[posn] = p0\n" ,
					"                    a2.sequences[1].sequence[posn] = p1\n" ,
					"                    a2.sequences[2].sequence[posn] = p2\n" ,
					"                    a2.sequences[3].sequence[posn] = p3\n" ,
					"                    posn += 1\n" ,
					"\n" ,
					"    for s in a2.sequences:\n" ,
					"        s.sequence = ''.join(s.sequence)\n" ,
					"\n" ,
					"    d2 = Data([a2])\n" ,
					"    t.data = d2\n" ,
					"    t.calcLogLike()\n" ,
					"    t.getSiteLikes()\n" ,
					"    print (sum(t.siteLikes))\n" ,
					"    sDict = {}\n" ,
					"    sList = []\n" ,
					"\n" ,
					"    for posn in range(", $p4_parameters{npattern}, "):\n" ,
					"        sl = ''.join(a2.sequenceSlice(posn))\n" ,
					"        sDict[sl] = [aLength * t.siteLikes[posn], 0]\n" ,
					"        sList.append(sl)\n" ,
					"\n" ,
					"    for posn in range(aLength):\n" ,
					"        sl = ''.join(a.sequenceSlice(posn))\n" ,
					"        sDict[sl][1] += 1\n" ,
					"\n" ,
					"    t.write()\n",
					"    print()\n",
					"    for sl in sList:\n" ,
					"        print (sl, \" %9i\" % sDict[sl][1], \" %11.2f\" % sDict[sl][0])\n" ,
					"sys.stdout.flush()\n" ,
	;

	close	OUTp4;
	##############################################
	
	######################################################################
}

sub write_spd_file{
	
	my $href_value_of_scr_par					= $_[0]  ; # key1: sealion user available parameter options; value: optional parameter setting				IN: defined; OUT: changed
	my $href_subfolder_q						= $_[1]	 ; # key: output-type (e.g. 'pdf'); value: pathdir													IN: defined; OUT: unchanged
	my $href_hoh_found_N_of_topo_of_pattern		= $_[2]  ; # key1: topology; key2: recoded site pattern; value: N site pattern in original data				IN: defined; OUT: unchanged
	my $href_hoh_expected_N_of_topo_of_pattern	= $_[3]  ; # key1: topology; key2: recoded site pattern; value: expected ML N site pattern					IN: defined; OUT: unchanged
	my $sref_outfile_p4							= $_[4]	 ; # P4 resultfile																					IN: defined; OUT: unchanged
	my $aref_qtrees_p4							= $_[5]	 ; # P4 trees in the same order as defined in P4 input skript										IN: defined; OUT: unchanged
	my $sref_tquartet							= $_[6]	 ; # 4 taxon-quartet names in alphab. order															IN: defined; OUT: unchanged
	my $sref_clan4								= $_[7]	 ; # actual 4clan combination																		IN: defined; OUT: unchanged
	my $sref_qcounter							= $_[8]	 ; # quartet internal number																		IN: defined; OUT: unchanged
	my $sref_seq_length							= $_[9]	 ; # common length of reduced quartet-sequences														IN: defined; OUT: unchanged
	my $href_clan_of_sequence					= $_[10] ; # key1: sequence-name; value: clan definition													IN: defined; OUT: unchanged
	
	######################################################################
	# called from subroutine(s)...
	# &quartet_building_and_p4_main
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	# Print P4 analysed split-pattern distribution to new file           #
	######################################################################
	
	##############################################
	# define and make new outfile pathdir
	my		@clan4_uncoded	= split ":", $$sref_clan4   ;
			@clan4_uncoded	= sort 		 @clan4_uncoded ;
	my		$clan4_uncoded	= join  "_", @clan4_uncoded ;
	
	my		$spdpath		= $href_subfolder_q->{spd}."/".$clan4_uncoded ;
	mkdir 	$spdpath 		;
	##############################################
	
	##############################################
	# open new spd-file
	my		$outfile = "SeaLion_detailed_split_calc_q".$$sref_qcounter.".txt" ;
	my	 	$spdfile = $spdpath."/".$outfile ;
	open 	OUTspd, ">$spdfile" || die "\nSPDFILE-ERROR:Cannot write spd file ", $spdfile, "!\n\n" ;
	##############################################
	
	##############################################
	# write header lines
	print 	OUTspd	"Q", $$sref_qcounter, ":\t", $$sref_tquartet ,
					"\n---------------------------------------------------------------"	,
					"\nAlignment observed and P4 expected pattern frequencies" ,
					"\n\nPattern\t\t\tObserved\tExp. ", $aref_qtrees_p4->[0], "\tExp. ", $aref_qtrees_p4->[1], "\tExp. ", $aref_qtrees_p4->[2] ;
	##############################################
	
	##############################################
	## print observed frequencies
	my @otus = split ":", $$sref_tquartet ;
	
	for my $pattern ( qw/A B C D E F G H J K L M N I V/ ){
		
		print OUTspd	"\nPattern\t"	, $pattern, 
						"\t:\t"			, $href_hoh_found_N_of_topo_of_pattern	 ->{$aref_qtrees_p4->[0]}{$pattern}	,
						"\t\t"			, $href_hoh_expected_N_of_topo_of_pattern->{$aref_qtrees_p4->[0]}{$pattern}	,
						"\t\t\t"		, $href_hoh_expected_N_of_topo_of_pattern->{$aref_qtrees_p4->[1]}{$pattern}	,
						"\t\t\t"		, $href_hoh_expected_N_of_topo_of_pattern->{$aref_qtrees_p4->[2]}{$pattern} ;
	}
	
	print	OUTspd	"\n---------------------------------------------------------------\n" 				,
					"\nSequence length of split informative site pattern (bp):\t", $$sref_seq_length 	,
					"\nClan assignment\t", $otus[0],"\t", $href_clan_of_sequence->{$otus[0]}			,
					"\nClan assignment\t", $otus[1],"\t", $href_clan_of_sequence->{$otus[1]}			,
					"\nClan assignment\t", $otus[2],"\t", $href_clan_of_sequence->{$otus[2]}			,
					"\nClan assignment\t", $otus[3],"\t", $href_clan_of_sequence->{$otus[3]} 			;
	
	close OUTspd ;
	##############################################
	
	##############################################
	# set path of input P4 results to main P4 print dir
	$href_value_of_scr_par->{sys_parameter}{i_main} = $href_subfolder_q->{spd} ;
	##############################################
	#exit;
	######################################################################
}

sub help{

	my	$sref_command			= $_[0]	; #	help print of specific command. if command = 'h' print everything		IN: defined; OUT: unchanged
	my	$href_value_of_scr_par	= $_[1] ; # key1: method_code; key2: method_option; value: float or integer			IN: defined; OUT: unchanged
	
	######################################################################
	# called from subroutine(s)...
	# &argv_handling
	# &check_msa_infile
	# &readin_fasta
	# &readin_phylip
	# &check_clan_infile
	# &evaluation_quartet_ratio_NcNa
	# &generate_clan_quartets
	# &clan4_extraction_from_filepaths
	# 
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	# print HELP menu & EXIT                                             #
	######################################################################
	my	$script_name = $href_value_of_scr_par->{sys_parameter}{pen_bin} ;
	my	$website_P4	 = $href_value_of_scr_par->{author_info}{p4} ;
	
	if		( $$sref_command eq 'h' ){
		
		print	"\n\t------------------------------------------------------------",
				"\n\tSeaLion Help Menu",
				"\n\t------------------------------------------------------------",
				"\n\tUsage (Linux):" ,
				"\n\tperl ", $script_name, " -i [msa infile] -p [clade infile] -o [outgroup] -[optional_command(s)] -s",
				"\n",
				"\n\t [i] <string>\t\tMultiple sequence alignment '*.phy' or '*.fas' or '*.fasta' (opt. with path DIR)",
				"\n\t [p] <string>\t\tClade-definition infile '*.txt' (opt. with path DIR)",
				"\n\t [o] <string>\t\tOutgroup defined clade-code of the clade-definition infile",
				"\n\t [s] \t\t\tStart-command, executing the process run",
				"\n",
				"\n\t[optional_command]",
				"\n\t",
				"\n\t [M] <integer>\t\tMinimum number of character-complete site positions in species-quartets (default: ", $href_value_of_scr_par->{menu_parameter}{minlengths}, " bp)",
				"\n\t [l] <integer>\t\tMaximum number of single species-quartets for each clade-quartet (default: ", $href_value_of_scr_par->{menu_parameter}{maxquartets}, " quartets)",
				"\n\t",
				"\n\t [r] <x>\t\tSpecies-quartet FILTER 'RISK' deactivation (default: ", $href_value_of_scr_par->{menu_parameter}{risk_default2} ,")",
				"\n\t [r] <integer>\t\t-Selection of 'RISK'-filter analyses (default: ", $href_value_of_scr_par->{menu_parameter}{risk_default} ,")",
				"\n\t [tlrisk] <float>\t-Lower (internally optimized) limit of species-quartet filter 'RISK' (default: ", $href_value_of_scr_par->{menu_parameter}{lowerlimit} ,")",
				"\n\t [turisk] <float>\t-Upper (fix) limit of species-quartet filter 'RISK' (default: ", $href_value_of_scr_par->{menu_parameter}{upperlimit}, ")",
				"\n\t [tsrisk] <float>\t-Scale steps from lower to upper limit of species-quartet filter 'RISK' (default: +", $href_value_of_scr_par->{menu_parameter}{rscale} ,")",
				"\n\t",
				"\n\t [d]\t\t\tSpecies-quartet FILTER 'DIST' - activation|deactivation (default: ", $href_value_of_scr_par->{menu_parameter}{qdist_default} ,")",
				"\n\t [tldist] <float>\t-Lower (fix) limit of species-quartet filter 'DIST' (default: ", $href_value_of_scr_par->{menu_parameter}{lowerlimit_d} ,")",
				"\n\t [tudist] <float>\t-Upper (internally optimized) limit of species-quartet filter 'DIST' (default: ", $href_value_of_scr_par->{menu_parameter}{upperlimit_d} ,")",
				"\n\t [tsdist] <float>\t-Scale steps from upper to lower limit of species-quartet filter 'DIST' (default: -", $href_value_of_scr_par->{menu_parameter}{dscale} ,")",
				"\n\t",
				"\n\t [m] <string>\t\tML substitution model NUC (default: ", $href_value_of_scr_par->{menu_parameter}{model_nu} ,") or AA (default: ", $href_value_of_scr_par->{menu_parameter}{model_aa} ,")",
				"\n\t [a] <float>\t\t-Start alpha-shape value for ML estimation (default: ", $href_value_of_scr_par->{menu_parameter}{start_alpha} ,")",
				"\n\t [I] <float>\t\t-Start proportion invariable sites for ML estimation (default: ", $href_value_of_scr_par->{menu_parameter}{start_pinv} ,")",
				"\n\t",
				"\n\t [restart]\t\tActivation|deactivation of a split-pattern re-analysis of existing SPD files (default: ", $href_value_of_scr_par->{menu_parameter}{restart_setting}, ")",
				"\n\t [imain] <string>\tMain SPD infilefolder for '-restart' option, optionally with path DIR (default: ", $href_value_of_scr_par->{menu_parameter}{spd_folder} ,"/*)",
				"\n\t [omain] <string>\tMain (new) resultfolder, optionally with path DIR (default: ", $href_value_of_scr_par->{menu_parameter}{output_folder} ,"/*)",
				"\n\t",
				"\n\t [average] <integer>\tAgrregate measure for single clade-quartet-, and rooted-clade tree support (default: ", $href_value_of_scr_par->{menu_parameter}{average}, ")",
				"\n\t [ry]\t\t\tActivation|deactivation of RY coding of character states (default: ", $href_value_of_scr_par->{menu_parameter}{ry_setting}, ")",
				"\n\t [prt] <integer>\tActivation|deactivation of additional grafic and table prints (default: ", $href_value_of_scr_par->{menu_parameter}{prt_setting}, ")",
				"\n\t [u]\t\t\tActivation|deactivation of script queries (default: ", $href_value_of_scr_par->{menu_parameter}{query_setting}, ")",
				"\n\t",
				"\n\tFor a detailed help menu about [optional_command] parameters type:",
				"\n\tperl ", $script_name, " -h [optional_command], e.g.:",
				"\n\tperl ", $script_name, " -h i",
				"\n\t",
				"\n\tFor license information type:",
				"\n\tperl ", $script_name, " -P",
				"\n\t",
				"\n\t", $script_name, " requires the P4 Python Packages",
				"\n\tFor download and P4 install instructions visit: ",
				"\n\t", $website_P4,
				"\n\t---------------------------------------------------\n\n",
		;
	}
	
	elsif	( $$sref_command eq 'a' ){

		print	"\n\t------------------------------------------------------------",
				"\n\tML Alpha Shape Start Parameter:" ,
				"\n\t-----------------------------------",
				"\n\tAdditional parameter:",
				"\n\t- Float number 0.1 to 100.0",
				"\n\t(default: ", $href_value_of_scr_par->{menu_parameter}{start_alpha},")",
				"\n\t",
				"\n\tMaximum-Likelihood (ML) alpha shape parameter for rate heterogeneity.",
				"\n\tStart parameter for  ML estimation of potentially convergently",
				"\n\tevolved split pattern frequencies.",
				"\n\t",
				"\n\tNOTE: For ML pattern estimation without ASRV type:",
				"\n\t\t\"-a 100\"",
				"\n\t",
				"\n\tSpecified via '-a' option:",
				"\n\t\t\"e.g. -a 1.0\"",
				"\n\t---------------------------------------------------\n\n",
		;
	}
	
	elsif	( $$sref_command eq 'average' ){

		print	"\n\t------------------------------------------------------------",
				"\n\tAverage Measure of Tree Support:" ,
				"\n\t-----------------------------------",
				"\n\tAdditional parameter:",
				"\n\t- Positive integer number: 1 to 3",
				"\n\t(default: '", $href_value_of_scr_par->{menu_parameter}{average}, "')",
				"\n\t",
				"\n\tSelectable parameter:",
				"\n\t\t1\tMEAN\t\t('mean')",
				"\n\t\t2\tMEDIAN\t\t('median')",
				"\n\t\t3\tMEAN+MEDIAN\t('mean,median')",
				"\n\t",
				"\n\tThe '-average' option defines the average measure of single species-quartet",
				"\n\tsupport for individual clade-quartet trees. This measure is used as the basis", 
				"\n\tfor final rooted-clade tree support. If the combined option ('-average 3')",
				"\n\tis selected, both measures are separately used as a measure of support, with",
				"\n\tindividual results printed for each of the two measures.",
				"\n\t",
				"\n\tSpecified via '-average' option:",
				"\n\t\t\"e.g. -average 1\"",
				"\n\t---------------------------------------------------\n\n",
		;
	}
	
	elsif	( $$sref_command eq 'd' ){

		print	"\n\t------------------------------------------------------------",
				"\n\tSpecies-Quartet Filter 'DIST' Activation:" ,
				"\n\t-------------------------------------------",
				"\n\tAdditional parameter:",
				"\n\t- No extra parameter necessary",
				"\n\t(default: ", $href_value_of_scr_par->{menu_parameter}{qdist_default}, ")",
				"\n\t",
				"\n\tThe '-d' option activates or deactivates (if already",
				"\n\tactivated) the species-quartet filter 'DIST' approach.", 
				"\n\tThis approach rejects quartets in single clade-quartets if",
				"\n\tcorresponding tree support difference between the best and",
				"\n\tsecond-best (species-quartet related) tree is below",
				"\n\tthe upper (optimized) 'DIST' threshold.",
				"\n\t",
				"\n\tSpecified via '-d' option:",
				"\n\t\t\"-d\"",
				"\n\t---------------------------------------------------\n\n",
		;
	}
	
	elsif	( $$sref_command eq 'i' ){

		print	"\n\t------------------------------------------------------------",
				"\n\tMultiple Sequence Alignment Infile:" ,
				"\n\t-----------------------------------",
				"\n\tAdditional parameter:",
				"\n\t- Name of the alignment infile, optionally with path DIR",
				"\n\t\t- FASTA  ('*.fas' or '*.fasta')",
				"\n\t\t- PHYLIP ('*.phy') (strict or relaxed)",
				"\n\t",
				"\n\tAllowed site conditions:",
				"\n\t\t- Sequences of equal length",
				"\n\t\t- Nucleotide states",
				"\n\t\t- Amino-acid states",
				"\n\t\t- Ambiguity states",
				"\n\t\t- Indel/GAP states ('-')",
				"\n\t",
				"\n\tAllowed sequence names:",
				"\n\t\t- Alphanumeric signs",
				"\n\t\t- Underscores ('_')",
				"\n\t\t- Only unique sequence names",
				"\n\t",
				"\n\tSpecified via '-i' option:",
				"\n\t\t\"e.g. -i MSA_infile_name.fas\"",
				"\n\t---------------------------------------------------\n\n",
		;
	}
	
	elsif	( $$sref_command eq 'imain' ){

		print	"\n\t------------------------------------------------------------",
				"\n\tMain SPD Infile PathDIR:" ,
				"\n\t-----------------------------------",
				"\n\tAdditional parameter:",
				"\n\t- String: Alphanumeric+Underscore ",
				"\n\t(default: '",$href_value_of_scr_par->{menu_parameter}{spd_folder} ,"'/*)",
				"\n\t",
				"\n\tPath to the main folder of already analyzed split-pattern",
				"\n\tdistribution (SPD) infiles. This parameter is used for starting",
				"\n\ta re-analysis (if '-restart' option is activated) without",
				"\n\tre-processing P4 of MSA observed and ML expected site pattern",
				"\n\tfrequencies, which is one of the most time-consuming",
				"\n\tprocesses in the overall analysis.",
				"\n\t",
				"\n\tSpecified via '-imain' Option:",
				"\n\t\t\"e.g. -imain MySPDfile_resultFolder_1\"",
				"\n\t---------------------------------------------------\n\n",
		;
	}
	
	elsif	( $$sref_command eq 'I' ){

		print	"\n\t------------------------------------------------------------",
				"\n\tML pINV Start Parameter:" ,
				"\n\t-----------------------------------",
				"\n\tAdditional parameter:",
				"\n\t- Float number between 0 and 1",
				"\n\t(default: ", $href_value_of_scr_par->{menu_parameter}{start_pinv} ,")",
				"\n\t",
				"\n\tML proportion of invariable site estimation. Start parameter", 
				"\n\tfor ML estimation of potentially convergently",
				"\n\tevolved split pattern frequencies.",
				"\n\t",
				"\n\tSpecified via '-I' option:",
				"\n\t\t\"e.g. -I 0.15\"",
				"\n\t---------------------------------------------------\n\n",
		;
	}
	
	elsif	( $$sref_command eq 'l' ){

		print	"\n\t------------------------------------------------------------",
				"\n\tMaximum Number of analysed Species-Quartets per Clade-Quartet:" ,
				"\n\t------------------------------------------------------------",
				"\n\tAdditional parameter:",
				"\n\t- Positive integer number",
				"\n\t(default: '", $href_value_of_scr_par->{menu_parameter}{maxquartets} ,"' bp)",
				"\n\t",
				"\n\tMaximum number of single analyzed species-quartets in each",
				"\n\tclade-quartet. The total computation time can get very", 
				"\n\tlong if the number of speies-quartets to be analyzed is large.",
				"\n\tDefining a maximum of single species-quartet analyses in each",
				"\n\tclade-quartet can help to reduce computation time.",
				"\n\t",
				"\n\tNOTE: If the number of possible species-quartets in a clade-",
				"\n\tquartet exceeds the maximum number of allowed species-quartet", 
				"\n\tanalyses, the process run will stop with a command line request.",
				"\n\tThe request can be suppressed by the '-u' option. In that case,",
				"\n\tspecies-quartets are drawn randomly from the overall pool of available,",
				"\n\tclade-quartet corresponding species-quartets.",
				"\n\t",
				"\n\tSpecified via '-l' option:",
				"\n\t\t\"e.g. -l 5000\"",
				"\n\t---------------------------------------------------\n\n",
		;
	}
	
	elsif	( $$sref_command eq 'm' ){

		print	"\n\t------------------------------------------------------------",
				"\n\tML Substitution Model of P4 Analysis:" ,
				"\n\t-----------------------------------",
				"\n\tAdditional parameter:",
				"\n\t- Positive integer number: 1 to 11",
				"\n\t",
				"\n\t\tdefault NUC: "	, $href_value_of_scr_par->{menu_parameter}{model_nu} ,
				"\n\t\tdefault AA: "	, $href_value_of_scr_par->{menu_parameter}{model_aa} ,
				"\n\t",
				"\n\tSelectable parameter:",
				"\n\t\t1\tGTR\t\t(NUC)",
				"\n\t\t2\tHKY\t\t(NUC)",
				"\n\t\t3\tK2P\t\t(NUC)",
				"\n\t\t4\tF81\t\t(NUC)",
				"\n\t\t5\tJC\t\t(NUC)" ,
				"\n\t\t6\twag\t\t(AA)",
				"\n\t\t7\tjtt\t\t(AA)",
				"\n\t\t8\td78\t\t(AA)",
				"\n\t\t9\tmtrev24\t\t(AA)",
				"\n\t\t10\tlg\t\t(AA)",
				"\n\t\t11\tblosum62\t(AA)",
				"\n\t",
				"\n\tSubstitution model used for Maximum Likelihood (ML)",
				"\n\testimation of potentially convergently evolved split",
				"\n\tpattern frequencies.",
				"\n\t",
				"\n\tSpecified via '-m' option:",
				"\n\t\t\"e.g. -m 1\"",
				"\n\t---------------------------------------------------\n\n",
		;
	}
	
	elsif	( $$sref_command eq 'M' ){

		print	"\n\t------------------------------------------------------------",
				"\n\tMinimum Sequence Length for Single Species-Quartets:" ,
				"\n\t------------------------------------------------------------",
				"\n\tAdditional parameter:",
				"\n\t- Positive integer number",
				"\n\t(default: ", $href_value_of_scr_par->{menu_parameter}{minlengths} ,")",
				"\n\t",
				"\n\tOnly character-complete site positions in single species-quartets are",
				"\n\tanalyzed. Character-incomplete sites, including e.g. gaps",
				"\n\tor ambiguities, are rejected from each species-quartet analysis.",
				"\n\tSpecies-Quartets are completely rejected if the number of remaining sites is",
				"\n\tlower than specified by the '-l' parameter.",
				"\n\t",
				"\n\tSpecified via -M' option:",
				"\n\t\t\"e.g. -M 2000\"",
				"\n\t---------------------------------------------------\n\n",
		;
	}
	
	elsif	( $$sref_command eq 'o' ){

		print	"\n\t------------------------------------------------------------",
				"\n\tOutgroup Definition:" ,
				"\n\t-----------------------------------",
				"\n\tAdditional parameter:",
				"\n\t- String: Alphanumeric+Underscore",
				"\n\t",
				"\n\tClade-code of the outgroup as defined in the analysis",
				"\n\tcorresponding clade-definition infile ('-p' option).",
				"\n\tThe clade-code of the outgroup must be identical to one of",
				"\n\tthe defined clades of the clade-definition infile (case",
				"\n\tsensitive).",
				"\n\t",
				"\n\tSpecified via '-o' option:",
				"\n\t\t\"e.g. -o O\"",
				"\n\t---------------------------------------------------\n\n",
		;
	}
	
	elsif	( $$sref_command eq 'omain' ){

		print	"\n\t------------------------------------------------------------",
				"\n\tMain Outfile PathDIR:" ,
				"\n\t-----------------------------------",
				"\n\tAdditional parameter:",
				"\n\t- String: Alphanumeric+Underscore",
				"\n\t(default: '", $href_value_of_scr_par->{menu_parameter}{output_folder} ,"'/*)",
				"\n\t",
				"\n\tName of the new result output folder, with path DIR if",
				"\n\tthe defined output folder is not located in the script DIR.",
				"\n\tWith '-restart', SPD files of 'imain' are copied to 'omain'.", 
				"\n\tFurthermore, already pre-existing result files in 'omain'",
				"\n\tare always deleted in advance of a new analysis.",
				"\n\t",
				"\n\tSpecified via '-omain' option:",
				"\n\t\t\"e.g. -omain NewResultFolder\"",
				"\n\t---------------------------------------------------\n\n",
		;
	}
	
	elsif	( $$sref_command eq 'p' ){

		print	"\n\t------------------------------------------------------------",
				"\n\tClade-Definition Infile:" ,
				"\n\t-----------------------------------",
				"\n\tAdditional parameter:",
				"\n\t- Name of the clade-definition infile, optionally with path DIR",
				"\n\t- Plain TEXT: ('*.txt')",
				"\n\t",
				"\n\tClade definition:",
				"\n\t\t- At least four clades must be defined in separate lines",
				"\n\t\t- Clades defined by first code in each Line",
				"\n\t\t- Only uniquely defined clade-codes allowed",
				"\n\t\t- Only alphanumeric signs and underscore(s) allowed" ,
				"\n\t",
				"\n\tSequence assignment to clades:",
				"\n\t\t- Defined clade follows assigned sequence names (comma separated)",
				"\n\t\t- Must be in the same line as assigned clade-code",
				"\n\t\t- No whitespace allowed",
				"\n\t\t- Only unique sequence names allowed",
				"\n\t\t- Sequence names must defined as in the original alignment (case sensitive)",
				"\n\t\t- Only alphanumeric signs and underscore(s) are allowed for sequence names" ,
				"\n\t",
				"\n\tExample format:",
				"\n\t\t--------",
				"\n\t\tline_1: CladeCode1,SeqName1,SeqName2,SeqName3...<linebreak>",
				"\n\t\tline_2: CladeCode2,SeqNameA,SeqNameB,SeqNameC...<linebreak>",
				"\n\t\tline_3: CladeCode3,SeqName4,SeqName5,SeqName6...<linebreak>",
				"\n\t\tline_4: CladeCode4,SeqNameD,SeqNameE,SeqNameF...<linebreak>",
				"\n\t\t...",
				"\n\t\tline_n: CladeCode4,SeqNameD,SeqNameE,SeqNameF...<linebreak>",
				"\n\t\t--------",
				"\n\t",
				"\n\tSpecified via '-p' option:",
				"\n\t\t\"e.g. -p clade_infile_name.txt\"",
				"\n\t---------------------------------------------------\n\n",
		;
	}
	
	elsif	( $$sref_command eq 'r' ){

		print	"\n\t------------------------------------------------------------",
				"\n\tSpecies-Quartet Filter 'RISK' Analysis:" ,
				"\n\t-------------------------------------------",
				"\n\tAdditional parameter:",
				"\n\t- Positive integer number: 1 to 15 or 'x' for deactivation",
				"\n\t(default: '", $href_value_of_scr_par->{menu_parameter}{risk_default},"')",
				"\n\t",
				"\n\tSelectable parameter:",
				"\n\t\tx\tDeactivate 'RISK' Filter\t('disabled')",
				"\n\t\t1\tRISK1\t\t\t\t('risk1')",
				"\n\t\t2\tRISK2\t\t\t\t('risk2')",
				"\n\t\t3\tRISK1+RISK2\t\t\t('risk1,risk2')",
				"\n\t",
				"\n\tActivates or deactivates parameter corresponding",
				"\n\t(depending on the parameter) 'RISK' filter settings to",
				"\n\tfilter single species-quartets based on their ratio of",
				"\n\tpotentially convergent- (Nc), and apomorphic (Na) evolved split",
				"\n\tsignal among species-quartet related site positions (ratio Nc/Na).",
				"\n\tRISK1 and RISK2 are based on the original, uncorrected",
				"\n\tratio of Nc/Na.",
				"\n\t",
				"\n\tSpecified via '-r' option:",
				"\n\t\t\"e.g. -r 3\"",
				"\n\t---------------------------------------------------\n\n",
		;
	}
	
	elsif	( $$sref_command eq 'ry' ){

		print	"\n\t------------------------------------------------------------",
				"\n\tRY-Coding of Site Characters:",
				"\n\t-----------------------------------",
				"\n\tAdditional parameter:",
				"\n\t- No extra parameter necessary,",
				"\n\t(default: ", $href_value_of_scr_par->{menu_parameter}{ry_setting},")",
				"\n\t",
				"\n\tWith the '-ry' option activated, site characters are", 
				"\n\ttranslated into R or Y states, depending on their",
				"\n\tchemical behaviour. This reduces the number of possible", 
				"\n\tcharacter states from four to two. An activated '-ry'",
				"\n\toption can be deactivated by typing '-ry' again.",
				"\n\t",
				"\n\tNucleotide: R-> Purine, Y-> Pyrimidin",
				"\n\tAmino-acid: R-> Hydrophobic, Y-> Hydrophilic",
				"\n\t",
				"\n\tIf RY coding is disabled, type '-ry' to activate RY coding",
				"\n\tIf RY coding is activated, type '-ry' to disable RY coding",
				"\n\t",
				"\n\tSpecified via '-ry' option:",
				"\n\t\t\"-ry\"",
				"\n\t---------------------------------------------------\n\n",
		;
	}
	
	elsif	( $$sref_command eq 'restart' ){

		print	"\n\t------------------------------------------------------------",
				"\n\tRestart without SPD Re-Calculation:" ,
				"\n\t-----------------------------------",
				"\n\tAdditional parameter:",
				"\n\t- No extra parameter necessary",
				"\n\t(default: ", $href_value_of_scr_par->{menu_parameter}{restart_setting}, ")",
				"\n\t",
				"\n\tThe '-restart' option allows a re-analysis of a previously",
				"\n\tconducted SeaLion analysis by already using existing split-pattern",
				"\n\tdistribution files (SPD) of a former analysis (defined via 'imain').",
				"\n\tSPD files of 'imain' are copied to the new result output path",
				"\n\t('omain') after deleting pre-existing 'omain' result files if", 
				"\n\tSPD infile path ('imain') and output path ('omain') are not",
				"\n\tidentical. An activated '-restart' option can be deactivated",
				"\n\tby typing '-restart' again.",
				"\n\t",
				"\n\tSpecified via '-restart' Option:",
				"\n\t\t\"-restart\"",
				"\n\t---------------------------------------------------\n\n",
		;
	}
	
	elsif	( $$sref_command eq 'tudist' ){

		print	"\n\t------------------------------------------------------------",
				"\n\tSpecies-Quartet Filter 'DIST' (upper) Threshold:" ,
				"\n\t------------------------------------------------------------",
				"\n\tAdditional parameter:",
				"\n\t- Positive float number: 0.0 to 1.0",
				"\n\t(default: ", $href_value_of_scr_par->{menu_parameter}{upperlimit_d} ,")",
				"\n\t",
				"\n\tStart value for the upper 'DIST' threshold optimization.",
				"\n\tThe upper 'DIST' threshold is stepwise decreased (following the", 
				"\n\t'tsdist' defined scaling) until either an optimum of rejected",
				"\n\tspecies-quartets is found, or the lower 'DIST' threshold is reached.",
				"\n\tThe optimized upper threshold scale is subsequently used as",
				"\n\tthe lower limit of allowed clade-quartet respective tree",
				"\n\tsupport distance between the best and second-best quartet-tree.",
				"\n\t",
				"\n\tSpecified via '-tudist' option:",
				"\n\t\t\"e.g. -tudist 0.3\"",
				"\n\t---------------------------------------------------\n\n",
		;
	}
	
	elsif	( $$sref_command eq 'tldist' ){

		print	"\n\t------------------------------------------------------------",
				"\n\tSpecies-Quartet Filter 'DIST' (lower) Threshold:" ,
				"\n\t------------------------------------------------------------",
				"\n\tAdditional parameter:",
				"\n\t- Positive float number: 0.0 to 1.0",
				"\n\t(default: ", $href_value_of_scr_par->{menu_parameter}{lowerlimit_d} ,")",
				"\n\t",
				"\n\tLower (fix) threshold limit for 'DIST' species-quartet filtering,",
				"\n\tused for upper 'DIST' threshold ('tudist') optimization. During",
				"\n\tthe threshold optimization, the upper 'DIST' threshold is",
				"\n\tstepwise decreased (following the 'tsdist' defined scaling)",
				"\n\tuntil either an optimum of rejected species-quartets is found, or the",
				"\n\tlower 'DIST' threshold is reached. The optimized upper threshold", 
				"\n\tscale. The adjustment of the lower threshold limit enables a", 
				"\n\tmore flexible range for threshold optimization. Instead of",
				"\n\tsearching for an optimum between 0 and the upper threshold,",
				"\n\tthe lower threshold limit can be shifted too, e.g., from 0 to",
				"\n\t0.05, excluding an optimization analysis of quartet tree",
				"\n\tdistances lower than 0.05.",
				"\n\t",
				"\n\tNote: The lower threshold",
				"\n\thas to be less/equal the upper threshold",
				"\n\tlimit ('tldist' <= 'tudist').",
				"\n\t",
				"\n\tSpecified via '-tldist' option:",
				"\n\t\t\"e.g. -tldist 0.01\"",
				"\n\t---------------------------------------------------\n\n",
		;
	}
	
	elsif	( $$sref_command eq 'tsdist' ){

		print	"\n\t------------------------------------------------------------",
				"\n\tSpecies-Quartet Filter 'DIST' Threshold Scaling:" ,
				"\n\t------------------------------------------------------------",
				"\n\tAdditional parameter:",
				"\n\t- Positive float number: 0.0 to 1.0",
				"\n\t(default: -", $href_value_of_scr_par->{menu_parameter}{dscale} ,")",
				"\n\t",
				"\n\tNegative scale value for upper DIST threshold optimization.",
				"\n\tStarting from the defined upper threshold value ('tudist'),",
				"\n\tthe number of 'DIST' rejected species-quartets is evaluated",
				"\n\tfor each negative scale step until either an optimum of",
				"\n\tfiltered/unfiltered quartets is found, or the lower threshold",
				"\n\tvalue ('tldist') is reached (whose value is then used",
				"\n\tas threshold criteria).",
				"\n\t",
				"\n\tSpecified via '-tsdist' option:",
				"\n\t\t\"e.g. -tsdist 0.05\"",
				"\n\t---------------------------------------------------\n\n",
		;
	}
	
	elsif	( $$sref_command eq 'tlrisk' ){

		print	"\n\t------------------------------------------------------------",
				"\n\tSpecies-Quartet Filter 'RISK' (lower) Threshold:" ,
				"\n\t------------------------------------------------------------",
				"\n\tAdditional parameter:",
				"\n\t- Positive float number: 0.0 to 1.0",
				"\n\t(default: ", $href_value_of_scr_par->{menu_parameter}{lowerlimit} ,")",
				"\n\t",
				"\n\tStart value of lower 'RISK' threshold optimization, regarding",
				"\n\t'RISK' filtering of potentially rNC/Na biased species-quartets",
				"\n\t(see -r option). The lower 'RISK' threshold is stepwise increased",
				"\n\t(following the 'tsrisk' defined scaling) until either an optimum",
				"\n\tof rejected species-quartets or the upper 'RISK' threshold is reached.",
				"\n\t",
				"\n\tNote: The lower (start) threshold",
				"\n\thas to be less/equal the upper threshold",
				"\n\tlimit ('tlrisk' <= 'turisk').",
				"\n\t",
				"\n\tSpecified via '-tlrisk' option:",
				"\n\t\t\"e.g. -tlrisk 0.92\"",
				"\n\t---------------------------------------------------\n\n",
		;
	}
	
	elsif	( $$sref_command eq 'turisk' ){

		print	"\n\t------------------------------------------------------------",
				"\n\tSpecies-Quartet Filter 'RISK' (upper) Threshold:" ,
				"\n\t------------------------------------------------------------",
				"\n\tAdditional parameter:",
				"\n\t- Positive float number: 0.0 to 1.0",
				"\n\t(default: ", $href_value_of_scr_par->{menu_parameter}{upperlimit} ,")",
				"\n\t",
				"\n\tUpper (fix) threshold limit for 'RISK' quartet-filtering.",
				"\n\tBesides its role as upper 'RISK' limit during the threshold",
				"\n\toptimization step, this value has also an impact on 'RISK1'",
				"\n\tand 'RISK3' filter since the rejection of quartets is based on",
				"\n\tthe proportion of how many of the three possible species-quartet trees",
				"\n\tare above the upper threshold, which is not the case for",
				"\n\t'RISK2' and 'RISK4' analyses.",
				"\n\t",
				"\n\tNote: The upper (fix) threshold has to be",
				"\n\tgreater/equal the lower threshold",
				"\n\tlimit ('tlrisk' <= 'turisk').",
				"\n\t",
				"\n\tSpecified via '-tlrisk' option:",
				"\n\t\t\"e.g. -turisk 1.01\"",
				"\n\t---------------------------------------------------\n\n",
		;
	}
	
	elsif	( $$sref_command eq 'tsrisk' ){

		print	"\n\t------------------------------------------------------------",
				"\n\tSpecies-Quartet-Filter 'RISK' Threshold Scaling:" ,
				"\n\t------------------------------------------------------------",
				"\n\tAdditional parameter:",
				"\n\t- Positive float number: 0.0 to 1.0",
				"\n\t(default: +", $href_value_of_scr_par->{menu_parameter}{rscale} ,")",
				"\n\t",
				"\n\tPositive scale factor for lower RISK threshold optimization.",
				"\n\tStarting from the defined lower threshold value ('tlrisk'),",
				"\n\tthe number of 'RISK' rejected species-quartets is evaluated",
				"\n\tfor each scale step until either an optimum of filtered species-quartets",
				"\n\tor the upper threshold value ('turisk') is reached (whose value",
				"\n\tis then used as threshold criteria).",
				"\n\t",
				"\n\tSpecified via '-tsrisk' option:",
				"\n\t\t\"e.g. -tsrisk 0.05\"",
				"\n\t---------------------------------------------------\n\n",
		;
	}
	
	elsif	( $$sref_command eq 'u' ){ # done

		print	"\n\t------------------------------------------------------------",
				"\n\tSupression of Script Queries:" ,
				"\n\t-----------------------------------",
				"\n\tAllowed format:",
				"\n\t- No extra parameter necessary",
				"\n\t(default: ", $href_value_of_scr_par->{menu_parameter}{query_setting}, ")",
				"\n\t",
				"\n\tThe '-u' option opposes possible script queries during",
				"\n\tthe process run, which would stop automatic pipeline processes", 
				"\n\tif the maximum number of allowed species-quartets is lower than the",
				"\n\tactual number of quartets ('-M' option).",
				"\n\tWith the '-u option', script queries are suppressed, meaning",
				"\n\tthat species-quartets are drawn randomly until the number of",
				"\n\tspecies-quartets exceeds the number of allowed quartets.",
				"\n\t",
				"\n\tAn activated '-u' option can be deactivated",
				"\n\tby typing '-u' again.",
				"\n\t",
				"\n\tSpecified via '-u' option:",
				"\n\t\t\"-u\"",
				"\n\t---------------------------------------------------\n\n",
		;
	}

	elsif	( $$sref_command eq 'prt' ){ # done

		print	"\n\t------------------------------------------------------------",
				"\n\tAdditional Print Option:" ,
				"\n\t-----------------------------------",
				"\n\tAllowed format:",
				"\n\t- Positive integer number: 1 to 3 or 'x' for deactivation",
				"\n\t(default: '", $href_value_of_scr_par->{menu_parameter}{prt_setting}, "')",
				"\n\t",
				"\n\tSelectable parameter:",
				"\n\t\tx\tDeactivate additional prints\t('none')",
				"\n\t\t1\tR-PLOTTING\t\t\t('Rplot')",
				"\n\t\t2\tLATEX-TABLE\t\t\t('LatexTable')",
				"\n\t\t3\tR-PLOTTING+LATEX-TABLE\t\t('LatexTable,Rplot')",
				"\n\t",
				"\n\tThe '-prt' option allows the output of additional",
				"\n\tresult-graphics (generated with R) and LATEX table",
				"\n\tprints. To use both print paramters succesfully,",
				"\n\tadditional R and Latex softwarepackages are needed",
				"\n\t(see manual).",
				"\n\t",
				"\n\twith the 'x' parameter activated extra prints",
				"\n\tare deactivated.",
				"\n\t",
				"\n\tSpecified via '-prt' option:",
				"\n\t\t\"e.g. -prt 3\"",
				"\n\t---------------------------------------------------\n\n",
		;
	}
	
	elsif	( $$sref_command eq 'P' ){

		&preface( \%$href_value_of_scr_par )
	}

	else{ die "\nBUG-ERROR: Cannot assign command ", $$sref_command, " in subroutine &help!\nPlease, report BUG to system developer!\n\n" }
	
	##############################
	#
	print	"\tFor main help menu\ttype\t-h <enter>"		,
			"\n\tFor config menu\t\ttype\t-b <enter>"		,
			"\n\tTo quit script process\ttype\t-q <enter>"	,
			"\n\t------------------------------------------------------------\n",
			"\n\tcommand:\t " ;
	
	chomp ( my $help_menu_answer = <STDIN> );
	
	if 	( $help_menu_answer =~ /^-b$/ ){ &argv_handling( \%$href_value_of_scr_par	) }
	else{ &check_user_parameter( \%$href_value_of_scr_par, \$help_menu_answer ) }
	##############################
	
	######################################################################
}

sub preface{

	my	$href_value_of_scr_par	= $_[0] ; # key1: method_code; key2: method_option; value: float or integer			IN: defined; OUT: unchanged
	
	print	"-------------------\n" , 
			":&preface\n";
	
	######################################################################
	# called from subroutine(s)...
	# &argv_handling
	# &help
	#
	############
	# using subroutine(s)...
	# ....
	######################################################################
	
	######################################################################
	# print skript information & EXIT                                    #
	######################################################################
	print	"\n\tVersion     : ", 	$href_value_of_scr_par->{author_info}{version}		,
			"\n\tLanguage    : ",	$href_value_of_scr_par->{author_info}{language}		,
			"\n\tLast Update : ",	$href_value_of_scr_par->{author_info}{lastUpdate}	,
			"\n\tAuthor      : ",	$href_value_of_scr_par->{author_info}{name}, ", "	, 
									$href_value_of_scr_par->{author_info}{institute}	,
			"\n\te-mail      : ",	$href_value_of_scr_par->{author_info}{email}		,
			"\n\tHomepage    : ",	$href_value_of_scr_par->{author_info}{homepage}		,
			"\n\t",
			"\n\tThis program is free software; you can distribute it ",
			"\n\tand/or modify it under the terms of the GNU General Public ",
			"\n\tLicense as published by the Free Software Foundation ; ",
			"\n\teither version 2 of the License, or (at your option) any ",
			"\n\tlater version.",
			"\n\t",
			"\n\tThis program is distributed in the hope that it will be",
			"\n\tuseful, but WITHOUT ANY WARRANTY; without even the",
			"\n\timplied warranty of MERCHANTABILITY or FITNESS FOR A",
			"\n\tPARTICULAR PURPOSE. See the GNU General Public License for",
			"\n\tmore details.",
			"\n\t",
			"\n\tYou should have received a copy of the GNU General Public",
			"\n\tLicense along with this program; if not, write to the Free",
			"\n\tSoftware Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139,",
			"\n\tUSA.",
			"\n\t---------------------------------------------------\n\n"
	;
	
	exit;
	######################################################################
}







 

















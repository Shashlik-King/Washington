#!/usr/bin/perl -w
use strict;

##############################################################################
#
# File   :  create-appendix.pl
# History:  2015-11-11 (Martin Kelm)    first implementation
#           2016-03-22 (Thorben Hamann) geotechnical input
#           2016-03-23 (Martin Kelm)    structural input, output path changed
#           2016-03-29 (Martin Kelm)    adapted to Hohe See
#           2016-03-24 (Thorben Hamann) geotechnical input adapted to Hohe See
#           2016-04-24 (Martin Kelm)	aligned
#
##############################################################################
#
#  
#
##############################################################################

=head1 NAME

perl create-appendix.pl

=cut

##############################################################################

our $script_version = "2016-04-24";


use lib "../commons/";
use lib "../../../lib";

use Math::Trig;
use POSIX;

use Debug;
use CAD;
use Defs;
use DBOps;
use Soil;
use CanProperties;


#******************************************************************************
# general parameters
#******************************************************************************
our $scale    =  2.0;		# scale of drawing 1:(100x value)
our $LAT_pos  = 34.0;		# vertical position of LAT on page in cm, measured from bottom

our $xcl   	=   7.8;		# center line of monopile
our $xpos1 	=  $xcl - 3.0;	# horizontal position of elevation markers

our $tb_scale   = 1.0;
our $tb_width   = 12.0*$tb_scale; 
our $tb_minipage_width = $tb_width - 0.5;


our $MP_ID       = "";
our $MP_ID_STR   = "";
our $CONF        = "00";
our $version     = "0.1";


our $prj_no = "A071540";
our $prj_name = "EnBW Hohe See OWF";

our $scriptname = __FILE__;

our $approved_by = "KIML";

our $FORMAT      = "A3";
our $ORIENTATION = "portrait"; # "landscape";"portrait"

# Hohe See specific modifications:
our $Alstom_load_factor = 1.00;		# Not used at Hohe See: load factor (1%) applied on characteristic ULS shear force to account for Alstom load factors
our $sec_order_factor 	= 1.06;		# factor (6%) applied on characteristic ULS overturning moment to account for Alstom load factor and 2nd order effect
our $psf_axial_load		= 1.35; 	# partial safety factor on axial load
our $psf_axial_res		= 1.40;		# partial safety factor on axial resistance
our $psf_lat_load		= 1.35; 	# partial safety factor on lateral load (DIN1054 approach)
our $psf_lat_res		= 1.40; 	# partial safety factor on lateral resistance (DIN1054 approach) 
our $perm_rot_allow		= 0.25;		# [°] allowed permanent rotation
our $crit_len_criterion	= 10.0;		# [%] applied percentage for critical pile length criterion

our $target_eigenfreq_shallow = 0.247;
our $target_eigenfreq_deep    = 0.245;

our $tower = "T83.80-01, Rev. 2"; 


MAIN:{
	Debug::set_debug_level( Debug::DEBUG_LEVEL_ALL );
	
	if( $#ARGV == 2 ){
		$MP_ID   = $ARGV[0];
		$CONF    = $ARGV[1];
		$version = $ARGV[2];
	}
	else{
		usage();
	}
	
	$MP_ID_STR = $MP_ID;
	$MP_ID_STR =~ s/_/\\_/g;
	
	# get data out of database
	
	DB_connect( Defs::DB );
	
	my @cfg = DBOps::DB_configuration( $MP_ID, $CONF );
	
	my $MP_STRUCT_REV = $cfg[2];
	my $TP_STRUCT_REV = $cfg[3];
	my $TP_TYPE       = $cfg[8];
	my $LOC           = $cfg[7];
	my $REV_ATTM      = $cfg[9];
	
	
	my $file_attm = DBOps::DB_file( $MP_ID, $CONF, "attm_file", "structural_pdf", "APP/input/", "attm.pdf" );
	my $file_fls  = DBOps::DB_file( $MP_ID, $CONF, "fls_file", "structural_pdf", "APP/input/", "fls.pdf" );
	my $file_uls  = DBOps::DB_file( $MP_ID, $CONF, "uls_file", "structural_pdf", "APP/input/", "uls.pdf" );
	
	my $file_soil_profile = DBOps::DB_file( $MP_ID, $CONF, "soil_profile_file", "soil_pdf", "APP/input/", "soil.pdf" );
	my $file_py_curves    = DBOps::DB_file( $MP_ID, $CONF, "py_curves_file", "soil_pdf", "APP/input/", "py.pdf" );
	
	
	my $code = "";
	
	# page format
	use constant WIDTH_PAGE  => 13.0; # cm
	use constant HEIGHT_PAGE => 24.0; # cm
	
	use constant WIDTH_TITLEPAGE  => 18.0; # cm
	use constant HEIGHT_TITLEPAGE => 26.2; # cm
	
	our @margins = ("4cm","0cm", WIDTH_PAGE."cm", HEIGHT_PAGE."cm","0cm");
	
	Latex::add_package( "pdfpages" );
	#Latex::add_package( "draftwatermark" );
	Latex::add_package( "float" );
	
	$code .= Latex::build_header( "article", "11pt,a4paper", "\\large", \@margins);
	
	
	my @doc_revs = DBOps::DB_doc_revisions( $MP_ID );
	my @rev_data = ();
	my( $doc_rev_cur, $doc_main_rev_cur, $doc_sub_rev_cur, $date_rev0, 
		$drawn_by_rev0, $checked_by_rev0, $approved_by_rev0 ) 
		= CAD::rev_data( \@doc_revs, \@rev_data, $approved_by );
		
		
	my $doc_rev = "";
    my $descr = "";
    my $issue_date = "";
    my $drawn_by = "";
    my $checked_by = "";
    my $approved = "";
    # last entry in @rev_data is taken
    if( @rev_data ){
        ( $doc_rev, $descr, $issue_date, $drawn_by, $checked_by, $approved ) = @{$rev_data[$#rev_data]};
    }
    else{
        $drawn_by   = $drawn_by_rev0;
        $checked_by = $checked_by_rev0;
    }


	# last entry in @rev_data is taken
	if( @rev_data ){
		( $doc_rev, $descr, $issue_date, $drawn_by, $checked_by, $approved ) = @{$rev_data[$#rev_data]};
	}
	else{
		$drawn_by   = $drawn_by_rev0;
		$checked_by = $checked_by_rev0;
	}
	
	
	$code .= build_titlepage( $version, $drawn_by, $checked_by, $approved_by );
	
	
	my $indentspace = -2.0;
	$code .=  "\\pagestyle{fancy}\n";
	$code .=  "\\renewcommand{\\headrulewidth}{0pt}\n";
#	print "\\lhead{\\scriptsize \\sf\\hspace{-4cm}  - TP and MP - Detailed Design Report / Appendix $INPUT_LOC}\n";
    $code .=  "\\rhead{\\scriptsize\\sf\\hspace{$indentspace cm} 
    			\\includegraphics[height=0.3cm]{../resources/cowi-logo.pdf} \\hspace*{0.9cm} \\\\ 
    			ENBW HOHE SEE OWF \\hspace{0.7cm} \\arabic{page} \\\\
    			DESIGN REPORT PRIMARY STEEL - APPENDIX $MP_ID_STR \\hspace*{0.9cm}
    		}\n";	
	
	$code .=  "\\lhead{}\n";

	$code .=  "\\rfoot{\\sf\\hspace{$indentspace cm}\\tiny TODO-PATH}\n";
#	print "\\lfoot{\\sf\\hspace{-4cm}\\tiny G/DAT/Prj/214_PNE_OWP_GWII_Ausfuerung/05 Berichte/...  }\n";#  ...\\tiny $path}\n; # ".$indentspace."
	$code .=  "\\cfoot{}\n";
	$code .=  "\\rfoot{}\n";
	
	
	
	#$code .= "\\tableofcontents\n";
	
	my %sections = ( );
	
	my $sum_code = build_summary( \%sections, $MP_ID, $MP_STRUCT_REV, $LOC, $TP_TYPE, $TP_STRUCT_REV, $tower, $CONF );
	my $prf_code = build_soil_profile( \%sections, $file_soil_profile );
	my $py_code  = build_py_curves( \%sections, $file_py_curves );
	
	my $geo_code = build_geo_design( \%sections );
	my $str_code = build_str_design( \%sections );
	
	
	$code .= build_contents( \%sections, $MP_STRUCT_REV, $TP_STRUCT_REV );
	$code .= "\n\\pagebreak\n";
	
	$code .= $sum_code;
	$code .= "\n\\pagebreak\n";
	
	$code .= $prf_code;
	$code .= $py_code;
	
	$code .= $geo_code;
	$code .= "\n\\pagebreak\n";
	
	$code .= $str_code;
	
	
	$code .= "\\includepdf[pages={1},fitpaper]{$file_attm}\n";
	$code .= "\\includepdf[pages={1},fitpaper]{$file_fls}\n";
	$code .= "\\includepdf[pages={1},fitpaper]{$file_uls}\n";
	
	$code .= Latex::build_footer();
	
	DB_disconnect();
	
	# generate PDF
	our $output_dir = "../../../../output/". Defs::DB . "/appendix/APP/$MP_ID/";
	
	#my $output_dir = "APP/$MP_ID";
	unless( -d $output_dir ){
		unless ( mkdir( $output_dir ) ){
			my $msg = "$output_dir not created";
			Debug::fatal( $msg, "MAIN", __LINE__ );
		}
	}
	
	
	my $filename   = "$prj_no-RP-DD-12-APP-$MP_ID-V$version"; 
#	my $str_main_rev = $doc_main_rev_cur;
#	if( $doc_main_rev_cur < 10 ){ $str_main_rev = "0" . $doc_main_rev_cur }
#	my $filename   = "$dwg_id-$str_main_rev"; # "overview-dwg"; #"$dwg_no-$rev_dwg";
#	if( $doc_sub_rev_cur > 0 ){
#		$filename   = "$dwg_id-$str_main_rev-prelim-$doc_rev_cur";
#	}
	
	Latex::create_pdf_from_str( $code, $output_dir, $filename, "pdflatex" );
	
}

##############################################################################
# subroutines
##############################################################################

sub usage
{
	print STDERR "USAGE: $scriptname <loc> <rev_conf> <version, e.g 1.0>\n";
	exit 0;
}


sub build_titlepage
{	
	Debug::check_num_params( 4, \@_, (caller(0))[3], __LINE__ );
	
	my $version  = shift;
	
	my $preparer = shift;
	my $checker  = shift;
	my $approver = shift;
	
	my $code = "";
	
	$code .= Latex::build_begin_titlepage();
	
	my $xpos = -4.0;
	
	my $xpos_logo_left = $xpos;
	my $ypos_logo_left = 25.0;
	
	my $xpos_logo_right = $xpos + 12.0;
	my $ypos_logo_right = 25.5;
	
	$code .= "\\put( $xpos_logo_left,$ypos_logo_left ){\\includegraphics[height=2.0cm]{../resources/gs-logo.pdf}}\n";
	$code .= "\\put( $xpos_logo_right,$ypos_logo_right ){\\includegraphics[height=1.5cm]{../resources/cowi-logo.pdf}}\n";
	
	my $t = "{\\small MAY 2016}\\\\
			{\\small GEOSEA}\\vspace*{3ex}\\\\
			{\\Huge ENBW HOHE SEE OWF}\\vspace*{3ex}\\\\
			{\\small DETAILED FOUNDATION DESIGN}\\vspace*{2ex}\\\\
			\\textbf{\\small DESIGN REPORT PRIMARY STEEL - APPENDIX $MP_ID_STR}";
	$code .= Latex::build_parbox( $xpos, 16.0, $t, WIDTH_TITLEPAGE );
	
	
	$code .= "\\put( $xpos, 0.0 ){\\begin{minipage}[b]{12cm}\\tiny\n";
	
	$code .= "\\begin{tabbing}\n";
	$code .= "XXXXXXXXXXXXXXXXXXXX \\=\\kill\\\\\n";
	$code .= "\\textcolor{red}{PROJECT NO.}          \\> $prj_no\\\\\n";
	$code .= " \\> \\\\\n";
	$code .= "\\textcolor{red}{DOCUMNENT NO.}        \\> $prj_no-RP-DD-12-APP-$MP_ID_STR\\\\\n";
	$code .= " \\> \\\\\n";
	$code .= "\\textcolor{red}{VERSION}              \\> $version\\\\\n";
	$code .= " \\> \\\\\n";
	$code .= "\\textcolor{red}{DATE OF ISSUE}        \\> 2016-05-14\\\\\n";
	$code .= " \\> \\\\\n";
	$code .= "\\textcolor{red}{PREPARED}             \\> $preparer\\\\\n";
	$code .= " \\> \\\\\n";
	$code .= "\\textcolor{red}{CHECKED}              \\> $checker\\\\\n";
	$code .= " \\> \\\\\n";
	$code .= "\\textcolor{red}{APPROVED}             \\> $approver\\\\\n";
	$code .= "\\end{tabbing}\n";
	
	$code .= "\\end{minipage}}\n";
	
	
	$code .= Latex::build_end_titlepage();
	
}

sub build_contents
{
	Debug::check_num_params( 3, \@_, (caller(0))[3], __LINE__ );
	
	my $ref_sections  = shift;
	my $MP_STRUCT_REV = shift;
	my $TP_STRUCT_REV = shift;
	
	my $code = "";
	
	$code .= "\\setcounter{page}{1}\n";

	# table of contents
	$code .= "\\vspace*{5cm}\n";
	$code .= "{\\LARGE\\color{red} CONTENTS}\n\n";
	#print "\\vspace*{0.5cm}\n";
	
	$code .= "\\begin{tabbing}\n";
	$code .= "\\hspace*{1.5cm}\\=\\hspace*{11.5cm}\\=\\kill\n";
	
	foreach my $secnum (sort keys %{$ref_sections}){
		my $sectext = ${$ref_sections}{$secnum};
		my $rounded = ceil($secnum );
		# print STDERR "CEIL IS: $rounded \n";
		
		if( $rounded == $secnum ){
			$code .=  "\\\\\n";
			$code .=  "\\bfseries{$secnum} \\>\\bfseries{$sectext}\\>\\makebox[0cm][r]{\\bfseries{\\pageref{sec_$secnum}}}\\\\\n";
		}
		else{
			$code .=  "$secnum \\>$sectext\\>\\makebox[0cm][r]{\\pageref{sec_$secnum}}\\\\\n";
		}
	}
	
	$code .=  "\\end{tabbing}\n";
	$code .=  "\\normalsize\n";
	
	$code .=  "\\vfill\n";
	$code .=  "{\\scriptsize Appendix generated on \\today~for configuration $CONF
				(structural revision TP $TP_STRUCT_REV and structural revision MP $MP_STRUCT_REV)}\n";
	
	#print "hello";
	return $code;
}

sub build_summary
{
	Debug::check_num_params( 8, \@_, (caller(0))[3], __LINE__ );
	
	my $ref_sections  = shift;
	my $MP_ID         = shift;
	my $MP_STRUCT_REV = shift;
	my $LOC           = shift;
	my $TP_TYPE       = shift;
	my $TP_STRUCT_REV = shift;
	my $tower         = shift;
	my $CONF          = shift;
	
	my $code = "";
	
	my $sec_count = 1;
	my $sec_count_old = 0;
	
	my $title = "Summary";
	
	$code .= Latex::build_section( $title, $ref_sections );
	
	$code .= build_table_sum( $MP_ID, $LOC, $TP_TYPE, $tower );
	
	$code .= Latex::build_par("In general, the fatigue limit state is the governing design situation.");
	
	$code .= Latex::build_par("The main dimensions and weights for location $MP_ID_STR are listed in the following table.");
	
	$code .= build_table_dim( $MP_ID, $MP_STRUCT_REV, $TP_TYPE, $TP_STRUCT_REV, $CONF );
	
	return $code;
}


# ----------------------------------------------------------------------------
# $code = build_table_sum( $INPUT_LOC );
# ----------------------------------------------------------------------------
# table with summarized data
# ----------------------------------------------------------------------------

sub build_table_sum 
{
	Debug::check_num_params( 4, \@_, (caller(0))[3], __LINE__ );
	
	my $MP_ID   = shift;
	my $LOC     = shift;
	my $TP_TYPE = shift;
	my $tower   = shift;
	
	my $code = "";
	
	my @data = ();
	my @row0 = ( "\\textbf{$MP_ID_STR}" );
	
	my $target_eigenfreq = 0.0;
	my $cluster = DBOps::DB_cluster( $MP_ID );
	if( $cluster eq "DEEP" ){ $target_eigenfreq = $target_eigenfreq_deep; }
	elsif( $cluster eq "SHALLOW" ){ $target_eigenfreq = $target_eigenfreq_shallow; }
	else{
		die;
	}
	
	my @row1 = ( "\\textbf{$cluster} / \\textbf{$target_eigenfreq}" );
	
	
	
	my $tp_type = $TP_TYPE;
	
	my $seabed = DBOps::DB_sea_floor_level( $LOC );
	my @row2 = ( $tower );
	my @row3 = ( "$tp_type" );
	
	my $str_seabed = sprintf( "%.1f",$seabed );
	my @row4 = ( $str_seabed );
	
	push( @data, \@row0 );
	push( @data, \@row1 );
	push( @data, \@row2 );
	push( @data, \@row3 );
	push( @data, \@row4 );
	
	my @col_names = ();
	my @row_names = ( "\\cellcolor{dgrey}Location", 
		"\\cellcolor{dgrey}Cluster / target eigenfrequency [Hz]",
		"\\cellcolor{dgrey}Tower structure",
		"\\cellcolor{dgrey}TP type", 
		"\\cellcolor{dgrey}Sea floor level [m LAT]" );
	
	$code .= Latex::build_tabular(\@data,undef,\@row_names, undef);
	
	$code .= "\\vspace{1ex}\n\n";
	
	return $code;
}


# ----------------------------------------------------------------------------
# $code = build_table_dim( $INPUT_LOC );
# ----------------------------------------------------------------------------
# table with summarized data
# ----------------------------------------------------------------------------

sub build_table_dim 
{
		Debug::check_num_params( 5, \@_, (caller(0))[3], __LINE__ );
	
		my $MP_ID         = shift;
		my $MP_STRUCT_REV = shift;
		my $TP_TYPE       = shift;
		my $TP_STRUCT_REV = shift;
		my $CONF          = shift;
		
		my $mcode = "";
		
		
		my @mdata = ( );
		
		my $top_tp = DBOps::DB_tp_top( $TP_TYPE, $TP_STRUCT_REV );
		my $bot_tp = DBOps::DB_tp_bottom( $TP_TYPE, $TP_STRUCT_REV );
		
		my $len_tp = $top_tp - $bot_tp;
		my $dia_mp = 0.0; # get_diameter($INPUT_LOC, MP_CANS);
		my $diam;
		
		our @mp_cans = DBOps::DB_pile( $MP_ID, $MP_STRUCT_REV );
		our @tp_cans = DBOps::DB_tp( $TP_TYPE, $TP_STRUCT_REV );
		
		my $top_dia_mp = CanProperties::top_outer_diameter( \@mp_cans );
		$top_dia_mp *= 1000.0;
		my $top_dia_tp = CanProperties::top_outer_diameter( \@tp_cans );
		$top_dia_tp *= 1000.0;
		
		my $bot_dia_mp = CanProperties::bot_outer_diameter( \@mp_cans );
		$bot_dia_mp *= 1000.0;
		my $bot_dia_tp = CanProperties::bot_outer_diameter( \@tp_cans );
		$bot_dia_tp *= 1000.0;
		
		my $tp_weight = CanProperties::calc_weight( \@tp_cans );
		my $mp_weight = CanProperties::calc_weight( \@mp_cans );
		
		my $w_tp   = sprintf( "%.1f", $tp_weight ); # get_tp_weight( $tp_type, TP_CANS ) );
		my $w_mp   = sprintf( "%.1f", $mp_weight ); # get_mp_weight( $INPUT_LOC, MP_CANS ) );
		
		my @mrow0 = ( $top_tp, $bot_tp, $len_tp, "$top_dia_tp/$bot_dia_tp", $w_tp );
		
		our ($pile_top, $pile_tip) = DBOps::DB_pile_top_bottom( $MP_ID, $MP_STRUCT_REV );
		
		my $top_mp = sprintf( "%.2f", $pile_top ); # get_top_mp( $INPUT_LOC ) );
		my $bot_mp = sprintf( "%.2f", $pile_tip ); # get_bot_mp( $INPUT_LOC ) );
		my $len_mp = abs($top_mp - $bot_mp);
		my $str_top_mp = $top_mp;
		unless( $top_mp < 0.0 ){ $str_top_mp = "+".$top_mp; }
		my @mrow1 = ( $str_top_mp, $bot_mp, $len_mp, "$top_dia_mp/$bot_dia_mp",  $w_mp);
		
		my $total_w = $w_tp + $w_mp;
		my @mrow2 = ( " ", " ", " ", "Total weight", $total_w);
	
		push( @mdata, \@mrow0 );
		push( @mdata, \@mrow1 );
		push( @mdata, \@mrow2 );

		my @mtitle = ( "Main dimension and masses","l" );
		my $mcol0 = "\\parbox{1.5cm}{Top [m~LAT]}";
		my $mcol1 = "\\parbox{1.5cm}{Bottom [m LAT]}";
		my $mcol2 = "\\parbox{1.5cm}{Length [m]}";
		my $mcol3 = "\\parbox{1.5cm}{Diameter Top / Bottom [mm]}";
		my $mcol4 = "\\parbox{1.7cm}{Steel weight [t]\$^{1)}\$}";
		
		my @mcol_names = ( "",$mcol0,$mcol1,$mcol2,$mcol3,$mcol4 );
		my @mrow_names = ( "Transition Piece", "Monopile", " " );
		
		$mcode .= Latex::build_tabular( \@mdata,\@mcol_names,\@mrow_names, \@mtitle );
		
		$mcode .= "\n\n";
		
		$mcode .= "\\textrm{\\textit{ \$^{1)}\$ Steel weight does not include appurtenances and part of the tower flange.}}\n\n";
		$mcode .= "\\vspace{1ex}\n\n";
		
		
		# table with eigenfrequencies
		
		my @results = ();
		DB_select( \@results, "eigenfreq_target,eigenfreq_lower_bound,eigenfreq_upper_bound", "structural_results", "id=\"$MP_ID\" and rev=\"$CONF\"" );
		my( $ef_target,$ef_lower_bound,$ef_upper_bound ) = @{$results[0]};
		
		$mcode .= Latex::build_par("The calculated natural frequencies are documented in the following table.");
		
		my @mdata_eq = ();
		my @mrow0_eq = ( "X.XXX" );
		my @mrow1_eq = ( "X.XXX" );
		
		push( @mdata_eq, [$ef_target] );
		push( @mdata_eq, [$ef_lower_bound] );
		push( @mdata_eq, [$ef_upper_bound] );
		
		my @mtitle_eq = ( "Calculated natural frequencies","l" );
		my $mcol0_eq = "\\parbox{1.5cm}{[Hz]}";
		
		my @mcol_names_eq = ( "", $mcol0_eq );
		my @mrow_names_eq = ( "Characteristic","Lower bound", "Upper bound" );
		
		$mcode .= Latex::build_tabular( \@mdata_eq,\@mcol_names_eq,\@mrow_names_eq, \@mtitle_eq );
		
		
		return $mcode;
}



sub build_soil_profile
{
	my $ref_sections = shift;
	my $file_soil_profile = shift;
	
	# ----------------------------------------------------------------------------
	# general definitions:
	# ----------------------------------------------------------------------------
	my $code 			= "";
	my $sec_count 		= 2; 		# needed? ->MNKM
	my $sec_count_old 	= 1; 	# needed? ->MNKM
	my $title = "Soil profile";
	
	
	# ----------------------------------------------------------------------------
	# generation of text:
	# ----------------------------------------------------------------------------
	#$code .= "\\newpage\n";
	$code .= Latex::build_section( $title, $ref_sections );
	#$code .= "\\includepdf[pages={1},fitpaper]{$file_soil_profile}\n";
	
	$code .= Latex::build_par("The soil profile applied in the design is based on the design soil profiles including
								characteristic soil parameters provided by the Geotechnical Expert in the Soil and
								Foundation expertise. The applied soil profile at location $MP_ID_STR is given on the following page.");
	
	$code .= "\\includepdf[pages={1},fitpaper]{$file_soil_profile}\n";
	
	return $code;
}


sub build_py_curves
{
	my $ref_sections = shift;
	my $file_py_curves = shift;
	
	# ----------------------------------------------------------------------------
	# general definitions:
	# ----------------------------------------------------------------------------
	my $code 			= "";
	my $sec_count 		= 2; 		# needed? ->MNKM
	my $sec_count_old 	= 1; 	# needed? ->MNKM
	my $title = "P-y curves";
	
	
	# ----------------------------------------------------------------------------
	# generation of text:
	# ----------------------------------------------------------------------------
	#$code .= "\\newpage\n";
	$code .= Latex::build_section( $title, $ref_sections );
	#$code .= "\\includepdf[pages={1},fitpaper]{$file_soil_profile}\n";
	
	$code .= Latex::build_par("The p-y-curves applied for the design at location $MP_ID_STR are based on characteristic soil properties and given on the following page.");
	
	$code .= "\\includepdf[pages={1},fitpaper]{$file_py_curves}\n";
	
	return $code;
}



sub build_geo_design
{
	my $ref_sections = shift;
	
	# ----------------------------------------------------------------------------
	# general definitions:
	# ----------------------------------------------------------------------------
	my $code 			= "";
	my $sec_count 		= 3; 		# needed? ->MNKM
	my $sec_count_old 	= 2; 	# needed? ->MNKM
	my $title 			= "Geotechnical design";		# title of section
	
	
	$code .= Latex::build_section( $title, $ref_sections );
	# ----------------------------------------------------------------------------
	# getting code other values from subsection subroutines
	# ----------------------------------------------------------------------------
	my @return_loads  			= build_geo_loads( $ref_sections );		# getting code and load from build_geo_loads subroutine
	my $axial_char_force		= $return_loads[1];						# [kN] characteristic axial force at mudline
	my $torsional_design_moment	= $return_loads[2];						# [kNm] design torsional moment at mudline
	
	# code of the individual subsections
	my $geo_sub_loads_code 		= $return_loads[0];						# code of the load subsection
	my $geo_pile_length_code 	= build_geo_pile_length( $ref_sections );#code of the pile length subsection
	my $geo_axial_check 		= build_geo_axial( $ref_sections, $axial_char_force, $torsional_design_moment);	# code of the axial capacity subsection
	
	
	
	# ----------------------------------------------------------------------------
	# preparing a table with the main results
	# ----------------------------------------------------------------------------
	
	# getting the relevant numbers from the database
	my $crit_pile_length 			= DBOps::DB_soil_results( $MP_ID, $CONF, "crit_pile_length");		# [m] critical pile embedment length
	my $deflec_mudline 				= DBOps::DB_soil_results( $MP_ID, $CONF, "hori_defl_mudline");		# [m] deflection at mudline
	my $deflec_tip		 			= DBOps::DB_soil_results( $MP_ID, $CONF, "hori_defl_pile_tip");		# [m] deflection at pile tip
	my $perm_rot_mudline_ULS 		= DBOps::DB_soil_results( $MP_ID, $CONF, "perm_rot_mudline");		# [°] permanent rotation at mudline under ULS loads
	my $perm_rot_mudline_oper 		= DBOps::DB_soil_results( $MP_ID, $CONF, "perm_rot_mudline_oper");	# [°] permanent rotation at mudline under SLS loads
	my $lat_utilisation 			= DBOps::DB_soil_results( $MP_ID, $CONF, "lat_util_ratio");		# [-] lateral utilisation ratio
	my $lat_res_char	 			= DBOps::DB_soil_results( $MP_ID, $CONF, "lat_capa");			# [kN] characteristic lateral resistance
	my $lat_force_char	 			= DBOps::DB_soil_results( $MP_ID, $CONF, "lat_force");		# [kN] characteristic mobilised lateral force
	my $axial_util_ratio 			= DBOps::DB_soil_results( $MP_ID, $CONF, "axial_util_ratio");	# [-] axial utilisation ratio
	my $axial_capa		 			= DBOps::DB_soil_results( $MP_ID, $CONF, "axial_capa_comp");			# [kN] characteristic axial pile capacity
	
	#formatting the values:
	my $crit_pile_length_str 		= sprintf( "%.1f", $crit_pile_length);								# [m] critical pile embedment length, formatted text
	my $deflec_mudline_str 			= sprintf( "%.1f", $deflec_mudline*1000);							# [m -> mm] deflection at mudline, formatted text
	my $deflec_tip_str 				= sprintf( "%.1f", $deflec_tip*1000);								# [m -> mm] deflection at pile tip, formatted text
	my $perm_rot_mudline_ULS_str 	= sprintf( "%.3f", abs($perm_rot_mudline_ULS));						# [°] permanent rotation at mudline under ULS loads, formatted text
	my $perm_rot_mudline_oper_str	= sprintf( "%.3f", abs($perm_rot_mudline_oper));					# [°] permanent rotation at mudline under SLS loads, formatted text
	my $perm_rot_mudline_str	 	= sprintf( "%.3f", abs($perm_rot_mudline_ULS + $perm_rot_mudline_oper)); # [°] total permanent rotation at mudline, formatted text
	my $perm_rot_allow_str 			= sprintf( "%.2f", $perm_rot_allow) ;
	my $lat_utilisation_str 		= sprintf( "%.2f", $lat_utilisation);							# [-] lateral utilisation ratio, formatted text
	my $lat_res_char_str 			= sprintf( "%.0f", $lat_res_char);								# [kN] characteristic lateral resistance, formatted text
	my $lat_force_char_str 			= sprintf( "%.0f", $lat_force_char);							# [kN] characteristic mobilised lateral force, formatted text
	my $axial_util_ratio_str 		= sprintf( "%.2f", $axial_util_ratio);							# [-] axial utilisation ratio, formatted text
	my $axial_capa_str 				= sprintf( "%.0f", $axial_capa);								# [kN] characteristic axial pile capacity, formatted text
	my $axial_char_force_str  		= sprintf( "%.0f", $axial_char_force);							# [kN] characteristic axial load, formatted text
	
	my @mdata 	= ( );			# initialise array containing table data
	my @mrow0 	= ( $crit_pile_length_str, "-");
	push( @mdata, \@mrow0 );
	my @mrow1 	= ( $deflec_mudline_str, "-");
	push( @mdata, \@mrow1 );
	my @mrow2 	= ( $deflec_tip_str, "-");
	push( @mdata, \@mrow2 );
	my @mrow3 	= ( $perm_rot_mudline_ULS_str, "-");
	push( @mdata, \@mrow3 );
	my @mrow4 	= ( $perm_rot_mudline_oper_str, "-");
	push( @mdata, \@mrow4 );
	my @mrow5 	= ( $perm_rot_mudline_str, $perm_rot_allow_str);
	push( @mdata, \@mrow5 );
	my @mrow6 	= ( $lat_force_char_str, "-");
	push( @mdata, \@mrow6 );
	my @mrow7 	= ( $lat_res_char_str, "-");
	push( @mdata, \@mrow7 );
	my @mrow8 	= ( $lat_utilisation_str, "1.0");
	push( @mdata, \@mrow8 );
	my @mrow9 	= ( $axial_char_force_str, "-");
	push( @mdata, \@mrow9 );
	my @mrow10 	= ( $axial_capa_str, "-");
	push( @mdata, \@mrow10 );
	my @mrow11 	= ( $axial_util_ratio_str, "1.0");
	push( @mdata, \@mrow11 );
		
	# defining table headings
	my @mtitle = ( "Summary of geotechnical design at location $MP_ID_STR","l" );
	my $mcol0 = "\\parbox\[l\]{1.5cm}{ }";
	my $mcol1 = "\\parbox{3.0cm}{calculated value}";
	my $mcol2 = "\\parbox{1.2cm}{limit}";
	my @mcol_names = ( $mcol0,$mcol1,$mcol2 );
	
	my @mrow_names = ( "Critical pile embedment length L \[m\]",
	 	 "Pile deflection at mudline \$u\_\\mathrm\{h,mud\}\$ \[mm\]",
	 	 "Pile tip deflection \$u\_\\mathrm\{h,tip\}\$ \[mm\]",
		 "Permanent rotation ULS loads \[\$\^\\circ\$\]",
		 "Permanent rotation SLS loads \[\$\^\\circ\$\]",
		 "Total permanent rotation \[\$\^\\circ\$\]",
		 "Characteristic mobilised lateral force \[kN\]",
		 "Characteristic lateral resistance \[kN\]",
		 "Lateral utilisation ratio \[-\]",
		 "Characteristic axial loading \[kN\]",
		 "Characteristic axial pile capacity \[kN\]",
		 "Axial utilisation ratio \[-\]",
		 );
	
	
	# ----------------------------------------------------------------------------
	# generation of text / adding the code of individual subsections
	# ----------------------------------------------------------------------------
	$code .= Latex::build_par("The geotechnical design is based on pile geometry given in Section 5 .
			 All calculations are performed in accordance with the methods described in the Design Brief. 
			 A summary of the geotechnical design is given in the following table.");
	# generation of table
	$code .= Latex::build_tabular( \@mdata,\@mcol_names,\@mrow_names, \@mtitle );
	$code .= "\\vspace{1ex}\n\n";;
	#generation of subsections:
	$code .= $geo_sub_loads_code;
	$code .= $geo_pile_length_code;
	$code .= $geo_axial_check;

	return $code;
}


sub build_geo_loads
{
	my $ref_subsections = shift;
	
	# ----------------------------------------------------------------------------
	# general definitions:
	# ----------------------------------------------------------------------------
	my $code 			= "";
	my $sec_count 		= 2; 		# needed? ->MNKM
	my $sec_count_old 	= 0; 	# needed? ->MNKM
	my $title 			= "Applied loads";	# title of subsection
	
	
	# ----------------------------------------------------------------------------
	# import and interpolation of the ULS loads at mudline:
	# ----------------------------------------------------------------------------
		
	# getting data from database
	my @cfg 			= DBOps::DB_configuration( $MP_ID, $CONF );
	my $ULS_load_rev  	= $cfg[4];									# current revision number for ULS loads
	my $loc_water_depth = DBOps::DB_sea_floor_level( $MP_ID );		# water depth at current location
	my @ULS_loads 		= DBOps::DB_loads( $MP_ID, $ULS_load_rev, "ULS" );	# array containing the ULS loads at current location
		
	# sort ULS loads by elevation (descending = beginning at pile head):
	my @sorted_ULS_loads 	= sort { $b->[1] <=> $a->[1] } @ULS_loads;
	$loc_water_depth 		= -abs($loc_water_depth); 					# ensure that water depth is negative
		
	# looking for first row number below mudline
	my $first_row_below_mud = 0;  # first row of array below mudline (for interpolation)
	foreach my $l (@sorted_ULS_loads ){
		my $elev = @{$l}[1];
		if ($elev > $loc_water_depth){
			$first_row_below_mud += 1;
		};
	}
		
	# interpolation of applied loads at mudline, if necessary:
	my $shear				= 0;	# characteristic shear force at mudline
	my $overturn_moment 	= 0;	# characteristic overturning moment at mudline
	my $axial 				= 0;	# characteristic axial load at mudline
	my $torsion_moment 		= 0;	# design torsinal moment at mudline
		
	if($sorted_ULS_loads[$first_row_below_mud][1] == $loc_water_depth){	# no interpolation necessary, if loads are specified at mudline
		
		$shear 				= ${$sorted_ULS_loads[$first_row_below_mud]}[12];
		$overturn_moment 	= ${$sorted_ULS_loads[$first_row_below_mud]}[14];
		$axial 				= ${$sorted_ULS_loads[$first_row_below_mud]}[6];
		$torsion_moment 	= ${$sorted_ULS_loads[$first_row_below_mud]}[9];
	}
	elsif($sorted_ULS_loads[$first_row_below_mud][1] < $loc_water_depth){ # otherwise interpolation is necessary
		
		my $delta_elev_DB 	= ${$sorted_ULS_loads[$first_row_below_mud]}[1] - ${$sorted_ULS_loads[$first_row_below_mud-1]}[1];
		my $delta_elev_wd 	= $loc_water_depth - ${$sorted_ULS_loads[$first_row_below_mud-1]}[1] ;
		$shear 				= ${$sorted_ULS_loads[$first_row_below_mud-1]}[12] + (${$sorted_ULS_loads[$first_row_below_mud]}[12] - 
								${$sorted_ULS_loads[$first_row_below_mud-1]}[12]) / $delta_elev_DB * $delta_elev_wd;
		$overturn_moment 	= ${$sorted_ULS_loads[$first_row_below_mud-1]}[14] + (${$sorted_ULS_loads[$first_row_below_mud]}[14] - 
								${$sorted_ULS_loads[$first_row_below_mud-1]}[14]) / $delta_elev_DB * $delta_elev_wd;
		$axial 				= ${$sorted_ULS_loads[$first_row_below_mud-1]}[6] + (${$sorted_ULS_loads[$first_row_below_mud]}[6] - 
								${$sorted_ULS_loads[$first_row_below_mud-1]}[6]) / $delta_elev_DB * $delta_elev_wd;
		$torsion_moment 	= ${$sorted_ULS_loads[$first_row_below_mud-1]}[9] + (${$sorted_ULS_loads[$first_row_below_mud]}[9] - 
								${$sorted_ULS_loads[$first_row_below_mud-1]}[9]) / $delta_elev_DB * $delta_elev_wd;
	}
		
		
	# ----------------------------------------------------------------------------
	# increase of loads due to 2nd order effect and Alstom-Loadfactor
	# ----------------------------------------------------------------------------
	$shear 				= $Alstom_load_factor * $shear; 			# increase of 1 % (only Alstom-Loadfactor)
	$overturn_moment 	= $sec_order_factor * $overturn_moment;		# increase of 10 % (2nd order effect and Alstom-Loadfactor)
		
		
	# ----------------------------------------------------------------------------
	# formatting and preparing values for table
	# ----------------------------------------------------------------------------
	
	#formatting the values:
	my $shear_str 			= sprintf( "%.1f", $shear ); # get_top_mp( $INPUT_LOC ) );
	my $overturn_moment_str = sprintf( "%.1f", $overturn_moment ); # get_bot_mp( $INPUT_LOC ) );
	my $axial_str 			= sprintf( "%.1f", $axial ); # get_top_mp( $INPUT_LOC ) );
	my $torsion_moment_str  = sprintf( "%.1f", $torsion_moment ); # get_bot_mp( $INPUT_LOC ) );

	my @mrow0 	= ( $shear_str, $overturn_moment_str, $torsion_moment_str, $axial_str);
	my @mdata 	= ( );			# initialise array containing table data
	push( @mdata, \@mrow0 );
		
	# defining table headings
	my @mtitle = ( "Applied loads at mudline","l" );
	my $mcol0 = "\\parbox{1.5cm}{Location }";
	my $mcol1 = "\\parbox{1.5cm}{Shear [kN]}";
	my $mcol2 = "\\parbox{2.5cm}{Overturning moment\$^{1)}\$ [kNm]}";
	my $mcol3 = "\\parbox{1.8cm}{Torsional moment\$^{2)}\$ [kNm]}";
	my $mcol4 = "\\parbox{1.5cm}{Axial load [kN]}";

	my @mcol_names = ( $mcol0,$mcol1,$mcol2,$mcol3,$mcol4 );
	my @mrow_names = ( $MP_ID_STR );
		
		
	# ----------------------------------------------------------------------------
	# generating text of section:
	# ----------------------------------------------------------------------------
	$code .= Latex::build_subsection( $title, $ref_subsections );	#title 
	$code .= Latex::build_par("The loads applied in the determination of pile embedment length and axial capacity check are shown in the following table.");
	# generation of table
	$code .= Latex::build_tabular( \@mdata,\@mcol_names,\@mrow_names, \@mtitle );
	$code .= "\n\n";
	# generation of table foot notes:
	#$code .= "\\textrm{\\textit{ \$^{1)}\$ the shear force includes an increase of 1.0\\% due to load factors provided by Alstom.}}\n\n";
	$code .= "\\textrm{\\textit{ \$^{1)}\$ the overturning moment includes an increase of 6.0\\% to account for the effect of theory 
				of second order.}}\n\n";
	$code .= "\\textrm{\\textit{ \$^{2)}\$ the torsional moment includes a psf according to the Siemens load document.}}\n\n";
	$code .= "\\vspace{1ex}\n\n";

	
	# ----------------------------------------------------------------------------
	# returning variables
	# ----------------------------------------------------------------------------
	return( $code, $axial, $torsion_moment);
}


# ----------------------------------------------------------------------------
# $geo_pile_length_code 	= build_geo_pile_length( $ref_sections );
# ----------------------------------------------------------------------------
	
sub build_geo_pile_length
{
	my $ref_subsections = shift;
	
	# ----------------------------------------------------------------------------
	# general definitions:
	# ----------------------------------------------------------------------------
	my $code = "";
	my $sec_count = 2; 		# needed? ->MNKM
	my $sec_count_old = 0; 	# needed? ->MNKM
	my $title = "Pile embedment length";
	
	
	# ----------------------------------------------------------------------------
	# formatting used values:
	# ----------------------------------------------------------------------------
	my $perm_rot_allow_str = sprintf( "%.2f", $perm_rot_allow) ;
	
	
	# ----------------------------------------------------------------------------
	# text of subsection:
	# ----------------------------------------------------------------------------
	$code .= Latex::build_subsection( $title, $ref_subsections );	
	$code .= Latex::build_par("Pile embedment length of the monopile foundations will be calculated with respect to 
			the following requirements as described in the Design Brief in detail:");
	$code .= Latex::build_par("\\begin\{itemize\} \\item Critical pile length \\item Accumulated rotation 
			\$\\leq$perm_rot_allow_str\^\\circ\$ \\item Full plastification of soil \\end\{itemize\}");
	
	
	# ----------------------------------------------------------------------------
	# critical pile length:
	# ----------------------------------------------------------------------------
	
	# importing and formatting used values:
	my $crit_pile_length 		= DBOps::DB_soil_results( $MP_ID, $CONF, "crit_pile_length");	# [m] critical pile embedment length
	my $crit_pile_length_str 	= sprintf( "%.1f", $crit_pile_length);							# [m] critical pile embedment length, formatted text
	my $crit_len_criterion_str 	= sprintf( "%.0f", $crit_len_criterion);						# [%] critical pile embedment length criterion, formatted text
	
	my $file_name_crit_len 		= $MP_ID . "-crit_pile_length.png"; # Critical pile length plot under ULS loads
	my $file_path_crit_len 		= DBOps::DB_file( $MP_ID, $CONF, "crit_pile_length_plot", "soil_results_geo", "APP/input/", $file_name_crit_len );
	my $caption_crit_len 		= "Relative change in pile head rotation under characteristic ULS loads as function of the embedment length at location $MP_ID_STR, best estimate soil properties.";
	my $label_crit_len 			= "label_crit_length";
	
	my $file_name_deflec 		= $MP_ID . "-hori_defl.png"; # deflection plot under ULS loads
	my $file_path_deflec  		= DBOps::DB_file( $MP_ID, $CONF, "hori_defl_plot", "soil_results_geo", "APP/input/", $file_name_deflec );
	my $caption_deflec 			= "Pile deflection over depth under characteristic ULS loads at location $MP_ID_STR, best estimate soil properties.";
	my $label_deflec 			= "label_deflection";
	
	
	# generating the text
	$code .= Latex::build_par("\\textbf\{Critical pile length\}");	# title
	$code .= Latex::build_par("Figure \\ref\{$label_crit_len\} illustrates that an embedment length of $crit_pile_length_str m fulfils 
			the critical pile length criteria, which requires that the change in pile head rotation relative to the minimum pile head 
			rotation derived for an excessively long pile is less than $crit_len_criterion_str\\%.");
	$code .= Latex::build_par("The corresponding deflection curve of a pile with an embedment length of $crit_pile_length_str m at location
			 $MP_ID_STR is given in Figure \\ref\{$label_deflec\}.");
	$code .= Latex::insert_graphic($file_path_crit_len, $caption_crit_len, $label_crit_len);	# insert plot of relative pile head rotation
	$code .= Latex::insert_graphic($file_path_deflec, $caption_deflec, $label_deflec);			# insert pile deflection plot
	
	
	
	# ----------------------------------------------------------------------------
	# Accumulated rotation:
	# ----------------------------------------------------------------------------
	
	# importing and formatting used values:
	my $perm_rot_mudline_ULS 		= DBOps::DB_soil_results( $MP_ID, $CONF, "perm_rot_mudline");		# [°] permanent rotation at mudline under ULS loads
	my $perm_rot_mudline_oper 		= DBOps::DB_soil_results( $MP_ID, $CONF, "perm_rot_mudline_oper");	# [°] permanent rotation at mudline under SLS loads
	my $perm_rot_mudline_ULS_str 	= sprintf( "%.3f", abs($perm_rot_mudline_ULS));						# [°] permanent rotation at mudline under ULS loads, formatted text
	my $perm_rot_mudline_oper_str	= sprintf( "%.3f", abs($perm_rot_mudline_oper));					# [°] permanent rotation at mudline under SLS loads, formatted text
	my $perm_rot_mudline_str	 	= sprintf( "%.3f", abs($perm_rot_mudline_ULS + $perm_rot_mudline_oper)); # [°] total permanent rotation at mudline, formatted text
	
	my $file_name_perm_rot_ULS 		= $MP_ID . "-perm_rot_mudline.png"; # Plot of permanent rotation at mudline under ULS loads
	my $file_path_perm_rot_ULS  	= DBOps::DB_file( $MP_ID, $CONF, "perm_rot_mudline_plot", "soil_results_geo", "APP/input/", $file_name_perm_rot_ULS );
	my $caption_perm_rot_ULS 		= "Estimate of permanent pile head rotation at seabed under characteristic ULS loads at location $MP_ID_STR, best estimate soil properties.";
	my $label_perm_rot_ULS 			= "label_perm_rot_ULS";
	
	my $file_name_perm_rot_oper 	= $MP_ID . "-perm_rot_mudline_oper.png"; # Plot of permanent rotation at mudline under SLS loads
	my $file_path_perm_rot_oper  	= DBOps::DB_file( $MP_ID, $CONF, "perm_rot_mudline_oper_plot", "soil_results_geo", "APP/input/", $file_name_perm_rot_oper );
	my $caption_perm_rot_oper 		= "Estimate of permanent pile head rotation at seabed under operational loads at location $MP_ID_STR, best estimate soil properties.";
	my $label_perm_rot_oper 		= "label_perm_rot_oper";
	
	
	# generating the text
	$code .= Latex::build_par("\\textbf\{Accumulated rotation\}");	# title
	$code .= Latex::build_par("The permanent pile head rotations at seabed consist of two contributions, the first from ULS 
			loading ($perm_rot_mudline_ULS_str\$\^\\circ\$) and the second from operational loading ($perm_rot_mudline_oper_str\$\^\\circ\$) 
			as shown in Figure \\ref\{$label_perm_rot_ULS\} and \\ref\{$label_perm_rot_oper\}. In case of ULS loading standard cyclic p-y curves 
			are applied which are valid up to 100 loading cycles. In case of operational loading the number of cycles is much larger than 100. In order to account for this larger number of load cycles, the cyclic
			the A-factor is modified as described in the Design Brief. The loads under operational conditions are 
			estimated to 30 \\% of the characteristic ULS loads. The total permanent pile head 
			rotations are estimated to $perm_rot_mudline_str\$\^\\circ\$. Hence, the predicted permanent pile head rotation 
			is less than the allowance of $perm_rot_allow_str\$\^\\circ\$.");
	$code .= Latex::insert_graphic($file_path_perm_rot_ULS, $caption_perm_rot_ULS, $label_perm_rot_ULS);		# insert permanent rotation plot under ULS loads
	$code .= Latex::insert_graphic($file_path_perm_rot_oper, $caption_perm_rot_oper, $label_perm_rot_oper);		# insert permanent rotation plot under SLS loads
	
	
	# check, if total permanent rotation is within the allowed range:
	if (abs($perm_rot_mudline_ULS + $perm_rot_mudline_oper)>$perm_rot_allow){
		$code .= Latex::build_par("\\textbf\{WARNING: permanent rotation criterion not fulfilled!\}");
	}
	
	
	
	
	# ----------------------------------------------------------------------------
	# Full plastification of soil
	# ----------------------------------------------------------------------------
	
	# importing and formatting used values:
	my $lat_utilisation 	= DBOps::DB_soil_results( $MP_ID, $CONF, "lat_util_ratio");		# [-] lateral utilisation ratio
	my $lat_res_char	 	= DBOps::DB_soil_results( $MP_ID, $CONF, "lat_capa");			# [kN] characteristic lateral resistance
	my $lat_force_char	 	= DBOps::DB_soil_results( $MP_ID, $CONF, "lat_force");		# [kN] characteristic mobilised lateral force
	my $lat_utilisation_str = sprintf( "%.2f", $lat_utilisation);							# [-] lateral utilisation ratio, formatted text
	my $lat_res_char_str 	= sprintf( "%.1f", $lat_res_char);								# [kN] characteristic lateral resistance, formatted text
	my $lat_force_char_str 	= sprintf( "%.1f", $lat_force_char);							# [kN] characteristic mobilised lateral force, formatted text
	
	# generating the text
	$code .= Latex::build_par("\\textbf\{Full plastification of soil\}");	# title
	$code .= Latex::build_par("In order to satisfy full plastification of soil resistance against local failure and 
					against global failure has to be verified. According to EA Piles (13.8.1) a check against local failure 
					is not necessary, when applying the p-y curves as described in the Design Brief: 
					\\textit{''If a modulus of subgrade reaction method with non-linear p-y curves (p-y method) is adopted, this condition is 
					normally met automatically, if the p-y curves include a strength limitation.''} 
					Global failure meaning verification of lateral pile capacity is achieved by 
					satisfying the following equation as described in the Design Brief:");
	$code .= Latex::build_par("\\[ B\_\\mathrm\{h,d\}=\\gamma\_Q B\_\\mathrm\{h,c\} =\\gamma\_Q \\int\\limits_0\\limits^\{z\_\{0\}\} p\_\\mathrm\{mob\}\\,
					\\mathrm\{d\}z \\leq R\_\\mathrm\{ult,d,cyc\} = \\frac\{R\_\\mathrm\{ult,c,cyc\} \}\{\\gamma\_\\mathrm\{R,e\}\}  = \\frac\{1\}\{\\gamma\_\\mathrm\{R,e\}\}\\int\\limits_0
					\\limits^\{z\_\{0\}\}e\_\\mathrm\{ph,c,cyc\}\\,\\mathrm\{d\}z\\]");
	$code .= Latex::build_par("with: \\begin\{itemize\} ");
	$code .= Latex::build_par("\\item \$p\_\\mathrm\{mob\}	 \$ characteristic mobilized soil pressure  ");
	$code .= Latex::build_par("\\item \$e\_\\mathrm\{ph,c,cyc\}\$ characteristic maximum mobilisable passive pressure 
					accounting for a reduction due to cyclic loading");
	$code .= Latex::build_par("\\end\{itemize\}");
	$code .= Latex::build_par("and \$\\gamma\_\\mathrm\{Q\} = $psf_lat_load\$ and \$\\gamma\_\\mathrm\{R,e\} = $psf_lat_res\$ 
					as well as \$z\_0\$ as the point with zero horizontal displacement (rotation point). Applying the equation 
					above for the monopile at location $MP_ID_STR (best estimate soil properties) an integrated mobilised lateral force of \$B\_\\mathrm\{h,c\}=$lat_force_char_str\~\\mathrm\{kN\}\$ 
					and a integrated lateral resistance of \$R\_\\mathrm\{ult,c,cyc\}=$lat_res_char_str\~\\mathrm\{kN\}\$ is determined. This results in a lateral 
					utilisation ratio of \$UR\_\\mathrm\{lat\} = $lat_utilisation_str\$.");
	
	
	# check, if lateral utilisation ratio is < 1.0:
	if ($lat_utilisation>1.00){
		$code .= Latex::build_par("\\textbf\{WARNING: Lateral capacity not sufficient!\}");
	}
	
	return $code;
}


# ----------------------------------------------------------------------------
# $geo_axial_check 		= build_geo_axial( $ref_sections, $axial_char_force, $torsional_design_moment);
# ----------------------------------------------------------------------------

sub build_geo_axial
{
	my $ref_subsections = shift;
	my $axial_char_force = shift;
	my $torsional_design_moment = shift;
	
	# ----------------------------------------------------------------------------
	# general definitions:
	# ----------------------------------------------------------------------------
	my $code 			= "";
	my $sec_count 		= 2; 		# needed? ->MNKM
	my $sec_count_old 	= 0; 	# needed? ->MNKM
	my $title 			= "Axial capacity";	# title of subsection
	
	
	# ----------------------------------------------------------------------------
	# importing and formatting used values:
	# ----------------------------------------------------------------------------
	my $axial_util_ratio 			= DBOps::DB_soil_results( $MP_ID, $CONF, "axial_util_ratio");	# [-] axial utilisation ratio
	my $axial_capa		 			= DBOps::DB_soil_results( $MP_ID, $CONF, "axial_capa_comp");			# [kN] characteristic axial pile capacity
	my $axial_util_ratio_str 		= sprintf( "%.2f", $axial_util_ratio);							# [-] axial utilisation ratio, formatted text
	my $axial_capa_str 				= sprintf( "%.1f", $axial_capa);								# [kN] characteristic axial pile capacity, formatted text
	my $axial_char_force_str  		= sprintf( "%.1f", $axial_char_force);							# [kN] characteristic axial load, formatted text
	my $torsional_design_moment_str	= sprintf( "%.1f", $torsional_design_moment);					# [kNm] design torsional moment, formatted text
	
	
	# ----------------------------------------------------------------------------
	# generating the text of subsection:
	# ----------------------------------------------------------------------------
	$code .= Latex::build_subsection( $title, $ref_subsections );	# title
	$code .= Latex::build_par("The ULS axial capacity check according to design approach 2 (GEO-2) of DIN1054 as described in the 
				Design Brief will be based on the unit base resistance and unit shaft friction given in the design soil profiles in the Soil and Foundation
				 Expertise. Verification of axial pile capacity is achieved by satisfying the following equation:");
	# equation:			
	$code .= Latex::build_par("\\[ E\_\\mathrm\{v,d\}=\\gamma\_Q E\_\\mathrm\{v,c\} \\leq R\_\\mathrm\{v,d\}= \\frac\{R\_\\mathrm\{v,c\}\}\{\\gamma\_t\}  \\]");
	$code .= Latex::build_par("with: \\begin\{itemize\} ");
	$code .= Latex::build_par("\\item \$E\_\\mathrm\{v,c\}	 \$ characteristic axial force  ");
	$code .= Latex::build_par("\\item \$R\_\\mathrm\{v,c\}   \$ characteristic axial resistance ");
	$code .= Latex::build_par("\\end\{itemize\}");
	$code .= Latex::build_par("and \$\\gamma\_\\mathrm\{Q\} = $psf_axial_load\$ and \$\\gamma\_\\mathrm\{t\} = $psf_axial_res\$. 
				For the monopile at location $MP_ID_STR a characteristic axial loading at mudline of 
				\$E\_\\mathrm\{v,c\}=$axial_char_force_str\~\\mathrm\{kN\}\$ is applied. Considering the soil profile at $MP_ID_STR an 
				axial ultimate resistance of \$R\_\\mathrm\{v,c\}=$axial_capa_str\~\\mathrm\{kN\}\$ and an axial utilisation ratio 
				of \$UR\_\\mathrm\{axial\} = $axial_util_ratio_str\$ is determined.");
	
	
	# ----------------------------------------------------------------------------
	# check, if axial utilisation ratio is < 1.0:
	# ----------------------------------------------------------------------------
	if ($axial_util_ratio>1.00){
		$code .= Latex::build_par("\\textbf\{WARNING: Axial capacity not sufficient!\}");
	}
	
	return $code;
}





sub build_str_design
{
	my $ref_sections = shift;
	
	my $code = "";
	
	my $sec_count = 2;
	my $sec_count_old = 0;
	
	my $title = "Structural design";
	
	$code .= Latex::build_section( $title, $ref_sections );
	
	$code .= Latex::build_par("The following pages show the results of the structural 
				design checks.");
	$code .= Latex::build_par("The first result sheet shows the fatigue limit state (FLS)
				design checks for the attachments to the MP and TP. Each attachment 
				has a unique identifier (Id). 
				In the table LIST OF ATTACHMENTS all checked attachments 
				are listed. A detailed description of the provided values can be found 
				in the Design Report." );
	$code .= Latex::build_par("The second result sheet shows the fatigue limit state (FLS) 
				design checks for the circumferential welds between the can sections 
				of the TP and MP.
				The design check is done separately for the inside and for the outside. 
				A detailed description of the provided values can be found in the 
				Design Report." );
	$code .= Latex::build_par("The third result sheet shows the ultimate limit state (ULS) 
				design checks. A detailed description of the provided values can be found 
				in the Design Report." );
	
	
	return $code;
}




# end of file create-appendix.pl

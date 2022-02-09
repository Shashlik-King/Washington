#!/usr/bin/perl -w

#*******************************************************************************
#    	Script for inserting  data of given Excel Sheet into LADB
#    	2010/04/29
#    	(c) Christian Ross, Martin Kelm
#
#*******************************************************************************
#		Has this script been fully adapted for use in 273_RaceBank?
#		(yes/no)
#		yes
#
#		Date: 120710
#		User: ros
#*******************************************************************************
#
#		Add. comments:
#
#		User: Ros
#		Date: 120710
#		Text: inserted default die option so that script cannot be run 
#			accidentally. Adjust file paths, names and tables to the .xls and
#			location you want to insert to. Then turn off die option.
# 		
#		
#*******************************************************************************

use strict;
use DBI;
use lib "/home/perl/lib";
use Spreadsheet::ParseExcel;
#use SMB_access;
#use program_utilities;
#use Database;
#use date;

#******************************************************************************
# general parameters
#******************************************************************************

# turn off die option when setting is done 
die;

use constant DB    => "rbdb";

use constant TABLE => "water_depths";

use constant XLSNAME => 
				"RB-D-EN-127-0000-000000-033_water_depths_import_2012-07-10";

use constant SHEETNAME=>"Table1";

my $user 					= getpwuid($<);

#+--------------------------------------------------------------------------------+#
#INPUT DATA

#my $filename = "101126_Bodendaten_aus_Baugrundgutachten.xls";
my $filename 		= XLSNAME.".xls";

my $win_path		= "/DAT/Prj/273_Race_Bank_OWF/04 Berechnungen/40-Datenbank/";

my $input_dir       = "input";

my $lin_path		= "$input_dir/";

my $target_path 	= $lin_path.$filename;



#+------------------------------------------------------------------------------+#

my $cp_direction = "windows2linux";

SMB_copy($win_path, $lin_path, $filename, $cp_direction);

#Definition of Excel Access
my $parser		= Spreadsheet::ParseExcel->new();
my $workbook	= $parser->Parse("$target_path");
my $worksheet	= $workbook->worksheet(SHEETNAME);


#Execution of Subroutines 
my ($row_start, $row_end, $col_min, $col_max,$row_title) = Position($parser, $workbook, $worksheet);
my ($ref_set_tbl, $ref_colnames) = Set_Table($row_start, $row_end, $col_min, $col_max, $row_title);
return_data($ref_set_tbl, $user, $ref_colnames);

###############Subroutine Defintion###################

######Determine Startposition#########
#my ($row_start, $row_end, $col_min, $col_max) = Position($parser, $workbook, $worksheet);
sub Position
{
	my $parser		= shift;
	my $workbook	= shift;
	my $worksheet	= shift;
	my ($row_min, $row_max ) = $worksheet->row_range();
	my ($col_min, $col_max ) = $worksheet->col_range();
	my $row_start;
	my $col_start;
	my $row_end;
	my $col_end;
	my $row_title;

	for my $row ($row_min .. $row_max)
	{
		for my $col ($col_min .. $col_max) 
		{
			my $cell = $worksheet->get_cell($row, $col);
			next unless $cell;
			my $val= $cell->value();

			if ($val=~ /id\s*/)
			{	
				$row_title = $row;
				$row_start = $row+1;
			}
				$row_end = $row_max;
		}
	}
return ($row_start, $row_end, $col_min, $col_max, $row_title);
}


#################Save Soil Data in Table#####################
#my ($ref_set_tbl) = Set_Table($row_start, $row_end, $col_min, $col_max);
sub Set_Table
{
	my $row_start	= shift;
	my $row_end	= shift;
	my $col_min	= shift;
	my $col_max	= shift;
	my $row_title;
	if (@_){$row_title = shift;}
	
	my @table=();
	my @col_names;
	my @zeile;
	# set column names
	if ($row_title){
		for my $col ($col_min .. $col_max) 
		{
			my $cell = $worksheet->get_cell($row_title, $col);
			next unless $cell;
			my $val= $cell->value();
			push (@col_names, $val);
		}
	}
	
	# set data table
	for my $mono_row ($row_start .. $row_end)
	{
		@zeile=();
		for my $col ($col_min .. $col_max) 
		{
			my $cell = $worksheet->get_cell($mono_row, $col);
			next unless $cell;
			my $val= $cell->value();
			push (@zeile, $val);
		}
		push (@table, [@zeile]);
	}


	#####die ersten beiden und letzten beiden Zeilen des Arrays löschen#####
#	pop @table;
#	pop @table;
#	pop @table;
#	shift @table;
#	shift @table;
	#####die ersten beiden und letzten beiden Zeilen des Arrays löschen#####
	return (\@table, \@col_names);
}



#################Return of Data##################
#Return_data($ref_set_tbl,$user);
sub return_data
{
	my $ref_set_tbl	= shift;
	my $user		= shift;
	my $ref_colnames;
	if(@_){$ref_colnames = shift;}

	my @col_names = ( "id",	"diameter",	"uls_shear_force",	"uls_moment",	"fls_shear_force","fls_moment","status","responsible", "inserted_by", "timestamp" );
	
	if ($ref_colnames){
		@col_names = @{$ref_colnames};
	}	
	

	foreach my $row (@{$ref_set_tbl})
	{
# 		print "@{${$ref_set_tbl}[$n]}\n";
		my (@data_row)
			= @{$row};
		
#		my $resp	= "COWI";
		my $ins		= "$user";
#		my $stat 	= "prelim";
		my $timest	= now( "%Y-%m-%d %H:%M:%S" );
		
		push(@data_row, ($ins,$timest));

		print STDERR "@data_row\n";
#		print "$id	$diameter	$uls_shear_force	$uls_moment	$fls_shear_force	$fls_moment $stat, $resp, $ins, $timest\n";
		
		# check whether given values have 3 decimals
		#check_number_of_decimals( $height, 3, "can height" );

		my @data_tbl = ();
#		my @data_row = ( $id,$diameter,$uls_shear_force,$uls_moment,$fls_shear_force,$fls_moment, $stat, $resp, $ins, $timest);
		push( @data_tbl, \@data_row );			
		
		my $table = DB . "." . TABLE;
	
		
		my $result = Database::Write( \@data_tbl, $table,"force" );
		
	}
#	print "\n\n mp cans inserted for status $qa_status \n\n "  #inserted by hof
}

#check_number_of_decimals( $height, 3, "can height" );
#sub check_number_of_decimals
#{
#	my $value        = shift;
#	my $num_decimals = shift;
#	my $text         = shift;
#	
#	# after the dot
#	my $pos = index( $value, ".",  );
#	my $decimals = substr( $value, $pos+1 );
#	
#	unless( length($decimals) == $num_decimals ){
#		die "FATAL: $text of $value has not enough digits! It should be $num_decimals!";
#	}
#}
#
#
#sub print_stars
#{
#	my $stars .= "*"x80;
#	print "$stars\n";
#}

sub usage
{
	my $prog = __FILE__;
	print STDERR "$prog <location>\n";
	exit 0;
}




#!/usr/bin/perl -w

#*******************************************************************************
#    	Script for inserting  data of given Excel Sheet into DB
#    	2010/04/29
#    	(c) Christian Ross, Martin Kelm
#
#*******************************************************************************
#*******************************************************************************
# HISTORY
# Date: 110106
# User: Ros
# Text: QM Block inserted 
#
#
# Date: 110127
# User: Ros
# Text: Script has been adapted for Riffgatt purposes. Edit data in INPUT 
# section if files are to be inserted from elsewhere  
#
#
# User: Ros
# Date: 110127
# Text: DONE : lib changed to standard lib
#
#
# User: Ros
# Date: 110228
# Text: Script adjusted to current naming convention of project. BE AWARE THAT 
# DATANAMES NEED TO STAY THE SAME! Otherwise INPUT DATA has to be adjusted.
#
# User: Ros
# Date: 120621
# Text: Script adjusted to 273_race_bank
#
#*******************************************************************************
# Has this script been adapted to current project? If so, DIE OPTION maybe 
# deleted!
# Adapted (yes/no):
# no
#
#print STDERR "SCRIPT NEEDS TO BE ADAPTED TO CURRENT PROJECT!", die;
#
#*******************************************************************************

use strict;
use DBI;
use lib "/home/perl/lib";
use Spreadsheet::ParseExcel;
#use SMB_access;
#use Database;
#use program_utilities;
#use date;

#*******************************************************************************
# general parameters
#*******************************************************************************

use constant DB    			=> "rbdb";
use constant TP_CANS		=> "tp_cans";
use constant GROUPING_LOAD 	=> "grouping_load";

#our $rev       = "rev00";
# TODO: RACE BANK: get qa status from status table
our $qa_status = "Q1";

our $status      = "testing"; # prelim or final

our $Access_file = $ENV{HOME} . "/.smb.cnf";

#+----------------------------------------------------------------------------+#
#INPUT DATA

if( @ARGV < 1 ){ usage(); }
my @LOCS = @ARGV;
foreach my $INPUT_LOC (@LOCS){
	
	$INPUT_LOC = uc( $INPUT_LOC );
	
	# Extract location number from location name 
	my $locnumber = $INPUT_LOC; 
	$locnumber =~ s/T([0-9]+)/$1/i;
	
	# get user	
	my ($user) = getpwuid($<);
	
	# getting cluster
	# TODO: RACE BANK: adjust to project environment
#	my $what 	= "cluster";
#	my $from 	= DB.".".GROUPING_LOAD;
#	my $where 	= "id=\"$INPUT_LOC\"";
#	my (@result) = Select($what,$from,$where);
#	my $cluster = ${$result[0]}[0];
	
	# TODO: RACE BANK: adjust to project environment
	# find actual structural rev from status 
#	my ( $rev,$struct_rev)= find_struct_rev($INPUT_LOC);
	my ( $rev,$struct_rev) = ("00","F3");
	
	# paths
	my $win_path = "G:/DAT/Prj/273_Race_Bank_OWF/04 Berechnungen/40-Datenbank/01 TESTING";
	
	my $filename	= "WTG$locnumber-DD-weights-rev$rev-$struct_rev.xls";
	
	my $lin_path   	= "input";
	
	my $target_path	= "$lin_path/$filename";
	
	my $cp_direction= "windows2linux";
	
	# create input if not existent
	unless (-e $lin_path){
		mkdir ($lin_path) || die $!;
		my $mode = 0770;
		chmod $mode, "$lin_path";
		print STDERR "\nNew location \"/$lin_path\" has been created!\n";		
	}
	
#+----------------------------------------------------------------------------+#
	
	my $res = SMB_copy($win_path, $lin_path, $filename, $cp_direction);
	
	unless ($res == 0){die "FILE $win_path/$filename NOT FOUND!\n"};
	
	#Definition of Excel Access
	my $parser	= Spreadsheet::ParseExcel->new();
	my $workbook	= $parser->Parse("$target_path");
	# TODO: RACE BANK adjust to project environment
	my $worksheet	= $workbook->worksheet("LC WTG$locnumber");
	
	
	#Execution of Subroutines 
	my ($row_start, $row_end, $col_min, $col_max) 
		= Position($parser, $workbook, $worksheet);
	my ($ref_set_tbl) 
		= Set_Table($row_start, $row_end, $col_min, $col_max, $worksheet);
	Return_data($ref_set_tbl, $INPUT_LOC, $user,$rev, $worksheet);
	
	print STDERR "\n"; 
	print STDERR 	"DB STATUS $struct_rev, REVISION $rev OF $INPUT_LOC HAS ".
					"BEEN INSERTED\n\n";
}

###############Subroutine Defintion###################

######Determine Startposition#########
#my ($row_start, $row_end, $col_min, $col_max) = Position($parser, $workbook, 
#$worksheet);
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

	for my $row ($row_min .. $row_max)
	{
		for my $col ($col_min .. $col_max) 
		{
			my $cell = $worksheet->get_cell($row, $col);
			next unless $cell;
			my $val= $cell->value();

			if ($val=~ /Transition Piece\s*\xd8\s*.*/)
			{
				$row_start = $row;
			}
			elsif ($val=~ /Length of Transition Piece\s*/)
			{
				$row_end = $row;
			}
		}
	}
return ($row_start, $row_end, $col_min, $col_max);
}


#################Save Monopile Block in table#####################
#my ($ref_set_tbl) = Set_Table($row_start, $row_end, $col_min, $col_max);
sub Set_Table
{
	my $row_start	= shift;
	my $row_end	= shift;
	my $col_min	= shift;
	my $col_max	= shift;
	my $worksheet=shift;
	my @table=();
	my @zeile;
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
	pop @table;
	pop @table;
	pop @table;
	shift @table;
	shift @table;
	#####die ersten beiden und letzten beiden Zeilen des Arrays löschen#####
	return (\@table);
}



#################Return of Data##################
#Return_data($ref_set_tbl, $INPUT_LOC, $user);
sub Return_data
{
	my $ref_set_tbl	= shift;
	my $id		= shift;
	my $user	= shift;
	my $rev		= shift;
	my $worksheet= shift;
	$id		= "$id";
	my $n=0;
	
	print STDERR "DELETING TP CANS FOR $id AND REVISION $rev\n";
	my $com = 	"DELETE FROM ".DB.".".TP_CANS.
				" WHERE id=\"$id\" AND revision=\"$rev\";";
	print STDERR "$com\n";
	Query($com);
	
	foreach my $wert (@{$ref_set_tbl})
	{
# 		print "${${$ref_set_tbl}[$n]}[1]\n";
		my $can_id	= (${${$ref_set_tbl}[$n]}[0]);
		my $top	    = (${${$ref_set_tbl}[$n]}[4])/1000.0;
 		my $bottom	= (${${$ref_set_tbl}[$n]}[5])/1000.0;
		my $height	= (${${$ref_set_tbl}[$n]}[7]);
		my $ca      = (${${$ref_set_tbl}[$n]}[9]);
		my $thick	= (${${$ref_set_tbl}[$n]}[10]);
		my $gr		= (${${$ref_set_tbl}[$n]}[13]);
		my $steel 	= (${${$ref_set_tbl}[$n]}[14]);
		my $SCF		= 0.0;
		my $stat	= $status;
		my $resp	= "IMS";
		my $ins		= "$user";
		my $timest	= now( "%Y-%m-%d %H:%M:%S" );

		print 	"$id, $can_id, $rev, $top, $bottom, $height, $thick, $ca, $gr,".
				" $steel, $stat, $resp, $ins, $timest\n"; #$SCF
		
		# check whether given values have 3 decimals
		check_number_of_decimals( $height, 3, "can height" );
		
		my @data_tbl = ();
		my @data_row = ( $id, $can_id, $rev,$top, $bottom, $height, $thick, $ca, 
							$gr, $steel, $stat, $resp, $ins, $timest ); #$SCF
		push( @data_tbl, \@data_row );
			
		my $table = DB . "." . TP_CANS;
		
		my $result = Database::Write( \@data_tbl, $table, "Force" );
		
		$n++;
	}
	print "\n\n tp cans inserted for status $qa_status \n\n "
}

sub usage
{
	my $prog = __FILE__;
	print STDERR "$prog <location>\n";# <revision>\n";
	exit 0;
}

sub find_struct_rev {
	my $loc = shift;
	
	my $revision;	
	
	my $what = "*";
	my $from = DB.".status";
	my $where = "id = \'$loc\'";
	
	my @results = Select($what,$from,$where);
		
	my $count = 0;
	foreach my $entry (@{$results[0]}){
		last if ($entry=~/-/i); 
		$count++;
	}
		
	
	my @results2 = Query("SHOW COLUMNS FROM $from;");
	my $count2 	= 0;
	my @columns;
	foreach my $result2 ( @results2 ){	
		push ( @columns, ${$results2[$count2]}[0] );
		$count2++;
	}	
	
	my $status = $columns[$count-1];
	
	if ($status=~ /^B|^F/i){
		$revision = "00";
	}
	elsif($status=~ /^S/i){
		$revision = "01";
	}
	elsif($status=~ /^T/i){
		$revision = "02";
	}
	
	print STDERR "CURRENT DB STATUS OF $loc IS $status\n";
	
	return ($revision,$status);
}




#!/usr/bin/perl -w
use strict;
#******************************************************************************
# insert-file-geo.pl
#******************************************************************************
# script for inserting pdf/png files into database
#******************************************************************************
# Martin Kelm, COWI, 2015
# date: 2015-11-16
# History:  2016-05-19 (Thorben Hamann) addition for saving pda plots (4 primary keys)
#******************************************************************************

use lib "../../../lib/";
use lib "../commons/";
use DBOps;
use Defs;
use date;


# initialising global variables
our $MP_ID       = "";
our $CONF        = "";
our $hammer      = "";
our $soil        = "";
our $table       = "";
our $name    	 = "";
our $file        = "";

# start of main routine
MAIN: {
	
	# getting input arguments
	if( @ARGV == 5 ){ # in case we have 5 input arguments (2 primary keys)
		$MP_ID    = $ARGV[0]; # primary key 1
		$CONF     = $ARGV[1]; # primary key 2
		$table    = $ARGV[2];
		$name     = $ARGV[3];
		$file     = $ARGV[4];
	}
	elsif(@ARGV == 7){ # in case we have 7 input arguments (4 primary keys)
		$MP_ID    = $ARGV[0]; # primary key 1
		$CONF     = $ARGV[1]; # primary key 2
		$hammer   = $ARGV[2]; # primary key 3
		$soil     = $ARGV[3]; # primary key 4
		$table    = $ARGV[4];
		$name     = $ARGV[5];
		$file     = $ARGV[6];
	}
	else{ usage() } # different number of input arguments -> error message

	DB_connect( Defs::DB );
	
	
	# open specified file, switch to binary mode and save it to $file_data
	my $file_data;
	open( my $fh2, $file ) or die $!;
	binmode( $fh2 );
	read( $fh2, $file_data, -s $fh2 );
	
	
	# date of inserting file to database
	my $insert_datetime = now( "%Y-%m-%d %H:%M:%S" );
	
	
	# name of the column, where the file shall be inserted
	my $col_name_file = $name;
	
	
	# inserting file into database
	if( @ARGV == 5 ){ # in case we have 5 input argument (2 primary keys)
		DBOps::DB_insert_file( $MP_ID, $CONF, $col_name_file, $file_data, $table);
	}
	elsif(@ARGV == 7){ # in case we have 7 input argument (4 primary keys)
		DBOps::DB_insert_file_pda( $MP_ID, $CONF, $hammer, $soil, $col_name_file, $file_data, $table);
	}
	
	
	DB_disconnect();
	
}


###############################################################################
# subroutines
###############################################################################


sub file_date_time_of_modification
{
	my $filename = shift;
	
	# get status of given filename
	my  ($dev,$ino,$mode,$nlink,$uid,$gid,$rdev,$size,$atime,$mtime,$ctime,$blksize,$blocks) = stat($filename);
	
	# convert modification time given in epoch into ISO format: YYYY-MM-DD hh:mm:ss
	my  ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime($mtime);
	$year += 1900;  # year is given since 1900
	$mon  += 1;	# month is given in the range 0..11
	
	return "$year-$mon-$mday $hour:$min:$sec";
}


sub usage
{
	my $msg = __FILE__;
	print STDERR "usage: $msg <mp_id> <conf> <table> <col_name> <file>\n";
	exit 0;
}

# end of file insert-p-y-curves-pdf.pl

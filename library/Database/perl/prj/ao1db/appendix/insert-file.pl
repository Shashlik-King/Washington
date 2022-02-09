#!/usr/bin/perl -w
use strict;
#******************************************************************************
# insert-file.pl
#******************************************************************************
# script for inserting pdf files into database
#******************************************************************************
# Martin Kelm, COWI, 2015
# date: 2015-11-16
#******************************************************************************

use lib "../../../lib/";
use lib "../commons/";
use DBOps;
use Defs;
use date;


our $MP_ID       = "";
our $CONF        = "";
our $table       = "";
our $name    = "";
our $file        = "";

MAIN: {
	
	if( @ARGV != 5 ){ usage() }
	else{
		$MP_ID    = $ARGV[0];
		$CONF     = $ARGV[1];
		$table    = $ARGV[2];
		$name     = $ARGV[3];
		$file     = $ARGV[4];
	}
	
	DB_connect( Defs::DB );
	
	my $datetime_dwg = file_date_time_of_modification( $file );
	
	my $file_data;
	open( my $fh2, $file ) or die $!;
	binmode( $fh2 );
	read( $fh2, $file_data, -s $fh2 );
	
	#TODO: needed?
	my $insert_datetime = now( "%Y-%m-%d %H:%M:%S" );
	
	my $pos = rindex( $file, "/" );
	if( $pos == -1 ){ $pos = rindex( $file, "\\" ); }
	my $file_name = substr( $file, $pos+1 );
	
	my $col_name_file = $name . "_file";
	my $col_name_name = $name . "_name";
	my $col_name_date = $name . "_date";
	
	DBOps::DB_insert_file( $MP_ID, $CONF, $col_name_file, $file_data, $table, $file_name, $col_name_name, $datetime_dwg, $col_name_date );
	
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

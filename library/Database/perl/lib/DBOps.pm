#! /usr/bin/perl -w
use strict;
##############################################################################
#
# File   :  DBOps.pm
# History:  2013-11-29 (Martin Kelm) first implementation
#           2014-08-29 (Martin Kelm) updated documentation
#           2014-10-08 (Martin Kelm) introduced subroutine DB_insert_data_set
#           2014-10-09 (Martin Kelm) introduced various other subroutines
#           2015-02-02 (Martin Kelm) Removed various TODOs
#           2015-06-12 (Martin Kelm) subroutine DB_splashzone added
#           2015-09-25 (Martin Kelm) subroutine DB_sea_floor_level added
#           2015-10-19 (Martin Kelm) new subroutines
#           2015-11-16 (Martin Kelm) new subroutines (DB_file, DB_insert_file)
#           2015-11-23 (Martin Kelm) new subroutines (DB_loads)
#           2016-02-24 (Martin Kelm) more documentation
#           2016-03-17 (Martin Kelm) new subroutines (DB_cluster)
#           2016-03-22 (Thorben Hamann) new subroutines (DB_soil_data, 
#                      DB_soil_results, DB_soil_p_mob, DB_soil_deflec, DB_py_data)
#           2016-03-23 (Martin Kelm) minor typos corrected
#           2016-04-06 (Martin Kelm) new subroutines (DB_locations)
#           2016-04-30 (Martin Kelm) insert_file subroutine updated
#			2016-05-19 (Thorben Hamann) new subroutines for pile driveability
#					   (DB_pda_results, DB_insert_file_pda, DB_file_pda)
#           2016-05-30 (Martin Kelm) added obsolete warning to DB_type_tp
#           2016-06-28 (Martin Kelm) dummy
#           2017-08-01 (Martin Kelm) closing function renamed to end
#           2018-09-12 (Martin Kelm) DB_cluster(ing) subroutine, rev added
#           2019-07-29 (Martin Kelm) DB_loads_input added.
#           2019-08-08 (Martin Kelm) DB_result_sheet_value added.
##############################################################################
#
#  Provides database operations to MP/TP database
#  uses DBConnect module.
#
##############################################################################
# (c) 2013 - 2017 COWI A/S - All Rights Reserved.
# This source file is part of the COPILOD software platform.
# 
# NOTICE:
# All information contained herein is, and remains the property of COWI A/S. 
# The intellectual and technical concepts contained herein are proprietary to 
# COWI A/S.
# Dissemination of this information or reproduction of this material
# is strictly forbidden unless prior written permission is obtained
# from COWI A/S.
##############################################################################

=head1 NAME

DBOps - interact with MySQL database

=cut

##############################################################################

package DBOps;

require Exporter;
use vars qw( $VERSION );

$VERSION = 0.5;

use lib ".";

use DBConnect;
use Debug;

# table names of database
use constant TBL_MP_MAIN_DIMENSIONS	=> "mp_main_dimensions";
use constant TBL_TP_MAIN_DIMENSIONS	=> "tp_main_dimensions";

use constant TBL_ADD_MASSES          => "add_masses";
use constant TBL_MP_CANS			 => "mp_cans";
use constant TBL_MARINE_GROWTH       => "marine_growth";
use constant TBL_TP_CANS			 => "tp_cans";
use constant TBL_WTG_LOCS			 => "wtg_locations";
use constant TBL_WATER_DEPTHS 		 => "water_depths";
use constant TBL_DOC_REV_HIST		 => "document_revision_history";
use constant TBL_DOC_NO              => "doc_no";
use constant TBL_STATUS				 => "status";
use constant TBL_LOADS_FLS           => "loads_FLS";
use constant TBL_LOADS_FLS_INPUT     => "loads_FLS_input";
use constant TBL_LOADS_FLS_FACTORS   => "loads_FLS_factors";
use constant TBL_LOADS_ULS           => "loads_ULS";
use constant TBL_LOADS_ULS_INPUT     => "loads_ULS_input";
use constant TBL_LOADS_ULS_FACTORS   => "loads_ULS_factors";
use constant TBL_GENERAL_DATA        => "general_data";
use constant TBL_RESULT_SHEET_VALUES => "result_sheet_values";
use constant TBL_CONFIGURATIONS      => "configurations";

use constant COL_SEABED_LEVEL        => "water_depth";
use constant COL_DSOIL_LEVEL         => "design_soil_profile_wd";


use constant COL_ALL                 => "*";
use constant COL_PILE_TOP            => "pile_top";
use constant COL_PILE_TIP            => "pile_tip";
use constant COL_BOTTOM_TP           => "bottom_tp";
use constant COL_TYPE_TP             => "type_tp";
use constant COL_INSERTED_BY         => "inserted_by";
use constant COL_STATUS				 => "status";
use constant COL_PREPARER			 => "preparer";
use constant COL_INTERFACE_LEVEL     => "interface_level";
use constant COL_VALUE               => "value";
use constant COL_STR_VALUE           => "str_value";


our (@ISA, @EXPORT, $VERSION);
@ISA 	= qw(Exporter);
@EXPORT = qw(DB_connect DB_disconnect DB_select);

our $gl_dbcon;
our $db_name = "";
our $gl_connected = 0;


# ----------------------------------------------------------------------------
# $connected = DB_connect( $db_name );
# ----------------------------------------------------------------------------
# Connects to the given database and returns a value of 1 if sucessful.
# 
# ----------------------------------------------------------------------------

sub DB_connect
{
	Debug::check_num_params( 1, \@_, (caller(0))[3], __LINE__ );
	
	$db_name = shift;
	
	my $OS = $^O;
	
	my $home_path = "";
	if( $OS eq "MSWin32" ){
		$home_path = $ENV{HOMEDRIVE} . $ENV{HOMEPATH} . "\\";
	}
	else{
		$home_path = $ENV{HOME} . "/";
	}
	
	my $access_file = "$home_path.my.cnf";
	$gl_dbcon = DBConnect->new( $db_name, $access_file, $db_name );
	my $err = $gl_dbcon->open();
	if($err){
		my $msg = $gl_dbcon->error();
		Debug::fatal( $msg, (caller(0))[3], __LINE__ );
	}
	$gl_connected = 1;
	
	return $gl_connected;
}


# ----------------------------------------------------------------------------
# DB_disconnect( );
# ----------------------------------------------------------------------------

sub DB_disconnect
{
	Debug::check_num_params( 0, \@_, (caller(0))[3], __LINE__ );
	
	$gl_dbcon->end();
	$gl_connected = 0;
}


# ----------------------------------------------------------------------------
# my @data = DB_added_masses_data( $id, $rev );
# ----------------------------------------------------------------------------
# Returns added masses data
# ----------------------------------------------------------------------------
# Params:	$id          Tower id or MP id
#           $rev         Revision
#
# Returns:	@data        Double array
# ----------------------------------------------------------------------------

sub DB_added_masses_data
{
	# ------------------------------------------------------------------------
	# Params
	Debug::check_num_params( 2, \@_, (caller(0))[3], __LINE__ );
	
	my $id    = shift;
	my $rev   = shift;
	
	# End of params
	# ------------------------------------------------------------------------
	
	
	#local variables:
	my @result = ();
	my $table  = TBL_ADD_MASSES;
	
	DBOps::DB_select( \@result, "*", $db_name.".".$table, "id = \"$id\" AND rev = \"$rev\"" );
	
	unless( @result){
		my $msg = "Result set empty. Requested value(s) not in database.table [$table] for revision [$rev].";
		Debug::fatal( $msg, (caller(0))[3], __LINE__ );
	}
	
	my $count = 0;
	our @am_data = ();
	foreach (@result){
		my @t_am_data = @{$result[$count]};
		#remove first 3 columns
		shift @t_am_data; # remove id
		shift @t_am_data; # remove mass_id
		shift @t_am_data; # remove rev
		# remove last 5 columns
		pop @t_am_data; # description
		pop @t_am_data; # timestamp
		pop @t_am_data; # inserted_by
		pop @t_am_data; # responsible
		pop @t_am_data; # status
		
		$am_data[$count] = \@t_am_data;
		$count++;
	}
	
	return @am_data;
}


# ----------------------------------------------------------------------------
# $value = DBOps::DB_general_data( $id, $type );
# ----------------------------------------------------------------------------
# This subroutine extracts values out of table "general_data".
# It returns.
# ----------------------------------------------------------------------------
# Params:
#
# Returns:  $value			
# ----------------------------------------------------------------------------

sub DB_general_data
{
	# ------------------------------------------------------------------------
	# Params
	Debug::check_num_params( 2, \@_, (caller(0))[3], __LINE__ );
	
	my $id   = shift;
	my $type = shift;
	
	Debug::check_string_in_array( $type, ["DOUBLE","STRING"], "type", (caller(0))[3], __LINE__ );
	
	# End of params
	# ------------------------------------------------------------------------
	
	my $val = "NDEF";
	
	my @result = ();
	
	my $col_name = COL_VALUE;
	
	if( $type eq "STRING" ){
		$col_name = COL_STR_VALUE;
	}
	
	DB_select( \@result, $col_name, $db_name.".".TBL_GENERAL_DATA, "id = \"$id\"" );
	
	unless( @result == 1 ){
		my $msg = "Result set empty. Requested value(s) ". $id . " not in table " . TBL_GENERAL_DATA . ".";
		Debug::warn( $msg, (caller(0))[3], __LINE__ );
		
	}
	else{
		
		if( defined ${$result[0]}[0] ){
			($val) = @{$result[0]};
		}
	}
	
	return $val;
}


# ----------------------------------------------------------------------------
# DB_attachments( $type, $rev );
# ----------------------------------------------------------------------------

sub DB_attachments
{
	Debug::check_num_params( 2, \@_, (caller(0))[3], __LINE__ );
	
	my $type = shift;
	my $rev  = shift;
	
	my @result = ();
	
	my $table = "attachments";
	
	DB_select( \@result, "*", $db_name.".".$table, "type_tp = \"$type\" and ".
								"rev = \"$rev\"" );
	
	unless( @result != 0 ){
		my $msg = "Result set empty. Requested value(s) ". $type . " for Rev " . $rev . " not in table " . $table . ".";
		Debug::fatal( $msg, (caller(0))[3], __LINE__ );
	}
	
	my $count = 0;
	
	my @data = ();
	
	
	foreach (@result){
		my @tdata = @{$result[$count]};
		my $id = shift @tdata;
		#my $can_id	= shift @tdata;
		shift @tdata; # remove rev
		shift @tdata; # remove type
		# remove last columns
		pop @tdata; # timestamp
		pop @tdata; # inserted_by
		pop @tdata; # responsible
		pop @tdata; # status
		
		# put back id
		unshift( @tdata, $id );
		
		push(@data, \@tdata);
		$count++;
	}
	
	return @data;
}


# ----------------------------------------------------------------------------
# DB_configuration( $id, $rev );
# ----------------------------------------------------------------------------
# Returns configuration as an array
# $id	identifier
# $rev	revision number of configuration
# ----------------------------------------------------------------------------

sub DB_configuration
{
	Debug::check_num_params( 2, \@_, (caller(0))[3], __LINE__ );
	my $loc = shift;
	my $rev = shift;
	
	
	my @result = ();
	my $table = $db_name . "." . TBL_CONFIGURATIONS;
	my $columns = COL_ALL;
	
	DB_select( \@result, $columns, $table, "id = \"$loc\" and ".
								"rev = \"$rev\"" );
	
	unless( @result ){
		my $msg = "Result set empty. Requested value(s) ". $columns . " not in table $table for id [$loc] and configuration [$rev].";
		Debug::fatal( $msg, (caller(0))[3], __LINE__ ); 
	}
	
	return @{$result[0]};
}


# ----------------------------------------------------------------------------
# DB_configuration_latest_rev( $id );
# ----------------------------------------------------------------------------
# Returns revision number of latest configuration
# $id	identifier
# ----------------------------------------------------------------------------

sub DB_configuration_latest_rev
{
	Debug::check_num_params( 1, \@_, (caller(0))[3], __LINE__ );
	my $loc = shift;
	
	my @result = ();
	my $table = $db_name .  "." . TBL_CONFIGURATIONS;
	my $columns = "rev";
	
	DB_select( \@result, $columns, $table, "id = \"$loc\"" );
	
	unless( @result ){
		my $msg = "Result set empty. Requested value(s) ". $columns . " not in table $table for id [$loc].";
		Debug::fatal( $msg, (caller(0))[3], __LINE__ ); 
	}
	
	my $rev = -1;
	foreach my $res (@result){
		my( $trev ) = @{$res};
		
		if( $trev > $rev ){ $rev = $trev; }
	}
	
	return $rev;
}


# ----------------------------------------------------------------------------
# DB_clustering( );
# ----------------------------------------------------------------------------
# Returns double array of the clustering table in database
# ----------------------------------------------------------------------------

sub DB_clustering
{
	Debug::check_num_opt_params( [0,1], \@_, (caller(0))[3], __LINE__ );
	
	my $rev = "01";
	my $has_rev = 0;
	if( @_ ){ $rev = shift; $has_rev = 1; }
	
	my @result = ();
	my $table = "clustering";
	my $columns = COL_ALL;
	my $where   = "";
	
	if( $has_rev ){
		$where = "rev = \"$rev\"";
	}
	
	DB_select( \@result, $columns, $db_name.".".$table, $where );
	
	unless( @result ){
		my $msg = "Result set empty. Requested value(s) ". $columns . " not in table $table in database.";
		Debug::fatal( $msg, (caller(0))[3], __LINE__ ); 
	}
	
	return @result;
}


# ----------------------------------------------------------------------------
# DB_cluster( $MP_ID );
# ----------------------------------------------------------------------------
# Returns name of cluster for given MP_ID
# $MP_ID identifier
# ----------------------------------------------------------------------------

sub DB_cluster
{
	Debug::check_num_opt_params( [1,2], \@_, (caller(0))[3], __LINE__ );
	
	my $MP_ID = shift;
	
	my $rev = "01";
	my $has_rev = 0;
	if( @_ ){ $rev = shift; $has_rev = 1; }
	
	
	my @result = ();
	my $table = "clustering";
	my $columns = "cluster";
	my $where = "id = \"$MP_ID\"";
	
	if( $has_rev ){
		$where .= " AND rev = \"$rev\"";
	}
	
	DB_select( \@result, $columns, $db_name.".".$table, $where );
	
	unless( @result ){
		my $msg = "Result set empty. Requested value(s) ". $columns . " not in table $table in database $db_name for id $MP_ID.";
		Debug::fatal( $msg, (caller(0))[3], __LINE__ ); 
	}
	
	return $result[0][0];
}


# ----------------------------------------------------------------------------
# DB_insert_data_set( $ref_data, $table );
# ----------------------------------------------------------------------------
# Returns nothing. If error, script is stopped. Uses INSERT SQL-command.
# 
# ----------------------------------------------------------------------------

sub DB_insert_data_set
{
	Debug::check_num_params( 2, \@_, (caller(0))[3], __LINE__ );
	
	my $ref_data = shift;
	my $table    = shift;
	
	unless( $gl_connected ){
		my $msg = "Not connected to database. Call DB_connect( <db> ) first"; 
		Debug::fatal( $msg, (caller(0))[3], __LINE__ ); 
	}
	
	foreach my $dat (@{$ref_data}){
		
		my @data = @{$dat};
		
		# put each datum into quotation marks
		my @tmp = ();
		foreach my $d (@data){
			push( @tmp, "\"$d\"" );
		}
		
		my $err = $gl_dbcon->insert( \@tmp, $table );
		if($err){
			my $msg = $gl_dbcon->error();
			Debug::fatal( $msg, (caller(0))[3], __LINE__ );
		}
	}
}

# ----------------------------------------------------------------------------
# DB_replace_data_set( $ref_data, $table );
# ----------------------------------------------------------------------------
# Returns nothing. If error, script is stopped. Uses REPLACE SQL-command.
# 
# ----------------------------------------------------------------------------

sub DB_replace_data_set
{
	Debug::check_num_params( 2, \@_, (caller(0))[3], __LINE__ );
	
	my $ref_data = shift;
	my $table    = shift;
	
	unless( $gl_connected ){
		my $msg = "Not connected to database. Call DB_connect( <db> ) first"; 
		Debug::fatal( $msg, (caller(0))[3], __LINE__ ); 
	}
	
	foreach my $dat (@{$ref_data}){
		
		my @data = @{$dat};
		
		# put each datum into quotation marks
		my @tmp = ();
		foreach my $d (@data){
			push( @tmp, "\"$d\"" );
		}
		
		# $what, $table, $where, $data
		my $err = $gl_dbcon->replace( \@tmp, $table );
		if($err){
			my $msg = $gl_dbcon->error();
			Debug::fatal( $msg, (caller(0))[3], __LINE__ );
		}
	}
}


# ----------------------------------------------------------------------------
# $val = DBOps::DB_result_sheet_value( $id, $type );
# ----------------------------------------------------------------------------
# This subroutine returns value for given $id. If not defined it returns "NDEF".
# ----------------------------------------------------------------------------
# Params:	$id         Name of item. Must be one in $gl_values.
#           $type       Either "DOUBLE" or "STRING"
#
# Returns:	$ret        Value.
# ----------------------------------------------------------------------------

sub DB_result_sheet_value
{
	# ------------------------------------------------------------------------
	# Params
	Debug::check_num_params( 2, \@_, (caller(0))[3], __LINE__ );
	
	my $id   = shift;
	my $type = shift;
	
	Debug::check_string_in_array( $type, ["DOUBLE","STRING"], "type", (caller(0))[3], __LINE__ );
	
	# End of params
	# ------------------------------------------------------------------------
	
	my $val = "NDEF";
	my @result = ();
	
	my $col_name = COL_VALUE;
	
	if( $type eq "STRING" ){
		$col_name = COL_STR_VALUE;
	}
	
	DB_select( \@result, $col_name, $db_name.".".TBL_RESULT_SHEET_VALUES, "id = \"$id\"" );
	
	unless( @result == 1 ){
		my $msg = "Result set empty. Requested value(s) ". $id . " not in table " . TBL_GENERAL_DATA . ".";
		Debug::warn( $msg, (caller(0))[3], __LINE__ );
	}
	else{
		($val) = @{$result[0]};
	}
	
	return $val;
}


# ----------------------------------------------------------------------------
# DB_insert_file( $id, $conf, $col_name, $col_file, $db_table );
# ----------------------------------------------------------------------------

sub DB_insert_file
{
	my $nparams = Debug::check_num_opt_params( [5,7,9], \@_, (caller(0))[3], __LINE__ );
	
	my $loc       = shift;
	my $conf      = shift;
	my $col_name  = shift;
	my $col_file  = shift;
	my $table     = shift;
	
	my $file_name = "";
	my $col_name_file = "";
	if( @_ ){ $file_name = shift; }
	if( @_ ){ $col_name_file = shift; }
		
	my $date_time = "";
	my $col_name_date = "";
	if( @_ ){ $date_time = shift; }
	if( @_ ){ $col_name_date = shift; }
	
	unless( $gl_connected ){
		my $msg = "Not connected to database. Call DB_connect( <db> ) first"; 
		Debug::fatal( $msg, (caller(0))[3], __LINE__ ); 
	}
	
	
	
	my $what  = $col_name;
	
	my $where = "id = \"$loc\" and rev = \"$conf\"";
	
	my $data  = $col_file;
	
	if( $nparams == 5 ){
		my $err = $gl_dbcon->update( $what, $table, $where, $data );
		if($err){
			my $msg = "UPDATE FAILED ($nparams) " . $gl_dbcon->error();
			Debug::fatal( $msg, (caller(0))[3], __LINE__ );
		}
	}
	elsif( $nparams == 7 ){
		my $err = $gl_dbcon->update( $what, $table, $where, $data );
		if($err){
			my $msg = "UPDATE FAILED ($nparams) " . $gl_dbcon->error();
			Debug::fatal( $msg, (caller(0))[3], __LINE__ );
		}
		
		$err = $gl_dbcon->update( $col_name_file, $table, $where, $file_name );
		if($err){
			my $msg = "UPDATE FAILED ($nparams) " . $gl_dbcon->error();
			Debug::fatal( $msg, (caller(0))[3], __LINE__ );
		}
	}
	elsif( $nparams == 9 ){
		my $err = $gl_dbcon->update( $what, $table, $where, $data );
		if($err){
			my $msg = "UPDATE FAILED ($nparams) " . $gl_dbcon->error();
			Debug::fatal( $msg, (caller(0))[3], __LINE__ );
		}
		
		$err = $gl_dbcon->update( $col_name_file, $table, $where, $file_name );
		if($err){
			my $msg = "UPDATE FAILED ($nparams) " . $gl_dbcon->error();
			Debug::fatal( $msg, (caller(0))[3], __LINE__ );
		}
		
		$err = $gl_dbcon->update( $col_name_date, $table, $where, $date_time );
		if($err){
			my $msg = "UPDATE FAILED ($nparams) " . $gl_dbcon->error();
			Debug::fatal( $msg, (caller(0))[3], __LINE__ );
		}
	}
	else{
		die;
	}
}




# ----------------------------------------------------------------------------
# DB_insert_file( $id, $conf, $hammer, $soil_prop,  $col_name, $col_file, $db_table );
# ----------------------------------------------------------------------------
# by THHM 2016-05-19

sub DB_insert_file_pda
{
	# check if we have [7,9 or 11] input arguments
	my $nparams = Debug::check_num_opt_params( [7,9,11], \@_, (caller(0))[3], __LINE__ );
	
	# getting input first 7 arguments (always the case):
	my $loc       = shift;
	my $conf      = shift;
	my $hammer 	  = shift;
	my $soil_prop = shift;
	my $col_name  = shift;
	my $col_file  = shift;
	my $table     = shift;
	
	# get 8th and 9th input argument, if provided
	my $file_name = "";
	my $col_name_file = "";
	if( @_ ){ $file_name = shift; }
	if( @_ ){ $col_name_file = shift; }

	# get 10th and 11th input argument, if provided
	my $date_time = "";
	my $col_name_date = "";
	if( @_ ){ $date_time = shift; }
	if( @_ ){ $col_name_date = shift; }
	
	# make an error message, if there is no connection to database
	unless( $gl_connected ){
		my $msg = "Not connected to database. Call DB_connect( <db> ) first"; 
		Debug::fatal( $msg, (caller(0))[3], __LINE__ ); 
	}
	
# ------------------------------------------	
# upload the file to the database
# ------------------------------------------	

	# column name, where the file shall be inserted
	my $what  = $col_name;

	# filter (applying the primary keys) to select the row, where the file shall be inserted	
	my $where = "id = \"$loc\" and rev = \"$conf\" and hammer_config = \"$hammer\" and soil_prop = \"$soil_prop\"";
	
	# data of file to be inserted
	my $data  = $col_file;
	
	if( $nparams == 7 ){ # in case we have 7 input arguments
		my $err = $gl_dbcon->update( $what, $table, $where, $data );	# save file data
		if($err){
			my $msg = "UPDATE FAILED ($nparams) " . $gl_dbcon->error();
			Debug::fatal( $msg, (caller(0))[3], __LINE__ );
		}
	}
	elsif( $nparams == 9 ){ # in case we have 9 input arguments (additionally: name of file is saved in the column specified by $col_name_file)
		my $err = $gl_dbcon->update( $what, $table, $where, $data );	# save file data
		if($err){
			my $msg = "UPDATE FAILED ($nparams) " . $gl_dbcon->error();
			Debug::fatal( $msg, (caller(0))[3], __LINE__ );
		}
		
		$err = $gl_dbcon->update( $col_name_file, $table, $where, $file_name );	# save file name
		if($err){
			my $msg = "UPDATE FAILED ($nparams) " . $gl_dbcon->error();
			Debug::fatal( $msg, (caller(0))[3], __LINE__ );
		}
	}
	elsif( $nparams == 11 ){# in case we have 11 input arguments (additionally: date is saved in the column specified by $col_name_date)
		my $err = $gl_dbcon->update( $what, $table, $where, $data );	# save file data
		if($err){
			my $msg = "UPDATE FAILED ($nparams) " . $gl_dbcon->error();
			Debug::fatal( $msg, (caller(0))[3], __LINE__ );
		}
		
		$err = $gl_dbcon->update( $col_name_file, $table, $where, $file_name );		# save file name
		if($err){
			my $msg = "UPDATE FAILED ($nparams) " . $gl_dbcon->error();
			Debug::fatal( $msg, (caller(0))[3], __LINE__ );
		}
		
		$err = $gl_dbcon->update( $col_name_date, $table, $where, $date_time );		# save file date
		if($err){
			my $msg = "UPDATE FAILED ($nparams) " . $gl_dbcon->error();
			Debug::fatal( $msg, (caller(0))[3], __LINE__ );
		}
	}
	else{
		die;
	}
}


# ----------------------------------------------------------------------------
# DB_file( $loc, $conf, $col_name, $col_file, $db_table, $path );
# ----------------------------------------------------------------------------
# Returns file name with path
# ----------------------------------------------------------------------------

sub DB_file
{
	Debug::check_num_params( 6, \@_, (caller(0))[3], __LINE__ );
	
	my $loc       = shift;
	my $conf      = shift;
	my $col_file  = shift;
	my $table     = shift;
	my $path      = shift;
	my $file_name = shift;
	
	
	my @result  = ();
	my $columns = "$col_file";
	DB_select( \@result, $columns, $db_name.".".$table, "id = \"$loc\" and ".
								"rev = \"$conf\"" );
	
	unless( @result ){
		my $msg = "Result set empty. Requested value(s) [". $columns . "] not in table [$table] for [$loc], conf: [$conf] in database.";
		Debug::fatal( $msg, (caller(0))[3], __LINE__ );
	}
	
	my $file_name_with_path = $path . $file_name;
	foreach my $res (@result){
		
		my ($file) = @{$res};
				
		open ( OUT, ">$file_name_with_path") or die "$file_name_with_path: ", $!;
		binmode( OUT );
		print OUT $file;
		close( OUT );
	}

	return $file_name_with_path;
}




# ----------------------------------------------------------------------------
# DB_file_pda( $loc, $conf, $hammer, $soil, $col_name, $col_file, $db_table, $path );
# ----------------------------------------------------------------------------

sub DB_file_pda
{
	# by THHM
	Debug::check_num_params( 8, \@_, (caller(0))[3], __LINE__ );
	
	my $loc       = shift;
	my $conf      = shift;
	my $hammer    = shift;
	my $soil      = shift;
	my $col_file  = shift;
	my $table     = shift;
	my $path      = shift;
	my $file_name = shift;
	
	
	my @result  = ();
	my $columns = "$col_file";
	DB_select( \@result, $columns, $db_name.".".$table, "id = \"$loc\" and ".
								"rev = \"$conf\" and hammer_config = \"$hammer\" and soil_prop = \"$soil\"" );
	
	unless( @result ){
		my $msg = "Result set empty. Requested value(s) [". $columns . "] not in table [$table] for [$loc], conf: [$conf] in database.";
		Debug::fatal( $msg, (caller(0))[3], __LINE__ );
	}
	
	my $file_name_with_path = $path . $file_name;
	foreach my $res (@result){
		
		my ($file) = @{$res};
				
		open ( OUT, ">$file_name_with_path") or die "$file_name_with_path: ", $!;
		binmode( OUT );
		print OUT $file;
		close( OUT );
	}

	return $file_name_with_path;
}


# ----------------------------------------------------------------------------
# DB_select( $ref_results, $what, $from, $where );
# ----------------------------------------------------------------------------

sub DB_select
{
	Debug::check_num_params( 4, \@_, (caller(0))[3], __LINE__ );
	
	my $ref_results = shift;
	my $what 	    = shift;
	my $from 	    = shift;
	my $where	    = shift;
	
	unless( $gl_connected ){
		my $msg = "Not connected to database. Call DB_connect( <db> ) first"; 
		Debug::fatal( $msg, (caller(0))[3], __LINE__ ); 
	}
	
	my $err = $gl_dbcon->select(  $ref_results, $what, $from, $where );
	if($err){
		my $msg = $gl_dbcon->error();
		Debug::fatal( $msg, (caller(0))[3], __LINE__ );
	}
}


# ----------------------------------------------------------------------------
# my( $cable_dir_0, $cable_dir_1 ) = DBOps::DB_cable_dirs( $INPUT_LOC );
# ----------------------------------------------------------------------------

sub DB_cable_dirs
{
	Debug::check_num_params( 1, \@_, (caller(0))[3], __LINE__ );
	
	my $loc = shift;
	
	my @result = ();
	
	#TODO:
	my $columns = "cable_dir_0, cable_dir_1, cable_dir_2";
	DB_select( \@result, $columns, $db_name.".".TBL_WTG_LOCS, "id = \"$loc\"" );
	
	#TODO: check if empty
	return @{$result[0]};
}


# ----------------------------------------------------------------------------
# my @results = 	DBOps::DB_document_revisions( $ID, $REV );
# ----------------------------------------------------------------------------
# Returns double array of document revision info for given ID and 
# document revision. Double array should include only one array of info.
# ----------------------------------------------------------------------------

sub DB_document_revisions
{
	Debug::check_num_params( 2, \@_, (caller(0))[3], __LINE__ );
	
	my $loc = shift;
	my $rev = shift;
									
								
	my @result = ();
	my $table = TBL_DOC_REV_HIST;
	my $columns = COL_ALL;
	DB_select( \@result, $columns, $db_name.".".$table, "id = \"$loc\" and ".
								"document_revision = \"$rev\"" );
	
	unless( @result ){
		my $msg = "Result set empty. Requested value(s) [". $columns . "] not in table [$table] for [$loc], rev: [$rev] in database.";
		Debug::fatal( $msg, (caller(0))[3], __LINE__ );
	}
	
	return @result;
}


# ----------------------------------------------------------------------------
# my @results = 	DBOps::DB_doc_revisions( $ID );
# ----------------------------------------------------------------------------
# Returns double array of all document revisions for given ID.
# ----------------------------------------------------------------------------

sub DB_doc_revisions
{
	Debug::check_num_params( 1, \@_, (caller(0))[3], __LINE__ );
	
	my $loc = shift;
									
								
	my @result = ();
	
	my $columns = COL_ALL;
	my $table = $db_name.".".TBL_DOC_REV_HIST;
	DB_select( \@result, $columns, $table, "id = \"$loc\"");
	
	unless( @result ){
		my $msg = "Result set empty. Requested value(s) ". $columns . " for id [" . $loc . "] not in database table [$table].";
		Debug::fatal( $msg, (caller(0))[3], __LINE__ );
	}
	
	return @result;
}


# ----------------------------------------------------------------------------
# my @results = DBOps::DB_doc_data( $doc_id );
# ----------------------------------------------------------------------------
# @results = ($doc_id, $client_id, $title_0, $title_1, $title_2, $status);
# ----------------------------------------------------------------------------

sub DB_doc_data
{
	Debug::check_num_params( 1, \@_, (caller(0))[3], __LINE__ );
	
	my $id  = shift;
									
							
	my @result = ();
	
	my $columns = "doc_id, client_id, title_0, title_1, title_2, status";
	my $table = $db_name.".".TBL_DOC_NO;
	DB_select( \@result, $columns, $table, "id = \"$id\"" );
	
	unless( @result ){
		my $msg = "Result set empty. Requested value(s) ". $columns . " for id [" . $id . "] not in database table [$table].";
		Debug::fatal( $msg, (caller(0))[3], __LINE__ ); 
	}
	
	#TODO: check if empty
	return @{$result[0]};
}


# ----------------------------------------------------------------------------
# my ($preparer, $checker) = DBOps::DB_preparer_checker( $INPUT_LOC, $struct_rev );
# ----------------------------------------------------------------------------

sub DB_preparer_checker
{
	Debug::check_num_params( 2, \@_, (caller(0))[3], __LINE__ );
	
	my $loc        = shift;
	my $struct_rev = shift;
	
														
	my @result = ();
	my $table  = TBL_STATUS;
	my $columns = "preparer, checker";
	DB_select( \@result, $columns, $db_name.".". $table, "id = \"$loc\" and ".
								"structural_revision = \"$struct_rev\"" );
	
	unless( @result ){
		my $msg = "Result set empty. Requested value(s) [". $columns . "] not in table [$table] for [$loc], structural_rev: [$struct_rev] in database.";
		Debug::fatal( $msg, (caller(0))[3], __LINE__ );
	}
	#TODO: check if empty
	return @{$result[0]};
	
}


# ----------------------------------------------------------------------------
# $seabed_level = DB_seabed_level( $INPUT_LOC );
# ----------------------------------------------------------------------------

# obsolete, use DB_sea_floor_level instead
sub DB_seabed_level
{
	Debug::check_num_params( 1, \@_, (caller(0))[3], __LINE__ );
	Debug::warn( "Obsolete subroutine, use DB_sea_floor_level instead",(caller(0))[3], __LINE__ );
	
	my $loc = shift;
	
	return DB_sea_floor_level( $loc );
}


# ----------------------------------------------------------------------------
# $seabed_level = DB_sea_floor_level( $loc );
# ----------------------------------------------------------------------------

sub DB_sea_floor_level
{
	Debug::check_num_params( 1, \@_, (caller(0))[3], __LINE__ );
	
	my $loc = shift;
	
	my @result = ();
	
	DB_select( \@result, COL_SEABED_LEVEL, $db_name.".".TBL_WATER_DEPTHS, "id = \"$loc\"" );
	
	unless( @result == 1 ){
		my $msg = "Result set empty. Requested value(s) ". COL_SEABED_LEVEL . " not in database for location [$loc].";
		Debug::fatal( $msg, (caller(0))[3], __LINE__ );
	}
	my  ($sea_floor_level) = @{$result[0]};
	
	return $sea_floor_level;
}


# ----------------------------------------------------------------------------
# $seabed_level = DB_design_soil_profile_level( $INPUT_LOC );
# ----------------------------------------------------------------------------

sub DB_design_soil_profile_level
{
	Debug::check_num_params( 1, \@_, (caller(0))[3], __LINE__ );
	
	my $loc = shift;
	
	my @result = ();
	
	DB_select( \@result, COL_DSOIL_LEVEL, $db_name.".".TBL_WATER_DEPTHS, "id = \"$loc\"" );
	
	unless( @result == 1 ){
		my $msg = "Result set empty. Requested value(s) ". COL_DSOIL_LEVEL . " not in database for location [$loc].";
		Debug::fatal( $msg, (caller(0))[3], __LINE__ );
	}
	my  ($seabed_level) = @{$result[0]};
	
	return $seabed_level;
}


# ----------------------------------------------------------------------------
# DB_soil_data();
# ----------------------------------------------------------------------------
#TODO: subroutine from THHM

sub DB_soil_data
{
	Debug::check_num_params( 3, \@_, (caller(0))[3], __LINE__ );
	my $loc   = shift;
	my $rev   = shift;
	my $table = shift;
	$table = "soil_data_detail";
	
	#print STDERR "$loc $rev\n";
	
	my @result = ();
	
	DB_select( \@result, "*", $db_name.".".$table, "id = \"$loc\" and ".
								"soil_revision = \"$rev\"" );
	my $count = 0;
	our %tsoil = ();
	foreach (@result){
		my @tsoil_data = @{$result[$count]};
		shift @tsoil_data; # remove loc
		my $layer_id	= shift @tsoil_data;
		shift @tsoil_data; # remove rev
		# remove last 4 columns
		pop @tsoil_data; # timestamp
		pop @tsoil_data; # inserted_by
		pop @tsoil_data; # responsible
		pop @tsoil_data; # status
		
		$tsoil{$layer_id} = \@tsoil_data;
		$count++;
	}
	
	my @soil = ();
	foreach my $layer_id ( sort{$a <=> $b} keys %tsoil ){
		push( @soil, $tsoil{$layer_id} );
	}
	
	return @soil;
}



# ----------------------------------------------------------------------------
# $value = DB_soil_results( $location, $configuration, $requested_type_of_value);
# ----------------------------------------------------------------------------
#TODO: subroutine from THHM

sub DB_soil_results
{
	Debug::check_num_params( 3, \@_, (caller(0))[3], __LINE__ );
	my $loc   		= shift;
	my $rev   		= shift;
	my $variable 	= shift;
	my $table 		= "soil_results_geo";
	
	#print STDERR "$loc $rev\n";
	
	my @result = ();
	
	DB_select( \@result, $variable, $db_name.".".$table, "id = \"$loc\" and ".
								"rev = \"$rev\"" );
	
	unless( @result == 1 ){
		my $msg = "Result set empty. Requested value(s) ". $variable . " not in database for location [$loc] and [$rev].";
		Debug::fatal( $msg, (caller(0))[3], __LINE__ );
	}
	my  ($value) = @{$result[0]};
	
	return $value;
}



# ----------------------------------------------------------------------------
# DB_soil_p_mob( $location, $configuration);
# ----------------------------------------------------------------------------
#TODO: subroutine from THHM

sub DB_soil_p_mob
{
	Debug::check_num_params( 2, \@_, (caller(0))[3], __LINE__ );
	my $loc   		= shift;
	my $rev   		= shift;
	my $table 		= "soil_p_mobilised";
	
	#print STDERR "$loc $rev\n";
	
	my @result = ();
	
	DB_select( \@result, "*", $db_name.".".$table, "id = \"$loc\" and ".
								"rev = \"$rev\" ORDER BY depth_node DESC" );

	my $count = 0;
	our @p_mob = ();
	foreach (@result){
		my @temp_data = @{$result[$count]};
		#remove first two columns
		shift @temp_data; # remove loc
		shift @temp_data; # remove rev
		# remove last 4 columns
		pop @temp_data; # timestamp
		pop @temp_data; # inserted_by
		pop @temp_data; # responsible
		pop @temp_data; # status
		
		$p_mob[$count] = \@temp_data;
		$count++;
	}
	
		
	return \@p_mob;
}




# ----------------------------------------------------------------------------
# DB_soil_deflec( $location, $configuration);
# ----------------------------------------------------------------------------
#TODO: subroutine from THHM

sub DB_soil_deflec
{
	Debug::check_num_params( 2, \@_, (caller(0))[3], __LINE__ );
	my $loc   		= shift;
	my $rev   		= shift;
	my $table 		= "soil_deflec_curve";
	
	#print STDERR "$loc $rev\n";
	
	my @result = ();
	
	DB_select( \@result, "*", $db_name.".".$table, "id = \"$loc\" and ".
								"rev = \"$rev\" ORDER BY depth_node DESC" );

	my $count = 0;
	our @deflec = ();
	foreach (@result){
		my @temp_data = @{$result[$count]};
		#remove first two columns
		shift @temp_data; # remove loc
		shift @temp_data; # remove rev
		# remove last 4 columns
		pop @temp_data; # timestamp
		pop @temp_data; # inserted_by
		pop @temp_data; # responsible
		pop @temp_data; # status
		
		$deflec[$count] = \@temp_data;
		$count++;
	}
	
		
	return \@deflec;
}



# ----------------------------------------------------------------------------
# DB_py_data( $LOC, $CONF, PY_DATA)
# ----------------------------------------------------------------------------
# get p-y data out of database
# ----------------------------------------------------------------------------

sub  DB_py_data
{
	Debug::check_num_params( 3, \@_, (caller(0))[3], __LINE__ );
	my $loc   = shift;
	my $rev   = shift;
	my $table = shift;
	
	
	#local variables:
	my @result = ();
	
	DBOps::DB_select( \@result, "*", $db_name.".".$table, "id = \"$loc\" and ".
								"rev = \"$rev\"" );
	
	unless( @result){
		my $msg = "Result set empty. Requested value(s) not in database.table [$table] for location [$loc], geo res revision [$rev].";
		Debug::fatal( $msg, (caller(0))[3], __LINE__ );
	}
	
	my $count = 0;
	our @py_data = ();
	foreach (@result){
		my @t_py_data = @{$result[$count]};
		#remove first two columns
		shift @t_py_data; # remove loc
		shift @t_py_data; # remove rev
		# remove last 6 columns
		pop @t_py_data; # timestamp
		pop @t_py_data; # inserted_by
		pop @t_py_data; # responsible
		pop @t_py_data; # status
		pop @t_py_data; # stat_cycl
		pop @t_py_data; # model_py
		
		$py_data[$count] = \@t_py_data;
		$count++;
	}
	
	return @py_data;
}



# ----------------------------------------------------------------------------
# DB_pda_results( $LOC, $CONF, pda_results_sum, $soil_prop, $hammer_cfg, $soil_prop, $what)
# ----------------------------------------------------------------------------
# get p-y data out of database
# ----------------------------------------------------------------------------

sub  DB_pda_results
{
	Debug::check_num_opt_params( [4,5,6], \@_, (caller(0))[3], __LINE__ );
	my $loc   		= shift;
	my $rev   		= shift;
	my $table 		= shift;
	my $soil_prop 	= shift;
	
	
	#local variables:
	my @result 		= ();
	my $flag 		= 0; # flag, how/what data shall be taken from database
	my $what		= "*";
	my $hammer_cfg 	= "1";
	
	if( @_ ){ 
		$hammer_cfg    	= shift; 
		$flag			= 1;
	}
	
	if( @_ ){ 
		$what    		= shift; 
	}
	
	if($flag == 0){ # select data, without specifying $hammer_config
		DBOps::DB_select( \@result, $what, $db_name.".".$table, "id = \"$loc\" and ".
								"rev = \"$rev\" and ".
								"soil_prop = \"$soil_prop\" ");
	}
	elsif($flag == 1){ # select data, with specifying $hammer_config
		DBOps::DB_select( \@result, $what, $db_name.".".$table, "id = \"$loc\" and ".
								"rev = \"$rev\" and hammer_config = \"$hammer_cfg\" and ".
								"soil_prop = \"$soil_prop\" ");
			#print STDERR "$what, $db_name.$table, id = \"$loc\" and rev = \"$rev\" and hammer_config = \"$hammer_cfg\" and soil_prop = \"$soil_prop\" \n";
	}
	else{
		Debug::info( "Error in MySQL query for id $loc, revision $rev, hammer configuration 
						$hammer_cfg and $soil_prop soil properties when trying to import data from $table", "MAIN" );
		exit 0;
	}
	
	
	unless( @result){
		my $msg = "Result set empty. Requested value(s) not in database.table [$table] for location [$loc], revision [$rev].";
		Debug::fatal( $msg, (caller(0))[3], __LINE__ );
	}
	
	
	my @pda_results 	= ();
	if ($what eq "*"){
		my $count 		= 0;
		foreach (@result){
			my @t_pda_results = @{$result[$count]};
			#remove first two columns
			shift @t_pda_results; # remove loc
			shift @t_pda_results; # remove rev
			# remove last 4 columns
			pop @t_pda_results; # timestamp
			pop @t_pda_results; # inserted_by
			pop @t_pda_results; # responsible
			pop @t_pda_results; # status
			
			$pda_results[$count] = \@t_pda_results;
			$count++;
		}
	}
	else{
		@pda_results	= @result
	}
	
	
	return @pda_results;
}



# ----------------------------------------------------------------------------
# DB_cpt_data( $LOC, $REV, soil_cpt_data, $type)
# ----------------------------------------------------------------------------
# get cpt data out of database
# ----------------------------------------------------------------------------

sub  DB_cpt_data
{
	Debug::check_num_params( 4, \@_, (caller(0))[3], __LINE__ );
	my $loc   		= shift;
	my $rev   		= shift;
	my $table 		= shift;
	my $type 		= shift;
	
	
	#local variables:
	my @result = ();
	DBOps::DB_select( \@result, "*", $db_name.".".$table, "id = \"$loc\" and ".
								"rev = \"$rev\" and type = \"$type\"");
	my $count = 0;
	our @cpt_data = ();
	foreach (@result){
		my @t_cpt_data = @{$result[$count]};
		#remove first two columns
		shift @t_cpt_data; # remove loc
		shift @t_cpt_data; # remove rev
		# remove last 4 columns
		pop @t_cpt_data; # timestamp
		pop @t_cpt_data; # inserted_by
		pop @t_cpt_data; # responsible
		pop @t_cpt_data; # status
		
		$cpt_data[$count] = \@t_cpt_data;
		$count++;
	}
	
	
	return @cpt_data;
}


# ----------------------------------------------------------------------------
# DB_delete_2( $table, $key0_name, $key0_val, $key1_name, $key1_val );
# ----------------------------------------------------------------------------

sub DB_delete_2{ 
	
	Debug::check_num_params( 5, \@_, (caller(0))[3], __LINE__ );
	
	my $table     = shift;
	my $key0_name = shift;
	my $key0_val  = shift;
	my $key1_name = shift;
	my $key1_val  = shift;
	
	unless( $gl_connected ){
		my $msg = "Not connected to database. Call DB_connect( <db> ) first"; 
		Debug::fatal( $msg, (caller(0))[3], __LINE__ ); 
	}

	my $where = $key0_name . "=\"" . $key0_val . "\" AND " . $key1_name . "=\"" . $key1_val . "\"";
	
	my $err = $gl_dbcon->delete( $table, $where );
	if($err < 0){
		my $msg = $gl_dbcon->error();
		Debug::fatal( $msg, (caller(0))[3], __LINE__ );
	}

}

#$dbcon->delete( $table, $where_condition );

# ----------------------------------------------------------------------------
# ($top, bottom) = DB_splashzone();
# ----------------------------------------------------------------------------

sub DB_splashzone
{
	Debug::check_num_params( 0, \@_, (caller(0))[3], __LINE__ );
	
	my @result0 = ();
	my @result1 = ();
	
	DB_select( \@result0, COL_VALUE, $db_name.".".TBL_GENERAL_DATA, "id = \"SPLASH_ZONE_BOTTOM\"" );
	
	unless( @result0 == 1 ){
		my $msg = "Result set empty. Requested value(s) ". COL_VALUE . " not in database for [SPLASH_ZONE_BOTTOM].";
		Debug::fatal( $msg, (caller(0))[3], __LINE__ );
	}
	
	DB_select( \@result1, COL_VALUE, $db_name.".".TBL_GENERAL_DATA, "id = \"SPLASH_ZONE_TOP\"" );
	
	unless( @result1 == 1 ){
		my $msg = "Result set empty. Requested value(s) ". COL_VALUE . " not in database for [SPLASH_ZONE_TOP].";
		Debug::fatal( $msg, (caller(0))[3], __LINE__ );
	}
	
	my  ($bot) = @{$result0[0]};
	my  ($top) = @{$result1[0]};
	
	return ( $top, $bot );
}


# ----------------------------------------------------------------------------
# my( $pile_top, $pile_bot ) = DB_pile_top_bottom( $INPUT_LOC, $rev );
# ----------------------------------------------------------------------------

sub DB_pile_top_bottom
{
	Debug::check_num_params( 2, \@_, (caller(0))[3], __LINE__ );
	
	my $loc = shift;
	my $rev = shift;
	
	my @result  = ();
	my $columns = COL_PILE_TOP . "," . COL_PILE_TIP;
	my $table   = $db_name.".". TBL_MP_MAIN_DIMENSIONS;
	DB_select( \@result, $columns, $table, "id = \"$loc\" and ".
								"structural_revision = \"$rev\"" );
	
	unless( @result == 1 ){
		my $msg = "Result set empty. Requested value(s) ". $columns . " not in database.table [$table] for location [$loc], structural revision [$rev].";
		Debug::fatal( $msg, (caller(0))[3], __LINE__ );
	}
	
	return @{$result[0]};
}


sub DB_interface_level
{
	Debug::check_num_params( 2, \@_, (caller(0))[3], __LINE__ );
	
	my $loc = shift;
	my $rev = shift;
	
	my @result = ();
	my $columns = COL_INTERFACE_LEVEL;
	my $table   = $db_name.".". TBL_MP_MAIN_DIMENSIONS;
	DB_select( \@result, $columns, $table, "id = \"$loc\" and ".
								"structural_revision = \"$rev\"" );
	
	unless( @result == 1 ){
		my $msg = "Result set empty. Requested value(s) ". $columns . " not in database.table [$table] for location [$loc], structural revision [$rev].";
		Debug::fatal( $msg, (caller(0))[3], __LINE__ );
	}
	
	my( $interface_level ) = @{$result[0]};
	return $interface_level;
}


# ----------------------------------------------------------------------------
# our $bottom_tp = DBOps::DB_tp_bottom( $INPUT_LOC, $REV_IN_DB );
# ----------------------------------------------------------------------------

sub DB_tp_bottom
{
	Debug::check_num_params( 2, \@_, (caller(0))[3], __LINE__ );
	
	my $tp_type = shift;
	my $rev = shift;
	
	my @result = ();
	my $columns = COL_BOTTOM_TP;
	my $table   = $db_name.".". TBL_TP_MAIN_DIMENSIONS;
	DB_select( \@result, $columns, $table, "id = \"$tp_type\" and ".
								"structural_revision = \"$rev\"" );
	
	unless( @result == 1 ){
		my $msg = "Result set empty. Requested value(s) ". $columns . " not in database.table [$table] for location [$tp_type], structural revision [$rev].";
		Debug::fatal( $msg, (caller(0))[3], __LINE__ );
	}
	
	my( $bottom_tp ) = @{$result[0]};
	
	return $bottom_tp;
}

sub DB_tp_top
{
	Debug::check_num_params( 2, \@_, (caller(0))[3], __LINE__ );
	
	my $tp_type = shift;
	my $rev = shift;
	
	my @result = ();
	my $columns = COL_INTERFACE_LEVEL;
	my $table   = $db_name.".". TBL_TP_MAIN_DIMENSIONS;
	DB_select( \@result, $columns, $table, "id = \"$tp_type\" and ".
								"structural_revision = \"$rev\"" );
	
	unless( @result == 1 ){
		my $msg = "Result set empty. Requested value(s) ". $columns . " not in database.table [$table] for location [$tp_type], structural revision [$rev].";
		Debug::fatal( $msg, (caller(0))[3], __LINE__ );
	}
	
	my( $bottom_tp ) = @{$result[0]};
	
	return $bottom_tp;
}


# ----------------------------------------------------------------------------
# our $bottom_tp = DBOps::DB_bottom_tp( $INPUT_LOC, $REV_IN_DB );
# ----------------------------------------------------------------------------

sub DB_bottom_tp
{
	Debug::check_num_params( 2, \@_, (caller(0))[3], __LINE__ );
	
	my $loc = shift;
	my $rev = shift;
	
	my @result = ();
	my $columns = COL_BOTTOM_TP;
	my $table   = $db_name.".". TBL_MP_MAIN_DIMENSIONS;
	DB_select( \@result, $columns, $table, "id = \"$loc\" and ".
								"structural_revision = \"$rev\"" );
	
	unless( @result == 1 ){
		my $msg = "Result set empty. Requested value(s) ". $columns . " not in database.table [$table] for location [$loc], structural revision [$rev].";
		Debug::fatal( $msg, (caller(0))[3], __LINE__ );
	}
	
	my( $bottom_tp ) = @{$result[0]};
	
	return $bottom_tp;
}


# ----------------------------------------------------------------------------
# our $type_tp = DBOps::DB_type_tp( $INPUT_LOC, $REV_IN_DB );
# ----------------------------------------------------------------------------

sub DB_type_tp
{
	Debug::check_num_params( 2, \@_, (caller(0))[3], __LINE__ );
	Debug::warn( "Obsolete subroutine, use table configuration in database instead",(caller(0))[3], __LINE__ );
	
	my $loc = shift;
	my $rev = shift;
	
	my @result = ();
	my $columns = COL_TYPE_TP;
	my $table   = $db_name.".". TBL_MP_MAIN_DIMENSIONS;
	DB_select( \@result, $columns, $table, "id = \"$loc\" and ".
								"structural_revision = \"$rev\"" );
	
	unless( @result == 1 ){
		my $msg = "Result set empty. Requested value(s) ". $columns . " not in database.table [$table] for location [$loc], structural revision [$rev].";
		Debug::fatal( $msg, (caller(0))[3], __LINE__ );
	}
	
	my( $type_tp ) = @{$result[0]};
	
	return $type_tp;
}


# ----------------------------------------------------------------------------
# my @data = DB_marine_growth_data( $id );
# ----------------------------------------------------------------------------
# Returns marine growth
# ----------------------------------------------------------------------------
# Params:	$id          MP id
#
# Returns:	@data        Double array
# ----------------------------------------------------------------------------

sub  DB_marine_growth_data
{
	# ------------------------------------------------------------------------
	# Params
	
	Debug::check_num_params( 1, \@_, (caller(0))[3], __LINE__ );
	my $id    = shift;
	
	# End of params
	# ------------------------------------------------------------------------
	
	my $table = TBL_MARINE_GROWTH;
	
	
	#local variables:
	my @result = ();
	
	DBOps::DB_select( \@result, "*", $db_name.".".$table, "id = \"$id\"" );
	
	unless( @result){
		my $msg = "Result set empty. Requested value(s) not in database.table [$table].";
		Debug::fatal( $msg, (caller(0))[3], __LINE__ );
	}
	
	
	my $count = 0;
	our @mg_data = ();
	foreach (@result){
		my @t_mg_data = @{$result[$count]};
		#remove first two columns
		shift @t_mg_data; # remove loc
		#shift @t_py_data; # remove rev
		# remove last 4 columns
		pop @t_mg_data; # description
		pop @t_mg_data; # timestamp
		pop @t_mg_data; # inserted_by
		pop @t_mg_data; # responsible
		pop @t_mg_data; # status
	
		
		$mg_data[$count] = \@t_mg_data;
		$count++;
	}
	
	return @mg_data;
}



# ----------------------------------------------------------------------------
# our($inserted_by, $status) = DBOps::DB_mp_status( $INPUT_LOC, $REV_IN_DB );
# ----------------------------------------------------------------------------

sub DB_mp_status
{
	Debug::check_num_params( 2, \@_, (caller(0))[3], __LINE__ );
	my $loc = shift;
	my $rev = shift;
	
	my @result = ();
	my $columns = COL_PREPARER . "," . COL_STATUS;
	my $table   = $db_name.".". TBL_STATUS;
	DB_select( \@result, $columns, $table, "id = \"$loc\" and ".
								"structural_revision = \"$rev\"" );
	
	unless( @result == 1 ){
		my $msg = "Result set empty. Requested value(s) ". $columns . " not in database.table [$table] for location [$loc], structural revision [$rev].";
		Debug::fatal( $msg, (caller(0))[3], __LINE__ );
	}
	
	return @{$result[0]};
}


# ----------------------------------------------------------------------------
# @cans = DBOps::DB_pile( $INPUT_LOC, $REV_IN_DB );
# ----------------------------------------------------------------------------
#TODO: merge with DB_tp
sub DB_pile
{
	Debug::check_num_params( 2, \@_, (caller(0))[3], __LINE__ );
	my $loc = shift;
	my $rev = shift;
	
	my @cans = DB_can_data( $loc, $rev, TBL_MP_CANS );
	return @cans;
}


# ----------------------------------------------------------------------------
# @cans = DBOps::DB_tp( $INPUT_LOC, $REV_IN_DB );
# ----------------------------------------------------------------------------
#TODO: merge with DB_pile
sub DB_tp
{
	Debug::check_num_params( 2, \@_, (caller(0))[3], __LINE__ );
	my $loc = shift;
	my $rev = shift;
	
	my @cans = DB_can_data( $loc, $rev, TBL_TP_CANS );
	return @cans;
}

sub DB_can_data
{
	Debug::check_num_params( 3, \@_, (caller(0))[3], __LINE__ );
	my $loc   = shift;
	my $rev   = shift;
	my $table = shift;
	
	#print STDERR "$loc $rev\n";
	
	my @result = ();
	
	DB_select( \@result, "*", $db_name.".".$table, "id = \"$loc\" and ".
								"structural_revision = \"$rev\"" );
	
	unless( @result){
		my $msg = "Result set empty. Requested value(s) not in database.table [$table] for location [$loc], structural revision [$rev].";
		Debug::fatal( $msg, (caller(0))[3], __LINE__ );
	}
	
	my $count = 0;
	our %tcans = ();
	foreach (@result){
		my @tdata = @{$result[$count]};
		shift @tdata; # remove loc
		my $can_id	= shift @tdata;
		shift @tdata; # remove rev
		# remove last columns
		pop @tdata; # timestamp
		pop @tdata; # inserted_by
		pop @tdata; # responsible
		pop @tdata; # status
		
		$tcans{$can_id} = \@tdata;
		$count++;
	}
	
	my @cans = ();
	foreach my $can_id ( sort{$a <=> $b} keys %tcans ){
		push( @cans, $tcans{$can_id} );
	}
	
	return @cans;
}


# ----------------------------------------------------------------------------
# my %piles = DB_all_piles();
# ----------------------------------------------------------------------------
#TODO: introduce structural revision
sub DB_all_piles
{
	Debug::check_num_params( 0, \@_, (caller(0))[3], __LINE__ );
	
	my @results = ();
	my %piles = ();
	DB_select( \@results, "*", TBL_MP_CANS );
	
	foreach my $d (@results){
		my( $loc, $can_id, $struct_rev, $top_od, $bottom_od, $height, $wall_thickness, $corrosion_allowance, $ground, $steel_grade, @tmp ) = @{$d};
		#print "@tmp\n";
		
		#TODO: include corrosion_allowance
		my @geom = ( $top_od, $bottom_od, $height, $wall_thickness/1000.0, $ground, $steel_grade );
		
		if( $piles{$loc} ){
			my %tcans = %{$piles{$loc}};
			$tcans{$can_id} = \@geom;
			$piles{$loc} = \%tcans;
		}
		else{
			my %tcans = ();
			$tcans{$can_id} = \@geom;
			$piles{$loc} = \%tcans;
		}
		
	}
	
	return %piles;
}


# ----------------------------------------------------------------------------
# DB_loads( $loc, $rev, $type );
# ----------------------------------------------------------------------------

sub DB_loads
{
	Debug::check_num_params( 3, \@_, (caller(0))[3], __LINE__ );
	
	my $loc  = shift;
	my $rev  = shift;
	my $type = shift;
	
	Debug::check_string_in_array( $type, ["FLS","ULS"], "Type", (caller(0))[3], __LINE__ );
	
	my $table = "";
	if( $type eq "FLS" ){
		$table = DBOps::TBL_LOADS_FLS;
	}
	elsif( $type eq "ULS" ){
		$table = DBOps::TBL_LOADS_ULS;
	}
	else{
		my $msg = "Unknown type = [$type], please use either [ULS] or [FLS].";
		Debug::fatal( $msg, "MAIN", __LINE__ );
	}
	
	my @results = ();
	DB_select( \@results, "*", $table, "id = \"$loc\" and ".
								"load_rev = \"$rev\"" );
	
	unless( @results ){
		my $msg = "Result set empty. Requested value(s) not in database.table [$table] for location [$loc], load revision [$rev].";
		Debug::fatal( $msg, (caller(0))[3], __LINE__ );
	}
	
	return @results;
}

# ----------------------------------------------------------------------------
# DB_loads_input( $loc, $rev, $type );
# ----------------------------------------------------------------------------

sub DB_loads_input
{
	Debug::check_num_params( 3, \@_, (caller(0))[3], __LINE__ );
	
	my $loc  = shift;
	my $rev  = shift;
	my $type = shift;
	
	Debug::check_string_in_array( $type, ["FLS","ULS"], "Type", (caller(0))[3], __LINE__ );
	
	my $table = "";
	if( $type eq "FLS" ){
		$table = DBOps::TBL_LOADS_FLS_INPUT;
	}
	elsif( $type eq "ULS" ){
		$table = DBOps::TBL_LOADS_ULS_INPUT;
	}
	else{
		my $msg = "Unknown type = [$type], please use either [ULS] or [FLS].";
		Debug::fatal( $msg, "MAIN", __LINE__ );
	}
	
	my @results = ();
	DB_select( \@results, "*", $table, "id = \"$loc\" and ".
								"load_rev = \"$rev\"" );
	
	unless( @results ){
		my $msg = "Result set empty. Requested value(s) not in database.table [$table] for location [$loc], load revision [$rev].";
		Debug::fatal( $msg, (caller(0))[3], __LINE__ );
	}
	
	return @results;	
}



#TODO: put main part into separate subroutine
sub DB_all_loads_FLS
{
	Debug::check_num_params( 1, \@_, (caller(0))[3], __LINE__ );
	my $rev = shift;
	
	my @results = ();
	my %loads = ();
	my $table = TBL_LOADS_FLS;
	DB_select( \@results, "*", $table, "" );
	
	if( @results == 0 ){
		my $msg = "Result set empty. Requested value(s) not in table [$table] for revision [$rev].";
		Debug::fatal( $msg, (caller(0))[3], __LINE__ );
	}
	
	foreach my $d (@results){
		#remove last four items
		my @res = @{$d};
		pop @res;
		pop @res;
		pop @res;
		pop @res;
		
		my( $loc, $elev, $load_rev, @loads ) = @res;
		
		#TODO: 
		unless( $rev eq $load_rev ){ next; }
		
		if( $loads{$loc} ){
			my %loads_loc = %{$loads{$loc}};
			$loads_loc{$elev} = \@loads;
			$loads{$loc} = \%loads_loc;
		}
		else{
			my %loads_loc = ();
			$loads_loc{$elev} = \@loads;
			$loads{$loc} = \%loads_loc;
		}
		
	}
	
	return %loads;
}


# ----------------------------------------------------------------------------
# DB_load_factors_latest_rev( $id, $type );
# ----------------------------------------------------------------------------
# Returns revision number of latest load factor revision
# $id	identifier
# ----------------------------------------------------------------------------

sub DB_load_factors_latest_rev
{
	Debug::check_num_params( 2, \@_, (caller(0))[3], __LINE__ );
	my $loc = shift;
	my $type = shift;
	
	Debug::check_string_in_array( $type, ["FLS","ULS"], "Type", (caller(0))[3], __LINE__ );
	
	my $table = "";
	if( $type eq "FLS" ){
		$table = DBOps::TBL_LOADS_FLS_FACTORS;
	}
	elsif( $type eq "ULS" ){
		$table = DBOps::TBL_LOADS_ULS_FACTORS;
	}
	else{
		my $msg = "Unknown type = [$type], please use either [ULS] or [FLS].";
		Debug::fatal( $msg, "MAIN", __LINE__ );
	}
	
	my @result = ();
	$table = $db_name .  "." . $table;
	my $columns = "load_rev";
	
	DB_select( \@result, $columns, $table, "id = \"$loc\"" );
	
	unless( @result ){
		my $msg = "Result set empty. Requested value(s) ". $columns . " not in table $table for id [$loc].";
		Debug::fatal( $msg, (caller(0))[3], __LINE__ ); 
	}
	
	my $rev = -1;
	foreach my $res (@result){
		my( $trev ) = @{$res};
		
		if( $trev > $rev ){ $rev = $trev; }
	}
	
	return $rev;
}


#
sub DB_load_factors
{
	Debug::check_num_params( 3, \@_, (caller(0))[3], __LINE__ );
	
	my $loc  = shift;
	my $rev  = shift;
	my $type = shift;
	
	Debug::check_string_in_array( $type, ["FLS","ULS"], "Type", (caller(0))[3], __LINE__ );
	
	my $table = "";
	if( $type eq "FLS" ){
		$table = DBOps::TBL_LOADS_FLS_FACTORS;
	}
	elsif( $type eq "ULS" ){
		$table = DBOps::TBL_LOADS_ULS_FACTORS;
	}
	else{
		my $msg = "Unknown type = [$type], please use either [ULS] or [FLS].";
		Debug::fatal( $msg, "MAIN", __LINE__ );
	}
	
	my @results = ();
	$table = $db_name .  "." . $table;
	
	
	DB_select( \@results, "*", $table, "id = \"$loc\" and ".
								"load_rev = \"$rev\"" );
	
	unless( @results ){
		my $msg = "Result set empty. Requested value(s) not in database.table [$table] for location [$loc], load revision [$rev].";
		Debug::fatal( $msg, (caller(0))[3], __LINE__ );
	}
	
	my @factors = ();
	
	my @tmp = @{$results[0]};
	shift @tmp;
	shift @tmp;
	pop @tmp;
	pop @tmp;
	pop @tmp;
	pop @tmp;
	
	
	for( my $i = 0; $i < @tmp/3; $i++  ){
		push( @factors, $tmp[3*$i+1] );
	}
	
	return @factors;
}


#TODO: put main part into separate subroutine
sub DB_all_loads_ULS
{
	Debug::check_num_params( 1, \@_, (caller(0))[3], __LINE__ );
	my $rev = shift;
	
	my @results = ();
	my %loads = ();
	my $table = TBL_LOADS_ULS;
	DB_select( \@results, "*", $table, "" );
	
	if( @results == 0 ){
		my $msg = "Result set empty. Requested value(s) not in table [$table] for revision [$rev].";
		Debug::fatal( $msg, (caller(0))[3], __LINE__ );
	}
	
	foreach my $d (@results){
		#remove last four items
		my @res = @{$d};
		pop @res;
		pop @res;
		pop @res;
		pop @res;
		
		my( $loc, $elev, $load_rev, @loads ) = @res;
		
		#TODO:
		unless( $rev eq $load_rev ){ next; }
		
		
		if( $loads{$loc} ){
			my %loads_loc = %{$loads{$loc}};
			$loads_loc{$elev} = \@loads;
			$loads{$loc} = \%loads_loc;
		}
		else{
			my %loads_loc = ();
			$loads_loc{$elev} = \@loads;
			$loads{$loc} = \%loads_loc;
		}
		
	}
	
	return %loads;
}



# ----------------------------------------------------------------------------
# my( $E_UTM32, $N_UTM32 ) = DBOps::DB_UTM_position( $INPUT_LOC );
# ----------------------------------------------------------------------------

sub DB_UTM_position
{
	Debug::check_num_params( 1, \@_, (caller(0))[3], __LINE__ );
	my $loc = shift;
	
	my @result = ();
	#TODO:
	my $columns = "E_UTM31,N_UTM31";
	DB_select( \@result, $columns, $db_name.".".TBL_WTG_LOCS, "id = \"$loc\"" );
	#TODO: put into Debug-module?
	unless( @result ){
		my $err = "ERROR in [" . (caller(0))[3] . "]:\nResult set empty. Requested value(s) ". $columns . " not in database.";
		die $err; 
	}
	
	return @{$result[0]};
}


# ----------------------------------------------------------------------------
# my @locs = DBOps::DB_locations( 1 );
# ----------------------------------------------------------------------------

sub DB_locations
{
	Debug::check_num_params( 0, \@_, (caller(0))[3], __LINE__ );
	
	my @result = ();
	#TODO:
	my $columns = "id";
	DB_select( \@result, $columns, $db_name.".".TBL_WTG_LOCS, "" );
	#TODO: put into Debug-module?
	unless( @result ){
		my $err = "ERROR in [" . (caller(0))[3] . "]:\nResult set empty. Requested value(s) ". $columns . " not in database.";
		die $err; 
	}
	
	my @ids = ();
	foreach my $res ( @result ){
		my( $id, @tmp ) = @{$res};
		push( @ids, $id );
	}
	
	return @ids;
}


1;
# end of file DBOps.pm

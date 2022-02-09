#! /usr/bin/perl -w
use strict;
##############################################################################
#
# File   :  DBConnect.pm
# History:  2013-11-03 (Martin Kelm) first implementation
#           2016-04-30 (Martin Kelm) update subroutine updated
#           2017-07-21 (Martin Kelm) cleaned
#           2017-08-01 (Martin Kelm) closing function renamed to end
#           2019-07-29 (Martin Kelm) replace added
#
##############################################################################
#
#  Provides basic functions to connect to a MySQL database by using a 
#  object-oriented approach.
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

DBConnect - connect to a MySQL database

=head1 SYNOPSIS

 use DBConnect;

 $dbcon   = DBConnect->new( $dbname, $access_file, $access_group );
 $err_msg = $dbcon->error();
 $err     = $dbcon->insert( $ref_data, $table );
 $err     = $dbcon->delete( $table, $where_condition );
 $err     = $dbcon->open();
 

=head1 DESCRIPTION

The B<DBConnect> module defines methods that make it easier to connect to a 
MySQL database. It uses an access file to get the user account data and server
name.

=over 4

=item $dbcon = DBConnect->new( $dbname, $access_file, $access_group );

Creates a new DBConnect object and gets and checks all information needed 
to open connection to database I<dbname>. The user account data (user and 
password) and the server (host) is extracted from I<access_file> using 
I<access_group>.

=item $err_msg = $dbcon->error();

Returns error message. If a method encounters an error, the error message
is saved and can be retrieved by this method.

=item $err = $dbcon->insert( $ref_data, $table );

Inserts data into I<table>. Returns I<err = 0>.
In case of an invalid argument, the function returns I<err = 1>.

=item $err = $dbcon->delete( $table, $where_condition );

Deletes items in I<table> which fulfill the given I<where_condition>.

=item $err = $dbcon->open();

Opens connection to database as defined in constructor.



=back

=head1 EXAMPLES

 use DBConnect;
 
 my $access_file = "$ENV{HOME}/.my.cnf";
 my $dbcon = DBConnect->new( "stocks", $access_file, "addresses" );
 my $err = $dbcon->open();
 if($err){
 	my $msg = $dbcon->error();
 	print "$msg\n";
 	die;
 }

 my @results = ();
 my $id=1;
 my $isodate="2013-05-28";
 $err = $dbcon->select(  \@results, "*", "quotes", "id=$id and date=\"$isodate\"" );
 if($err){
 	my $msg = $dbcon->error();
 	print "$msg\n";
 	die;
 }
 
 foreach my $d (@results){
 	my @row = @{$d};
 	print "@row\n";
 }

=head1 REQUIRES

Perl 5.8, DBI

=head1 AUTHOR

Martin Kelm, martin.kelm@gmx.de

=cut

##############################################################################

package DBConnect;

use DBI;


# ----------------------------------------------------------------------------
# $dbcon = DBConnect->new( $dbname, $access_file, $access_group );
# ----------------------------------------------------------------------------
# constructor
# ----------------------------------------------------------------------------

sub new
{
	my $classname    = shift;
	my $dbname       = shift;
	my $access_file  = shift;
	my $access_group = shift;
	
	# check if access file is present and readable
	open( FH, "<$access_file" ) || die "Cannot open file $access_file";
	my @lines = ();
	while(<FH>){
		chomp;
		if( $_ eq "" ){ next; }
		if( substr( $_, 0, 1 ) eq "#" ){ next; }
		push( @lines, $_ );
	}
	close( FH );
	# parse file and read data for each group
	my %groups = ();
	my $group = "";
	foreach my $line (@lines){
		if( substr($line, 0, 1) eq "[" ){
			my $num = index( $line, "]" );
			$group = substr( $line, 1, $num-1 );
			next;
		}
		if( $group ){
			my %vars;
			if( $groups{$group} ){
				%vars = %{$groups{$group}};
			}
			my( $key, $val ) = split( "=", $line );
			$vars{$key} = $val;
			$groups{$group} = \%vars;
		}
	}
	
	# check whether access group is defined in access file
	my $found = 0;
	foreach my $g (keys %groups){
		if( $g eq $access_group ){ $found = 1; }
	}
	unless($found){ die "Access group: $access_group not found in $access_file"; }
	
	my $self = {};
	bless( $self, $classname );
	
	# check whether user, host and password is defined for access group
	if( $found ){
		my %vars = %{$groups{$access_group}};
		unless( $vars{"user"} ){ die "User not defined in access file $access_file"; }
		unless( $vars{"host"} ){ die "Host not defined in access file $access_file"; }
		unless( $vars{"password"} ){ die "Password not defined in access file $access_file"; }
		$self->{USER} = $vars{"user"};
		$self->{HOST} = $vars{"host"};
		$self->{PASSWORD} = $vars{"password"};
	}
	
	
	$self->{CONNECTED}    = 0;
	$self->{DB_NAME}      = $dbname;
	$self->{ACCESS_FILE}  = $access_file;
	$self->{ACCESS_GROUP} = $access_group;
	
	$self->{ERROR_MSG}    = "No error";
	
	return $self;
}


# ----------------------------------------------------------------------------
# destructor
# ----------------------------------------------------------------------------

sub DESTROY
{
	my $self = shift;
	end();
}


# ----------------------------------------------------------------------------
# $err_msg = $dbcon->error();
# ----------------------------------------------------------------------------

sub error
{
	my $self = shift;
	
	return $self->{ERROR_MSG};
}


# ----------------------------------------------------------------------------
# $err = $dbcon->insert( $ref_data, $table );
# ----------------------------------------------------------------------------

sub insert
{
	my $self = shift;
	
	my $ref_data = shift;
	my $table = shift;
	
	my $str_vals = "";
	foreach my $data (@{$ref_data}){
		$str_vals .= $data;
		$str_vals .= ",";
	}
	$str_vals = substr( $str_vals, 0, length($str_vals)-1 );
		
	my $string = "INSERT INTO " . $table . " VALUES(" . $str_vals . ")";
	
	#print STDERR "$string\n";
	#return 0;
	
	my $dbh = $self->{DBH};
	my $sth = $dbh->prepare("$string");
	unless( $sth ){ 
		$self->{ERROR_MSG} = $DBI::errstr;
		return 1; 
	}
	
	my $numrows = 0;
	unless( $numrows = $sth->execute() ){
		$self->{ERROR_MSG} = $DBI::errstr;
		return 1;
	}
	
	return 0;
}

sub replace
{
	my $self = shift;
	
	my $ref_data = shift;
	my $table = shift;
	
	my $str_vals = "";
	foreach my $data (@{$ref_data}){
		$str_vals .= $data;
		$str_vals .= ",";
	}
	$str_vals = substr( $str_vals, 0, length($str_vals)-1 );
		
	my $string = "REPLACE INTO " . $table . " VALUES(" . $str_vals . ")";
	
	#print STDERR "$string\n";
	#return 0;
	
	my $dbh = $self->{DBH};
	my $sth = $dbh->prepare("$string");
	unless( $sth ){ 
		$self->{ERROR_MSG} = $DBI::errstr;
		return 1; 
	}
	
	my $numrows = 0;
	unless( $numrows = $sth->execute() ){
		$self->{ERROR_MSG} = $DBI::errstr;
		return 1;
	}
	
	return 0;
}


# ----------------------------------------------------------------------------
# $err = $dbcon->delete( $table, $where_condition );
# ----------------------------------------------------------------------------

sub delete
{
	my $self = shift;
	
	my $table = shift;
	my $where = shift;
	
	
	my $string = "DELETE FROM " . $table;
	if( $where ){
		$string .= " WHERE " . $where;
	}
	
	#print "$string\n";
	
	my $dbh = $self->{DBH};
	my $sth = $dbh->prepare("$string");
	unless( $sth ){ 
		$self->{ERROR_MSG} = $DBI::errstr;
		return -1; 
	}
	
	unless( $sth->execute() ){
		$self->{ERROR_MSG} = $DBI::errstr;
		return -1;
	}
	
	my $numrows = 0;
	unless( $numrows = $sth->execute() ){
		$self->{ERROR_MSG} = $DBI::errstr;
		return 1;
	}
	# TODO: check and define better return code
	return $numrows;
}


# ----------------------------------------------------------------------------
# $err = $dbcon->open();
# ----------------------------------------------------------------------------

sub open
{
	my $self = shift;
	
	if( $self->{CONNECTED} ){ return; }
	
	my $db_name      = $self->{DB_NAME};
	my $access_file  = $self->{ACCESS_FILE};
	my $access_group = $self->{ACCESS_GROUP};
	
	my $user         = $self->{USER};
	my $host         = $self->{HOST};
	my $password     = $self->{PASSWORD};
	
	# under windows it is not possible to read a config file
	my $datasource = "DBI:mysql:database=$db_name;host=$host";
	my $dbh = DBI->connect($datasource, $user, $password, {'PrintError' => 0});
	
	#TODO: REMOVE
	#my $datasource = "DBI:mysql:database=$db_name;".
	#"mysql_read_default_file=$access_file;".
	#"mysql_read_default_group=$access_group";
	#my $dbh = DBI->connect($datasource, undef, undef, {'PrintError' => 0});
	
	if( $dbh ){
		$self->{DBH} = $dbh;
		$self->{CONNECTED} = 1;
		#TODO: get tables and column names for checking
	}
	else{
		$self->{ERROR_MSG} = $DBI::errstr;
		return 1;
	}
	
	return 0;
}


# ----------------------------------------------------------------------------
# $err = $dbcon->select(  );
# ----------------------------------------------------------------------------

sub select
{
	my $self  = shift;
	
	my $ref_res = shift;
	
	my $what  = shift;
	my $table = shift;
	my $where = shift;
	
	
	my $string = "SELECT " . $what . " FROM " . $table;
	if( $where ){
		$string .= " WHERE " . $where;
	}
	
	#TODO: Debug module
	#print STDERR "$string\n";
	
	my $dbh = $self->{DBH};
	my $sth = $dbh->prepare("$string");
	unless( $sth ){ 
		$self->{ERROR_MSG} = $DBI::errstr;
		return 1; 
	}
	
	unless( $sth->execute() ){
		$self->{ERROR_MSG} = $DBI::errstr;
		return 1;
	}
		
	while (my @row = $sth->fetchrow_array())
	{
		push (@{$ref_res},[@row]);
	}
	
	return 0;
}

sub columns
{
	my $self = shift;
	
	my $ref_res = shift;
	
	my $table = shift;
	
	my $string = "SHOW COLUMNS FROM " . $table;
}

sub update
{
	my $self = shift;
	
	my $what  = shift;
	my $table = shift;
	my $where = shift;
	
	my $data  = shift;
	
	my $string = "UPDATE " . $table . " SET " . $what . "=?" . " WHERE " . $where;
	#print STDERR "$string\n";
	
	my $dbh = $self->{DBH};
	my $sth = $dbh->prepare("$string");
	unless( $sth ){ 
		$self->{ERROR_MSG} = $DBI::errstr;
		return 1; 
	}
	
	my $OK = $sth->execute($data);
	#print STDERR "*** return: $OK\n";
	unless( $OK == 1 ){
		$self->{ERROR_MSG} = "Update failed: $string";
		return 1;
	}
	
	return 0;
}

sub end
{
	my $self = shift;
	
	#print "closing\n";
	
	if( $self->{CONNECTED} ){
		$self->{DBH}->disconnect();
	}
}


sub handle
{
	my $self = shift;
	return $self->{DBH};
}


1;
# end of file DBConnect.pm

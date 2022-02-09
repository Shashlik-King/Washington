#! /usr/bin/perl -w
use strict;
##############################################################################
#
# File   :  Debug.pm
# History:  2014-12-04 (Martin Kelm) first implementation
#           2015-06-20 (Martin Kelm) update documentation
#           2015-11-16 (Martin Kelm) check_num[_opt]_params with return value
#
##############################################################################
#
#  Provides basic subroutines to support debugging
#
##############################################################################

=head1 NAME

Debug - support subroutines for debugging

=head1 SYNOPSIS

 use Debug;

 set_debug_level( $level );
 info( $msg, $sub );
 debug( $msg, $sub );
 warn( $msg, $sub );
 error( $msg, $sub );
 fatal( $msg, $sub );
 check_num_params( $req_nparams, $ref_params, $sub );
 check_num_opt_params( $ref_req_nparams, $ref_params, $sub, $line );

=head1 DESCRIPTION

The B<Debug> module defines subroutines that make it easier to debug a script.

=over 4

=item set_debug_level( $level );

Sets debug level to I<level>. The I<level> defines whether a certain message
is printed to STDERR, see other subroutines. A high debug I<level> means 
that more is printed.

=item info( $msg, $sub );

Prints info message I<msg>. I<sub> is the subroutine out of which the message
shall be printed.

=item debug( $msg, $sub );



=item warn( $msg, $sub );



=item error( $msg, $sub );


=item fatal( $msg, $sub );


=back

=head1 EXAMPLES

=head1 REQUIRES

Perl 5.8

=head1 AUTHOR

Martin Kelm, martin.kelm@gmx.de

=cut

package Debug;

require Exporter;
use vars qw( $VERSION );

$VERSION = 0.2;


use constant DEBUG_LEVEL_ZERO         => 0;
use constant DEBUG_LEVEL_LOW          => 1;
use constant DEBUG_LEVEL_INTERMEDIATE => 2;
use constant DEBUG_LEVEL_HIGH         => 3;
use constant DEBUG_LEVEL_ALL          => 4;


our $gl_debug_level = DEBUG_LEVEL_ZERO;
our $gl_output      = "STDERR";


# ----------------------------------------------------------------------------
# set_debug_level( $level );
# ----------------------------------------------------------------------------

sub set_debug_level
{
	my $level = shift;
	
	$gl_debug_level = $level;
}


# ----------------------------------------------------------------------------
# set_output_file( $name );
# ----------------------------------------------------------------------------

sub set_output_file
{
	my $file = shift;
	
	$gl_output = $file;
	open( OUT, ">$file" );
}


# ----------------------------------------------------------------------------
# info( $msg, $sub );
# ----------------------------------------------------------------------------

sub info
{
	my $msg = shift;
	my $sub = shift;
	
	if( $gl_debug_level >= DEBUG_LEVEL_HIGH  ){
		#print STDERR "INFO in [$sub]: $msg\n";
		print_msg( "INFO in [$sub]:\n$msg\n" );
	}
}


# ----------------------------------------------------------------------------
# debug( $msg, $sub, $line );
# ----------------------------------------------------------------------------

sub debug
{
	my $msg  = shift;
	my $sub  = shift;
	my $line = shift;
	
	if( $gl_debug_level >= DEBUG_LEVEL_ALL  ){
		print_msg( "DEBUG in subroutine [$sub], line [$line]:\n$msg\n" );
	}
}


# ----------------------------------------------------------------------------
# warn( $msg, $sub, $line );
# ----------------------------------------------------------------------------

sub warn
{
	my $msg  = shift;
	my $sub  = shift;
	my $line = shift;
	
	if( $gl_debug_level >= DEBUG_LEVEL_INTERMEDIATE  ){
		print_msg( "WARNING in [$sub], line [$line]:\n$msg\n" );
	}
}


# ----------------------------------------------------------------------------
# error( $msg, $sub, $line );
# ----------------------------------------------------------------------------

sub error
{
	my $msg  = shift;
	my $sub  = shift;
	my $line = shift;
	
	if( $gl_debug_level >= DEBUG_LEVEL_LOW  ){
		print_msg( "ERROR in subroutine [$sub], line [$line]:\n$msg\n" );
	}
}


# ----------------------------------------------------------------------------
# fatal( $msg, $sub, $line );
# ----------------------------------------------------------------------------

sub fatal
{
	my $msg  = shift;
	my $sub  = shift;
	my $line = shift;
	
	if( $gl_debug_level >= DEBUG_LEVEL_ZERO  ){
		print_msg( "FATAL in subroutine [$sub], line [$line]:\n$msg\n" );
		die;
	}
}


# ----------------------------------------------------------------------------
# check_num_params( $req_nparams, $ref_params, $sub, $line );
# ----------------------------------------------------------------------------
# Checks whether the given reference to array of parameters contains the 
# required number of parameters $req_nparams.
# $sub shall be the name of the checked subroutine.
# ----------------------------------------------------------------------------

sub check_num_params
{
	my $req_nparams = shift;
	my $ref_params  = shift;
	my $sub         = shift;
	my $line        = shift;
	
	unless( @{$ref_params} == $req_nparams ){
		my $prov_nparams = @{$ref_params};
		my $msg = "Wrong number of parameters [$prov_nparams]. Must be [$req_nparams].";
		fatal( $msg, $sub, $line );
	}
	
	return $req_nparams;
}


# ----------------------------------------------------------------------------
# check_num_opt_params( $ref_req_nparams, $ref_params, $sub, $line );
# ----------------------------------------------------------------------------
# Checks whether the given reference to array of parameters contains the 
# required number of parameters $req_nparams.
# $sub shall be the name of the checked subroutine.
# ----------------------------------------------------------------------------

sub check_num_opt_params
{
	my $ref_req_nparams = shift;
	my $ref_params      = shift;
	my $sub             = shift;
	my $line            = shift;
	
	my $prov_nparams = @{$ref_params};
	
	my $found = 0;
	
	foreach my $val (@{$ref_req_nparams}){
		if( $val == $prov_nparams ){ $found = 1; last; }
	}
	
	unless( $found ){
		my @tmp = @{$ref_req_nparams};
		my $msg = "Wrong number of parameters [$prov_nparams]. Must be [@tmp].";
		fatal( $msg, $sub, $line );
	}
	
	return $prov_nparams;
}


# ----------------------------------------------------------------------------
# check_limits( $value, $lower_limit, $upper_limit, $name );
# ----------------------------------------------------------------------------
# Checks whether the given $value is within the given interval. The given 
# limits are part of the interval, i.e. [$lower_limit, $upper_limit].
# ----------------------------------------------------------------------------

sub check_limits
{
	my $value       = shift;
	my $lower_limit = shift;
	my $upper_limit = shift;
	my $name        = shift;
	my $sub         = shift;
	my $line        = shift;
	
	unless( ($value >= $lower_limit) && ($value <= $upper_limit) ){ 
		fatal( "Given $name [$value] out of range [$lower_limit .. $upper_limit].", $sub, $line );
	}
}


# ----------------------------------------------------------------------------
# check_num_positive( $value, $name );
# ----------------------------------------------------------------------------

sub check_num_positive
{
	my $value = shift;
	my $name  = shift;
	my $sub   = shift;
	my $line  = shift;
	
	unless( is_numeric( $value ) ){
		fatal( "Given $name [$value] is not a number.", $sub, $line );
	}
	
	unless( $value > 0.0 ){
		fatal( "Given $name [$value] is not a positive number.", $sub, $line );
	}
}


# ----------------------------------------------------------------------------
# check_string_in_array( $value, $ref_array, $name, $sub, $line );
# ----------------------------------------------------------------------------

sub check_string_in_array
{
	my $value     = shift;
	my $ref_array = shift; 
	my $name      = shift;
	my $sub       = shift;
	my $line      = shift;
	
	my $found = 0;
	
	foreach my $val (@{$ref_array}){
		if( $val eq $value ){ $found = 1; last; }
	}
	
	unless( $found ){
		my @tmp = @{$ref_array};
		fatal( "Given $name [$value] is not in allowed array of [@tmp].", $sub, $line );
	}
	
	return $found;
}

sub check_value_in_array
{
	my $value     = shift;
	my $ref_array = shift; 
	my $name      = shift;
	my $sub       = shift;
	my $line      = shift;
	
	my $found = 0;
	
	foreach my $val (@{$ref_array}){
		if( $val == $value ){ $found = 1; last; }
	}
	
	unless( $found ){
		fatal("Given $name [$value] is not in given array of values.", $sub, $line );
	}
	
	return $found;
}

# test if string is numeric
# $res = is_numeric($string)
sub is_numeric { defined getnum($_[0]) }

# convert string to numeric
sub getnum
{
	use POSIX qw(strtod);
	use POSIX qw(locale_h);

	# --> use this if set mistakenly to wrong setting
	# setlocale(LC_NUMERIC, "de_DE.UTF-8"); 
	
	# query and save the old locale
	my $old_locale = setlocale(LC_NUMERIC);
	
	# set new locale
	my $nlocale = setlocale(LC_NUMERIC, "C");

	
	my $str = shift;
	if(defined $str){
		$str =~ s/^\s+//;
		$str =~ s/\s+$//;
		$! = 0;
		my( $num, $unparsed ) = strtod($str);

		if( ($str eq '') || ($unparsed != 0) || $! ){
			return undef;
		}
		else{
#			print "getnum: $num $str $unparsed\n";		
			return $num;
		}			
	}
	else{
		return undef;
	}
	
	# return to old locale
	setlocale(LC_NUMERIC, $old_locale);
}

sub print_msg
{
	my $msg = shift;
	
	unless( $gl_output eq "STDERR" ){ print OUT $msg; }
	else{ print STDERR $msg; }
	
}


1;
# end of file Debug.pm

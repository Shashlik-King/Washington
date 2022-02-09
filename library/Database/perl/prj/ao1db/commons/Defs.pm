#!/usr/bin/perl -w
use strict;
##############################################################################
#
# File   :  Defs.pm
# History:  2013-11-17 (Martin Kelm) first implementation
#           2014-10-20 (Martin Kelm) documentation updated 
#
##############################################################################
#
#  Includes definitions used in various perl scripts 
#
##############################################################################

=head1 NAME

Defs - includes definitions used in various perl scripts

=head1 SYNOPSIS

 use Defs;

=head1 DESCRIPTION

The B<Defs> module provides definitions of constants etc.

=over 4

=cut

package Defs;


# ----------------------------------------------------------------------------
# constants which change with every project
# ----------------------------------------------------------------------------

use constant DB 				=> "ao1db";
use constant APPROVER			=> "SIA";


# ----------------------------------------------------------------------------
# constants which do not change with every projects
# ----------------------------------------------------------------------------





1;
# end of file Defs.pm

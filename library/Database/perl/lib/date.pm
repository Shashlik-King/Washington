#!/usr/bin/perl
# -*- cperl -*-
# $Id: date.pm,v 1.1 2001/08/31 13:26:23 uli Exp $

# Erstellt am 18.01.01 von  (<Ulrich.Herbst@gmx.de>)
#
# Letzte Aenderung: <date.pm: 19.09.02, 19:25 (uli)>
# 
# Editiert zur Ausgabe von Monatsnamen durch C. Ross
# Letzte Aenderung: <date.pm: 08.03.10, 10:11 (Ros)>

package        date;
require        Exporter;
@ISA        = qw(Exporter);
@EXPORT     = qw(yesterday now);    # Diese Funktionen werden
									# defaultmaessig exportiert
#@EXPORT_OK = qw()                  # Diese Funktionen werden
                                    # bei Bedarf exportiert

use strict;
use English;
my $Version='$Id: date.pm,v 1.1 2001/08/31 13:26:23 uli Exp $$';

#----------------------------------------------------------------------

=head1 NAME

date.pm - Datumsangaben  in verschiedenen Formaten ausgeben

=head1 VERSION

$Id: date.pm,v 1.1 2001/08/31 13:26:23 sai1731 Exp $

=head1 SYNTAX

  use date;
  my $DATEFORMAT="%Y%m%d";
  my $HEUTE=now($DATEFORMAT);
  my $GESTERN=yesterday($DATEFORMAT);

=head1 BESCHREIBUNG

$DATEFORMAT hat folgenden Aufbau:

=over 4

=item %y

Jahr 2-stellig

=item %Y

Jahr 4-stellig

=item %m

Monat 2-stellig

=item %B

Monat ausgeschrieben

=item %d

Tag 2-stellig

=item %H

Stunden 2-stellig (24h)

=item %M

Minuten 2-stellig

=item %S

Sekunden 2-stellig

=back 

alles andere wird buchstaeblich wiedergegeben:

-> %Y-%m-%d -> 2001-07-09

now() liefert das aktuelle Datum bzw. Uhrzeit, yesterday() das Datum/Uhrzeit
vor 24h.

=head1 AUTHORS

Ulrich Herbst <ulrich.herbst@gmx.de>

=cut

#----------------------------------------------------------------------

sub yesterday {

  use Time::Local;
  my $FORMAT=$_[0];
  my @month_names = ("January","February","March","April","May","June","July","August","September","October","November","December");

  my $Gestern=timelocal(localtime)-24*60*60;
  # Das ist jetzt nicht Y2K-Sicher....
  my $y=sprintf("%02d",(localtime($Gestern))[5]-100);
  my $Y=sprintf("%04d",(localtime($Gestern))[5]+1900);
  my $m=sprintf("%02d",(localtime($Gestern))[4]+1);
  my $B=sprintf(@month_names[(localtime($Gestern))[4]]);
  my $d=sprintf("%02d",(localtime($Gestern))[3]);
  my $H=sprintf("%02d",(localtime($Gestern))[2]);
  my $M=sprintf("%02d",(localtime($Gestern))[1]);
  my $S=sprintf("%02d",(localtime($Gestern))[0]);

  $FORMAT =~ s/%y/$y/;
  $FORMAT =~ s/%Y/$Y/;
  $FORMAT =~ s/%m/$m/;
  $FORMAT =~ s/%B/$B/;
  $FORMAT =~ s/%d/$d/;
  $FORMAT =~ s/%H/$H/;
  $FORMAT =~ s/%M/$M/;
  $FORMAT =~ s/%S/$S/;

  return $FORMAT;
}


# ----------------------------------------------------------------------------
# my $today = now("%Y-%m-%d");
# ----------------------------------------------------------------------------

sub now {

  my $FORMAT=$_[0];
  my @month_names = ("January","February","March","April","May","June","July","August","September","October","November","December");

  my $NOW=timelocal(localtime);
  # Das ist jetzt nicht Y2K-Sicher....
  my $y=sprintf("%02d",(localtime($NOW))[5]-100);
  my $Y=sprintf("%04d",(localtime($NOW))[5]+1900);
  my $m=sprintf("%02d",(localtime($NOW))[4]+1);
  my $B=sprintf(@month_names[(localtime($NOW))[4]]);
  my $d=sprintf("%02d",(localtime($NOW))[3]);
  my $H=sprintf("%02d",(localtime($NOW))[2]);
  my $M=sprintf("%02d",(localtime($NOW))[1]);
  my $S=sprintf("%02d",(localtime($NOW))[0]);

  $FORMAT =~ s/%y/$y/;
  $FORMAT =~ s/%Y/$Y/;
  $FORMAT =~ s/%m/$m/;
  $FORMAT =~ s/%B/$B/;
  $FORMAT =~ s/%d/$d/;
  $FORMAT =~ s/%H/$H/;
  $FORMAT =~ s/%M/$M/;
  $FORMAT =~ s/%S/$S/;

  return $FORMAT;
}

1;

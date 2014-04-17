#!/usr/bin/perl

package Connect;

use assembler::Alignment;
use strict;

our @ISA = qw(Alignment);

#######################################################################################################################
sub new {
	my $class = shift;

	(@_ == 11) || die "\nError: Connect::new expected 11 arguments but recieved only " . scalar(@_) . ":\n@_\n\n";
	my ($seq1, $side1, $seq2, $side2, $s1, $e1, $s2, $e2, $len1, $len2, $pidentity) = @_;
	
	my $self = $class->SUPER::new($seq1, $seq2, $s1, $e1, $s2, $e2, $len1, $len2, $pidentity);
	$self->{_side1} = $side1;
	$self->{_side2} = $side2;
	$self->{_type}  = undef;

	if(($self->side1 eq 'Middle') || ($self->side2 eq 'Middle')) {
		# If both regions are in the middle of sequences - middle-middle shared.
		# Else - edge-middle conneection
		$self->{_type} = (($self->side1 eq 'Middle') && ($self->side2 eq 'Middle'))? 'MM_SHARED' : 'EM_CONNECT';
	}
	elsif((($self->side1 == $self->side2) && ($self->direction eq 'forward')) || (($self->side1 != $self->side2) && ($self->direction eq 'reverse'))) {
		$self->{_type} = 'EE_SHARED'		# Edge-edge shared
	}
	else {
		$self->{_type} = 'EE_CONNECT';		# Edge-edge connection
	}
 
	bless $self, $class;
	return $self;
}

#######################################################################################################################
sub type {
        my ( $self ) = @_;

	return $self->{_type};
}

#######################################################################################################################
sub other_end {
        my ( $self, $seq, $side) = @_;

	if(($seq eq $self->seq1) && ($side eq $self->side1)) {
		return ($self->seq2, $self->side2);
	}
	if(($seq eq $self->seq2) && ($side eq $self->side2)) {
		return ($self->seq1, $self->side1);
	}

	die "\nConnect::other_end: end parameters ($seq:$side) does not match neither end1 (" . $self->seq1 . ":" . $self->side1 . 
	    ") nor end2 (" . $self->seq2 . ":" . $self->side2 . ")\n\n"; 
}

#######################################################################################################################
sub side1 {
        my ( $self ) = @_;
	return $self->{_side1};
}

#######################################################################################################################
sub side2 {
        my ( $self ) = @_;
	return $self->{_side2};
}

#######################################################################################################################
sub edge1 {
        my ( $self ) = @_;
	return $self->seq1 . ":" . $self->side1;
}

#######################################################################################################################
sub edge2 {
        my ( $self ) = @_;
	return $self->seq2 . ":" . $self->side2;
}

1;

#!/usr/bin/perl

use strict;
use Bio::SeqIO;

package ScafLink;

#######################################################################################################################
sub new {
	my $class = shift;
	my $self = {
		_seq_obj     => shift,
		_start       => shift,
		_end         => shift,
		_orientation => shift,
		_pos_on_scaf => undef,
	};

	bless $self, $class;
	return $self;
}

#######################################################################################################################
sub copy {
	my ( $self, $link ) = @_;

	$self->{_seq_obj} = $link->seq_obj;
	$self->{_start} = $link->start;
	$self->{_end} = $link->end;
	$self->{_orientation} = $link->orientation;
	$self->{_pos_on_scaf} = $link->get_pos_on_scaffold;
}

#######################################################################################################################
sub clone {
	my ( $self ) = @_;
	my $c = new ScafLink($self->{_seq_obj}, $self->{_start}, $self->{_end}, $self->{_orientation});

	$c->set_pos_on_scaffold($self->{_pos_on_scaf});

	return $c;
}

#######################################################################################################################
sub contribute_seq {
	my ( $self ) = @_;

	($self->start >= 1) || die "\nError, ScafLink::contribute_seq: start coordinate, " . $self->start . " is illegal (" . $self->{_seq_obj}->display_id . ", " .
				   $self->{_seq_obj}->length . " bps)\n\n";
	($self->end <= $self->{_seq_obj}->length) || die "\nError, ScafLink::contribute_seq: end coordinate, " . $self->end . " is illegal (" . $self->{_seq_obj}->display_id . ", " .
				   $self->{_seq_obj}->length . " bps)\n\n";
	($self->end > $self->start) || die "\nError, ScafLink::contribute_seq: start and end coordinates (" . $self->start . ", " . $self->end . ") are illegal (" . 
					   $self->{_seq_obj}->display_id . ", " . $self->{_seq_obj}->length . " bps)\n\n";

	my $seq = $self->seq_obj->subseq($self->start, $self->end);

	# Just to be sure
	$seq =~ tr/acgtn/ACGTN/;
	if($self->orientation eq 'reverse') {
		$seq =~ tr/ACGT/TGCA/;
		$seq = reverse($seq);
	}

	return $seq;
}

#######################################################################################################################
sub set_pos_on_scaffold {
	my ( $self, $pos ) = @_;

	$self->{_pos_on_scaf} = $pos;
}

#######################################################################################################################
sub get_pos_on_scaffold {
	my ( $self, $pos ) = @_;

	return $self->{_pos_on_scaf};
}

#######################################################################################################################
sub seq_obj {
	my ( $self ) = @_;

	return $self->{_seq_obj};
}

#######################################################################################################################
sub start {
	my ( $self ) = @_;

	return $self->{_start};
}

#######################################################################################################################
sub end {
	my ( $self ) = @_;

	return $self->{_end};
}

#######################################################################################################################
sub orientation {
	my ( $self ) = @_;

	return $self->{_orientation};
}

1;

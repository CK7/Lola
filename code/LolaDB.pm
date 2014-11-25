#!/usr/bin/perl

package LolaDB;

use strict;
use assembler::AlignmentDB;
use assembler::LolaAssembler;
use Bio::SeqIO;

#######################################################################################################################
sub new {
        my $class = shift;
        (@_ == 3) || die "\nError: LolaDB::new expected 2 arguments but recieved only " . scalar(@_) . "\n\n";
        my $self = {
                _sequences_db => shift,
                _alignment_db => shift,
		_pidentity_threshold => shift,
                _clusters     => {},
                _seq2cluster  => {},
                _unclustered  => {},
        };

        bless $self, $class;

	$self->init_clusters();

        return $self;
}

#######################################################################################################################
sub num_clusters {
	my( $self ) = @_;

	return scalar(keys %{$self->{_clusters}});
}

#######################################################################################################################
sub get_cluster {
	my( $self, $cluster_id ) = @_;

	return $self->{_clusters}{$cluster_id};
}

#######################################################################################################################
sub get_unclustered {
	my( $self ) = @_;

	return $self->{_unclustered};
}

#######################################################################################################################
sub get_seq_cluster {
	my( $self, $seq_id ) = @_;

	return $self->{_seq2cluster}{$seq_id};
}

#######################################################################################################################
sub init_clusters {
	my( $self ) = @_;

        # We cluster together all sequences that are either CONNECTED, CONTAINS or IDENTICAL
	my $next_cluster_id = 1;
        foreach my $seq_id (keys %{$self->{_sequences_db}}) {
                next if(exists($self->{_seq2cluster}{$seq_id}));
		my $next_cluster_sequences = {};
		my @check = ($seq_id);
		$next_cluster_sequences->{$seq_id} = $self->{_sequences_db}{$seq_id};

		while(@check > 0) {
			my $check_seq_id = shift(@check);

			# Sanity check
			!exists($self->{_seq2cluster}{$check_seq_id}) || die "\nError, Lola::init_clusters: $check_seq_id already belongs to a cluster but is re-clustered\n\n"; 

			my $all_connections5 = $self->{_alignment_db}->get_connections($check_seq_id, 5);
			my $all_connections3 = $self->{_alignment_db}->get_connections($check_seq_id, 3);
			my $all_connectionsMiddle = $self->{_alignment_db}->get_connections($check_seq_id, 'Middle');
			my $all_containers = $self->{_alignment_db}->get_containers($check_seq_id);
			my $all_contained = $self->{_alignment_db}->get_contained($check_seq_id);
			my $all_unknowns = $self->{_alignment_db}->get_unknowns($check_seq_id);

			foreach my $db ($all_connections5, $all_connections3, $all_connectionsMiddle, $all_containers, $all_contained, $all_unknowns) {
				foreach my $c (values %{$db}) {
					my $new_seq_id = ($c->seq1 eq $check_seq_id)? $c->seq2 : $c->seq1;
					if(!exists($next_cluster_sequences->{$new_seq_id})) {
						$next_cluster_sequences->{$new_seq_id} = $self->{_sequences_db}{$new_seq_id};
						push(@check, $new_seq_id);
					}
				}
			}
		}

		# If the following is true then this is a single sequence cluster - only the query cluster
		if(scalar(keys %{$next_cluster_sequences}) == 1) {
			$self->{_unclustered}{$seq_id} = $self->{_sequences_db}{$seq_id};
			next;
		}

		my $lola = new LolaAssembler($next_cluster_id, $next_cluster_sequences, $self->{_alignment_db}, $self->{_pidentity_threshold});
		$self->{_clusters}{$next_cluster_id} = $lola;
		foreach my $seq_id (keys %{$next_cluster_sequences}) {

			$self->{_seq2cluster}{$seq_id} = $lola;
		}
		$next_cluster_id++;
        }
}

1;

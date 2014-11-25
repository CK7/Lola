#!/usr/bin/perl

use strict;
use assembler::ScafLink;
use assembler::Connect;
use assembler::AlignmentDB;

package Scaffold;

#######################################################################################################################
sub new {
	my $class = shift;
	my $self = {
		_scaf_links  => [],
                _participating_seqs => {},
                _associated_seqs    => {},
                _inconsistent_seqs  => {},
		_size => 0,
        };

        bless $self, $class;
        return $self;
}

#######################################################################################################################
sub clone {
	my ( $self ) = @_;
	my $c = new Scaffold();

	$c->copy($self);
	return $c;
}

#######################################################################################################################
sub copy {
	my ( $self, $c ) = @_;

	$self->{_scaf_links} = [];
	$self->{_associated_seqs} = {};
	$self->{_participating_seqs} = {};
	$self->{_size} = $c->length();

	foreach my $i (1 .. $c->num_links) {
		my $link = $c->get_link($i)->clone();
		push(@{$self->{_scaf_links}}, $link);
		$self->{_participating_seqs}{$link->seq_obj->display_id} = $link;
	}

	foreach my $clink ($c->get_associated) {
		my $link = $clink->clone();
		$self->{_associated_seqs}{$link->seq_obj->display_id} = $link;
	}

	foreach my $clink ($c->get_inconsistent) {
		my $link = $clink->clone();
		$self->{_inconsistent_seqs}{$link->seq_obj->display_id} = $link;
	}
}

#######################################################################################################################
sub seq {
	my ( $self ) = @_;
	my $seq = '';

	foreach my $i (1 .. $self->num_links) {
		$seq .= $self->get_link($i)->contribute_seq();
	}

	return $seq;
}

#######################################################################################################################
sub get_associated {
	my ( $self, $seq_id ) = @_;

	return defined($seq_id)? $self->{_associated_seqs}{$seq_id} : values(%{$self->{_associated_seqs}});
}

#######################################################################################################################
sub get_inconsistent {
	my ( $self, $seq_id ) = @_;

	return defined($seq_id)? $self->{_inconsistent_seqs}{$seq_id} : values(%{$self->{_inconsistent_seqs}});
}


#######################################################################################################################
sub clear {
	my ( $self ) = @_;

	$self->{_scaf_links} = [];
	$self->{_participating_seqs} = {};
	$self->{_associated_seqs} = {};
	$self->{_inconsistent_seqs} = {};
	$self->{_size} = 0;

	return $self;
}

#######################################################################################################################
sub length {
	my ( $self ) = @_;

	return $self->{_size};
}

#######################################################################################################################
sub push_link {
	my ( $self, $link ) = @_;

	$self->{_size} += ($link->end - $link->start + 1);
	my $pos = $self->{_size} - $link->seq_obj->length + 1;

	$link->set_pos_on_scaffold($pos);
	$self->{_size} = $pos+$link->seq_obj->length-1;
	push(@{$self->{_scaf_links}}, $link);
	$self->{_participating_seqs}{$link->seq_obj->display_id} = $link;
}

#######################################################################################################################
sub pop_link {
	my ( $self ) = @_;

	(@{$self->{_scaf_links}} > 0) || return undef;

	my $link = pop(@{$self->{_scaf_links}});
	$self->{_size} -= ($link->end - $link->start + 1);
	$link->set_pos_on_scaffold(undef);
	delete($self->{_participating_seqs}{$link->seq_obj->display_id});

	return $link;
}

#######################################################################################################################
sub get_link {
	my ( $self, $i ) = @_;

	return $self->{_scaf_links}[$i-1];
}

#######################################################################################################################
sub num_links {
	my ( $self, $link ) = @_;

	return scalar(@{$self->{_scaf_links}});
}

#######################################################################################################################
sub num_associated {
	my ( $self, $link ) = @_;

	return scalar(keys %{$self->{_associated_seqs}});
}

#######################################################################################################################
sub num_inconsistent {
	my ( $self, $link ) = @_;

	return scalar(keys %{$self->{_inconsistent_seqs}});
}

#######################################################################################################################
sub add_associate {
	my ( $self, $seq_obj, $alignment_db, $OUT_ref ) = @_;

	# The following cases are not errors - it is ok to try and check if something is associated by simply adding it
	(!exists($self->{_participating_seqs}{$seq_obj->display_id}) && !exists($self->{_associated_seqs}{$seq_obj->display_id})) || return 1;

	# I don't expect this to happen and would like to raise arrention if it does
	!exists($self->{_inconsistent_seqs}{$seq_obj->display_id}) || die "\nScaffold::add_associate: trying to add " . $seq_obj->display_id .
		" that is already known to be inconsistent\n\n";

	# Ok - see comment above
	$self->is_consistent($seq_obj, $alignment_db, $OUT_ref) || return 0;

	my ($pos_on_scaf, $orientation) = $self->find_position($seq_obj, $alignment_db);
	my $link = new ScafLink($seq_obj, 1, $seq_obj->length, $orientation, $pos_on_scaf);
	$link->set_pos_on_scaffold($pos_on_scaf);

	$self->{_associated_seqs}{$seq_obj->display_id} = $link;

	return 1;
}

#######################################################################################################################
sub add_inconsistent {
	my ( $self, $seq_obj, $alignment_db, $OUT_ref ) = @_;

	# Not an error... probably
	return 1 if(exists($self->{_inconsistent_seqs}{$seq_obj->display_id}));

	# Sanity check
	if(exists($self->{_participating_seqs}{$seq_obj->display_id}) || exists($self->{_associated_seqs}{$seq_obj->display_id})) {
		die "Scaffold::add_inconsistent: trying to add " . $seq_obj->display_id . " but this sequence already appears to be" .
		    " either part of or associated with the scaffold\n\n";
	}

	if($self->is_consistent($seq_obj, $alignment_db, $OUT_ref)) {
		die "Scaffold::add_inconsistent: trying to add " . $seq_obj->display_id . " but this sequence appears to be consistent with the scaffold\n\n";
	}

	my ($pos_on_scaf, $orientation) = $self->find_position($seq_obj, $alignment_db);

	if(!defined($pos_on_scaf)) {
		die "Scaffold::add_inconsistent: trying to add " . $seq_obj->display_id . " but could not find its estimated position on the scaffold\n\n";
	}

	my $link = new ScafLink($seq_obj, 1, $seq_obj->length, $orientation, $pos_on_scaf);
	$link->set_pos_on_scaffold($pos_on_scaf);

	$self->{_inconsistent_seqs}{$seq_obj->display_id} = $link;
	return 1;
}

#######################################################################################################################
sub find_position {
	my ( $self, $seq_obj, $alignment_db ) = @_;

	# The trivial cases
	if(exists($self->{_participating_seqs}{$seq_obj->display_id})) {
		return ($self->{_participating_seqs}{$seq_obj->display_id}->get_pos_on_scaffold, $self->{_participating_seqs}{$seq_obj->display_id}->orientation);
	}
	elsif(exists($self->{_associated_seqs}{$seq_obj->display_id})) {
		return ($self->{_associated_seqs}{$seq_obj->display_id}->get_pos_on_scaffold, $self->{_associated_seqs}{$seq_obj->display_id}->orientation);
	}

        # If this one is contained in a consistent sequence - calculate the coordinate accordingly
        my $containers = $alignment_db->get_containers($seq_obj->display_id);
        foreach my $c (values %{$containers}) {
		my $cseq = ($seq_obj->display_id eq $c->seq1)? $c->seq2 : $c->seq1;
                if(exists($self->{_participating_seqs}{$cseq}) || exists($self->{_associated_seqs}{$cseq})) {
			my $overlapping_seq = exists($self->{_participating_seqs}{$cseq})? $self->{_participating_seqs}{$cseq} : $self->{_associated_seqs}{$cseq};
			my $pos = $overlapping_seq->get_pos_on_scaffold;
			$pos += ($overlapping_seq->orientation eq 'forward')? ($c->start($cseq)-1) : ($overlapping_seq->seq_obj->length - $c->end($cseq));
			my $orientation = ($overlapping_seq->orientation eq $c->direction)? 'forward' : 'reverse';
			return ($pos, $orientation);
                }
        }

	# Else - look for overlapping sequences
	foreach my $side(5, 3) {
		my $connections = $alignment_db->get_connections($seq_obj->display_id, $side, 'EE_CONNECT');
	        foreach my $c (values %{$connections}) {
			my ($seq2, $side2) = $c->other_end($seq_obj->display_id, $side);

			next if (!exists($self->{_participating_seqs}{$seq2}) && !exists($self->{_associated_seqs}{$seq2}));
			my $overlapping_seq = exists($self->{_participating_seqs}{$seq2})? $self->{_participating_seqs}{$seq2} : $self->{_associated_seqs}{$seq2};
			my $pos = $overlapping_seq->get_pos_on_scaffold;
			if($overlapping_seq->orientation eq 'forward') {
				$pos += ($side2 == 3)? ($c->start($seq2)-1) : ($c->end($seq2)-$seq_obj->length);
			}
			else {
				$pos += ($side2 == 5)? ($c->length($seq2)-$c->end($seq2)) : ($c->length($seq2)-$c->start($seq2)-$seq_obj->length+1);
			}
			my $orientation = ($overlapping_seq->orientation eq $c->direction)? 'forward' : 'reverse';
			return ($pos, $orientation);
                }
        }

	return (undef, undef);
}

#######################################################################################################################
sub is_consistent {
	my ( $self, $seq_obj, $alignment_db, $OUT_ref ) = @_;

        # Nothing to do in this case - sequence is already known to be consistent
        if(exists($self->{_participating_seqs}{$seq_obj->display_id}) || exists($self->{_associated_seqs}{$seq_obj->display_id})) {
                print $OUT_ref $seq_obj->display_id, " is already recognized as consistent\n";
                return 1;
        }

        # At the current time we do not deal with sequences that are connected to the middle of something
#        my $connections5 = $alignment_db->get_connections($seq_obj->display_id, 5, 'EM_CONNECT');
#        my $connections3 = $alignment_db->get_connections($seq_obj->display_id, 3, 'EM_CONNECT');
#        if(((keys %{$connections5}) != 0) || ((keys %{$connections3}) != 0)) {
#                print $OUT_ref $seq_obj->display_id, " is inconsistent - has a middle connection\n";
#                return 0;
#        }

        # And also share connections with other sequences
#        my $connections5 = $alignment_db->get_connections($seq_obj->display_id, 5, 'EE_SHARED');
#        my $connections3 = $alignment_db->get_connections($seq_obj->display_id, 3, 'EE_SHARED');
#        if(((keys %{$connections5}) != 0) || ((keys %{$connections3}) != 0)) {
#                print $OUT_ref $seq_obj->display_id, " is inconsistent - has a shared connection\n";
#                return 0;
#        }

        # If this one is contained in a consistent sequence - we are done
        my $containers = $alignment_db->get_containers($seq_obj->display_id);
        foreach my $container (keys %{$containers}) {
                if(exists($self->{_participating_seqs}{$container}) || exists($self->{_associated_seqs}{$container})) {
                        print $OUT_ref $seq_obj->display_id, " is consistent - contained within $container which is consistent\n";
                        return 1;
                }
        }

        # Try to find consistent sequences on both sides and middle. Go through all sequences that are either connected to
	# or contained by this sequence
        my $orientation = undef;
        my $connections_info = ''; 
	my %regions_covered = ();

        my $contained = $alignment_db->get_contained($seq_obj->display_id);
        foreach my $c (values %{$contained}) {
		my $cseq = $c->contained;
                if(exists($self->{_participating_seqs}{$cseq}) || exists($self->{_associated_seqs}{$cseq})) {
			my ($s, $e) = $c->get_coordinates($seq_obj->display_id);
			my $overlapping_seq = exists($self->{_participating_seqs}{$cseq})? $self->{_participating_seqs}{$cseq} : $self->{_associated_seqs}{$cseq};
			$regions_covered{$s} = $e if(!exists($regions_covered{$s}) || ($regions_covered{$s} < $e));
                        my $alignment_orientation = ($overlapping_seq->orientation eq $c->direction)? 'forward' : 'reverse';

			$connections_info .= ", " if(length($connections_info) > 0);
			$connections_info .= "$cseq (contained, $s-$e, $alignment_orientation)";

                        if(!defined($orientation)) {
                                $orientation = $alignment_orientation;
                        }
                        elsif($orientation ne $alignment_orientation) {
				print STDERR "\nError: orientation of " . $seq_obj->display_id . " based on $cseq differs from a previous one.\nInfo: $connections_info\n\n";
				print STDERR "\nScaffold (", $self->length, " bps):\n";
				foreach my $i (1 .. $self->num_links()) {
					my $link = $self->get_link($i);
					print STDERR "$i\t", $link->seq_obj->display_id, "\t", $link->seq_obj->length, "\t|\t", $link->get_pos_on_scaffold(), "\t", 
						     $link->orientation, "\n";
				}
				print STDERR "\nAssociated:\n";

				my $i = $self->num_links()+1;
				foreach my $link ($self->get_associated()) {
					print STDERR "$i\t", $link->seq_obj->display_id, "\t", $link->seq_obj->length, "\t|\t", $link->get_pos_on_scaffold(), "\t", 
						     $link->orientation, "\n";
					$i++;
				}

                                die "\nAborting (1)\n\n";
                        }
                }
        }

        foreach my $side (5, 3) {
                my $connections = $alignment_db->get_connections($seq_obj->display_id, $side, 'EE_CONNECT');
                foreach my $c (values %{$connections}) {
                        my ($seq2, $side2) = $c->other_end($seq_obj->display_id, $side);
			next if(!exists($self->{_participating_seqs}{$seq2}) && !exists($self->{_associated_seqs}{$seq2}));  

			my $overlapping_seq = exists($self->{_participating_seqs}{$seq2})? $self->{_participating_seqs}{$seq2} : $self->{_associated_seqs}{$seq2};
			my ($s, $e) = $c->get_coordinates($seq_obj->display_id, $side);
			$regions_covered{$s} = $e if(!exists($regions_covered{$s}) || ($regions_covered{$s} < $e));
			
                        my $alignment_orientation = ($overlapping_seq->orientation eq $c->direction)? 'forward' : 'reverse';
			$connections_info .= ", " if(length($connections_info) > 0);
			$connections_info .= "$seq2:$side2 ($side\', $s-$e, $alignment_orientation)";

                        if(!defined($orientation)) {
                                $orientation = $alignment_orientation;
                        }
                        elsif($orientation ne $alignment_orientation) {
				print STDERR "\nError: orientation of " . $seq_obj->display_id . " based on $seq2:$side2 differs from a previous one.\nInfo: $connections_info\n\n";
				print STDERR "\nScaffold (", $self->length, " bps):\n";
				foreach my $i (1 .. $self->num_links()) {
					my $link = $self->get_link($i);
					print STDERR "$i\t", $link->seq_obj->display_id, "\t", $link->seq_obj->length, "\t|\t", $link->get_pos_on_scaffold(), "\t", 
						     $link->orientation, "\n";
				}

				print STDERR "\nAssociated:\n";
				my $i = $self->num_links()+1;
				foreach my $link ($self->get_associated()) {
					print STDERR "$i\t", $link->seq_obj->display_id, "\t", $link->seq_obj->length, "\t|\t", $link->get_pos_on_scaffold(), "\t", 
						     $link->orientation, "\n";
					$i++;
				}

                                die "\nAborting (2)\n\n";
                        }
                }
        }

	if(!defined($orientation)) {
		print $OUT_ref $seq_obj->display_id, " is inconsistent - no consistent sequence aligned to it\n";
		return 0; 
	}

	my $last = 0;
	foreach my $s (sort {$a <=> $b} keys %regions_covered) {
		if($s > ($last+1)) {
			print $OUT_ref $seq_obj->display_id, " is inconsistent - not all regions were covered, e.g. (", ($last+1), ", ", ($s-1), ") ($connections_info)\n";
			return 0; 
		}
		$last = $regions_covered{$s} if($last < $regions_covered{$s});
	}

	if($last != $seq_obj->length) {
		print $OUT_ref $seq_obj->display_id, " is inconsistent - not all regions were covered ($connections_info)\n";
		return 0; 
	}
	
        # If we've reached this point then we are good
        print $OUT_ref $seq_obj->display_id, " is consistent ($connections_info)\n";

        return 1;
}

1;

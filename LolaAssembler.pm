
#!/usr/bin/perl

package LolaAssembler;

use strict;

use assembler::AlignmentDB;
use assembler::Scaffold;
use assembler::ScafLink;
use Bio::SeqIO;

#######################################################################################################################
sub new {
        my $class = shift;
        (@_ == 4) || die "\nError: LolaAssembler::new expected 3 arguments but recieved " . scalar(@_) . "\n\n";
        my $self = {
		_id   => shift,
                _sequences_db => shift,
                _alignment_db => shift,
		_pidentity_threshold => shift,
		_used_seqs    => {},
		_next_seq_id  => 1,
        };

        bless $self, $class;
        return $self;
}

#######################################################################################################################
sub num_sequences {
	my( $self ) = @_;

	return scalar(keys %{$self->{_sequences_db}});
}

#######################################################################################################################
sub id {
	my( $self ) = @_;

	return $self->{_id};
}

#######################################################################################################################
sub sequence_db {
	my( $self ) = @_;

	return $self->{_sequences_db};
}

#######################################################################################################################
# This subroutine will return one of the following:
# - A new contig that consists of several reads
# - the next read (chosen randomly) if no scaffold could be constructed
# - undef, if no unused read was left
#######################################################################################################################
sub next_scaffold {
	my( $self, $PROGRESS_ref, $LOG_ref, $GFF_ref ) = @_;

	print $PROGRESS_ref "\@LolaAssembler::next_scaffold/staring\n";
	print $PROGRESS_ref "\n---------------------------------------------------------------------------------------------------\n";
	# If all sequences were used - just return undef
	if(scalar(keys %{$self->{_used_seqs}}) == scalar(keys %{$self->{_sequences_db}})) {
		print $LOG_ref "No more sequeces left to check\n\n";
		print $PROGRESS_ref "\@LolaAssembler::next_scaffold/leaving\n";
		return undef;
	}

	# Look for the next scaffold
	my $new_scaffold = undef;
	my %links_addition = (); # Will store, for each sequence, the round in which it was added to the contig
        foreach my $seq_id (keys %{$self->{_sequences_db}}) {
		next if(exists($self->{_used_seqs}{$seq_id}));

                my ($can_continue5, $can_continue3) = ($self->can_continue($seq_id, 5, $self->{_used_seqs}), $self->can_continue($seq_id, 3, $self->{_used_seqs}));
                next if(($can_continue5 && $can_continue3) || (!$can_continue5 && !$can_continue3));

                # We've got a read to begin with
		$new_scaffold = new Scaffold();
                my $first_link = new ScafLink($self->{_sequences_db}{$seq_id}, 1, $self->{_sequences_db}{$seq_id}->length, ($can_continue3? 'forward' : 'reverse'));
                $new_scaffold->push_link($first_link);
                $self->{_used_seqs}{$seq_id} = 1;
		$links_addition{$seq_id} = 1;
		print $PROGRESS_ref "New contig\n";
		print $PROGRESS_ref "Link 1: $seq_id, ", $first_link->orientation, "\n";

                $self->elongate_scaffold($new_scaffold, \%links_addition, $PROGRESS_ref);
		last;
        }
	# If a contig was assembled - generate and return
	if(defined($new_scaffold)) {
		# If the new scaffold is only one sequence long - simply return the sequence
		if($new_scaffold->num_links() == 1) {
			print $PROGRESS_ref "\nWas not able to elongate the contig beyond the first link, will return ", $new_scaffold->get_link(1)->seq_obj->display_id, "\n";
			print $PROGRESS_ref "\@LolaAssembler::next_scaffold/leaving\n";
			return $new_scaffold->get_link(1)->seq_obj;
		}
		# Else - generate the scaffold
		print $PROGRESS_ref "\@LolaAssembler::next_scaffold/leaving\n";
		return $self->generate_scaffold($new_scaffold, $LOG_ref, $GFF_ref);
	}

	# If we are here then there are sequences, but none of them could be elongated. Just pick one and return it
	my $seq_id = undef;
	foreach my $s (keys %{$self->{_sequences_db}}) {
		if(!exists($self->{_used_seqs}{$s})) {
			$self->{_used_seqs}{$s} = 1;
			print $PROGRESS_ref "Returning $s (no read to assemble)\n\n";
			print $PROGRESS_ref "\@LolaAssembler::next_scaffold/leaving\n";
			return $self->{_sequences_db}{$s};
		}
	}

	# Houston, we have a problem: this was not supposed to happen!!!
	print $PROGRESS_ref "\nLolaAssembler::next_scaffold: something very weird happened, could not find an unused sequence\n\n";
	die "\nLolaAssembler::next_scaffold: something very weird happened, could not find an unused sequence\n\n";
}

##########################################################################################################################
# Recursively elongate the current contig ($curr_scaffold) and make sure that the elongated sequence is consistent with
# the data. Stop when no more elongation could be done, or when additional elongation results with inconsistencies.
##########################################################################################################################
sub elongate_scaffold {
        my ($self, $curr_scaffold, $links_addition_ref, $PROGRESS_ref) = @_;
        my $last_link = $curr_scaffold->get_link($curr_scaffold->num_links());
        my ($last_seq_obj, $last_side) = ($last_link->seq_obj, ($last_link->orientation eq 'forward')? 3 : 5);

	print $PROGRESS_ref "\@LolaAssembler::elongate_scaffold/starting\n";
        # We will pick the longest extension we can at this point so that we can check consistencies for the last two reads
        my $edge_connections = $self->{_alignment_db}->get_connections($last_seq_obj->display_id, $last_side, 'EE_CONNECT');
	my ($extending_link, $longest_extension) = (undef, 0);
        foreach my $c (values %{$edge_connections}) {

                my ($next_seq_id, $next_side) = $c->other_end($last_seq_obj->display_id, $last_side);
                my $orientation = ($c->direction eq $last_link->orientation)? 'forward' : 'reverse';

                # Sanity check
                exists($self->{_sequences_db}{$next_seq_id}) || 
			die "\nError, LolaAssembler::find_longest_scaffold: unable to find length for $next_seq_id (derived from connection " . $last_seq_obj->display_id . ":$last_side -> $next_seq_id:$next_side)\n\n";

                # If next sequence was used already - we cannot use it again
                next if(exists($self->{_used_seqs}{$next_seq_id}));

                # If the next sequence cannot be used - do not continue
                next if(!$self->can_continue($next_seq_id, $next_side, {}));

                # Add this scaffold and recursively call this subroutine
                my $overlap = $c->overlap($next_seq_id, $next_side);
                my ($start, $end) = ($orientation eq 'forward')? ($c->end($next_seq_id)+1, $self->{_sequences_db}{$next_seq_id}->length) : (1, $c->start($next_seq_id)-1);

                my $link = new ScafLink($self->{_sequences_db}{$next_seq_id}, $start, $end, $orientation);

                $curr_scaffold->push_link($link);
		if($curr_scaffold->length > $longest_extension) {
			$longest_extension = $curr_scaffold->length;
			$extending_link = $link;
		}
                $curr_scaffold->pop_link();
	}

	# Stop: no way to extend this one
	if(!defined($extending_link)) {
		print $PROGRESS_ref "No more sequences to add - finishing\n";
		print $PROGRESS_ref "\@LolaAssembler::elongate_scaffold/leaving\n";
		return;
	}

	print $PROGRESS_ref "Link ", ($curr_scaffold->num_links()+1), ": ", $extending_link->seq_obj->display_id, ", ", $extending_link->orientation, "\n";
	# Try to add the new sequence, recursively call this sub again if addition succeeded
	if($self->add_link_to_contig($curr_scaffold, $extending_link, $links_addition_ref, $PROGRESS_ref)) {
		$self->elongate_scaffold($curr_scaffold, $links_addition_ref, $PROGRESS_ref);
	}
	print $PROGRESS_ref "\@LolaAssembler::elongate_scaffold/leaving\n";
}

##########################################################################################################################
# This subroutine confirms that the proposed scaffolds is correct by looking at all sequences connected to any of the
# sequences in the proposed scaffold and validating that they are consistent with the reads in the scaffolds. At this 
# point we consider all connections, including those with % identity lower than threshold (assuming that all connections 
# provided by the user are ok to use).
#
# A sequence S is consitent with the scaffold if
# 1. It is one of the sequences of which the scaffold consists
# 2. It is contained within a sequence that is consistent with the scaffold
# 3. It has a 5-(5/3) connection with sequence A and a 3-(5/3) connection with sequence B were both A and B are consistent
#    with the scaffold and B comes after A in the proposed scaffold, in the proper orientation
#    and
#    Its edge-edge connections with A and B cover a longer area than any edge-middle or shared connections on both sides.
# 4. It does not have any edge-edge connections.
#
# As for (4) - note that middle-middle connections are ok, and shared connections with sequences that do have edge-edge
# connections with sequences that could be otherwise consistent with the scaffold will make those sequences inconsistent.
# As for 3 - I currently do not allow any shred/edge-middle connections for a sequence to be considered as consistent. This
# is just for making things simple, I plan to change this later.
##########################################################################################################################
sub add_link_to_contig {
        my ($self, $contig, $next_link, $links_addition_ref, $PROGRESS_ref) = @_;

	print $PROGRESS_ref "\@LolaAssembler::add_link_to_contig/staring\n";
	# Create a copy of $longest_scaffold and do everything on the copy. If all goes well then we will copy the information 
	# back to the original
	my $test_contig = new Scaffold();
	$test_contig->copy($contig);
        $test_contig->push_link($next_link);

	my %already_confirmed = ();
	my %new_confirmed = ($next_link->seq_obj->display_id => 1);		# All sequences that were confirmed in this round
	# Recruitment round for every sequence. Could be backdated if the sequence is actually connected to a previously 
	# added sequence but could be detected only once enough elongation was done
        my %recruited_from_round = ($next_link->seq_obj->display_id => $test_contig->num_links()); 
        my @sequences_to_check = ();						# A list of all sequences whose connections we still need to check
	my ($contig_first_seq_id, $contig_first_seq_side) = ($test_contig->get_link(1)->seq_obj->display_id, ($test_contig->get_link(1)->orientation eq 'forward')? 5 : 3);
	my ($contig_last_seq_id, $contig_last_seq_side) = ($next_link->seq_obj->display_id, ($next_link->orientation eq 'forward')? 3 : 5);

        # fetch all participating sequences and their orientation in the proposed scaffold
        for(my $i = 1; $i <= $test_contig->num_links(); $i++) {
                my $seq_id = $test_contig->get_link($i)->seq_obj->display_id;
                push(@sequences_to_check, $seq_id);
		if($i<$test_contig->num_links()) {
			exists($links_addition_ref->{$seq_id}) || die "\nError: something wrong (172633)\n\n";
	                $recruited_from_round{$seq_id} = $links_addition_ref->{$seq_id};
		}
		$already_confirmed{$seq_id} = 1;
        }

	# Fetch all consistent sequences
	foreach my $link ($test_contig->get_associated()) {
		my $seq_id = $link->seq_obj->display_id;
		push(@sequences_to_check, $seq_id);
		$already_confirmed{$seq_id} = 1;
		exists($links_addition_ref->{$seq_id}) || die "\nError: something wrong (172634)\n\n";
		$recruited_from_round{$seq_id} = $links_addition_ref->{$seq_id};
	}

        print $PROGRESS_ref "\n\nVerifying related sequences\n\n";

        # Go over all candidate sequences and check whether they are consistent. If they are not - push them back to the end of @sequences_to_check,
        # maybe when more consistent sequences are found they will become consistent. $left_to_check is set to the length of @sequences_to_check
        # every time a new consistent sequence was found since this can change the status of all sequences currently in line. Otherwise, if the
        # checked sequence was not found to be consistent then $left_to_check is decreased by 1. When we reach $left_to_check == 0 then either all
        # sequences in @sequences_to_check are inconsistent or no candidate sequences were left (all sequences are consistent with the scaffold).
        my $left_to_check = @sequences_to_check;
        while($left_to_check > 0) {
                my $seq = shift(@sequences_to_check);
                if(!$test_contig->add_associate($self->{_sequences_db}{$seq}, $self->{_alignment_db}, $PROGRESS_ref)) {
                        # We weren't able to verify this one - push again, maybe after more sequences are verfied it will be possible to verify this one as well
                        push(@sequences_to_check, $seq);
                        $left_to_check--;
                        next;
                }

		if(!exists($already_confirmed{$seq}) && !exists($new_confirmed{$seq})) {
			print $PROGRESS_ref "$seq: consistent, from link ", $recruited_from_round{$seq}, "\n";
			$new_confirmed{$seq} = 1 
		}

                # This one is good, add all sequences that were not checked yet and are connected to it
                foreach my $side (5, 3) {
			# Do not cehck first/last sequence ends
			next if(($seq eq $contig_first_seq_id) && ($side eq $contig_first_seq_side));
			next if(($seq eq $contig_last_seq_id) && ($side eq $contig_last_seq_side));

                        my $connections = $self->{_alignment_db}->get_connections($seq, $side, 'EE_CONNECT');
                        foreach my $c (values %{$connections}) {
                                my ($other_seq_id, $other_side) = $c->other_end($seq, $side);

				if(exists($self->{_used_seqs}{$other_seq_id})) {
					print $PROGRESS_ref "$other_seq_id ($seq:$side -> $other_seq_id:$other_side) ... skipping, already used\n";
					next;
				}

				if(exists($recruited_from_round{$other_seq_id})) {
					# In this case we need to replace the current round assignment if $seq was assigned in an earlier round.
					# $other_seq_id will be considered for sure.
					$recruited_from_round{$other_seq_id} = $recruited_from_round{$seq} if($recruited_from_round{$seq} < $recruited_from_round{$other_seq_id});
					next;
				}

                                
				print $PROGRESS_ref "$other_seq_id ($seq:$side -> $other_seq_id:$other_side) ... ";
                                # Check the position of this contig 
                                my ($pos_on_scaf, $orientation) = $test_contig->find_position($self->{_sequences_db}{$other_seq_id}, $self->{_alignment_db});

                                if(!defined($pos_on_scaf)) {
                                        print $PROGRESS_ref "\nError: could not find position for $other_seq_id even though it is linked to $seq that is associated/belongs to the scaffold\n\\n";
                                        die "\nError: could not find position for $other_seq_id even though it is linked to $seq that is associated/belongs to the scaffold\n\n";
                                }

                                # If sequence starts before the 5' end of the contig (negative position) then this implies 
				# inconsistency caused by the sequence we are trying to add (all other "backward" connections were already added).
				if($pos_on_scaf < 1) {
					# This one must be an inconsistency. In order to be correct the first sequence must be connected to other_seq
					# on its outwards (contig 5') end, but if this was the case we would have not been using it to start a contig!
					print $PROGRESS_ref "INCONSISTENT: expected beginning position of $other_seq_id is before the 5\' end of the contig\n\n";
					# We can simply return from here because up until now we modified only local data structures
					print $PROGRESS_ref "\@LolaAssembler::add_link_to_contig/leaving\n";
					return 0; 
				}
				# If sequence ends afer 3' end - skip it as it does not have to be consistent at this stage
				elsif($pos_on_scaf+$self->{_sequences_db}{$other_seq_id}->length > $test_contig->length) {
					print $PROGRESS_ref "ignoring: expected end position is beyond the contig\'s end\n";
                                        next;
                                }

                                # Add sequence to candidates list
				print $PROGRESS_ref "checking\n";
                                push(@sequences_to_check, $other_seq_id);
                                $recruited_from_round{$other_seq_id} = $recruited_from_round{$seq};
			}
		}

		# Also add the contained sequences: these all can go directly to the consistent list. Add them also the the candidate list
		# so that it will be possible to check sequences to which they are connected
		my $contained_seqs = $self->{_alignment_db}->get_contained($seq);

		foreach my $s (keys %{$contained_seqs}) {
			if(exists($self->{_used_seqs}{$s})) {
				print $PROGRESS_ref "$s ($s contained in $seq) ... skipping, already used\n";
				next;
			}
			if(!exists($recruited_from_round{$s})) {
				push(@sequences_to_check, $s);
				$recruited_from_round{$s} = $recruited_from_round{$seq};
				print $PROGRESS_ref "$s ($s contained in $seq) ... checking\n";
			}
		}
		# We do not add sequences that contain the newly added sequence since new sequences should be added only on the basis of 5'/3' connections.

                # Now start the counting again
                $left_to_check = @sequences_to_check;
        }

        # Finished: If we don't have inconsistent sequences then we can elongate
        if(@sequences_to_check == 0) {
                print $PROGRESS_ref "\nAll related sequences are consistent\n\n";
		$contig->copy($test_contig);
		foreach my $seq_id(keys %new_confirmed) {
			$self->{_used_seqs}{$seq_id} = 1;
			$links_addition_ref->{$seq_id} = $recruited_from_round{$seq_id};
		}

		print $PROGRESS_ref "\@LolaAssembler::add_link_to_contig/leaving\n";
		return 1;
	}

	# Otherwise - we need to check whether the inconsistency is the result of the current sequence that was added or
	# is it because of a sequence that was added in one of the previous rounds but could only be discovered now  

	print $PROGRESS_ref "\nInconsistent sequences exists\n\n";
	my $round = undef;
	foreach my $seq_id (@sequences_to_check) {
		print $PROGRESS_ref "$seq_id: INCONSISTENT, from link ", $recruited_from_round{$seq_id}, "\n";
                $round = $recruited_from_round{$seq_id} if(!defined($round) || ($round > $recruited_from_round{$seq_id}));
	}

	# No corrections are required in the following case. This means that either
	# - a sequence was found that is connected to the newly added sequence towards the 5' end of contig, which simply means
	#   that the new sequence cannot be added
	# - a sequence was found that is connected towards the 3' end of the contig of the last link that is inconsistent with
	#   the new one we are trying to add. This means that the contig should simply be left as it is but not elongated.  
	#
	if($round >= $contig->num_links()) {
		print $PROGRESS_ref "\nNo need to change the contig, finishing elongation\n";
		print $PROGRESS_ref "\@LolaAssembler::add_link_to_contig/leaving\n";
		return 0; 
	}

	# We have to remove all sequences that were added after round $round. 
	print $PROGRESS_ref "Inconsistencies were found that result from previously added sequences, will remove links after round $round\n";
	$self->remove_links_from_contig($contig, $round, $links_addition_ref, $PROGRESS_ref);

	print $PROGRESS_ref "\@LolaAssembler::add_link_to_contig/leaving\n";
	return 0;
}

##########################################################################################################################
# This subroutine will create a new scaffold with only reads that were added up to round $round
sub remove_links_from_contig {
	my ($self, $contig, $round, $links_addition_ref, $PROGRESS_ref) = @_;
	print $PROGRESS_ref "\@LolaAssembler::remove_links_from_contig/starting\n";

	return if($round == $contig->num_links);

	print $PROGRESS_ref "Removing sequences that were added between rounds $round-", $contig->num_links(), "\n";

	my $new_contig = new Scaffold();
	my %seq_obj_to_add = ();

	# Free the links we are removing, keep ids of those we are about to keep
	my @seq_ids = keys %{$links_addition_ref};
	foreach my $seq_id (@seq_ids) {
		if(($links_addition_ref->{$seq_id} > $round) || 
		   (($links_addition_ref->{$seq_id} == $round) && ($seq_id ne $contig->get_link($round)->seq_obj->display_id))) {
			print $PROGRESS_ref "$seq_id: removing (added on round ", $links_addition_ref->{$seq_id}, ")\n";
			delete($links_addition_ref->{$seq_id});
			delete($self->{_used_seqs}{$seq_id});
		}
		else {
			print $PROGRESS_ref "$seq_id: keeping (added on round ", $links_addition_ref->{$seq_id}, ")\n";
			$seq_obj_to_add{$seq_id} = 0;
		}
	}

	# Add the backbone of the contig
        for(my $i = 1; $i <= $round; $i++) {
		my $link = $contig->get_link($i);
                $new_contig->push_link($link);
		delete($seq_obj_to_add{$link->seq_obj->display_id});
	}

	# At this point we have to add the associated sequences in several iterations since maybe some of them need others
	# to be added. We expect all sequenced to be added, if not this is an error
	my $left_to_add = scalar(keys %seq_obj_to_add);
	while($left_to_add > 0) {
		my @seq_ids = sort {$seq_obj_to_add{$a} <=> $seq_obj_to_add{$b}} keys %seq_obj_to_add;
		foreach my $seq_id (@seq_ids) {
			if($new_contig->add_associate($self->{_sequences_db}{$seq_id}, $self->{_alignment_db}, $PROGRESS_ref)) {
				delete($seq_obj_to_add{$seq_id});
			}
		}
		# The following case means that we were not able to add any more sequences. This is an error
		if($left_to_add == scalar(keys %seq_obj_to_add)) {
			die "\nError: could not add all sequences that were supposed to be added\n\n";
		}
		$left_to_add = scalar(keys %seq_obj_to_add);
	}

	# If we are here then all is well. Just copy the information from the new contig to the old one
	$contig->copy($new_contig);
	print $PROGRESS_ref "\@LolaAssembler::remove_links_from_contig/leaving\n";
}


##########################################################################################################################
# The following sub checks whether sequence $seq_id could be elongated to side $side. For elongation there must not be 
# a connection of $seq:$side to a middle/shared region and also there must be an alignment that passes % identity 
# threshold ($self->{_pidentity_threshold}).
##########################################################################################################################
sub can_continue {
	my ($self, $seq_id, $side, $used_seqs_ref) = @_;
        return 0 if(!$self->can_continue_middle_and_shared($seq_id, $side));

        # Now check edge-edge
        my $connections = $self->{_alignment_db}->get_connections($seq_id, $side, 'EE_CONNECT');

        return 0 if(scalar(keys %{$connections}) == 0);

        my $can = 0;
        foreach my $c (values %{$connections}) {
		# New: for elongation we consider only connections at a given % identity
		next if($c->pidentity < $self->{_pidentity_threshold});
                my ($other_seq_id, $other_side) = $c->other_end($seq_id, $side);
		# If the other connection cannot be elongated towards this one - we cannot elongate this one
                return 0 if(!$self->can_continue_middle_and_shared($other_seq_id, $other_side));
		# If the other connection is not used - it can be used for elongating this one (if all other connections are ok)
                $can = 1 if(!exists($used_seqs_ref->{$other_seq_id}));
        }

        return $can;
}

##########################################################################################################################
sub can_continue_middle_and_shared {
        my ($self, $seq_id, $side) = @_;

	# If $seq_id has an unknown connection we cannot continue
	my $connections = $self->{_alignment_db}->get_unknowns($seq_id);
	return 0 if(scalar(keys %{$connections}) > 0);

        # Check edge-middle connections
        $connections = $self->{_alignment_db}->get_connections($seq_id, 'Middle', 'EM_CONNECT');
        foreach my $c (values %{$connections}) {
		# New: we consider only connections above threshold 
		next if($c->pidentity < $self->{_pidentity_threshold});
                my ($other_seq_id, $other_end) = $c->other_end($seq_id, 'Middle');

                if(($side == 5) && (($c->direction eq 'forward') && ($other_end == 5)) || (($c->direction eq 'reverse') && ($other_end == 3))) {
                        return 0;
                }
                if(($side == 3) && (($c->direction eq 'forward') && ($other_end == 3)) || (($c->direction eq 'reverse') && ($other_end == 5))) {
                        return 0;
                }
        }

        # Check shared connections
        $connections = $self->{_alignment_db}->get_connections($seq_id, $side, 'EE_SHARED');
        foreach my $c (values %{$connections}) {
		return 0 if($c->pidentity >= $self->{_pidentity_threshold});
	}

        return 1;
}

##########################################################################################################################
sub generate_scaffold {
        my ($self, $scaffold_assembly_obj, $LOG_ref, $GFF_ref) = @_;

        my $new_scaf = undef;

	my $display_id = "M" . $self->id() . "." . $self->{_next_seq_id};
	$display_id .= ".unsafe" if($scaffold_assembly_obj->num_inconsistent() > 0);
	$self->{_next_seq_id}++;   

        print $LOG_ref "\n*** $display_id ***\n\n";
	print $LOG_ref "Number of links: ", $scaffold_assembly_obj->num_links(), "\n";

        my $seq = $scaffold_assembly_obj->seq;
	if(length($seq) == 0) {
		print $LOG_ref "\nError (LolaAssembler::generate_scaffold): length of new scaffold is 0\n";
		die "\nError (LolaAssembler::generate_scaffold): length of new scaffold is 0\n\n";
	}

        my $gc = int(1000 * ($seq =~  tr/GC/GC/) / (length($seq) - ($seq =~ tr/N/N/))) / 10;
        my $desc = "/length=" . length($seq) . " /%G+C=" . $gc;

        $new_scaf = new Bio::Seq(-display_id => $display_id, -seq => $seq, -desc => $desc);

        print $LOG_ref "Size: ", $new_scaf->length, "\n";
        print $LOG_ref "\%G+C: $gc\n";
        print $LOG_ref "\nComposing sequences:\n\n";
        print $LOG_ref "\t                 \t      \t|\tSequence  \t        \t           \t|\tContributed region\n";
        print $LOG_ref "\tsequence         \tlength\t|\tscaf_start\tscaf_end\torientation\t|\tseq_start\tseq_end\tscaf_start\tscaf_end\n";

	my $msg = '';
        foreach my $i (1 .. $scaffold_assembly_obj->num_links()) {
                my $link = $scaffold_assembly_obj->get_link($i);
                if(!exists($self->{_used_seqs}{$link->seq_obj->display_id})) {
			$msg .= "ERROR: " . $link->seq_obj->display_id . " was not marked as used\n";
		}

                my $region_pos = $link->get_pos_on_scaffold;
                $region_pos += ($link->orientation eq 'forward')? ($link->start-1) : ($link->seq_obj->length - $link->end);

                print $LOG_ref "$i\t", $link->seq_obj->display_id, "\t", $link->seq_obj->length, "\t|\t",
                          $link->get_pos_on_scaffold, "\t", ($link->get_pos_on_scaffold+$link->seq_obj->length-1), "\t", $link->orientation, "\t|\t",
                          $link->start, "\t", $link->end, "\t", $region_pos, "\t", ($region_pos+($link->end-$link->start)), "\n";
		# M2.1	Lola	mol-32-1605-037016	1	9371	.	-	.	seq_length=16279; read_length=9371
		print $GFF_ref "$display_id\tLola\t", $link->seq_obj->display_id, "\t", $link->get_pos_on_scaffold, "\t", 
				($link->get_pos_on_scaffold+$link->seq_obj->length-1), "\t.\t", (($link->orientation eq 'forward')? '+' : '-'), "\t.\t",
				"seq_length=", $new_scaf->length, "; seq_gc=$gc; read_length=", $link->seq_obj->length, "\n"; 
        }

        my $i = $scaffold_assembly_obj->num_links()+1;
        if($scaffold_assembly_obj->num_associated() > 0) {
                print $LOG_ref "\nConsistent sequences:\n\n";
                print $LOG_ref "\tsequence          \tlength\t|\tscaf_start\tscaf_end\torientation\n";
                foreach my $link ($scaffold_assembly_obj->get_associated()) {
			if(!exists($self->{_used_seqs}{$link->seq_obj->display_id})) {
				$msg .= "ERROR: " . $link->seq_obj->display_id . " (consistent) was not marked as used\n";
			}

                        print $LOG_ref "$i\t", $link->seq_obj->display_id, "\t", $link->seq_obj->length, "\t|\t", $link->get_pos_on_scaffold, "\t",
                                                ($link->get_pos_on_scaffold+$link->seq_obj->length-1), "\t", $link->orientation, "\n";
                        $i++;
			print $GFF_ref "$display_id\tLola\t", $link->seq_obj->display_id, "\t", $link->get_pos_on_scaffold, "\t", 
					($link->get_pos_on_scaffold+$link->seq_obj->length-1), "\t.\t", (($link->orientation eq 'forward')? '+' : '-'), "\t.\t",
					"seq_length=", $new_scaf->length, "; seq_gc=$gc; read_length=", $link->seq_obj->length, "\n"; 
                }
        }
        else {
                print $LOG_ref "\nNo consistent sequences\n\n";
        }

        if($scaffold_assembly_obj->num_inconsistent() > 0) {
                print $LOG_ref "\n*** INCONSISTENT sequences: ***\n\n";
                print $LOG_ref "\tsequence          \tlength\t|\tscaf_start\tscaf_end\torientation\n";

                foreach my $link ($scaffold_assembly_obj->get_inconsistent()) {
                        # We don't put inconsistent sequences in the used sequence list - it's ok to use them
                        print $LOG_ref "$i\t", $link->seq_obj->display_id, "\t", $link->seq_obj->length, "\t|\t", $link->get_pos_on_scaffold, "\t",
                                                ($link->get_pos_on_scaffold+$link->seq_obj->length-1), "\t", $link->orientation, "\n";
                        $i++;
                }
        }
        else {
                print $LOG_ref "\nNo inconsistent sequences\n\n";
        }

	if($msg ne '') {
		print $LOG_ref "==> Error message(s) <==\n\n$msg\n\n";
	}
        return $new_scaf;
}

1;

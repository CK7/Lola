#!/usr/bin/perl

package AlignmentDB;

use assembler::Alignment;
use assembler::Connect;
use assembler::Container;
use assembler::Unknown;
use strict;

#######################################################################################################################
sub new {
	my $class = shift;

	my $self = {
		_files       => [],
		_connections => {},
		_containers  => {},
		_contained  => {},
		_unknowns  => {},
	};

	bless $self, $class;
	return $self;
}

#######################################################################################################################
sub get_sequences_with_connections {
	my( $self ) = @_;

	my @seqs = keys %{$self->{_connections}};

	return \@seqs; 
}

#######################################################################################################################
sub get_sequences_containers {
	my( $self ) = @_;

	my @seqs = keys %{$self->{_containers}};

	return \@seqs; 
}

#######################################################################################################################
sub get_sequences_contained {
	my( $self ) = @_;

	my @seqs = keys %{$self->{_contained}};

	return \@seqs; 
}

#######################################################################################################################
# get_connection
# 
# get_connection(seq1, side1, seq2, side2)
# get_connection(seq1, side1, edge2)
# get_connection(edge1, edge2)
#
# edge = "seq:side"
#######################################################################################################################
sub get_connection {
	my( $self, $seq1, $side1, $seq2, $side2 ) = @_;

	if(!defined($seq2)) {
		$side1 =~ /^(.+):([53]|Middle)$/ || die "\nError, AlignmentDB::get_connection: received 2 arguments but second was not in edge format: $side1\n\n";
		($seq2, $side2) = ($1, $2);
		$seq1 =~ /^(.+):([53]|Middle)$/ || die "\nError, AlignmentDB::get_connection: received 2 arguments but first was not in edge format: $seq1\n\n";
		($seq1, $side1) = ($1, $2);
	}
	elsif(!defined($side2)) {
		return $self->{_connections}->{$seq1}{$side1}{$seq2};
	}

	return $self->{_connections}->{$seq1}{$side1}{"$seq2:$side2"};
}

#######################################################################################################################
sub get_connections {
	my( $self, $seq, $side, $type ) = @_;

	$self->{_connections}->{$seq}{$side} = {} if(!exists($self->{_connections}->{$seq}{$side}));


	return $self->{_connections}->{$seq}{$side} if(!defined($type));

	my $ret_val = {};
	foreach my $connection (keys %{$self->{_connections}->{$seq}{$side}}) {
		$ret_val->{$connection} = $self->{_connections}->{$seq}{$side}{$connection} if($self->{_connections}->{$seq}{$side}{$connection}->type eq $type);
	}

	return $ret_val;
}

#######################################################################################################################
sub get_containers {
	my( $self, $contained_seq ) = @_;

#	$self->{_contained}->{$contained_seq} = {} if(!exists($self->{_contained}->{$contained_seq}));
	return {} if(!exists($self->{_contained}->{$contained_seq}));

	return $self->{_contained}->{$contained_seq};
}

#######################################################################################################################
sub get_contained {
	my( $self, $container_seq ) = @_;

#	$self->{_containers}->{$container_seq} = {} if(!exists($self->{_containers}->{$container_seq}));
	return {} if(!exists($self->{_containers}->{$container_seq}));

	return $self->{_containers}->{$container_seq};
}

#######################################################################################################################
sub get_unknowns {
	my( $self, $seq_id ) = @_;

#	$self->{_unknowns}->{$seq_id} = {} if(!exists($self->{_unknowns}->{$seq_id}));
	return {} if(!exists($self->{_unknowns}->{$seq_id}));

	return $self->{_unknowns}->{$seq_id};
}

#######################################################################################################################
sub add_connections {
	my( $self, $file ) = @_;
	
	open(IN, $file) || die "\nError: failed to open $file for reading\n\n";
	push(@{$self->{_files}}, $file);

	my ($ncontained, $nedge_connections, $nmiddle_connections, $nshared, $nidentical) = (0, 0, 0, 0, 0);
        while(<IN>) {
                chomp;
                my @fs = split(/\t/);
                my $op = shift(@fs);

                # CONNECTED     13ft.mol-32-15fa-049956 3       16ft.mol-32-1605-074376 3       6173    10387   8788    4574    10387   8788    ok
                # SHARED        16ft.mol-32-1605-053525 3       16ft.mol-32-1605-045204 3       4911    9146    3092    7327    9148    7327    Major
                if((($op eq 'CONNECTED') || ($op eq 'SHARED')) && (@fs == 11)) {
			my $c = new Connect(@fs);
			$self->{_connections}->{$c->seq1}{$c->side1}{$c->seq2 . ':' . $c->side2} = $c;
			$self->{_connections}->{$c->seq2}{$c->side2}{$c->seq1 . ':' . $c->side1} = $c;
                }
		# IDENTICAL     16ft.mol-32-1605-019817 16ft.mol-32-1605-023231 1       7009    7009    1       7055    7036    ok
                # CONTAINS      19ft.mol-32-15ef-029213 13ft.mol-32-15fa-068442 2049    7102    1       5053    7852    5053    ok
                elsif((($op eq 'CONTAINS') || ($op eq 'IDENTICAL')) && (@fs == 9)) {
			my $c = new Container(@fs);
			$self->{_containers}->{$c->container}{$c->contained} = $self->{_contained}->{$c->contained}{$c->container} = $c;

                        if($op eq 'IDENTICAL') {
				$self->{_containers}->{$c->contained}{$c->container} = $self->{_contained}->{$c->container}{$c->contained} = $c;
                        }
                }
		# UNKNOWN mol-32-1605-031898      mol-32-1605-084167      13488   9037    Major
		elsif(($op eq 'UNKNOWN') && (@fs == 5)) {
			my $c = new Unknown(@fs);
			$self->{_unknowns}->{$c->seq1}{$c->seq2} = $self->{_unknowns}->{$c->seq2}{$c->seq1} = $c;
		}
                else {
                        die "\nAlignmentDB::add($file): Unexpeected input line $_\n\n";
                }
        }
        close(IN);
}

1;

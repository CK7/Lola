#!/usr/bin/perl

package Alignment;

#######################################################################################################################
sub new {
	my $class = shift;
	(@_ == 9) || die "\nError: Alignment::new expected 9 arguments but recieved only " . scalar(@_) . ": " . (join(", ", @_)) . "\n\n"; 
	my $self = {
        	_seq1   => shift,
        	_seq2   => shift,
		_start1 => shift,
		_end1   => shift,
		_start2 => shift,
		_end2   => shift,
		_len1   => shift,
		_len2   => shift,
		_pidentity => shift, 
    	};

	$self->{_direction} = (($self->{_end1} - $self->{_start1})*($self->{_end2} - $self->{_start2}) > 0)? 'forward' : 'reverse';	
	($self->{_start1}, $self->{_end1}) = ($self->{_end1}, $self->{_start1}) if($self->{_start1} > $self->{_end1});
	($self->{_start2}, $self->{_end2}) = ($self->{_end2}, $self->{_start2}) if($self->{_start2} > $self->{_end2});

	bless $self, $class;
	return $self;
}

#######################################################################################################################
sub type {
	my( $self ) = @_;
	die "\nError: Alignment::type subroutine called, should be pure virtual\n\n";
}

#######################################################################################################################
sub seq1 {
	my( $self ) = @_;

	return $self->{_seq1};
}

#######################################################################################################################
sub seq2 {
	my( $self ) = @_;

	return $self->{_seq2};
}

#######################################################################################################################
sub start {
	my( $self, $seq ) = @_;
	return ($seq eq $self->{_seq1})? $self->{_start1} : $self->{_start2};
}

#######################################################################################################################
sub end {
	my( $self, $seq ) = @_;

	return ($seq eq $self->{_seq1})? $self->{_end1} : $self->{_end2};
}

#######################################################################################################################
sub get_coordinates {
	my ( $self, $seq ) = @_;

	return ($seq eq $self->{_seq1})? ($self->{_start1}, $self->{_end1}) : ($self->{_start2}, $self->{_end2});
}

#######################################################################################################################
sub length {
	my( $self, $seq ) = @_;
	return ($seq eq $self->{_seq1})?  $self->{_len1} : $self->{_len2};
}

#######################################################################################################################
sub overlap {
	my( $self, $seq ) = @_;

	return ($seq eq $self->{_seq1})? $self->{_overlap1} : $self->{_overlap2};
}

#######################################################################################################################
sub direction {
	my( $self ) = @_;
	return $self->{_direction}; 
}

#######################################################################################################################
sub pidentity {
	my( $self ) = @_;
	return $self->{_pidentity}; 
}

#######################################################################################################################
sub break_edge_name {
	my( $self, $edge ) = @_;

	($edge =~ /^(.+):([53]|Middle)$/) || die "\nAlignment::break_edge_name: unexpected edge name, $edge\n\n";

	return ($1, $2);	
}	

1;

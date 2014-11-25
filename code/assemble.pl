#!/usr/bin/perl
use strict;
use Bio::SeqIO;

BEGIN {
	unshift(@INC, '/data4/Rifle/Moleculo/2.Moleculo-assembly');
}

my $VERSION = 'v0.62';
my $pidentity_threshold = 99;

use assembler::Connect;
use assembler::Container;
use assembler::AlignmentDB;
use assembler::ScafLink;
use assembler::Scaffold;
use assembler::LolaAssembler;
use assembler::LolaDB;

my %sequences = ();

($#ARGV >= 2) || die "\nUsage: $0 [-p <\% identity threshold] <sequence-file> <connection-file1> [... <connection-filen>] <out-prefix>\n\n";

my $out_prefix = pop(@ARGV);
my @connection_files = pop(@ARGV);
my $sequence_file = pop(@ARGV);
if(@ARGV > 0) {
	if(($ARGV[0] ne '-p') || (@ARGV != 2) || ($ARGV[1] <= 0) || ($ARGV[1] > 100)) {
		die "\nError: unknown paraameters\n\n";
	}
	$pidentity_threshold = pop(@ARGV);
}

my $out_scaffolds_file = "$out_prefix.scaffolds.fna";
my $out_progress_file = "$out_prefix.progress.log";
my $out_log_file = "$out_prefix.scaffolds.log";
my $out_gff_file = "$out_prefix.scaffolds.gff";

open(LOG, ">$out_log_file") || die "\nCannot write to $out_log_file\n\n";
open(GFF, ">$out_gff_file") || die "\nCannot write to $out_gff_file\n\n";
open(PROGRESS, ">$out_progress_file") || die "\nCannot write to $out_progress_file\n\n";

print STDERR "Lola $VERSION, ", get_str_time(), "\n\n";
print LOG "Lola $VERSION, ", get_str_time(), "\n\n";
print PROGRESS "Lola $VERSION, ", get_str_time(), "\n\n";

print STDERR "Reading sequence lengths from $sequence_file\n";
print PROGRESS "Reading sequence lengths from $sequence_file\n";

my $in = new Bio::SeqIO(-file => $sequence_file);
while(my $seq = $in->next_seq) {
        $sequences{$seq->display_id} = $seq;
}
print STDERR "Finished, ", scalar(keys %sequences), " sequences read\n\n";
print PROGRESS "Finished, ", scalar(keys %sequences), " sequences read\n\n";

my $alignment_db = new AlignmentDB;

#### Loop over all connection files provided ####
while(@connection_files > 0) {
	my $connection_file = shift(@connection_files);
	print STDERR "Reading $connection_file ... ";
	print PROGRESS "Reading $connection_file ... ";
	$alignment_db->add_connections($connection_file);
	print STDERR "ok\n";
	print PROGRESS "ok\n";
}

#### Generate scaffolds ####
my $scaffolds = {};

my $asm_db = new LolaDB(\%sequences, $alignment_db, $pidentity_threshold);

print STDERR $asm_db->num_clusters(), " clusters\n\n";
foreach my $i (1 .. $asm_db->num_clusters()) {
	my $asm = $asm_db->get_cluster($i);
	print STDERR "$i (", scalar(keys %{$asm->sequence_db}), " sequences)\n";
	print STDERR join("\n", keys %{$asm->sequence_db}), "\n";

	print PROGRESS "\n########################################### CLUSTER $i ###########################################\n\n"; 
	print PROGRESS $asm->num_sequences(), " reads\n";
	while (my $scaf_seq_obj = $asm->next_scaffold(\*PROGRESS, \*LOG, \*GFF)) {
		$scaffolds->{$scaf_seq_obj->display_id} = $scaf_seq_obj;
		if($scaf_seq_obj->display_id =~ /unsafe/) {
			print STDERR "Check: ", $scaf_seq_obj->display_id, " (", $scaf_seq_obj->length, " bps)\n"; 
		}
		print PROGRESS "\nFINISHED: ", $scaf_seq_obj->display_id, "\n"; 
	}
}

########## NEW CODE v0.62 ###########
# The following code is supposed to remove reads that were not assembled and are contained within other reads that were not assembled. 
# Such reads must appear in clusters only as they are clustered together with their container.
my %contained_reads = ();
foreach my $seq_contained_id (@{$alignment_db->get_sequences_contained()}) {
	$contained_reads{$seq_contained_id} = 1;
}

foreach my $seq_contained_id (keys %contained_reads) {
	if(exists($scaffolds->{$seq_contained_id})) {
		my $containers_of_ref = $alignment_db->get_containers($seq_contained_id);
		foreach my $container_id (keys %{$containers_of_ref}) {
			# We want a container that will be reported 
			next if(exists($contained_reads{$container_id}));
			if(exists($scaffolds->{$container_id})) {
				my $container_obj = $containers_of_ref->{$container_id};
				# M2.8	Lola	mol-32-1605-051725	6172	11883	.	+	.	seq_length=11883; seq_gc=40.7; read_length=5712
				print GFF "$container_id\tMoleculo\t$seq_contained_id\t", $container_obj->start($container_id), "\t", $container_obj->end($container_id), "\t.\t",
						(($container_obj->direction eq 'forward')? '+' : '-'), "\t.\tseq_length=", $scaffolds->{$container_id}->length, 
						"; seq_gc=N/A; read_length=", $scaffolds->{$seq_contained_id}->length, "\n";
				delete($scaffolds->{$seq_contained_id});
				last;
			}
		}
	}
}
########## NEW CODE v0.62 ###########

my $out = new Bio::SeqIO(-file => ">$out_scaffolds_file", -format => 'fasta');

foreach my $scaf (reverse sort {$scaffolds->{$a}->length <=> $scaffolds->{$b}->length} keys %{$scaffolds}) {
	$out->write_seq($scaffolds->{$scaf});
}

foreach my $seq_obj (values %{$asm_db->get_unclustered()}) {
	$out->write_seq($seq_obj);
}

print STDERR "\nFinished, ", get_str_time(), "\n\n";
print LOG "\nFinished, ", get_str_time(), "\n\n";
print PROGRESS "\nFinished, ", get_str_time(), "\n\n";

close(LOG);
close(GFF);
close(PROGRESS);

##########################################################################################################################
sub get_str_time {
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);

	return sprintf("%02d/%02d/%04d %02d:%02d:%02d", 1+$mon, $mday, 1900+$year, $hour, $min, $sec);
}

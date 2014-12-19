#!/usr/bin/env perl

# Determine pseudogene candidates by tiling and filtering tblastn output, where
# queries are the amino acid sequences of functional/high-confidence genes, and
# subject is the DNA sequence of the genome in question.  It produces a GFF file
# describing the candidate pseudogenes found.
#
# For now this is explicitly tied to the Norway spruce genome project products
# and pipelines.  This will be loosened as pseudogene analyses progress.
#
# The name 'psari' is a hat-tip to 'Sari', María Rosario García-Gil, who has
# worked on pseudogenes in pine at Umeå Plant Science Centre for some time.
#
# Douglas G. Scofield, Evolutionary Biology Centre
# (formerly Umeå Plant Sciences Centre, Umeå Sweden)
# douglas.scofield@ebc.uu.se
# douglasgscofield@gmail.com
#
# Originally I'd intended to build tiling using BioPerl's mechanisms under
# Bio::Search::Tiling::*, but that didn't quite do what I wanted, so I wrote my
# own tiler from scratch.  I still use BioPerl to parse the tblastn output.
#
# TODO
# - Generalise interface
# - Better documentation
# - Incorporate inference re: pathway of pseudogene origin (could be separate)
# - Test against known set of pseudogenes, of course!
# - Comparisons against pseudopipe and other detectors
#
# CHANGELOG
#
# 2014-12-19
# - indentation via perltidy
# 2014-11-20
# - little updates, not yet generalised
# 2013-05-17
# - initial upload to github
# 2012-12-12
# - completion for genome paper

use strict;
use warnings;

use Data::Dumper;
$Data::Dumper::Sortkeys = 1;
use List::Util qw/sum/;

# to find out the methods for an object:
# print "Instance METHOD is " . Dumper(\%{ref ($this_hsp)."::"});

use Bio::SearchIO;
use Bio::Search::Result::BlastResult;

#use Bio::Search::Hit::PsiBlastHit;
#use Bio::Search::BlastUtils;
#use Bio::Search::Tiling::MapTiling;
#use Bio::Search::Tiling::MapTileUtils;

# ----------------------------

# pipeline utility variables

my $num_sources                   = 0;         # cumulative counter
my $num_hits                      = 0;         # cumulative counter
my $num_hsps                      = 0;         # cumulative counter
my $num_candidates                = 0;         # cumulative counter
my $num_blocks                    = 0;         # cumulative counter
my $num_stop_codons_in_candidates = 0;
my $stopchar                      = "\\\*";    # character used for stop codon

# ----------------------------

# pipeline identity and output options

my $this_pipeline                 = "psari";
my $query_source                  = "augustus_aa";           # will join algorithm and this_pipeline later
my $candidate_algorithm           = undef;
my $candidate_source              = undef;                   # col 2 of GFF; will join query_source, algorithm and this_pipeline
my $candidate_type                = "pseudogenic_region";    # what is or Sequence Ontology type of our candidates?
my $candidate_number_format       = "pr%03d";                # pseudogenic region
my $candidate_block_type          = "pseudogenic_block";     # what is a type of our candidate blocks?
my $candidate_block_number_format = "pb%02d";                # pseudogenic block
my $fmt                           = "%1.4f";                 # format used for printing fractional numbers
my $scorefmt                      = "%1.3f";                 # format used for printing scores

my $produce_bed = 0;

my $progress_mod   = 100000;                                 # print a message every this many hits, set to 0 for no message
my $print_progress = 0;                                      # set to 1 if we are to print a progress message this iteration

my $verbose = 0;

# ----------------------------

sub fix_query_name($) {
    my $ans = shift;
    $ans =~ s/^lcl\|//;
    return $ans;
}

sub fix_subject_name($) {
    my $ans = shift;
    $ans =~ s/^lcl\|//;
    return $ans;
}

# ----------------------------

# pipeline filtering and configuration options

my @trim_units      = ();    # part of dealing with overlap, defined depending on algorithm, e.g. TBLASTN is (1, 3)
my $query_overlap   = 1;     # 0 is non-overlapping, -1 is 1-unit gap
my $subject_overlap = 0;
die "subject overlap isn't supported" if $subject_overlap;
my @num_query_overlap_trim   = ( 0, 0 );
my @num_subject_overlap_trim = ( 0, 0 );

# filter suggestions mostly following Zou et al 2009 Plant Physiology

my $min_frac_covered_candidate            = 0.4;    # if a candidate covers less than this of the source, drop it
my $num_filter_min_frac_covered_candidate = 0;

my $min_frac_identical                              = 0.4;    # if an HSP or covered part of a candidate has less than this identity, drop it
my $num_filter_min_frac_identical_hsp               = 0;
my $min_covered_frac_identical_candidate            = 0.4;    # if an HSP or covered part of a candidate has less than this identity, drop it
my $num_filter_min_covered_frac_identical_candidate = 0;

my $min_single_hit_length        = 15;                        # if an HSP is shorter than this number of AA, drop it amino acids
my $num_filter_single_hit_length = 0;

my $overlap_filter                 = "bits";                  # if two HSPs in a subject overlap, keep the one with the higher ->{BITS}
my $num_filter_overlap_filter_bits = 0;

my $min_source_aligned_length_candidate            = 40;      # minimum length of source covered by candidate, Zou et al 2009 Plant Physiology
                                                              # if i make this be the length of the subject alignment, too many gaps and too many
                                                              # bad alignments, it seems
my $num_filter_min_source_aligned_length_candidate = 0;

my $max_gap_length            = 12338;                        # max gap between two subject hits for them to be joined.  Zou et al use 0.99
                                                              # quantile of intron length, which is reasonable.  For spruce,  from high confidence
                                                              # set of transcripts at http://sprucewiki.scilifelab.se/sprucetranscriptome/node/2803,
                                                              # quantiles of intron length are
                                                              # 0.99 = 12338 bp
                                                              # 0.95 = 4805 bp
                                                              # 0.90 = 2584 bp
my $num_filter_max_gap_length = 0;

sub final_report() {
    print STDERR "
Pseudogene pipeline                                         = $this_pipeline
Blast algorithm used to generate HSPs                       = $candidate_algorithm\n";
    print STDERR "
NOTE: I tried to avoid AA query/nt subject assumptions but it hasn't been tested with $candidate_algorithm...\n" if $candidate_algorithm ne "TBLASTN";
    print STDERR "
--

Total number of source/query sequences used in the search   = $num_sources
Total number of hits                                        = $num_hits
Total number of HSPs within hits                            = $num_hsps
Total number of candidates passing filtering                = $num_candidates
Total number of blocks within candidates                    = $num_blocks
Total number of stops '$stopchar' observed in candidates           = $num_stop_codons_in_candidates
Candidates per source                                       = " . sprintf( "%.3f\n", ( $num_candidates / $num_sources ) ) . "

--

Parameters and results for trimming of overlapping HSPs:

trim_units when overlap detected [query, subject]           = " . join( ", ", @trim_units ) . "
query_overlap allowed                                       = $query_overlap
num_query_overlap_trim [upstream block, downstream block]   = " . sum(@num_query_overlap_trim) . " [" . join( ", ", @num_query_overlap_trim ) . "]
subject_overlap allowed                                     = $subject_overlap
num_subject_overlap_trim [upstream block, downstream block] = " . sum(@num_subject_overlap_trim) . " [" . join( ", ", @num_subject_overlap_trim ) . "]

--

Filtering results:

num_filter_min_frac_identical_hsp                   [< $min_frac_identical] = $num_filter_min_frac_identical_hsp
num_filter_single_hit_length                         [< $min_single_hit_length] = $num_filter_single_hit_length
num_filter_overlap_filter_bits                              = $num_filter_overlap_filter_bits
num_filter_min_source_aligned_length_candidate       [< $min_source_aligned_length_candidate] = $num_filter_min_source_aligned_length_candidate
num_filter_max_gap_length                         [> $max_gap_length] = $num_filter_max_gap_length
num_filter_min_frac_covered_candidate               [< $min_frac_covered_candidate] = $num_filter_min_frac_covered_candidate
num_filter_min_covered_frac_identical_candidate     [< $min_covered_frac_identical_candidate] = $num_filter_min_covered_frac_identical_candidate\n";
    my $tot =
      $num_filter_min_frac_identical_hsp +
      $num_filter_single_hit_length +
      $num_filter_overlap_filter_bits +
      $num_filter_min_source_aligned_length_candidate +
      $num_filter_max_gap_length +
      $num_filter_min_frac_covered_candidate +
      $num_filter_min_covered_frac_identical_candidate;
    print STDERR "
Total number of filtering events                            = $tot\n";
}

# ----------------------------

sub candidate_score($)    # higher fractional value in the argument is a better score
{
    my $fraction    = shift;
    my $lower_bound = 1e-6;
    my $score       = ( -10 * ( log( 1 - $fraction + $lower_bound ) / log(10) ) );
    return $score;
}

# ----------------------------

die "single argument required, input blast results" if ( !$ARGV[0] );

my $input = $ARGV[0];

# ----------------------------

#  handle blocks of Ns, some day
#my %Ns;
#my %Ns = read_Ns($ARGV[1] or "t.fa.gaps");
#check_Ns();

# ----------------------------

# extra annotation for parent and children, tied to query name

open( ANNOT_1, "<sprucev01_full-TEConf-only.txt" ) or die "Could not open TE-Conf file.";
my %ANNOT_1 = ();
while (<ANNOT_1>) {
    chomp;
    my ( $key, $value ) = /(.*)\t(.*)/g;
    $ANNOT_1{$key} = $value;
}
close(ANNOT_1);

# ----------------------------

my $blast_homology = "GPAVLIMCFYWHKRQNEDST";
my $blast_mismatch = " +";

my $BED_intervals = "";
if ($produce_bed) {
    $BED_intervals = $input . ".bed";
    open( BED_intervals, ">$BED_intervals" ) or die "couldn't open $BED_intervals: $!" if $BED_intervals;
    print BED_intervals "# BED containing psari intervals\n";
}

my $report = new Bio::SearchIO( -file => "gzip -d -c -f $input |", -format => "blast" );

my $result;    # each query
my $hit;       # each hit of the query
my $tiling;    # built from hits

my $mode = "self";    # "self"
my $sep  = "\t";

print STDOUT "#gff-version 3\n";

#printf "%s\n", join($sep, qw/query qlen subject tlen hsp cov covident/);
#printf "%s\n", join($sep,
#qw/query q_len same_scaffold q_name1 q_name2 q_name3 subject t_len t_ident t_fraclen t_nblocks q_blocks t_blocks/);

sub hsp_qual($) {
    my $hsp = shift;
    my ( $hspl, $hspg, $id, $ql ) = ( $hsp->{HSP_LENGTH}, $hsp->{HSP_GAPS}, $hsp->{IDENTICAL}, $hsp->{QUERY_LENGTH} );
    my ( $frac_identical, $frac_covered ) = ( ( $id / $hspl ), ( $hspl - $hspg ) / $ql );
    my $query_frac_identical = $frac_identical * $frac_covered;
    return ( $frac_identical, $frac_covered, $query_frac_identical, $hsp->{BITS} );
}

sub sort_interval_q_sp {
    $a->[0][0] <=> $b->[0][0] || $a->[0][1] <=> $b->[0][1]    # query +
      || $a->[1][0] <=> $b->[1][0]
      || $a->[1][1] <=> $b->[1][1]                            # subject +
}

sub sort_interval_q_sm {
    $a->[0][0] <=> $b->[0][0] || $a->[0][1] <=> $b->[0][1]    # query +
      || $b->[1][1] <=> $a->[1][1]
      || $b->[1][0] <=> $a->[1][0]                            # subject -
}

sub filter_hsps(@) {
    my @hsps = @_;
    my $pre;

    # --- filters acting on one HSP at at time

    # hit length
    $pre = @hsps;
    @hsps = grep { $_->hsp_length >= $min_single_hit_length } @hsps;
    $num_filter_single_hit_length += ( $pre - @hsps );

    # frac identical
    $pre = @hsps;
    @hsps = grep { ( $_->{IDENTICAL} / $_->hsp_length ) >= $min_frac_identical } @hsps;
    $num_filter_min_frac_identical_hsp += ( $pre - @hsps );

    return @hsps if @hsps <= 1;

    # --- filters acting on more than 1 HSP at a time

    # the subject overlap filter: if two hits overlap in the subject, apply the $overlap_filter
    $pre = @hsps;
    if ( $overlap_filter eq "bits" ) {

        my @hsps_p;    # separate + and - strands
        my @hsps_m;
        foreach (@hsps) {
            if   ( $_->{HIT_FRAME} > 0 ) { push @hsps_p, scalar($_); }
            else                         { push @hsps_m, scalar($_); }
        }

        sub sort_interval_b_sp {    # subject +, use reverse bits and 0 out loser if overlapping
            my $aint = [ $a->start('subject'), $a->end('subject') ];
            my $bint = [ $b->start('subject'), $b->end('subject') ];
            if ( !are_disjoint( $aint, $bint, 0 ) && $a->{BITS} != $b->{BITS} ) {
                if   ( $a->{BITS} < $b->{BITS} ) { $a->{BITS} = 0; }
                else                             { $b->{BITS} = 0; }    # 0 bits if low overlap
                $b->{BITS} <=> $a->{BITS};
            }
            else {
                $$aint[0] <=> $$bint[0] || $$aint[1] <=> $$bint[1];
            }
        }

        sub sort_interval_b_sm {                                        # subject -, use reverse bits and 0 out loser if overlapping
            my $aint = [ $a->start('subject'), $a->end('subject') ];
            my $bint = [ $b->start('subject'), $b->end('subject') ];
            if ( !are_disjoint( $aint, $bint, 0 ) && $a->{BITS} != $b->{BITS} ) {
                if   ( $a->{BITS} < $b->{BITS} ) { $a->{BITS} = 0; }
                else                             { $b->{BITS} = 0; }    # 0 bits if low overlap
                $a->{BITS} <=> $b->{BITS};
            }
            else {
                $$bint[1] <=> $$aint[1] || $$bint[0] <=> $$aint[0];
            }
        }
        @hsps_p = grep { $_->{BITS} } sort sort_interval_b_sp @hsps_p;
        @hsps_m = grep { $_->{BITS} } sort sort_interval_b_sm @hsps_m;
        @hsps = ( @hsps_p, @hsps_m );
        $num_filter_overlap_filter_bits += ( $pre - @hsps );
    }
    elsif ( !$overlap_filter ) {

        # do nothing
    }
    else {
        die("unrecognized overlap_filter: $overlap_filter");
    }
    if ( $verbose > 1 ) {
        foreach ( get_intervals_from_hsps( 'subject', @hsps ) ) {
            printf "*** filter_hsps returning s [%d, %d]\n", @{$_};
        }
    }
    return @hsps;
}

sub are_neighbors_p($$$) {    # intervals $x and $y are neighbors on +
    my ( $x, $y, $z ) = @_;
    return ( ( $$x[1] < ( $$y[0] + $z ) ) ? 1 : 0 );
}

sub are_neighbors_m($$$) {    # intervals $x and $y are neighbors on -
    my ( $x, $y, $z ) = @_;
    return ( ( $$x[0] > ( $$y[1] - $z ) ) ? 1 : 0 );
}
my $are_hsp_neighbors_pp = sub($$) {    # hsps are neighbors on q + s +
    my ( $x, $y ) = @_;
    return ( ( are_neighbors_p( $x->[0], $y->[0], $query_overlap ) && are_neighbors_p( $x->[1], $y->[1], $subject_overlap ) ) ? 1 : 0 );
};
my $are_hsp_neighbors_pm = sub($$) {    # hsps are neighbors on q + s -
    my ( $x, $y ) = @_;
    return ( ( are_neighbors_p( $x->[0], $y->[0], $query_overlap ) && are_neighbors_m( $x->[1], $y->[1], $subject_overlap ) ) ? 1 : 0 );
};

while ( $result = $report->next_result ) {

    ++$num_sources;

    my $query_name             = fix_query_name( $result->query_name );
    my $query_length           = $result->query_length;
    my @query_name_fields      = split( /\./, $query_name );              # scaffold.gene.transcript
    my $query_hit_number       = 0;
    my $query_candidate_number = 0;

    $candidate_algorithm = $result->algorithm;

    # $candidate_source = $query_source . ":" . $result->algorithm . ":" . $this_pipeline;
    $candidate_source = $this_pipeline;                                   # above is too long

    if ( $candidate_algorithm eq "TBLASTN" ) {
        @trim_units = ( 1, 3 );                                           # for query (AA) and subject (nt)
    }
    elsif ( $candidate_algorithm eq "BLASTN" ) {
        @trim_units = ( 1, 1 );                                           # for query (nt) and subject (nt)
        die "likely to mistakes handling overlap with $candidate_algorithm" if $query_overlap;
    }

    print STDERR "\n=== query $query_name  len=$query_length\n" if $verbose;

    while ( $hit = $result->next_hit ) {

        ++$query_hit_number;

        ++$num_hits;

        $print_progress = ( $progress_mod and ( ( $num_hsps + $hit->hsps ) % $progress_mod < $num_hsps % $progress_mod ) ) ? 1 : 0;
        $num_hsps += $hit->hsps;

        die("exceeded hit number") if $query_hit_number > 500000000;    # 20;

        my $subject_name   = fix_subject_name( $hit->name );
        my $subject_length = $hit->length;
        my $same_scaffold  = $query_name_fields[0] eq $subject_name ? 1 : 0;

        print STDERR "\n====== subject $subject_name  len=$subject_length", ( $same_scaffold ? " * SAME_SCAFFOLD" : "" ), "\n" if $verbose;

        # DONE: sort out what i see when an hsp hits the - strand, do i need to munge the
        #       start and end?
        #       start and end are in numerical order, regardless of strand; this is how they
        #       should be in the output GFF too.
        # DONE: reject hsps here below an identity cutoff?
        # DONE: add gap tolerance
        # DONE: consider adding frac_identical?? $identical / $hsp_length;

        my $dump_hsp_intervals = sub {
            printf STDERR "hsp  %s  %d,%d => %d,%d\n", $_->[2], $_->[0][0], $_->[0][1], $_->[1][0], $_->[1][1] if $verbose;
        };
        my $dump_hsp_element = sub {
            printf STDERR "hsp 0:[ %d %d ] 1:[ %d %d] 2:%s 3:%s\n", $_->[0][0], $_->[0][1], $_->[1][0], $_->[1][1], $_->[2], $_->[3] if $verbose;
        };

        # -------------  preliminary HSP filtering: identity, length, overlap

        my @hsps = filter_hsps( $hit->hsps );

        if ( !@hsps ) {
            print STDERR "\n*** NO HSPS PASSED FILTERING\n" if $verbose;
            next;
        }

        # -------------  annotated lists of HSPs separated by subject + and -

        my @hsps_p;
        my @hsps_m;

        foreach (@hsps) {

            # split HSPs into subject + and - and annotate them a bit to help out processing
            my @h = (
                [ $_->start('query'),   $_->end('query') ],      # [ 0 ][ 0 , 1 ]
                [ $_->start('subject'), $_->end('subject') ],    # [ 1 ][ 0 , 1 ]
                ( $_->{HIT_FRAME} > 0 ? "+" : "-" ),             # [ 2 ]  frame
                $_                                               # [ 3 ]  hsp object
            );
            if ( $h[2] eq "+" ) {
                push @hsps_p, \@h;
            }
            elsif ( $h[2] eq "-" ) {
                push @hsps_m, \@h;
            }
            else { die "some silly error"; }
            if ( $verbose > 1 ) {
                $dump_hsp_intervals->($_) for ( \@h );
                printf STDERR "frac_id= %f  frac_cov= %f  query_frac_id= %f  bits= %d\n", hsp_qual($_);
            }
        }

        # -------------  interval-sorted annotated lists of HSPs

        @hsps_p = sort sort_interval_q_sp @hsps_p;
        @hsps_m = sort sort_interval_q_sm @hsps_m;

        if ( $verbose > 1 ) {
            if (@hsps_p) { print STDERR "*** + sorted hsps_p\n"; $dump_hsp_element->($_) foreach @hsps_p; }
            if (@hsps_m) { print STDERR "*** - sorted hsps_m\n"; $dump_hsp_element->($_) foreach @hsps_m; }
        }

        # -------------  build list of candidates

        # -------------  finalize_candidate() takes a set of 1+ blocks and finalizes it with
        #                everything needed for output.  It also does a little filtering
        #                so if it returns an undefined value, that candidate has failed a
        #                filtering step
        #
        #                finalize_candidate() is called by create_candidates()

        my $finalize_candidate = sub {

            # -------------  return candidate created from tiling candidate

            print STDERR "*** finalize_candidate received ", scalar(@_), " blocks\n" if $verbose > 1;
            my @blocks    = @_;
            my $strand    = $blocks[0]->[2];    # candidate strand
            my $candidate = {};                 # hashref
            if ( $verbose > 1 ) {
                foreach my $block (@blocks) {
                    print STDERR "*** finalize_candidate received block: ";
                    $dump_hsp_intervals->($_) for ($block);
                }
            }

            # statistics
            $candidate->{candidate_source} = $candidate_source;
            $candidate->{candidate_type}   = $candidate_type;
            $candidate->{source_name}      = $query_name;
            $candidate->{source_length}    = $query_length;
            $candidate->{source_scaffold}  = $query_name_fields[0];
            $candidate->{source_strand}    = "+";
            $candidate->{source_blocks}    = scalar(@blocks);         # n blocks the source divided into in this tiling
            $candidate->{source_gaps}      = scalar(@blocks) - 1;     # n gaps between source blocks

            $candidate->{candidate_scaffold}        = $subject_name;
            $candidate->{candidate_scaffold_length} = $subject_length;
            $candidate->{candidate_strand}          = $strand;           # all blocks the same
                                                                         #print STDERR "0 block : " . Dumper($blocks[0]);
                                                                         #print STDERR "1 block : " . Dumper($blocks[1]);

            # gaps between chunks of source that hit
            my $source_gaplength = 0;

            # gaps between and within chunks of candidate
            my ( $candidate_ngaps, $candidate_gaplength, $gapchar ) = ( 0, 0, $blocks[0][3]->{GAP_SYMBOL} );

            # trim blocks : if there is an overlap, then trip the overlap off the one with
            # the lower quality hit.  if there is a tie, trim off the downstream one.

            warn("finalize_candidate can't handle subject overlap") if $subject_overlap;
            for ( my ( $i, $j ) = ( 0, 1 ) ; $query_overlap and $j < @blocks ; ++$i, $j = $i + 1 ) {    # will not enter if one block
                my ( $bi, $bj ) = ( $blocks[$i], $blocks[$j] );
                if ( !are_disjoint( $bi->[0], $bj->[0], 0 ) ) {

                    # they overlap, so trim one or the other so they don't
                    die "can't handle - query" if $candidate->{source_strand} eq "-";
                    die "can't handle non-TBLASTN hits" if $bi->[3]->{ALGORITHM} ne "TBLASTN";
                    my $n_overlap = amount_overlap( $bi->[0], $bj->[0] );
                    die "violated overlap" if $n_overlap > $query_overlap;
                    my ( $q_trim, $s_trim ) = map { $_ * $n_overlap } @trim_units;                      # algorithm dependency
                    my ( $bic, $bjc ) = ( substr( $bi->[3]->{HOMOLOGY_SEQ}, -$q_trim, $q_trim ), substr( $bj->[3]->{HOMOLOGY_SEQ}, 0, $q_trim ) );

                    # modifying hsp directly, yuck TODO: rewrite to not do this
                    # don't need to adjust the entire candidate intervals since overlaps can only
                    # be internal to the candidate
                    if ( $bic eq $bjc or $bjc =~ /[$blast_mismatch]/ ) {

                        # keep upstream, trim start of downstream
                        $bj->[3]->{QUERY_SEQ}    = substr( $bj->[3]->{QUERY_SEQ},    $q_trim );
                        $bj->[3]->{HOMOLOGY_SEQ} = substr( $bj->[3]->{HOMOLOGY_SEQ}, $q_trim );
                        $bj->[3]->{HIT_SEQ}      = substr( $bj->[3]->{HIT_SEQ},      $q_trim );
                        $bj->[3]->{QUERY_START} += $q_trim;
                        $bj->[0][0] += $q_trim;
                        $bj->[3]->{IDENTICAL} -= $q_trim if $bjc =~ /[$blast_homology]/;
                        $bj->[3]->{HSP_LENGTH} -= $q_trim;
                        if ( $strand eq "+" ) {
                            $bi->[3]->{HIT_START} += $s_trim;
                            $bj->[1][0] += $s_trim;
                        }
                        elsif ( $strand eq "-" ) {
                            $bi->[3]->{HIT_END} -= $s_trim;
                            $bi->[1][1] -= $s_trim;
                        }
                        else { die "some mistake"; }
                        $bj->[3]->{trimmed} = [ +$q_trim, +$s_trim ];    # urg adding new key to hsp object...
                        ++$num_query_overlap_trim[1];
                    }
                    else {
                        # keep downstream, trim end of upstream
                        $bi->[3]->{QUERY_SEQ}    = substr( $bi->[3]->{QUERY_SEQ},    0, -$q_trim );
                        $bi->[3]->{HOMOLOGY_SEQ} = substr( $bi->[3]->{HOMOLOGY_SEQ}, 0, -$q_trim );
                        $bi->[3]->{HIT_SEQ}      = substr( $bi->[3]->{HIT_SEQ},      0, -$q_trim );
                        $bi->[3]->{QUERY_END}  -= $q_trim;
                        $bi->[0][1]            -= $q_trim;
                        $bi->[3]->{IDENTICAL}  -= $q_trim if $bic =~ /[$blast_homology]/;
                        $bi->[3]->{HSP_LENGTH} -= $q_trim;

                        if ( $strand eq "+" ) {
                            $bi->[3]->{HIT_END} -= $s_trim;
                            $bi->[1][1] -= $s_trim;
                        }
                        elsif ( $strand eq "-" ) {
                            $bi->[3]->{HIT_START} += $s_trim;
                            $bi->[1][0] += $s_trim;
                        }
                        else { die "some mistake"; }
                        $bi->[3]->{trimmed} = [ -$q_trim, -$s_trim ];    # urg adding new key to hsp object...
                        ++$num_query_overlap_trim[0];
                    }
                }

                # gaps between pieces of the query hit
                my $gap;
                if ( ( $gap = $bj->[0][0] - $bi->[0][1] - 1 ) ) {
                    $source_gaplength += $gap;
                    print STDERR "source interblock gap : $gap bp\n" if $verbose > 1;
                }

                # account for gap between blocks if there is one
                if ( $strand eq "+" and ( $gap = $bj->[1][0] - $bi->[1][1] - 1 ) ) {
                    ++$candidate_ngaps;
                    $candidate_gaplength += $gap;
                    print STDERR "candidate interblock gap $strand : $gap\n" if $verbose > 1;
                }
                elsif ( $strand eq "-" and ( $gap = $bi->[1][0] - $bj->[1][1] - 1 ) ) {
                    ++$candidate_ngaps;
                    $candidate_gaplength += $gap;
                    print STDERR "candidate interblock gap $strand : $gap\n" if $verbose > 1;
                }
            }
            if ( $verbose > 1 ) {
                print STDERR "interblock source gaps    : $source_gaplength bp\n";
                print STDERR "interblock candidate gaps : $candidate_gaplength bp in $candidate_ngaps gaps\n";
            }

            # final statistics for each block and for the entire candidate

            my ( $candidate_identical, $candidate_covered ) = ( 0, 0 );
            my $candidate_nstops = 0;
            my @block_stats;
            foreach (@blocks) {

                my $thisb = {};

                $thisb->{block_source}             = $candidate_source;                   # same as the candidate
                $thisb->{block_type}               = $candidate_block_type;               # not the same as the candidate
                $thisb->{block_source_strand}      = $candidate->{source_strand};
                $thisb->{block_source_interval}    = $_->[0];
                $thisb->{block_source_length}      = $_->[0][1] - $_->[0][0] + 1;
                $thisb->{block_candidate_scaffold} = $candidate->{candidate_scaffold};    # same as the candidate
                $thisb->{block_candidate_strand}   = $candidate->{candidate_strand};      # same as the candidate
                $thisb->{block_candidate_phase}    = ".";                                 # TODO: maybe later
                $thisb->{block_candidate_interval} = $_->[1];
                $thisb->{block_candidate_length}   = $_->[1][1] - $_->[1][0] + 1;
                $thisb->{block_candidate_trimmed}  = $_->[3]->{trimmed};                  # query, start; neg if right, pos if left
                my $ngaps = ( $_->[3]->{QUERY_SEQ} =~ s/$gapchar+//g );
                $thisb->{block_gaps}      = $ngaps ? $ngaps : 0;
                $thisb->{block_gaplength} = $_->[3]->{HSP_GAPS};
                $thisb->{block_frac_gap}  = $thisb->{block_gaplength} / $thisb->{block_candidate_length};
                my $nstops = ( $_->[3]->{HIT_SEQ} =~ s/$stopchar+//g );
                $thisb->{block_stops}                  = $nstops ? $nstops : 0;
                $thisb->{block_covered}                = ( $_->[3]->{HSP_LENGTH} - $_->[3]->{HSP_GAPS} );
                $thisb->{block_frac_covered}           = $thisb->{block_covered} / $thisb->{block_source_length};
                $thisb->{block_identical}              = $_->[3]->{IDENTICAL};
                $thisb->{block_covered_frac_identical} = $thisb->{block_identical} / $thisb->{block_covered};
                $thisb->{block_frac_identical}         = $thisb->{block_covered_frac_identical} * $thisb->{block_frac_covered};
                $thisb->{block_score}                  = candidate_score( $thisb->{block_frac_identical} );
                $thisb->{block_covered_score}          = candidate_score( $thisb->{block_covered_frac_identical} );

                $candidate_ngaps     += $thisb->{block_gaps};
                $candidate_gaplength += $thisb->{block_gaplength};
                $candidate_nstops    += $thisb->{block_stops};
                $candidate_identical += $thisb->{block_identical};
                $candidate_covered   += $thisb->{block_covered};

                push @block_stats, scalar($thisb);
            }
            if ( $verbose > 0 ) {
                print STDERR "total source gaps          : $source_gaplength bp\n";
                print STDERR "total candidate gaps       : $candidate_gaplength bp in $candidate_ngaps gaps\n";
                print STDERR "total candidate stops '*'  : $candidate_nstops\n";
            }

            $num_stop_codons_in_candidates += $candidate_nstops;

            # ---- set statistics
            # source statistics

            $candidate->{block_stats} = \@block_stats;

            $candidate->{source_blocks} = scalar(@blocks);        # n blocks the source divided into in this tiling
            $candidate->{source_gaps}   = scalar(@blocks) - 1;    # n gaps between source blocks

            ############################
            $candidate->{source_interval} = [ $blocks[0]->[0][0], $blocks[-1]->[0][1] ];
            $candidate->{source_aligned_length} = ( $candidate->{source_interval}->[1] - $candidate->{source_interval}->[0] + 1 );
            $candidate->{source_gaplength} = $source_gaplength;    # total length of gaps between blocks in the source

            # candidate statistics

            $candidate->{candidate_blocks} = scalar(@blocks);      # may be different from source when we start merging
                 #$candidate->{candidate_gaps} = scalar(@blocks) - 1; # may be different from source when we start merging
            $candidate->{candidate_scaffold}           = $subject_name;
            $candidate->{candidate_scaffold_length}    = $subject_length;
            $candidate->{candidate_on_source_scaffold} = $same_scaffold;
            $candidate->{candidate_strand}             = $strand;           # all blocks the same
            $candidate->{candidate_phase}              = ".";               # TODO: maybe later
                 # across the candidate, what is the start position and what is the end position?
            if ( $strand eq "+" ) {
                $candidate->{candidate_interval} = [ $blocks[0]->[1][0], $blocks[-1]->[1][1] ];
            }
            elsif ( $strand eq "-" ) {
                $candidate->{candidate_interval} = [ $blocks[-1]->[1][0], $blocks[0]->[1][1] ];
            }

            # across the candidate, how long is it from end to end?
            $candidate->{candidate_length} = ( $candidate->{candidate_interval}->[1] - $candidate->{candidate_interval}->[0] + 1 );

            # across the candidate, how many gaps are there, both within blocks and between blocks?
            $candidate->{candidate_gaps} = $candidate_ngaps;

            # across the candidate, how much of the candidate is gaps, both within blocks and between blocks?
            $candidate->{candidate_gaplength} = $candidate_gaplength;

            # how many stops across the blocks?
            $candidate->{candidate_stops} = $candidate_nstops;

            # across the candidate, what fraction of the length is gap?
            $candidate->{candidate_frac_gap} = $candidate->{candidate_gaplength} / $candidate->{candidate_length};

            # across all blocks, how many positions of the source block does each block cover?
            $candidate->{candidate_covered} = $candidate_covered;

            # across the candidate, what fraction of the source is covered by the candidate?
            $candidate->{candidate_frac_covered} = $candidate->{candidate_covered} / $candidate->{source_length};

            # across the candidate, how many positions are identical to the source block?
            $candidate->{candidate_identical} = $candidate_identical;

            # across the candidate, what fraction of covered positions are identical to the source?
            $candidate->{candidate_covered_frac_identical} = $candidate->{candidate_identical} / $candidate->{candidate_covered};

            # across the candidate, what fraction of the *covered* source is identical to the *entire* source?
            $candidate->{candidate_frac_identical} = $candidate->{candidate_covered_frac_identical} * $candidate->{candidate_frac_covered};

            # across the candidate, what fraction of the *entire* source is identical to the *entire* source?
            # Phred-scaled covered frac identical with maximum 60 (if covered frac identical = 1)
            $candidate->{candidate_score}         = candidate_score( $candidate->{candidate_frac_identical} );
            $candidate->{candidate_covered_score} = candidate_score( $candidate->{candidate_covered_frac_identical} );

            # gather stats across the whole subject and candidate

            print STDERR "candidate pre-filtering: ", Dumper($candidate) if $verbose > 1;

            # do some filtering
            if ( $candidate->{source_aligned_length} < $min_source_aligned_length_candidate ) {
                $candidate = undef;
                ++$num_filter_min_source_aligned_length_candidate;
            }
            elsif ( $candidate->{candidate_frac_covered} < $min_frac_covered_candidate ) {
                $candidate = undef;
                ++$num_filter_min_frac_covered_candidate;
            }
            elsif ( $candidate->{candidate_covered_frac_identical} < $min_covered_frac_identical_candidate ) {
                $candidate = undef;
                ++$num_filter_min_covered_frac_identical_candidate;
            }
            else {
                # candidate passed
            }

            print STDERR "CANDIDATE DID NOT PASS FINAL FILTERING\n" if !defined($candidate) and $verbose > 0;

            # return a hashref

            return $candidate;
        };

        # -------------  create_candidates() takes a set of HSPs for a single query-subject
        #                pair and scans the HSPs building up tilings of blocks that could be
        #                pseudogene candidates.  It does a little filtering, by refusing
        #                to add a block to a tiling if it is more than $max_gap_length
        #                away from the last-added block.
        #
        #                It then calls finalize_candidate() to finalize the candidate.
        #                If the value returned is undef, then the candidate failed some
        #                filtering step there.
        #
        #                After this function is completed, the candidates contain all
        #                information needed to create GFF output.

        my $create_candidates = sub {
            my @hsps = @_;
            return if !@hsps;
            my @ans    = ();
            my $strand = $hsps[0]->[2];
            if ( $verbose > 1 ) {
                print STDERR "*** create_candidates: strand = $strand\n";
                print STDERR "tiling with ", scalar(@hsps), " q + s $strand candidates:\n";
                $dump_hsp_intervals->($_) foreach (@hsps);
            }

            # subroutine pointer, strand-correct neighbor determination
            my $are_hsp_neighbors = ( $strand eq "+" ) ? $are_hsp_neighbors_pp : $are_hsp_neighbors_pm;

            for ( my $i = 0 ; $i < @hsps ; ++$i ) {
                next if ( !$hsps[$i] );    # we've already used it
                if ( $verbose > 1 ) {
                    print STDERR "tiling + $strand : starting with $i: ";
                    $dump_hsp_intervals->($_) for ( $hsps[$i] );
                }
                my ( $j, $k ) = ( $i, $i );
                my @itile = ($j);          # our tile list is at least this hsp
                while ( ++$k < @hsps ) {
                    next if ( !$hsps[$k] );    # we've already used it
                    ##### filter: add hsp as a tile only if the separation is <= $max_gap_length
                    if ( $are_hsp_neighbors->( $hsps[$j], $hsps[$k] ) ) {
                        if ( amount_separation( $hsps[$j]->[1], $hsps[$k]->[1] ) <= $max_gap_length ) {

                            # TODO: is the first neighbor i find the best?  it is definitely the closest
                            if ( $verbose > 1 ) {
                                print STDERR "tiling + $strand : adding neighbor $k: ";
                                $dump_hsp_intervals->($_) for ( $hsps[$k] );
                            }
                            push @itile, $k;
                            $j = $k;
                        }
                        else {
                            ++$num_filter_max_gap_length;
                        }
                    }
                }
                if ($verbose) {
                    print STDERR "*** final tiling + + : \n";
                    foreach my $tile (@itile) {
                        $dump_hsp_intervals->($_) for ( $hsps[$tile] );
                    }
                }

                my $cand = $finalize_candidate->( @hsps[@itile] );

                if ($cand) {
                    $cand->{source_candidate_number} = ++$query_candidate_number;
                    $cand->{candidate_name}          = sprintf "%s:$candidate_number_format", $cand->{source_name}, $cand->{source_candidate_number};
                    $cand->{candidate_id}            = $cand->{candidate_name};
                    my $block_num = 0;
                    foreach ( @{ $cand->{block_stats} } ) {
                        ++$block_num;
                        $_->{block_name} = sprintf "%s:$candidate_block_number_format", $cand->{candidate_name}, $block_num;
                        $_->{block_id} = $_->{block_name};
                    }
                    print STDERR "finalized candidate: ", Dumper($cand) if $verbose > 1;
                    push @ans, scalar($cand);
                }
                @hsps[@itile] = 0;    # get rid of the hsps we've swallowed
            }

            return (@ans);
        };

        my @psari_pp = $create_candidates->(@hsps_p);
        my @psari_pm = $create_candidates->(@hsps_m);

        # write the GFF

        # first create the parent
        # then create the children

        # -------------------

        my $write_candidate_gff = sub {
            ++$num_candidates;
            my $cand          = shift;
            my @parent_fields = ();
            push @parent_fields, $cand->{candidate_scaffold};                    # where we found the hit
            push @parent_fields, $cand->{candidate_source};                      # source of the GFF, this script
            push @parent_fields, $cand->{candidate_type};                        # type of data for the hit
            push @parent_fields, $cand->{candidate_interval}->[0];               # candidate start
            push @parent_fields, $cand->{candidate_interval}->[1];               # candidate end
            push @parent_fields, sprintf $scorefmt, $cand->{candidate_score};    # Phred-scaled score of frac identical
            push @parent_fields, $cand->{candidate_strand};                      # must be + or -
            push @parent_fields, $cand->{candidate_phase};                       # TODO: for now, is just "."

            my @parent_attr = ();
            push @parent_attr, sprintf "ID=%s", $cand->{candidate_id};

            # push @parent_attr, sprintf "Name=%s", $cand->{candidate_name};  # redundant with ID=
            push @parent_attr, sprintf "Parent=%s",          $cand->{source_name};
            push @parent_attr, sprintf "parentCov=%d-%d",    $cand->{source_interval}->[0], $cand->{source_interval}->[1];
            push @parent_attr, sprintf "onSourceScaff=%d",   $cand->{candidate_on_source_scaffold};
            push @parent_attr, sprintf "numBlocks=%d",       $cand->{candidate_blocks};
            push @parent_attr, sprintf "fracIdent=$fmt",     $cand->{candidate_frac_identical};
            push @parent_attr, sprintf "fracCov=$fmt",       $cand->{candidate_frac_covered};
            push @parent_attr, sprintf "covFracIdent=$fmt",  $cand->{candidate_covered_frac_identical};
            push @parent_attr, sprintf "covScore=$scorefmt", $cand->{candidate_covered_score};
            push @parent_attr, sprintf "numGaps=%d",         $cand->{candidate_gaps};
            push @parent_attr, sprintf "lenGaps=%d",         $cand->{candidate_gaplength};
            push @parent_attr, sprintf "fracGaps=$fmt",      $cand->{candidate_frac_gap};
            push @parent_attr, sprintf "numStops=%d",        $cand->{candidate_stops};
            push @parent_attr, $ANNOT_1{ $cand->{source_name} };

            push @parent_fields, join( ";", @parent_attr );
            print STDOUT join( $sep, @parent_fields ), "\n";

            print BED_intervals join( "\t", $cand->{candidate_scaffold}, $cand->{candidate_interval}->[0], $cand->{candidate_interval}->[1] ), "\n"
              if $BED_intervals;

            my $parent_id = $cand->{candidate_id};

            foreach ( @{ $cand->{block_stats} } ) {
                ++$num_blocks;
                my @child_fields = ();
                push @child_fields, $_->{block_candidate_scaffold};
                push @child_fields, $_->{block_source};                      # source of the GFF, this script
                push @child_fields, $_->{block_type};                        # type of data for the hit
                push @child_fields, $_->{block_candidate_interval}->[0];     # candidate start
                push @child_fields, $_->{block_candidate_interval}->[1];     # candidate end
                push @child_fields, sprintf $scorefmt, $_->{block_score};    # Phred-scaled score of frac identical
                push @child_fields, $_->{block_candidate_strand};            # must be + or -
                push @child_fields, $_->{block_candidate_phase};             # TODO: for now, is just "."

                my @child_attr = ();
                push @child_attr, sprintf "ID=%s",     $_->{block_id};
                push @child_attr, sprintf "Parent=%s", $parent_id;

                # push @child_attr, sprintf "Name=%s", $_->{block_name};  # redundant with ID=
                push @child_attr, sprintf "parentCov=%d-%d", $_->{block_source_interval}->[0], $_->{block_source_interval}->[1];
                push @child_attr, sprintf "fracIdent=$fmt",  $_->{block_frac_identical};
                push @child_attr, sprintf "fracCov=$fmt",    $_->{block_frac_covered};
                push @child_attr, sprintf "covFracIdent=$fmt",  $_->{block_covered_frac_identical};
                push @child_attr, sprintf "covScore=$scorefmt", $_->{block_covered_score};
                push @child_attr, sprintf "numGaps=%d",         $_->{block_gaps};
                push @child_attr, sprintf "lenGaps=%d",         $_->{block_gaplength};
                push @child_attr, sprintf "fracGaps=$fmt",      $_->{block_frac_gap};
                push @child_attr, sprintf "numStops=%d",        $_->{block_stops};
                die "data inconsistency, couldn't find ANNOT_1 entry for $cand->{source_name}" if !exists $ANNOT_1{ $cand->{source_name} };
                push @child_attr, $ANNOT_1{ $cand->{source_name} };

                push @child_fields, join( ";", @child_attr );
                print STDOUT join( $sep, @child_fields ), "\n";

                print BED_intervals
                  join( "\t", $_->{block_candidate_scaffold}, $_->{block_candidate_interval}->[0], $_->{block_candidate_interval}->[1] ), "\n"
                  if $BED_intervals;

            }
        };

        $write_candidate_gff->($_) foreach (@psari_pp);

        $write_candidate_gff->($_) foreach (@psari_pm);

        if ($print_progress) {
            print STDERR "$this_pipeline: examined $num_hsps HSPs in $num_hits hits from $num_sources sources,"
              . " and produced $num_candidates candidates containing $num_blocks blocks...\n";
        }

    }
}

final_report();

#sub read_Ns($) {
#    my $Ns_file = shift;
#    my @Ns;
#    open (Ns, "<$Ns_file") or die "can't open file of N gaps '$Ns_file': ", $!;
#    while (<Ns>) {
#        my @line = split /\t/;  # last field (with \n) is number of Ns, we ignore it here
#        push @Ns, $line[0];
#        my $Nintvls = [];
#        push @{$Nintvls}, split /-/ foreach (split (/,/, $line[1]));
#        push @Ns, $Nintvls;
#    }
#    return @Ns;
#}
#
#sub check_Ns {
#    foreach my $key (sort keys %Ns) {
#        printf "$key N interval: [%s, %s]\n", @{$_} foreach ($Ns{$key});
#    }
#}

sub dump_methods {
    my $obj = shift;
    print "Instance METHOD is " . Dumper( \%{ ref($obj) . "::" } );
}

# inspired by utility routines in Bio::Search::Tiling::*

sub are_disjoint($$$) {
    my ( $x, $y, $z ) = @_;
    return ( ( $$x[1] < ( $$y[0] + $z ) ) || ( $$y[1] < ( $$x[0] + $z ) ) ) ? 1 : 0;
}

sub amount_overlap($$) {
    my ( $x, $y ) = @_;
    if ( !are_disjoint( $x, $y, 0 ) ) {
        if ( $$x[1] >= $$y[0] ) {
            return ( $$x[1] - $$y[0] + 1 );
        }
        elsif ( $$y[1] >= $$x[0] ) {
            return ( $$y[1] - $$x[0] + 1 );
        }
    }
    else {
        return 0;
    }
}

sub amount_separation($$) {
    my ( $x, $y ) = @_;
    if ( my $overlap = amount_overlap( $x, $y ) ) {
        return -$overlap;
    }
    if ( $$x[1] < $$y[0] ) {
        return ( $$y[0] - $$x[1] );
    }
    elsif ( $$y[1] < $$x[0] ) {
        return ( $$x[0] - $$y[1] );
    }
    else {
        die "unknown situation in amount_separation";
    }
}


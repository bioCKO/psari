#!/bin/bash
set -e
set -x

# TODO: test with bed header, with GFF header
# TODO: describe how the VCF sorting must be done for bedtools

TmpID=${SLURM_JOB_ID:-$$}

Arg1=$1
Input=$1
Err=0
if [ ! -e "$Input" ] ; then
    echo "Input file $Input cannot be found"
    Err=1
fi
Het=${2:-hets_sorted.vcf}
Indel=${3:-indels_sorted.vcf}
if [ ! -e "$Het" -o ! -e "$Indel" ] ; then
    echo "Sorted file of heterozygotes $Het or sorted file of indels $Indel or both cannot be found"
    Err=1
fi
if [ "$Err" != "0" ] ; then
    echo "There were errors, I cannot drink your milkshake"
    exit 1
fi

IsGFF=0
if [ "$Input" == "${Input%.bed}" ] ; then
    if [ "$Input" != "${Input%.gff}" -o "$Input" != "${Input%.gff3}" ] ; then
    #  we received a GFF file
        echo -n "Creating unsorted BED from GFF...  "
        IsGFF=1
        NewBed="${Input%.gff*}.$TmpID.bed"
        cat "$Input" | cut -f1,4,5 > "$NewBed"
        Input="$NewBed"
        echo "Done!"
    elif [ "$Input" != "${Input%.gff.gz}" -o "$Input" != "${Input%.gff3.gz}" ] ; then
        echo -n "We cannot handle a compressed GFF, sorry"
        exit 1
    fi
fi



HetOutput=$Input.$TmpID.het.dat
IndelOutput=$Input.$TmpID.indel.dat
TmpOutput=$Input.$TmpID.$TmpID.dat
Output=$Input.$TmpID.final.txt

BEDSorted=${Input%.bed}.sorted.bed
echo sorting BED into $BEDSorted
sort -k1,1 -k2,2n $Input > $BEDSorted

echo gathering hets...
intersectBed -sorted -header -c -a $BEDSorted -b $Het > $HetOutput
echo gathering indels...
intersectBed -sorted -header -c -a $BEDSorted -b $Indel | cut -f4 > $IndelOutput
echo pasting output...
paste $HetOutput $IndelOutput > $TmpOutput

echo generating annotation...
cat > psari_produce-annotation.awk << '__end1__'
BEGIN{ OFS="\t" }
{
    if (NR <= 1) next;
    len = $3-$2+1;
    n_het = $4; 
    n_indel = $5;
    if (n_het == 0) {
        s1 = "hetPerBP=0";
        s3 = "hetCount=0";
    } else {
        s1 = sprintf("hetPerBP=%0.6f", n_het / len);
        s3 = "hetCount=" n_het;
    }
    if (n_indel == 0) {
        s2 = "indelPerBP=0";
        s4 = "indelCount=0";
    } else {
        s2 = sprintf("indelPerBP=%0.6f", n_indel / len);
        s4 = "indelCount=" n_het;
    }
    annot = s1 ";" s2 ";" s3 ";" s4;
    print $1, $2, $3, annot;
}
__end1__
awk -f psari_produce-annotation.awk < $TmpOutput > $Output

if [ "$IsGFF" != "1" ] ; then
    echo done, annotation in $Output
    exit
fi
Output2=${Arg1%.gff}.$TmpID.merged.gff
echo merging annotation onto $Input, producing merged output to $Output2
cat > psari_join-annotation.pl << '__end2__'
my $ann = $ARGV[0];
open(ANN,"<$ann") or die "couldn't open annotation file $ann";
my $key;
my %ANN = {};
while (<ANN>) {
    next if /^#/;
    chomp;
    @l = split /\t/;
    $key = $l[0] . ":" . $l[1] . ":" . $l[2];  # key describing interval
    $ANN{$key} = $l[3];
}
close(ANN);
my $gff = $ARGV[1];
open(GFF,"<$gff") or die "couldn't open input GFF file $gff";
while (<GFF>) {
    /^#/ && do { print; next; };
    chomp;
    @l = split /\t/;
    $key = $l[0] . ":" . $l[3] . ":" . $l[4];  # key describing interval
    die "line $., key $key not in annotation" if ! exists $ANN{$key};
    $l[8] .= (";" . $ANN{$key});
    print join("\t", @l), "\n";
}
close(GFF);
__end2__
perl psari_join-annotation.pl $Output $Arg1 > $Output2
echo Done!  Hope it worked...


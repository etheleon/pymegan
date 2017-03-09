#!/usr/bin/env perl

use Modern::Perl '2013';
use autodie;
use Getopt::Long;

die "$0 -l <kegg links directory> -o <output file> -x <gi2taxid_input> -t <gi2taxid.refseq>\n" unless $#ARGV == 3;

my ($links_dir, $outputfile, $gi2taxid_vanilla, $gi2taxid_refseq, $nrProtDB);
my (%genesko, %genesncbi, %ncbirefseq);
GetOptions(
            'l|links=s'       	=> \$links_dir,
            'o|output=s' 	    => \$outputfile,
            #'x|taxoninput=s'	=> \$gi2taxid_vanilla,
            #'t|taxoutput=s'	    => \$gi2taxid_refseq,
            'n|nr=s'            => \$nrProtDB
           );

#Unmapped
open my $unmapped, '>', 'unmapped';

#Mapping of WP -> GI
open my $errorlog, '>', 'errorlog';
say $errorlog "gi\tNRrefseq";

#gi <- refseq (NR) gi <- NP <- WP
    my @nr = map { if ($_ =~ /\d$/){$_}else{()} }  (split /\n/, `ls ~/db/nr/*`); #not sure whats the reason for reading each one by one
    giRefseq_NR($_) foreach @nr;    #note this only includes Refseq sequences from Bacteria & Archea;
#gi <- refseq (NP)
    giRefseq_prot($_) foreach split(/\n/, `ls /export2/home/uesu/db/refseq/arch_prot/* | grep -v nonredundant`);
    #note this only includes only Refseq sequences from Bacteria & Archea folder;
#ko <-> gi
    parseGeneGI("$links_dir".'/genes_ncbi-gi.list.gz');
#gi<->taxid (only refseq sequences)
    #giTaxon($gi2taxid_vanilla, $gi2taxid_refseq);
#KO <-> refseq
    linkRefseq2KO("$links_dir/genes_ko.list.gz", $outputfile);

####################################################################################################
sub giRefseq_NR {
my ($inputFile) = @_;
say "Reading ".$inputFile;
open my $input , "<", "$inputFile";
    while(<$input>){
    	if (m/^\>/){
		my @sequences = split /gi\|/;
		shift @sequences;	#removes the >
		#process headnode
    		my ($headgi, $wpref) = (split(/\|/, shift @sequences))[0,2];
    		$ncbirefseq{$headgi} = $wpref;
    		say $errorlog join "\t", $headgi, $wpref;

		#childnodes
		foreach (@sequences)
			{
    		my ($gi, $refstatus, $ref) = (split(/\|/, $_))[0,1,2];
    		$ncbirefseq{$gi} = $wpref if ($refstatus eq 'ref');	#only refseq matches
    		say $errorlog join "\t", $gi, $wpref;			#mapping
		    	}
	    	    }
    }
}

sub giRefseq_prot {
my ($inputFile) = @_;
say "Reading ".$inputFile;
open my $input , "<", "$inputFile";
    while(<$input>){
    	if (m/^\>/){
	my ($gi, $ref) = (split /\|/)[1,3];
	$ncbirefseq{$gi} = $ref if (!exists $ncbirefseq{$gi})
    	}
    }
}

sub parseGeneGI {
    my ($file)  = @_;
    open my $input, "-|","zcat $file";
    while(<$input>) {
	    chomp;
	    my ($kegggene, $gi)  = split(/\t/);
	    $gi =~ s/ncbi-gi://;
	    $genesncbi{$kegggene}= $gi;
    }
    say "Finished storing genes-gi";
}

sub giTaxon {
    my ($file, $outputfile) = @_;
    open my $input, '<', $file;
    open my $output, '>', $outputfile;
    while(<$input>) {
	chomp;
	my ($gi, $taxID) = split /\t/;
	say $output join "\t", 'gi|'.$gi, $taxID if exists $ncbirefseq{$gi};
    }
}

sub linkRefseq2KO {
	my ($file, $outputfile)  = @_;
	open my $input, "-|","zcat $file";
	open my $output, ">", $outputfile;
    	while(<$input>) {
	    chomp;
	    my ($kegggene,$ko) = split /\t/;
	    $ko =~ s/ko\:K//;
	    if (exists $genesncbi{$kegggene} && exists $ncbirefseq{$genesncbi{$kegggene}}) {
		#refseqID		ko	GI
	    	say $output join "\t", $ncbirefseq{$genesncbi{$kegggene}}, $ko, $genesncbi{$kegggene}, $kegggene;		#means that ko-(gene:gene no mapping to gi) but how can u
	    }else{
	    	if(exists $genesncbi{$kegggene}){
	    	say $unmapped "This ".$kegggene." is not associated with a refseq gene (KO:".$ko." is mapped to GI:".$genesncbi{$kegggene}."ie. gi is not mapped to a refseqID)"
	    	}else{
		say $unmapped "This ".$kegggene." was not recorded with a gi ";	#all genes were mapped to a gi
	    	}
	    	}
    	}
}

__END__
genes_ncbi-geneid.list.gz
hsa:1   ncbi-geneid:1
hsa:10  ncbi-geneid:10
hsa:100 ncbi-geneid:100
hsa:1000        ncbi-geneid:1000
hsa:10000       ncbi-geneid:10000

genes_ko.list.gz
hsa:10004       ko:K01301
hsa:100049587   ko:K06549
hsa:10005       ko:K11992
hsa:10007       ko:K02564
hsa:10008       ko:K04897
hsa:10009       ko:K10507
hsa:1001        ko:K06796

final output
YP_114235       08684
YP_115248       08684
YP_114234       08684
YP_115247       08684
YP_114236       08684
YP_115249       08684
YP_001940163    08684
YP_001940162    08684
YP_001940158    08684
YP_004511252    08684

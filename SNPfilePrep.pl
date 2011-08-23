#!/usr/bin/perl -w
# use this to prepare a SNP input file for phyloviz

use Bio::SeqIO;
use strict;
use Getopt::Long; 


my ($infile,$inref,$outfile);
&GetOptions(
	    'in:s'      => \$infile,#input alignment in fasta format
	    'ref:s'     => \$inref,#input reference in fasta format
	    'output:s'   => \$outfile,#output file for phyloviz 
           );

my $ref = Bio::SeqIO->new(-file => "$inref" , '-format' => 'fasta');
my $in  = Bio::SeqIO->new(-file => "$infile" , '-format' => 'fasta');
open (OUT, ">SNP$outfile")|| die "Can't open SNP$outfile\n";
open (OUT2, ">IsolateData$outfile")|| die "Can't open IsolateData$outfile\n";
my (@refbases, %bases, %mismatch, $i, %polysites);
while ( my $refseq = $ref->next_seq() ) {
    @refbases=split(//,uc($refseq->seq));
    
}
my $cnt=0;

while ( my $seq_obj = $in->next_seq() ) {
    my $id=$seq_obj->display_id;
    #print "$id \n".$seq_obj->seq."\n";
    $bases{$id}=uc($seq_obj->seq);
    #print $seq_obj->display_id, "\t", $seq_obj->length, "\n";
    
}

foreach my $key(keys %bases){ 
   my @splitseq=split(//,$bases{$key});
   #print "$key $splitseq[0] $splitseq[1]\n";
   for ($i=0; $i<scalar(@splitseq); $i++){
     if ($splitseq[$i]!~/$refbases[$i]/){
       print "mismatch at $i ref is $refbases[$i] and variant is $splitseq[$i]\n";
       $mismatch{$key}{$i}=$splitseq[$i];
       $polysites{$i}++;
     }else{
       $mismatch{$key}{$i}=$splitseq[$i];
     }
   }
}
# print OUT "Site\t# sequence with different aa\n";
# for ($i=0; $i<scalar(@refbases); $i++){
# #  print "$base\n";
# $site=$i+1;
#   if ($mismatch{$i}){
#     
#     print OUT "$site\t$mismatch{$i}\n";
#   }else{
#     print OUT "$site\t0\n";
#   }
# }
my (@SNPsites, %seenSNPprofile);
my $uniqid=1;
foreach my $polysite (keys %polysites){
  print "$polysite\n";
  push(@SNPsites,$polysite);
}
my @ascending = sort { $a <=> $b } @SNPsites;

print OUT "ID\t";
foreach my $SNP(@ascending){
   print OUT "SNP$refbases[$SNP]".($SNP+1)."\t";
}
print OUT "\n";
print OUT2 "ID\tSeqID\thost\tday\n";
for my $id (keys %mismatch){
   my $SNPprofile="";
#  print OUT "$id\t";
  for my $site (keys %{$mismatch{$id}}){
  #  print "$site=$mismatch{$id}{$site} ";
  }
  foreach my $SNP(@ascending){
    $SNPprofile.="$mismatch{$id}{$SNP}\t";
  }
  if (!$seenSNPprofile{$SNPprofile}){    
    print OUT "Profile$uniqid\t$SNPprofile\n";
    $seenSNPprofile{$SNPprofile}=$uniqid;
    print OUT2 "Profile$uniqid\t$id\t$1\t$2\n" if $id=~/^(\d+)D(\d+)S\d+/;
    $uniqid++;
  }elsif($seenSNPprofile{$SNPprofile}){  
    print OUT2 "Profile$seenSNPprofile{$SNPprofile}\t$id\t$1\t$2\n" if $id=~/^(\d+)D(\d+)S\d+/;
  }
}
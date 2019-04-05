#!/usr/bin/perl

#/***************************************************************************
# *   Copyright (C) 2018-2019 by Ari Loytynoja                              *
# *   ari.loytynoja@gmail.com                                               *
# *                                                                         *
# *   This program is free software; you can redistribute it and/or modify  *
# *   it under the terms of the GNU General Public License as published by  *
# *   the Free Software Foundation; either version 2 of the License, or     *
# *   (at your option) any later version.                                   *
# *                                                                         *
# *   This program is distributed in the hope that it will be useful,       *
# *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
# *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
# *   GNU General Public License for more details.                          *
# *                                                                         *
# *   You should have received a copy of the GNU General Public License     *
# *   along with this program; if not, write to the                         *
# *   Free Software Foundation, Inc.,                                       *
# *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
# ***************************************************************************/

#
# exonerate -m p2g -t region.nuc -q prot.pep --showalignment F --showtargetgff F --useaatla F > exonerate.vul
#
# perl vulgar2gff3.pl exonerate.vul
#
 
use strict;
use strict;
use Getopt::Long;
use Bio::SeqIO;
use List::Util qw(first);

my ($gnp,$prp);
my $genenum;
my %txs;

while(<>) {
    if(/^Command/) {
        ($gnp,$prp) = /.* (\w+.nuc).* (\w+.pep)/;
        $genenum=0;
    }
    if(/^vulgar/) {
        my ($tmp,$prf,$prb,$pre,$prs,$gnf,$gnb,$gne,$gns,$sco,@r) = split / /,$_;
        if($gnb>$gne) {$tmp=$gne;$gne=$gnb;$gnb=$tmp;}
        $prb++; 
        $gnb++;
        
        my ($gnseq,$prseq);        
        my $gnin = Bio::SeqIO->new(-file => $gnp,-format => 'fasta');
        while (my $seq = $gnin->next_seq) {
            if($seq->display_id =~ /$gnf/) { 
                $gnseq = Bio::Seq->new(-seq => $seq->subseq($gnb,$gne)); 
                if($gns=~/-/) {
                    $gnseq = $gnseq->revcom; 
                }
                last; 
            }
        }
 
        my $prin = Bio::SeqIO->new(-file => $prp,-format => 'fasta');
        while (my $seq = $prin->next_seq) {
            if($seq->display_id =~ /$prf/) { $prseq = Bio::Seq->new(-seq => $seq->subseq($prb,$pre)); last; }
        }

        #M      Match
        #C      Codon
        #G      Gap
        #N      Non-equivalenced region
        #5      5' splice site
        #3      3' splice site
        #I      Intron
        #S      Split codon
        #F      Frameshift

        my (@t,$mp,$mg,$nfs);
        
        for(my $i=0;$i<$#r; $i+=3) {
            if($r[$i]=~/[CN53IF]/) {
                if($mg>0) { push @t,"M",$mp,$mg; $mp=0; $mg=0; }
                push @t,$r[$i],$r[$i+1],$r[$i+2];
                ++$nfs if($r[$i]=~/F/);
            }
            elsif($r[$i]=~/G/) {
                if($r[$i+2]>0 && $r[$i+2]%3==0) {
                    $mp+=$r[$i+1];
                    $mg+=$r[$i+2];
                } else {
                    push @t,"G",$r[$i+1],$r[$i+2];
                }
            } elsif($r[$i]=~/[MS]/) {
                $mp+=$r[$i+1];
                $mg+=$r[$i+2];
            }
        }
        if($mg>0) { push @t,"M",$mp,$mg; $mp=0; $mg=0; }
        
        ##

        ++$genenum;
        my $gname=$prf.".".$genenum;
        my $tname=$prf.".".$genenum;

        my $prof = 1;
        my $gnof = 1;
        my $gnln = 0;
        
        my @out;
        my $cds;
        while(@t) {
            my $sg = shift @t; 
            my $pn = shift @t; 
            my $gn = shift @t; 
            #print join " ",$sg,$pn,$gn,"\n";
            if( $sg =~ /M/ ) {
                my $of = 0;
                if($gns =~ /\+/) {
                    if($gnln%3!=0) { $of = 3-$gnln%3; }
                    push @out, join "\t",$gnf,"exnrt","exon",$gnb+$gnof-1,$gnb+$gnof+$gn-2,".",$gns,".","Parent=transcript:$tname\n";
                    push @out, join "\t",$gnf,"exnrt","CDS",$gnb+$gnof-1,$gnb+$gnof+$gn-2,".",$gns,$of,"ID=CDS:$tname;Parent=transcript:$tname\n";
                    $cds .= $gnseq->subseq($gnof,$gnof+$gn-1);
                    $gnln += $gn;
                } else {
                    if($gnln%3!=0) { $of = 3-$gnln%3; }
                    unshift @out, join "\t",$gnf,"exnrt","exon",$gne-$gnof-$gn+2,$gne-$gnof+1,".",$gns,".","Parent=transcript:$tname\n";
                    unshift @out, join "\t",$gnf,"exnrt","CDS",$gne-$gnof-$gn+2,$gne-$gnof+1,".",$gns,$of,"ID=CDS:$tname;Parent=transcript:$tname\n";
                    $cds .= $gnseq->subseq($gnof,$gnof+$gn-1);
                    $gnln += $gn;
                }
            }
            $prof += $pn;
            $gnof += $gn;
        }

        my $errors;
        if($nfs>0) { $errors .= ";frameshifts=$nfs"; }
        my $pep = Bio::Seq->new(-seq => $cds, -alphabet   => "dna" )->translate->seq,"\n";
        if(substr($pep,0,1)!~/M/) {$errors .= ";nonmetstart"; }
        my $nsc = ($pep =~ tr/\*//);
        if($nsc>0) { $errors .= ";stopcodons=$nsc"; }
        
        if($gns =~ /\+/) {
            print join "\t",$gnf,"exnrt","gene",$gnb,$gne,".",$gns,".","ID=gene:$gname$errors\n";
            print join "\t",$gnf,"exnrt","mRNA",$gnb,$gne,".",$gns,".","ID=transcript:$tname;Parent=gene:$gname$errors\n";
        } else {
            print join "\t",$gnf,"exnrt","gene",$gnb,$gne,".",$gns,".","ID=gene:$gname$errors\n";
            print join "\t",$gnf,"exnrt","mRNA",$gnb,$gne,".",$gns,".","ID=transcript:$tname;Parent=gene:$gname$errors\n";
        }
        
        print @out;
    }
}

### gfa to fasta
```
for i in `cat ../00.ccs_reads/zall.143.ID `;do grep "^S" $i.best.gfa.nochlo.gfa | awk -v a=$i '{print ">"$2"#"a"\n"$3}' > $i.best.gfa.nochlo.gfa.fasta ;done
```
### index for blastn
```
for i in `ls | grep "fasta"`;do makeblastdb -in $i -dbtype nucl ;done
```
### obtain all the orf longer than 100aa
```
for i in `ls | grep "fasta$"`;do getorf -sequence $i -outseq $i.pep -minsize 150 -find 1 ;done
for i in `ls | grep "pep$"`;do perl ~/script/getorf2gff.pl $i > $i.getorf.gff;done
```
### gff already check

#use miniprot to obtain the annotation of the reference genes.
#for i in `cat ../00.ccs_reads/zall.143.ID`;do /tmp/global2/wxian/software/miniprot/miniprot -P $i.mito.miniprot $i.best.gfa.nochlo.gfa.fasta /tmp/global2/wxian/00.data/Col-0.M.uniq.pep.fa --gff | grep -v "##PAF" > $i.best.gfa.nochlo.gfa.fasta.miniprot.gff ;done
#extract from the miniprot pep
#for i in `cat ../00.ccs_reads/zall.143.ID`;do gffread -g $i.best.gfa.nochlo.gfa.fasta -x $i.best.gfa.nochlo.gfa.fasta.miniprot.cds.fa $i.best.gfa.nochlo.gfa.fasta.miniprot.gff ;done
#for i in `cat ../00.ccs_reads/zall.143.ID`;do transeq  -table 11 $i.best.gfa.nochlo.gfa.fasta.miniprot.cds.fa $i.best.gfa.nochlo.gfa.fasta.miniprot.pep.fa ;done

#for i in `cat ../00.ccs_reads/zall.143.ID`;do bedtools intersect -a $i.best.gfa.nochlo.gfa.fasta.pep.getorf.gff -b $i.best.gfa.nochlo.gfa.fasta.miniprot.gff -v > $i.best.gfa.nochlo.gfa.fasta.pep.getorf.uniq.gff ;done

#for i in `cat ../00.ccs_reads/zall.143.ID`;do grep "\sgene\s" $i.best.gfa.nochlo.gfa.fasta.pep.getorf.uniq.gff | awk '{print $9}' | sed "s/ID=//" > $i.best.gfa.nochlo.gfa.fasta.pep.getorf.uniq.ID ;done
#for i in `cat ../00.ccs_reads/zall.143.ID`;do seqtk subseq $i.best.gfa.nochlo.gfa.fasta.pep $i.best.gfa.nochlo.gfa.fasta.pep.getorf.uniq.ID > $i.best.gfa.nochlo.gfa.fasta.pep.getorf.uniq.pep ;done
#for i in `cat ../00.ccs_reads/zall.143.ID`;do cat $i.best.gfa.nochlo.gfa.fasta.miniprot.pep.fa $i.best.gfa.nochlo.gfa.fasta.pep.getorf.uniq.pep > $i.final.pep.fasta ;done

#for i in `cat ../00.ccs_reads/zall.143.ID`;do gffread -g $i.best.gfa.nochlo.gfa.fasta -x $i.best.gfa.nochlo.gfa.fasta.pep.getorf.uniq.cds $i.best.gfa.nochlo.gfa.fasta.pep.getorf.uniq.gff ;done

#for i in `cat ../00.ccs_reads/zall.143.ID`;do cat $i.best.gfa.nochlo.gfa.fasta.pep.getorf.uniq.gff $i.best.gfa.nochlo.gfa.fasta.miniprot.gff  > $i.final.gff ;done

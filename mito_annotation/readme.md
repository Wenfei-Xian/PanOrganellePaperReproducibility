### Convert GFA to FASTA
```
for i in `cat ../00.ccs_reads/zall.143.ID `;do grep "^S" $i.best.gfa.nochlo.gfa | awk -v a=$i '{print ">"$2"#"a"\n"$3}' > $i.best.gfa.nochlo.gfa.fasta ;done
```
### Extract ORFs Longer Than 50aa
```
for i in `ls | grep "fasta$"`;do getorf -sequence $i -outseq $i.pep -minsize 150 -find 1 ;done
for i in `ls | grep "pep$"`;do perl ~/script/getorf2gff.pl $i > $i.getorf.gff;done
```
### Annotate Reference Genes Using MiniProt
```
for i in `cat ../00.ccs_reads/zall.143.ID`;do /tmp/global2/wxian/software/miniprot/miniprot -P $i.mito.miniprot $i.best.gfa.nochlo.gfa.fasta /tmp/global2/wxian/00.data/Col-0.M.uniq.pep.fa --gff | grep -v "##PAF" > $i.best.gfa.nochlo.gfa.fasta.miniprot.gff ;done
```
### Extract CDS and Translate to Proteins from MiniProt Results
```
for i in `cat ../00.ccs_reads/zall.143.ID`;do gffread -g $i.best.gfa.nochlo.gfa.fasta -x $i.best.gfa.nochlo.gfa.fasta.miniprot.cds.fa $i.best.gfa.nochlo.gfa.fasta.miniprot.gff ;done
for i in `cat ../00.ccs_reads/zall.143.ID`;do transeq  -table 11 $i.best.gfa.nochlo.gfa.fasta.miniprot.cds.fa $i.best.gfa.nochlo.gfa.fasta.miniprot.pep.fa ;done
```
### Filter Out Known Protein-Coding Annotations from GetORF Results
```
for i in `cat ../00.ccs_reads/zall.143.ID`;do bedtools intersect -a $i.best.gfa.nochlo.gfa.fasta.pep.getorf.gff -b $i.best.gfa.nochlo.gfa.fasta.miniprot.gff -v > $i.best.gfa.nochlo.gfa.fasta.pep.getorf.uniq.gff ;done
for i in `cat ../00.ccs_reads/zall.143.ID`;do grep "\sgene\s" $i.best.gfa.nochlo.gfa.fasta.pep.getorf.uniq.gff | awk '{print $9}' | sed "s/ID=//" > $i.best.gfa.nochlo.gfa.fasta.pep.getorf.uniq.ID ;done
for i in `cat ../00.ccs_reads/zall.143.ID`;do seqtk subseq $i.best.gfa.nochlo.gfa.fasta.pep $i.best.gfa.nochlo.gfa.fasta.pep.getorf.uniq.ID > $i.best.gfa.nochlo.gfa.fasta.pep.getorf.uniq.pep ;done
for i in `cat ../00.ccs_reads/zall.143.ID`;do cat $i.best.gfa.nochlo.gfa.fasta.miniprot.pep.fa $i.best.gfa.nochlo.gfa.fasta.pep.getorf.uniq.pep > $i.final.pep.fasta ;done
for i in `cat ../00.ccs_reads/zall.143.ID`;do gffread -g $i.best.gfa.nochlo.gfa.fasta -x $i.best.gfa.nochlo.gfa.fasta.pep.getorf.uniq.cds $i.best.gfa.nochlo.gfa.fasta.pep.getorf.uniq.gff ;done
```
### Merge Final Annotations
```
for i in `cat ../00.ccs_reads/zall.143.ID`;do cat $i.best.gfa.nochlo.gfa.fasta.pep.getorf.uniq.gff $i.best.gfa.nochlo.gfa.fasta.miniprot.gff  > $i.final.gff ;done
```

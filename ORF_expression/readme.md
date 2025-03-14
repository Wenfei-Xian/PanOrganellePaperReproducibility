
### Alignment with STAR
```
for i in `cat ../01.Ref/ID`;do echo "/tmp/global2/wxian/software/anaconda3/envs/star/bin/STAR --twopassMode Basic --quantMode TranscriptomeSAM GeneCounts --runThreadN 24 --genomeDir ../01.Ref/$i.STAR.index --alignIntronMin 20 --alignIntronMax 50000 --outSAMtype BAM SortedByCoordinate --sjdbOverhang 99 --outSAMattrRGline ID:$i SM:$i PL:ILLUMINA --outFilterMismatchNmax 2 --outSJfilterReads Unique --outSAMmultNmax 1 --outFileNamePrefix $i --outSAMmapqUnique 60 --readFilesCommand gunzip -c --readFilesIn ../00.RNAseq/$i.fastq.gz" > STAR.$i.sh;done
```

### Extraction of Transcripts with RSEM
```
for i in `cat ../01.Ref/ID `;do echo "nohup rsem-prepare-reference --gtf ../01.Ref/$i.merged.exon.gtf ../01.Ref/$i.merged.fasta $i.rsem_reference &" ;done
```

### Expression Quantification with Salmon
```
for i in `cat ../01.Ref/ID`;do salmon quant -t $i.rsem_reference.transcripts.fa -l A -a $i\Aligned.toTranscriptome.out.bam -p 64 -o $i.salmon ;done
```

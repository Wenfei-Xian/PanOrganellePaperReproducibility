### dwnloadthe phenotype data from https://zenodo.org/records/3701176#.XmX9u5NKhhE (SupplementaryDateSet1/Arabidopsis)

```
for i in `ls -l | grep "drwxr" |awk 'NR>2{print $9}'`;do for j in `ls $i`;do echo "Rscript Spearman.MultipleLinearRegression4PCs.plot.r 1001G.mito.type.txt $i/$j/used_phenotype_values.txt 1001g.pca_result.eigenvec" ;done ;done
```

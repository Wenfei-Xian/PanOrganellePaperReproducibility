### Download the bioclimate data
```
you can download the data with this link
https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_2.5m_bio.zip
```
### Get the bioclimate data for each accession by using longtitude and latitude.
```
Rscript long.lat.get.bioclimate.R
#The input are 1001G.txt and 1001G_Eastern-Asia.txt.
#The output are 1001G.plus.bio1-19.txt and 1001G_Eastern-Asia.plus.bio1-19.txt
```
### Correlation between the mitochondrial genome structure and geograph/bioclimate
```
Rscript correlation_geograph_bioclimate.R
#The input are 1001G.plus.bio1-19.txt and 1001G_Eastern-Asia.plus.bio1-19.txt
#The output are pdfs and correlation_results.txt
```

# Genetic microarrays

The most common genetic variations are measured with microchips, which provide genotypes of millions of markers. Most genetic clinical studies use them to compare markers frequencies between populations. These studies are called [Genome Wide Association Studies](https://en.wikipedia.org/wiki/Genome-wide_association_study)

For around 150$, private companies measure around one million markers.

## Catalog of published associations

The [EBI catalog](https://www.ebi.ac.uk/gwas/) archives the associations of markers from published genetic studies. The studies are published in scientific journals and peer-reviewed. Many are reported with replication studies, where results from a first publication are verified in a second study.

Around 15,000 markers associations are included. Most concern clinical trials and several concern appearance traits. Markers are reported with odds ratios (effect sizes) and associations scores. The 20 Mb database can be freely [downloaded](https://www.ebi.ac.uk/gwas/api/search/downloads/full).

## SNPedia

[SNPedia](https://snpedia.com) is a wiki for genetic markers. They produce summaries for single markers, by referencing and comparing published studies. Details about diseases and genotypes frequencies in the general population are available.

Their complete catalog is not downloadable but the website can be queried by batches with the [SNPediaR R package](https://github.com/genometra/SNPediaR/).

# Comparing genomes to published studies

## Input

You can either input a microarray genome file from your computer or download a publicly available one from [OpenSNP](https://opensnp.org/genotypes). The input can be either a zip or a text file. Only 23andMe files were tested.

## Search EBI markers and query SNPedia

Markers published in the EBI catalog with high asscoiations scores are searched in the genome (pvalue_mlog > 10). 

The detected markers are then queried on SNPedia to provide summaries and further information.

## Table report

Markers informations are displayed in a table and ordered by [magnitude](https://www.snpedia.com/index.php/Magnitude), which indicates their importance. Summaries and details from SNPedia are included along with EBI catalog informations, as study details and links.

## Report generation

The report is generated using a downloadable [Docker container](https://hub.docker.com/r/thomaschln/genome).

```
wget https://opensnp.org/data/5582.23andme.4073
file_path=5582.23andme.4073
docker run -v `pwd`:`pwd` -w `pwd` --rm -t thomaschln/genome \
  R -e "rmarkdown::render('SNPdiagnostic.Rmd', params = list(genome = '$file_path'))"
```

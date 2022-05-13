# Ancestry estimation by SNPweights

There are many software tools available to estimate ancestry using genetic data. The method "SNPweights" 
[Chen et al. 2013 Bioinformatics](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3661048/) uses SNP weights 
precomputed from large external reference panels. It presented a unique approach to leverage the rich ancestry 
information that is available in WGS or WES data.

The software package "Ancestry-SNPweights" is a wrapper that integrated varies software tools to ease the usage
of the software 'SNPweights', and enables users to take genetic data in popular formats such as vcf, bam/cram, plink
and Illumina studio report format of genotyping arrays.

## Software dependencies

```
python 3.6
pyvcf (depends on python 3.6)
bcftools
vcfutils.pl from bcftools package
samtools
Picard tools
plink

Notes: if you don't take array data, pyvcf and python3.6 are not required.
```

## Software installation

```
git clone https://github.com/tgen/Ancestry-SNPweights

cd Ancestry-SNPweights/data

unzip snpwt.NA.zip

cat snpwt.AS.gz.0* > snpwt.AS.gz 

gunzip snpwt.AS.gz
```

## Usage examples

#### 1. VCF input in Hg38/GRCh38 for a single sample

```
python ~/compute/github/Ancestry-SNPweights/infer_ancestry_vcf.py \
    -v /path/to/vcf_file \
    -o /output/dirctory
```

#### 2. Bam or cram input in Hg38/GRCh38 for a single sample

```
python ~/compute/github/Ancestry-SNPweights/infer_ancestry_vcf.py \
    -b /path/to/bam_file \
    -o /output/dirctory
```

#### 3. Array genotype report input in GRch37 for a single sample

```
python ~/compute/github/Ancestry-SNPweights/infer_ancestry_vcf.py \
    -v /path/to/Illumina_genome_studio_report_file \
    -o /output/dirctory
```

### 4. Plink input data for a cohort of samples

```
To be implemented.

```

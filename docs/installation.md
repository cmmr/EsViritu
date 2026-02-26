## Installation

### Current Versions

Code: **v1.2.0**

Database: **v3.2.4**


### Stable release via Bioconda (recommended)

*NOTE: 2026-01-27 EsViritu v1.1.6 released and available on bioconda.*

**1)  Create conda environment. `mamba` is preferable to `conda` for environment creation.**

`mamba create -n EsViritu bioconda::esviritu`

**2)  Activate the environment with `conda`**

`conda activate EsViritu`

*you should be able to run help menu:*

`EsViritu -h`

**3)  Download the database (\~400 MB when decompressed). EsViritu v1.0.0 or higher requires DB v3.1.0 or higher!**

[![Database DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17716199.svg)](https://doi.org/10.5281/zenodo.17716199)


`cd` *where you want the database to reside*

`mkdir esviritu_DB && cd esviritu_DB`

Download the tarball of DB `v3.2.4` (most recent version) from Zenodo:

`wget https://zenodo.org/records/17716199/files/esviritu_db_v3.2.4.tar.gz`

Check that the download was successful:

`md5sum esviritu_db_v3.2.4.tar.gz`

should return `24d85c1ec3cbffff12e921d2f39c91b2  esviritu_db_v3.2.4.tar.gz`

Unpack and remove the tarball:

`tar -xvf esviritu_db_v3.2.4.tar.gz`

`rm esviritu_db_v3.2.4.tar.gz`

DB files should be in `v3.2.4/`

**4)  Set the database path (optional but recommended):**

`conda env config vars set ESVIRITU_DB=/path/to/esviritu_DB/v3.2.4`

**5)  (OPTIONAL BUT RECOMMENDED) Install the `R` package `dataui` manually in an R session. Without `dataui` reports won't show genome coverage sparklines.**

`R`

then:

`remotes::install_github("timelyportfolio/dataui")`

### Developmental verision
??? "Detailed Instructions"
    
    1)  Clone repo
    
    `git clone https://github.com/cmmr/EsViritu.git`
    
    2)  Go to `EsViritu` directory.
    
    `cd EsViritu`
    
    3)  use the file `environment/EsViritu.yml` with `mamba create` to generate the environment used with this tool
    
    `mamba env create --file environment/EsViritu.yml`
    
    4)  Activate the environment
    
    `conda activate EsViritu`
    
    5)  Make it command-line executable. From repo directory:
    
    `pip install .`
    
    *Now follow the database set up instructions above*


### Docker
??? "Basic Instructions"
    
    **Please note that, while I WAS able to get this to run using `Docker`/`Docker Desktop` on my Mac, I am not a `Docker` expert, and I may be unable to troubleshoot issues.**
    
    1)  Pull Docker image (v1.1.0 shown below)
    
    `docker pull quay.io/biocontainers/esviritu:1.1.0--pyhdfd78af_0`

    *Notes:* 
    
    * be sure to mount your volumes/directories with the `EsViritu` database as well as those with input read files 
    
    * I believe you can save environmental variables like ESVIRITU_DB in `Docker` containers


## (OPTIONAL) Database for filtering out host reads and spike-ins

You could filter unwanted sequences out upstream of this tool, but this will allow you to do it within `EsViritu` using `minimap2`. The pipeline script will look for a file at `filter_seqs/filter_seqs.fna` which could be any fasta-formatted sequence file you want to use to remove matching reads (e.g. from host or spike-in).

Here are instructions for downloading and formatting the human genome and phiX spike-in (3 GB decompressed).

**NOTE: When analyzing sequences from human tissues processed via hybrid capture virome sequencing, quantification may be more accurate if human reads are NOT removed**

```
cd EsViritu ### or `cd` where you want the filter_seqs to reside
mkdir filter_seqs && cd filter_seqs

## download phiX genome and gunzip
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/Sinsheimervirus_phiX174/latest_assembly_versions/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz
gunzip GCF_000819615.1_ViralProj14015_genomic.fna.gz

## download human genome and gunzip
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz
gunzip GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz

## concatenate files
cat GCF_000819615.1_ViralProj14015_genomic.fna GCF_009914755.1_T2T-CHM13v2.0_genomic.fna > filter_seqs.fna

## optionally delete separate files
rm GCF_000819615.1_ViralProj14015_genomic.fna GCF_009914755.1_T2T-CHM13v2.0_genomic.fna

## set the filter_seqs directory as an environmental variable
conda env config vars set ESVIRITU_FILTER=/path/to/filter_seqs
```

Remember to set `-f True` to run the filtering step.

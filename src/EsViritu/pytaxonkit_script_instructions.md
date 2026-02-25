# Database Download

1. Download taxdump from NCBI
```bash
mkdir taxdump && cd taxdump
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz

tar -xvf taxdump.tar.gz
```

2. Download and format acc2taxid (keep .gz for space considerations)
```bash
wget https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz

## need cols 2 (accession.version) and 3 (taxid)
zcat nucl_gb.accession2taxid.gz | cut -f 2,3 | gzip -c > acc2taxid.tsv.gz
```

3. (OPTIONAL) Delete NCBI tarballs

# Script dependencies

Available from conda/pixi:
- pytaxonkit
- polars
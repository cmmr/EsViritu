## Notable changes from `v0` to `v1`

From `v0.2.3` to `v1.0.0` the tool, while having the same goal, has undergone very extensive code refactoring, logic, and database updates. These include:
- The database is completely redone and increased in size by ~2.5X.

  - Subspecies taxonomy labels are added (not in GenBank records) for SARS-CoV-2, Norovirus, Mpox, Rotavirus.

  - Segmented virus genomes are treated as coherent assemblies by the tool rather than unaffiliated contigs.

  - Genbank records that appeared incomplete based on length and/or number of segments were removed.

  - Genbank records with sequences highly similar to human genomic DNA, vectors, or other common contaminants were removed.

  - Genbank records without complete taxonomy information (e.g. a virus called picornaviridae sp. lacks genus and species labels) were excluded.

- All code (other than the .hmtl report) is rewritten in python, making it more robust and reducing dependencies.

- Read ANI and nucleotide diversity (Pi) calculations are now reported.

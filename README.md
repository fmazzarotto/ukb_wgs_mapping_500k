# ukb_wgs_mapping_500k
This repo contains a map of the genomic coordinates for each WGS VCF block of the 500k WGS UK Biobank release, and the code used to create it.

The R script file "create_wgs_block_map.R" contains a single function which takes a table as input ("n_chunks_chr_size.csv"). The input table is structured in a 1-row-per-chromosome basis, with the content of each column described in the R script (see commented lines at the beginning of the script).

Each block contains 20kb of sequence. The R script creates a table ("WGS_500k_block_map.tsv") which is a map of exact coordinates contained by each WGS block.

**WARNING**: this is a quick and efficient way to derive the WGS block map without actually having to access each chunk file, but rather inferring their size based on the number of chunks on DNAnexus and the chromosome size. The various QC steps performed (e.g. comparing the number of expected chunks with the number of actual chunks) did not highlight any potential exception, as for example in the case of the 200k WGS release. However, this would not be robust to situations in which exceptions even out each other to the average of 20kb (e.g. if a given chunk is 25kb but it's "balanced" by another chunk 15kb in size). For what it's worth - I myself have used this map to retrieve >10k variants and all were in the chunk specified in the map, so I am fairly confident about its reliability.

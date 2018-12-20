## Requirements
python3

## Usage
```shell
cd scripts
python3 gs-process-dgv.py
python3 gs-process-svclassify.py
python3 gs-statistics.py
python3 gs-check-overlap.py
```

## Gold Standard sources
* [DGV](http://dgv.tcag.ca/dgv/app/downloads)
```shell
wget http://dgv.tcag.ca/dgv/docs/GRCh37_hg19_variants_2016-05-15.txt
```
* [SVClassify](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/technical/svclassify_Manuscript/Supplementary_Information/)
```shell
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/technical/svclassify_Manuscript/Supplementary_Information/Personalis_1000_Genomes_deduplicated_deletions.bed
```
Modify paths to GS sources in gs-process-* after downloading
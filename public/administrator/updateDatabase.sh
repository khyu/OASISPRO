#!/bin/bash
# updateDatabase.sh
#
# update the database: download the newest version from TCGA and format the files
#
#

## download the newest version
echo "Downloading files" > milestone.txt
bash ../../../downloaderV2/download.sh

## format the files
echo "Formatting files" > milestone.txt

# DNA methylation
bash ../../../format-DNA-methylation27
bash ../../../format-DNA-methylation450

# RNA-seq
bash ../../../rnaSeq.sh

# Proteomics
php -f parse_proteomics.php

echo "Done" > milestone.txt

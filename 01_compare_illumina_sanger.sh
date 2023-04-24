#!/bin/bash

#conda activate qiime2-2022.8

mkdir -p 01_compare_illumina_sanger

# First, copy and extract reads from illumina run
qiime feature-table filter-samples \
--i-table 00_ILLUMINASEQUENCES/dada2_output_16S_250-220/filtered_table_nomitochlorarcheuk.qza \
--m-metadata-file 00_ILLUMINASEQUENCES/meta/meta_tab.tsv \
--p-where "sampleid='Mls1-16S_S264'" \
--output-dir 01_compare_illumina_sanger/filtered_table

qiime feature-table filter-seqs \
--i-data 00_ILLUMINASEQUENCES/dada2_output_16S_250-220/representative_sequences.qza \
--i-table 01_compare_illumina_sanger/filtered_table/filtered_table.qza \
--output-dir  01_compare_illumina_sanger/filtered_seqs

# Next, copy sanger sequences and remove those that don't exist
mkdir -p 01_compare_illumina_sanger/sanger_seqs

qiime tools import \
--input-path 00_SANGERSEQUENCES/04_extract_projects/all_dna-sequences.fasta \
--output-path 01_compare_illumina_sanger/sanger_seqs/merged_data.qza \
--type FeatureData[Sequence]

# Then, extract reads from both illumina and sanger to be 515/806
qiime feature-classifier extract-reads \
--i-sequences 01_compare_illumina_sanger/sanger_seqs/merged_data.qza \
--p-f-primer GTGYCAGCMGCCGCGGTAA \
--p-r-primer GGACTACNVGGGTWTCTAAT \
--output-dir 01_compare_illumina_sanger/extract_sanger_reads



qiime feature-classifier extract-reads \
--i-sequences 01_compare_illumina_sanger/filtered_seqs/filtered_data.qza \
--p-f-primer GTGYCAGCMGCCGCGGTAA \
--p-r-primer GGACTACNVGGGTWTCTAAT \
--output-dir 01_compare_illumina_sanger/extract_illumina_reads


# Finally, cluster at 100% seq ID and 99% seq ID
qiime feature-classifier blast \
--i-query 01_compare_illumina_sanger/extract_sanger_reads/reads.qza \
--i-reference-reads 01_compare_illumina_sanger/extract_illumina_reads/reads.qza \
--p-perc-identity 1 \
--output-dir 01_compare_illumina_sanger/blast_100

qiime feature-classifier blast \
--i-query 01_compare_illumina_sanger/extract_sanger_reads/reads.qza \
--i-reference-reads 01_compare_illumina_sanger/extract_illumina_reads/reads.qza \
--p-perc-identity 0.99 \
--output-dir 01_compare_illumina_sanger/blast_99

qiime feature-classifier blast \
--i-query 01_compare_illumina_sanger/extract_sanger_reads/reads.qza \
--i-reference-reads 01_compare_illumina_sanger/extract_illumina_reads/reads.qza \
--p-perc-identity 0.97 \
--output-dir 01_compare_illumina_sanger/blast_97

# Reverse, full-length
qiime feature-classifier blast \
--i-query 01_compare_illumina_sanger/filtered_seqs/filtered_data.qza \
--i-reference-reads 01_compare_illumina_sanger/sanger_seqs/merged_data.qza \
--p-perc-identity 0.97 \
--p-maxaccepts 1 \
--output-dir 01_compare_illumina_sanger/blast_97_ill_against_san

# export
qiime tools export \
--input-path 01_compare_illumina_sanger/blast_100/search_results.qza \
--output-path 01_compare_illumina_sanger/blast_100/exported

qiime tools export \
--input-path 01_compare_illumina_sanger/blast_99/search_results.qza \
--output-path 01_compare_illumina_sanger/blast_99/exported

qiime tools export \
--input-path 01_compare_illumina_sanger/blast_97/search_results.qza \
--output-path 01_compare_illumina_sanger/blast_97/exported

qiime tools export \
--input-path 01_compare_illumina_sanger/blast_97_ill_against_san/search_results.qza \
--output-path 01_compare_illumina_sanger/blast_97_ill_against_san/exported

# Try SEPP insertion?? Can't because I don't have a SEPP tree for my query seqs. filter tree to include only target seq and assign by sepp insertion

######## Compare sanger to itself? See if it matches anything 100% seq identity

#qiime feature-classifier blast \
#--i-query 01_compare_illumina_sanger/sanger_seqs/myc_sangerseq.qza \
#--i-reference-reads 01_compare_illumina_sanger/sanger_seqs/myc_sangerseq.qza \
#--p-perc-identity 0.99 \
#--output-dir 01_compare_illumina_sanger/blast_99_sanger_sanger

#qiime tools export \
#--input-path 01_compare_illumina_sanger/blast_99_sanger_sanger/search_results.qza \
#--output-path 01_compare_illumina_sanger/blast_99_sanger_sanger/exported

# SCIE4002_Project_Orthogroups
Supplementary data and code for the project report

## Supplementary Tables S1 to S5. 
Provided in a single Excel file (Supplementary_Tables.xlsx), with one sheet per table:
•	Table S1. Orthogroup results: sequence counts per species, species category, single-copy status and number of duplication events for the 91 orthogroups.
•	Table S2. Isoforms: genes with more than one isoform, with isoform IDs, species and orthogroup.
•	Table S3. MAFFT results: pairwise protein identity, conserved positions and gaps for the 59 orthogroups containing all three species.
•	Table S4. Duplication events: the 108 events, classed as isoform-only or gene-level.
•	Table S5. Gene names and functions: names, WRKY groups and reported functions of the genes in gene-level events.

## Commands and scripts

### 1. Data preparation and OrthoFinder (by Chairunnisa Amanda)

```bash
#!/bin/bash

# List the species names in the file
grep ">" seq.fas | cut -d'|' -f2 | sort -u

# Split the file into one FASTA file per species, keeping only the gene ID in each header
mkdir -p project/WRKY_fastas_final

extract_species () {
    # $1 = species name in the header, $2 = output file
    awk -F'|' -v species="$1" '
        /^>/ {
            keep = ($2 == species)
            if (keep) {
                gene_id = $1
                sub(/^>/, "", gene_id)
                print ">" gene_id
            }
            next
        }
        keep { print }
    ' seq.fas > "$2"
}

extract_species "Oryza_sativa_subsp._indica" project/WRKY_fastas_final/Osativa_indica_WRKY.fasta
extract_species "Sorghum_bicolor"            project/WRKY_fastas_final/Sbicolor_WRKY.fasta
extract_species "Zea_mays"                   project/WRKY_fastas_final/Zmays_WRKY.fasta

# Check the number of sequences per species 
for f in project/WRKY_fastas_final/*.fasta; do
    echo "$f: $(grep -c '>' "$f") sequences"
done

# Check for duplicated identifiers
for f in project/WRKY_fastas_final/*.fasta; do
    grep ">" "$f" | sort | uniq -d
done

# Run OrthoFinder
# 8 threads for sequence searches (-t)
# 8 for parallel analysis steps (-a)
cd project/
orthofinder -f WRKY_fastas_final/ -t 8 -a 8 -o WRKY_results_1

# Inspect the main outputs
RES=$(echo WRKY_results_1/Results_*)
tail -50 "$RES/Log.txt"
ls "$RES/"
cat "$RES/Orthogroups/Orthogroups.GeneCount.tsv"
cat "$RES/Comparative_Genomics_Statistics/Duplications_per_Species_Tree_Node.tsv"
head -20 "$RES/Gene_Duplication_Events/Duplications.tsv"
cat "$RES/Comparative_Genomics_Statistics/Statistics_Overall.tsv"
cat "$RES/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv"
```

### 2. Alignment of orthogroups with MAFFT (by Xin Ho)

```bash
#!/bin/bash

INPUT="/mnt/s-ws/s-29/WRKY_project/analysis/shared_orthogroups"
OUTPUT="/mnt/s-ws/s-29/WRKY_project/analysis/mafft_alignments"

mkdir -p "$OUTPUT"

# Align each orthogroup with MAFFT
for file in "$INPUT"/OG*.fa; do
    [ -f "$file" ] || continue
    name=$(basename "$file" .fa)
    echo "Aligning $name"
    mafft --auto "$file" > "$OUTPUT/${name}_aligned.fa"
done

echo "All MAFFT alignments completed."
```

## 3. Phylogenetic tree with MAFFT and IQ-TREE2 (by Jeremy Veder)

```bash
#!/bin/bash

# Align all protein sequences with MAFFT v7.490
mafft --maxiterate 1000 --globalpair allprotein.protein.fa > allprotein.msa

# Infer the tree with IQ-TREE2: 8 threads, 1,000 SH-like approximate likelihood
# ratio test replicates (--alrt) and 1,000 ultrafast bootstrap replicates (-B)
cmd="iqtree2 \
    -T 8 \
    -s allprotein.msa \
    --alrt 1000 \
    -B 1000"

# Submit as a batch job
sbatch --cpus-per-task=8 --wrap="$cmd"
```

## R Script for Figure Making

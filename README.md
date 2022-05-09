# AFpredictoR
Sometimes sumstats report allele frequency (AF) information but is not clear which allele the AF belongs to. 
`AFpredictoR` predicts which allele the AF belongs to and outputs a new summary statistics file. 


## Required packages
- `data.table`
- `VarianceGamma`


## How to run
```bash
./AFpredictoR.R \
    /path/to/sumstats.txt \
    /path/to/cohort.bim \
    /path/to/cohort.frqx.gz \
    predicted-sumstats
```

### `sumstats.txt`
The `data.table::fread` function is used to read the data.  
Assumes that the sumstats has the following columns:
- SNP: Variant ID in the summary statistics file
- chrom: Chromosome position of the variant
- pos: Basepair position of the variant on the chromosome
- A1: Usually the effect allele. 
- A2: Usually the other allele
- AF: Allele frequency for the unknown allele (A1 or A2)
- N: Sample size of the variant

### `cohort.bim`
The `.bim` file from the Plink binary genotype dataset (See [link](https://www.cog-genomics.org/plink/1.9/input#bed)).

### `cohort.frqx.gz`
The `.frqx.gz` file created using Plink's `--freqx gz` command (See [link](https://www.cog-genomics.org/plink/1.9/basic_stats#freq))

### `/path/to/<prefix>`
The output prefix. `<prefix>.csv` and `<prefix>.debug.csv` file will be generated.

## Output
### `<prefix>.csv`
Same format as the input summary statistics file, except that the `AF` column now refers to the frequency for `A1`
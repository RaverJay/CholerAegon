# CholerAegon
CholerAegon is a nextflow pipeline for assembly and antimicrobial resistance gene detection with long and/or short read data.
It will assemble your data with a suitable method, run [RGI](https://github.com/arpcard/rgi) and [Abricate](https://github.com/tseemann/abricate) to detect AMR genes, and then try to predict the resistances conferred by the found genes.

## Quick start

To run CholerAegon, you need [Nextflow](https://www.nextflow.io/) and [Singularity](https://github.com/apptainer/singularity) (now called Apptainer).
If you use [Conda](https://docs.conda.io/en/latest/), you can create an environment for CholerAegon like so:
```
conda create -n choleraegon nextflow singularity
conda activate choleraegon
```

Required parameters:
* pick a release revision of this repository with `-r`, e.g. `-r 0.3.0`
* pick your executor (local, slurm, etc) and engine (singularity) with `-profile`, e.g. `-profile local,singularity`
* supply read data with `--samples`, `--longreads`, `--shortreads` or `--fasta`

Optional parameters:
* supply a genome reference with `--genome_reference` to get average nucleotide identity values for the assemblies
* specify the output folder with `--output` (default `results_CholerAegon`)

## Examples

* run on Nanopore reads (one .fq file per sample)
```
nextflow run RaverJay/CholerAegon -r 0.2.0 -profile local,singularity \
--longreads 'sample_lr_*.fq' --genome_reference my_pathogen.fa --output results_lr
```
* run on short paired-end reads (supply the filename pattern in quotes, with a '*' and '{1,2}')
```
nextflow run RaverJay/CholerAegon -r 0.2.0 -profile local,singularity \
--shortreads 'sample_sr_ID*_{1,2}.fq' --genome_reference my_pathogen.fa --output results_sr
```
* run hybrid assembly with both long and short reads
```
nextflow run RaverJay/CholerAegon -r 0.2.0 -profile local,singularity \
--samples list_of_samples.csv --genome_reference my_pathogen.fa --output results_samples
```
where `list_of_samples.csv` has the following structure:
```
sample1_name,s1_nanopore_reads.fq,s1_short_reads_pair1.fq.gz,s1_short_reads_pair2
sample2_name,s2_nanopore_reads.fq,s2_short_reads_pair1.fq.gz,s2_short_reads_pair2
...
```
## Additional parameters

```
--genome_reference <reference>      optional: supply a genome reference for the analyzed pathogen, will produce %ANI values
--do_all_assemblies                 do not skip short-read-only and long-read-only assemblies in hybrid mode
--fasta <assembly/assemblies>       supply pre-assembled genomes for AMR detection only
--output <folder>                   specify output folder
```

## Pipeline overview

![Pipeline diagram of CholerAegon](https://github.com/RaverJay/CholerAegon/blob/main/figures/pipeline_github.png)


## Further information

CholerAegon makes use of the bioconda docker containers hosted on quai.io.
Docker support coming soon.
Support for cluster/cloud execution coming soon.


## What's with the name?

It loosely emerged from **Cholera** **A**ntibiotic r**E**sistance **G**ene detecti**ON**. But Aegon is also the true name of *A Song of Ice and Fire*'s Jon Snow, meant as a reference to Cholera scientist [John Snow](https://en.wikipedia.org/wiki/John_Snow).

## Citation

If you found CholerAegon helpful for your work or research kindly cite:

Fuesslin, V., Krautwurst, S., Srivastava, A., Winter, D., Liedigk, B., Thye, T., Herrera-León, S., Wohl, S., May, J., Fobil, J.N. and Eibach, D., 2022. **Prediction of Antibiotic Susceptibility Profiles of Vibrio cholerae Isolates From Whole Genome Illumina and Nanopore Sequencing Data: CholerAegon**. *Frontiers in Microbiology*, 13.

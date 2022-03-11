# CholerAegon
CholerAegon is a nextflow pipeline for assembly and antimicrobial resistance gene detection with long and/or short read data.
It will assemble your data with a suitable method, run [RGI](https://github.com/arpcard/rgi) and [Abricate](https://github.com/tseemann/abricate) to detect AMR genes, and then try to predict the resistances conferred by the found genes.

## Quick start

CholerAegon makes use of the bioconda docker containers hosted on quai.io.
We recommend running it using singularity. Docker support coming soon.
Support for cluster/cloud execution coming soon.

Required:
* pick a release revision of this repository with `-r`, e.g. `-r 0.1.0`
* pick your executor and engine with `-profile`, e.g. `-profile local,singularity`
* supply data with `--samples`, `--longreads`, `--shortreads` or `--fasta`

Optional:
* supply a genome reference with `--genome_reference` to get average nucleotide identity values for the assemblies
* specify the output folder with `--output` (default `results_CholerAegon`)

### Examples

* run on Nanopore reads
```
nextflow run RaverJay/CholerAegon -r 0.1.0 -profile local,singularity \
--longreads 'sample_lr_*.fq' --genome_reference my_pathogen.fa --output results_lr
```
* run on short reads
```
nextflow run RaverJay/CholerAegon -r 0.1.0 -profile local,singularity \
--shortreads 'sample_sr_*.fq' --genome_reference my_pathogen.fa --output results_sr
```
* run hybrid assembly with both long and short reads
```
nextflow run RaverJay/CholerAegon -r 0.1.0 -profile local,singularity \
--samples list_of_samples.csv --genome_reference my_pathogen.fa --output results_samples
```
where `list_of_samples.csv` has the following structure:
```
sample1_name,s1_nanopore_reads.fq,s1_short_reads_pair1.fq.gz,s1_short_reads_pair2
sample2_name,s2_nanopore_reads.fq,s2_short_reads_pair1.fq.gz,s2_short_reads_pair2
...
```
## Full usage

```
--genome_reference <reference>      optional: supply a genome reference for the analyzed pathogen, will produce %ANI values
--do_all_assemblies                 do not skip short-read-only and long-read-only assemblies in hybrid mode
--fasta <assembly/assemblies>       supply pre-assembled genomes for AMR detection only
--output <folder>                   specify output folder
```

etc

## Pipeline overview

![Pipeline diagram of CholerAegon](https://github.com/RaverJay/CholerAegon/blob/main/figures/pipeline_github.png)



## What's with the name?

It loosely emerged from **Cholera** **A**ntibiotic r**E**sistance **G**ene detecti**ON**. But Aegon is also the true name of *A Song of Ice and Fire*'s Jon Snow, meant as a reference to Cholera scientist [John Snow](https://en.wikipedia.org/wiki/John_Snow).

## Citation

Preprint coming soon

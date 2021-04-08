snakemake singularity
================
Sylvain Schmitt
April 8, 2021

Test of the `snakemake`
[tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/setup.html#running-the-tutorial-on-your-local-machine)
on my local machine. **Next step is to combine it with `singularity` as
[here](https://forgemia.inra.fr/umr-1202-biogeco/snakemake_singularity_hpc#test-example-).**

# Setup

## Requirements

  - [x] Python ≥3.5 `sudo apt-get install python3.5`
  - [x] Snakemake ≥5.24.1 `sudo apt install snakemake`
  - [x] BWA 0.7 `sudo apt-get install bwa`
  - [x] SAMtools 1.9 `sudo apt-get install samtools`
  - [x] Pysam 0.15 `pip install pysam`
  - [x] BCFtools 1.9 `sudo apt-get install bcftools`
  - [x] Graphviz 2.42 `pip install Graphviz`
  - [x] Jinja2 2.11 `pip install Jinja2`
  - [x] NetworkX 2.5 `pip install NetworkX`
  - [x] Matplotlib 3.3 `pip install Matplotlib`

## Installing Mambaforge

``` bash
cd ~/Téléchargements
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
bash Mambaforge-Linux-x86_64.sh
```

## Preparing data

``` bash
wget https://github.com/snakemake/snakemake-tutorial-data/archive/v5.24.1.tar.gz
tar --wildcards -xf v5.24.1.tar.gz --strip 1 "*/data" "*/environment.yaml"
```

Vizualise obtained data.

``` bash
tree data
```

    [01;34mdata[00m
    ├── genome.fa
    ├── genome.fa.amb
    ├── genome.fa.ann
    ├── genome.fa.bwt
    ├── genome.fa.fai
    ├── genome.fa.pac
    ├── genome.fa.sa
    └── [01;34msamples[00m
        ├── A.fastq
        ├── B.fastq
        └── C.fastq
    
    1 directory, 10 files

Visualize environment file.

``` bash
less environment.yaml
```

    channels:
      - bioconda
      - conda-forge
    dependencies:
      - snakemake-minimal >=5.24.1
      - jinja2 =2.11
      - networkx =2.5
      - matplotlib =3.3
      - graphviz =2.42
      - bcftools =1.9
      - samtools =1.9
      - bwa =0.7
      - pysam =0.15

## Creating an environment with the required software

**Slow with bad internet connection\!**

``` bash
conda activate base
mamba env create --name snakemake-tutorial --file environment.yaml
```

## Activating the environment

``` bash
conda activate snakemake-tutorial
snakemake --help
```

To stop it (not now):

``` bash
conda deactivate
```

# Basics

## Mapping reads

`Snakefile`:

``` python
rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/A.fastq"
    output:
        "mapped_reads/A.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"
```

dry run:

``` bash
snakemake -np mapped_reads/A.bam
```

test:

``` bash
snakemake --cores 1 mapped_reads/A.bam
```

## Generalizing the read mapping rule

``` python
rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/{sample}.fastq"
    output:
        "mapped_reads/{sample}.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"
```

dry run:

``` bash
snakemake -np mapped_reads/B.bam
snakemake -np mapped_reads/A.bam mapped_reads/B.bam
snakemake -np mapped_reads/{A,B}.bam
```

`snakemake` is not producing `A.bam` as already done. Test the
follwongin and rerun dry run:

``` bash
touch data/samples/A.fastq
```

## Sorting read alignments

Add to `Snakefile`:

``` python
rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"
```

dry run:

``` bash
snakemake -np sorted_reads/B.bam
```

## Indexing read alignments and visualizing the DAG of jobs

Add to `Snakefile`:

``` python
rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"
```

see dag:

``` bash
snakemake --dag sorted_reads/{A,B}.bam.bai | dot -Tsvg > dag/dag1.svg
```

``` r
knitr::include_graphics("dag/dag1.svg")
```

![](dag/dag1.svg)<!-- -->

*Dashed frames are already done.*

## Calling genomic variants

Add to top of `Snakefile`:

``` python
SAMPLES = ["A", "B"]
```

Add to `Snakefile`:

``` python
rule bcftools_call:
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES)
    output:
        "calls/all.vcf"
    shell:
        "samtools mpileup -g -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"
```

``` bash
snakemake --dag calls/all.vcf | dot -Tsvg > dag/dag2.svg
```

``` r
knitr::include_graphics("dag/dag2.svg")
```

![](dag/dag2.svg)<!-- -->

## Using custom scripts

Add to `Snakefile`:

``` python
rule plot_quals:
    input:
        "calls/all.vcf"
    output:
        "plots/quals.svg"
    script:
        "scripts/plot-quals.py"
```

Write `scripts/plots-qual.py`:

``` python
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pysam import VariantFile

quals = [record.qual for record in VariantFile(snakemake.input[0])]
plt.hist(quals)

plt.savefig(snakemake.output[0])
```

It is also possible with R using `snakemake@input[[1]]`. See this
\[page\]’<https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#snakefiles-external-scripts>)
for automatic Rmd report.

## Adding a target rule

Add to top of `Snakefile`:

``` python
rule all:
    input:
        "plots/quals.svg"
```

dry run:

``` bash
snakemake -n
```

dag:

``` bash
snakemake --dag | dot -Tsvg > dag/dagAll.svg
```

``` r
knitr::include_graphics("dag/dagAll.svg")
```

![](dag/dagAll.svg)<!-- -->

run:

``` bash
snakemake --cores 4
```

``` r
knitr::include_graphics("plots/quals.svg")
```

![](plots/quals.svg)<!-- -->

# Advanced: Decorating the example workflow

## Specifying the number of used threads

Change `Snakefile`:

``` python
rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/{sample}.fastq"
    output:
        "mapped_reads/{sample}.bam"
    threads: 8
    shell:
        "bwa mem -t {threads} {input} | samtools view -Sb - > {output}"
```

## Config files

Add to top of `Snakefile`:

``` python
configfile: "config.yaml"
```

Write `config.yml`:

``` python
samples:
    A: data/samples/A.fastq
    B: data/samples/B.fastq
```

Change `Snakefile`:

``` python
rule bcftools_call:
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=config["samples"]),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=config["samples"])
    output:
        "calls/all.vcf"
    shell:
        "samtools mpileup -g -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"
```

## Input functions

Change `Snakefile`:

``` python
def get_bwa_map_input_fastqs(wildcards):
    return config["samples"][wildcards.sample]
    
rule bwa_map:
    input:
        "data/genome.fa",
        get_bwa_map_input_fastqs
    output:
        "mapped_reads/{sample}.bam"
    threads: 8
    shell:
        "bwa mem -t {threads} {input} | samtools view -Sb - > {output}"
```

## Rule parameters

Change `Snakefile`:

``` python
rule bwa_map:
    input:
        "data/genome.fa",
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        "mapped_reads/{sample}.bam"
    params:
        rg=r"@RG\tID:{sample}\tSM:{sample}"
    threads: 8
    shell:
        "bwa mem -R '{params.rg}' -t {threads} {input} | samtools view -Sb - > {output}"
```

## Logging

Change `Snakefile`:

``` python
rule bwa_map:
    input:
        "data/genome.fa",
        get_bwa_map_input_fastqs
    output:
        "mapped_reads/{sample}.bam"
    params:
        rg=r"@RG\tID:{sample}\tSM:{sample}"
    log:
        "logs/bwa_mem/{sample}.log"
    threads: 8
    shell:
        "(bwa mem -R '{params.rg}' -t {threads} {input} | "
        "samtools view -Sb - > {output}) 2> {log}"
```

## Temporary and protected files

Change `Snakefile`:

``` python
rule bwa_map:
    input:
        "data/genome.fa",
        get_bwa_map_input_fastqs
    output:
        temp("mapped_reads/{sample}.bam")
    params:
        rg=r"@RG\tID:{sample}\tSM:{sample}"
    log:
        "logs/bwa_mem/{sample}.log"
    threads: 8
    shell:
        "(bwa mem -R '{params.rg}' -t {threads} {input} | "
        "samtools view -Sb - > {output}) 2> {log}"
```

Change `Snakefile`:

``` python
rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        protected("sorted_reads/{sample}.bam")
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"
```

## Benchmarking

Change `Snakefile`:

``` python
rule bwa_map:
    input:
        "data/genome.fa",
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        temp("mapped_reads/{sample}.bam")
    params:
        rg="@RG\tID:{sample}\tSM:{sample}"
    log:
        "logs/bwa_mem/{sample}.log"
    benchmark:
        "benchmarks/{sample}.bwa.benchmark.txt"
    threads: 8
    shell:
        "(bwa mem -R '{params.rg}' -t {threads} {input} | "
        "samtools view -Sb - > {output}) 2> {log}"
```

## Automatic deployment of software dependencies

*Not tested.*

Change `Snakefile`:

``` python
rule bwa_map:
    input:
        "data/genome.fa",
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        temp("mapped_reads/{sample}.bam")
    params:
        rg="@RG\tID:{sample}\tSM:{sample}"
    log:
        "logs/bwa_mem/{sample}.log"
    benchmark:
        "benchmarks/{sample}.bwa.benchmark.txt"
    threads: 8
    shell:
        "(bwa mem -R '{params.rg}' -t {threads} {input} | "
        "samtools view -Sb - > {output}) 2> {log}"
```

Add `envs/samtools.yaml`:

``` python
channels:
  - bioconda
  - conda-forge
dependencies:
  - samtools =1.9
```

run:

``` bash
snakemake --use-conda --cores 1
```

## Tool wrappers

Some tools have wrapper to avoid writing the shell command, *e.g.*:

``` python
rule bwa_mem:
  input:
      ref="data/genome.fa",
      sample=lambda wildcards: config["samples"][wildcards.sample]
  output:
      temp("mapped_reads/{sample}.bam")
  log:
      "logs/bwa_mem/{sample}.log"
  params:
      "-R '@RG\tID:{sample}\tSM:{sample}'"
  threads: 8
  wrapper:
      "0.15.3/bio/bwa/mem"
```

## Cluster

[Newt version for
genotoul](https://forgemia.inra.fr/umr-1202-biogeco/snakemake_singularity_hpc#bash)

## Constraining wildcards

Regular expressions can be used, *e.g.*:
`"sorted_reads/{sample,[A-Za-z0-9]+}.bam"`.

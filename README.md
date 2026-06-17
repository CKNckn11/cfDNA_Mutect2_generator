# cfDNA Mutect2 Pipeline Generator

**cfDNA Mutect2 Pipeline Generator** is a lightweight pipeline generator for tumor-only cfDNA somatic mutation calling. It creates standalone shell scripts for Mutect2 calling, contamination estimation, filtering, ANNOVAR annotation, and final mutation table aggregation from a YAML configuration file.

## Introduction

![cfDNA Mutect2 workflow](./Mutect2.png)

This tool does not run all analysis jobs automatically inside Python. Instead, it generates an independent shell-based project that can be submitted and monitored manually on a local server or computing cluster.

```bash
$ python generate_pipeline.py -h
usage: cfdna-mutect2-generate [-h] [-v] --config CONFIG [--outdir OUTDIR]

cfDNA Mutect2 Pipeline Generator (Version = 0.1.0): Generate tumor-only cfDNA Mutect2 shell scripts from a YAML config.

optional arguments:
  -h, --help       show this help message and exit
  -v, --version    show the version of cfdna-mutect2-generate and exit.
  --config CONFIG  YAML config file.
  --outdir OUTDIR  Override project.outdir from the YAML config.
```

### Main processes

**cfDNA Mutect2 Pipeline Generator** creates one tumor-only cfDNA somatic mutation calling workflow from a YAML configuration file. The generated project is organized into independent shell scripts so each step can be inspected, submitted, rerun, or monitored separately.

**1. Metadata and BAM index preparation**

The generator first builds a sample manifest from either `input.sample_table` or BAM files discovered under `input.bam_dir`. It writes the normalized sample metadata to `00_metadata/` and prepares BAM index checks before downstream analysis.

This stage includes:

- Copying the input sample information into `00_metadata/input_samples.tsv`.
- Creating `00_metadata/sample_manifest.tsv` and `00_metadata/bam_list.tsv`.
- Checking whether each BAM file exists.
- Reusing existing `.bai` files or generating missing BAM indexes with `samtools index`.

**2. Mutect2 tumor-only variant calling**

The main mutation calling module uses GATK Mutect2 for tumor-only cfDNA samples. It calls somatic variants for each sample with a reference genome, panel of normals, germline resource, and allele-frequency prior from the configuration file.

This stage includes:

- Running `GATK Mutect2` for each BAM file.
- Using `--panel-of-normals` to suppress recurrent technical artifacts.
- Using `--germline-resource` and `--af-of-alleles-not-in-resource` for tumor-only filtering support.
- Creating compressed and indexed per-sample somatic VCF files in `02_mutect2/`.
- Running samples in parallel according to `resources.max_jobs`.

**3. Contamination estimation and Mutect2 filtering**

After raw Mutect2 calling, the workflow estimates contamination from common variants and uses the contamination table to filter Mutect2 calls.

This stage includes:

- Running `GATK GetPileupSummaries` with the configured interval list and common variant resource.
- Running `GATK CalculateContamination` for each sample.
- Running `GATK FilterMutectCalls` with the corresponding contamination table.
- Writing filtered, indexed VCF files to `05_filter/`.

**4. ANNOVAR annotation and final aggregation**

The annotation module converts filtered VCF files into ANNOVAR input format, annotates variants with configured databases, and merges per-sample annotation results into a final mutation table.

This stage includes:

- Converting filtered VCF files with `convert2annovar.pl`.
- Annotating each sample with `table_annovar.pl`.
- Using `annovar.genome_build`, `annovar.protocols`, and `annovar.operations` from the YAML configuration.
- Aggregating annotated variants into the file specified by `outputs.aggregate_file`.

---------------------

**cfDNA Mutect2 Pipeline Generator** does not directly execute all analysis tasks inside Python for two main reasons:

- First, it keeps task scheduling simple and transparent. Different local servers and computing clusters use different submission systems, so the generator only creates shell scripts and lets users decide how to submit, monitor, and rerun them.

- Second, it keeps the generated workflow independent from the Python package. After scripts are generated, the analysis can be executed without calling the generator again. Each step is a standalone shell script, which makes it easier to split large jobs, rerun failed steps, and adapt execution to available computing resources.

The generated project includes `00.shell/run_pipeline.sh`, which can run individual modes such as `mutect2`, `pileup`, `contamination`, `filter`, `annovar_convert`, `annovar_table`, and `aggregate`, or run the complete workflow with `all`.

## Installation

Run directly from the source directory:

```bash
cd /400T/ckn/cfDNA_Mutect2_generator
python generate_pipeline.py -h
```

Optional local installation:

```bash
python -m pip install -e .
```

After installation, the command-line entry point is:

```bash
cfdna-mutect2-generate -h
```

## Quick Start

Generate the example bladder cfDNA Mutect2 project:

```bash
cd /400T/ckn/cfDNA_Mutect2_generator
python generate_pipeline.py --config examples/config.bladder.yaml
```

Or after installation:

```bash
cfdna-mutect2-generate --config examples/config.bladder.yaml
```

Run the generated project:

```bash
cd /400T/ckn/database_code/bladder_mutect2_project
bash 00.shell/run_pipeline.sh help
bash 00.shell/run_pipeline.sh all
```

For background execution:

```bash
nohup bash 00.shell/run_pipeline.sh all > loginfo/pipeline.nohup.log 2>&1 &
```

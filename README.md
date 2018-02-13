# MitoSort 
MitoSort is a comprehensive python package designed to rapidly perform mitochondrial (mtDNA) analysis of Next Generation Sequence (NGS) data. 

![mito outline](./src/images/mitoflow.png)

## Dependencies
 MitoSort supports Python 2.7.
Installation requires [samtools](http://samtools.sourceforge.net/), [bwa](http://bio-bwa.sourceforge.net/), [qsub]() and python packages [numpy](http://www.numpy.org/), [pandas](http://pandas.pydata.org/), [matplotlib](https://matplotlib.org/),  [pysam](http://pysam.readthedocs.io/en/latest/api.html), and [viennaRNA](https://www.tbi.univie.ac.at/RNA/).
## Installation
To install:
```bash
git clone https://github.com/sjrowley8/MitoSort
```
Once installed:
```bash	
python setup.py
```
This will check for any necessary python dependencies and install them.  You will need to have samtools, bwa, and qsub installed prior to this, however.
## Config.txt	
Once setup.py has successfully run, you will likely need to edit some additional information into your config.txt.  To run MitoSort on your cluster, you will need to fill out the cluster parameters within this file.

## Usage
Inputs can be in either FASTQ or BAM format from Whole Genome Sequencing (WGS), targeted mtDNA sequencing, or even some exome sequencing protocols that include mtDNA coverage. MitoSort can run on a single machine or optionally submit all samples in parallel on a Linux cluster.
```bash
EXAMPLE USAGE:
python MitoSort.py -s sample_info.txt -o /path/to/output_dir/

OPTIONS:
-s,--sample_info        A tab-separated text file defining sample information (required)
-o,--output             Full path to output directory (required)
-c,--config             Config file (required)
-p,--parallel           Parallel mode. Submits jobs to an SGE cluster using qsub. Default: False
-i,--ipsc               Add iPSC annotations. Default: False
-w,--wgs                Samples are from Whole Genome Sequencing. Default: False
-f,--pedigree           Pedigree file for family analysis.
-t,--somatic            Perform a Somatic analysis with Normal and Tumor samples. Default: False
-h,--help               help message

```

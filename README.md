# BamQVFilter
## Installation
### Option1. Pre-built binary
Pre-built binaries for the earlier v0.1.0 release are available
[here](https://github.com/JMencius/BamQVFilter/releases/tag/0.2.0).

You may have to change the file permissions to execute it with `chmod +x bamqvfilter`.

### Option2. Build from source
```bash
git clone https://github.com/JMencius/BamQVFilter
cd BamQVFilter
cargo build --release
./target/release/bamqvfilter --help
```

## Performance
In internal testing, `BamQVFilter` processed a **119 GB BAM** file in under **25 minutes** using 24 threads on a system equipped with an AMD EPYC 7K62 CPU and a SATA SSD.

Actual performance may vary depending on CPU performance, memory bandwidth, and storage I/O speed.

## Usage
```
Usage: bamqvfilter [OPTIONS] --quality <QUALITY> --input <INPUT> --output <OUTPUT>

Options:
  -q, --quality <QUALITY>
          Sets a minimum Phred average quality score
  -t, --threads <THREADS>
          Worker-pool budget: QV workers plus extra BAM compression workers [default: 4]
      --compression-thread-ratio <RATIO>
          Maximum fraction of --threads assigned to extra BAM compression workers [default: 0.25]
      --batch-records <BATCH_RECORDS>
          Maximum number of BAM records held in one processing batch [default: 4096]
      --batch-mib <BATCH_MIB>
          Approximate maximum allocated BAM record memory per batch, in MiB [default: 64]
      --primary
          Keep only mapped primary alignments [aliases: primary-only]
  -i, --input <INPUT>
          Input BAM filename
  -o, --output <OUTPUT>
          Output BAM filename
  -h, --help
          Print help
  -V, --version
          Print version
```

## Example
A simple example:
```
bamqvfilter -i input.bam -t 24 -q 10 -o output.bam;
```
**If the `input.bam` is sorted then `output.bam` is also sorted.**


## Known limitation
1. Only tested on ONT data, but in theory compatible with other sequencing platforms, such as PacBio sequencing.


## Validation script
This tool is validated by a single-thread Python script in [here](./check_output.py) using `pysam`. The validation script calculates the minimum read QV of a given BAM file and if the output is sorted.
```
# build environment
conda create -n valid-env python=3.7;
conda activate valid-env;
pip install pysam tqdm;

# validate BamQVFilter
python check_output output.bam;
```

## Citation
`BamQVFilter` is a tool developed for my own convenience during my research [LongBow](https://github.com/JMencius/LongBow), so if you want to cite `BamQVFilter`, please cite:


Mencius, J., Chen, W., Zheng, Y. et al. Restoring flowcell type and basecaller configuration from FASTQ files of nanopore sequencing data. Nat Commun 16, 4102 (2025). 

<https://doi.org/10.1038/s41467-025-59378-x>
```
@article{mencius_restoring_2025,
	title = {Restoring flowcell type and basecaller configuration from {FASTQ} files of nanopore sequencing data},
	volume = {16},
	issn = {2041-1723},
	url = {https://doi.org/10.1038/s41467-025-59378-x},
	doi = {10.1038/s41467-025-59378-x},
	abstract = {As nanopore sequencing has been widely adopted, data accumulation has surged, resulting in over 700,000 public datasets. While these data hold immense potential for advancing genomic research, their utility is compromised by the absence of flowcell type and basecaller configuration in about 85\% of the data and associated publications. These parameters are essential for many analysis algorithms, and their misapplication can lead to significant drops in performance. To address this issue, we present LongBow, designed to infer flowcell type and basecaller configuration directly from the base quality value patterns of FASTQ files. LongBow has been tested on 66 in-house basecalled FAST5/POD5 datasets and 1989 public FASTQ datasets, achieving accuracies of 95.33\% and 91.45\%, respectively. We demonstrate its utility by reanalyzing nanopore sequencing data from the COVID-19 Genomics UK (COG-UK) project. The results show that LongBow is essential for reproducing reported genomic variants and, through a LongBow-based analysis pipeline, we discovered substantially more functionally important variants while improving accuracy in lineage assignment. Overall, LongBow is poised to play a critical role in maximizing the utility of public nanopore sequencing data, while significantly enhancing the reproducibility of related research.},
	number = {1},
	journal = {Nature Communications},
	author = {Mencius, Jun and Chen, Wenjun and Zheng, Youqi and An, Tingyi and Yu, Yongguo and Sun, Kun and Feng, Huijuan and Feng, Zhixing},
	month = may,
	year = {2025},
	pages = {4102},
}
```

## Issues & Contributions
If you encounter any problems, bugs, or unexpected results while using MethQC, please open an issue in this repository.

We welcome all forms of contributions — whether it’s reporting bugs, suggesting new features, improving documentation, or submitting pull requests.


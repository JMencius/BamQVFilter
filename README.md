# BamQVFilter
## Installation
Download a ready-to-use binary from the release [here](https://github.com/JMencius/BamQVFilter/releases/tag/0.1.0)


You may have to change the file permissions to execute it with `chmod +x bamqvfilter`.

## Usage
```
Filters BAM files based on read quality values.
Options:
  -q, --quality <QUALITY>  Sets a minimum Phred average quality score
  -t, --threads <THREADS>  Use N parallel threads [default: 4]
  -i, --input <INPUT>      Input filename
  -o, --output <OUTPUT>    Output filename
  -h, --help               Print help
  -V, --version            Print version
```

example:
```
bamqvfilter -i input.bam -t 24 -q 10 -o output.bam;
```

If you encounter `bamqvfilter: command not found` on Ubuntu or Debian, try:
```
./bamqvfilter -i input.bam -t 24 -q 10 -o output.bam;
```


## Known limitation
1. Non-primary reads (example below, without sequence and QV) will also be filtered.
```
07cfd828-eadc-4776-807b-86539c404dc9    272     chr1    5289755 0       18510S11M7D2M1D26M3D13M4D20M8D9M1D22M3D26M5D61M4D10M1I73M1I51M588S      *       0       0       *       *  qs:f:19.4802     du:f:37.6886    ns:i:188443     ts:i:10 mx:i:2  ch:i:2935       st:Z:2023-04-24T12:31:05.240+00:00      rn:i:39063      fn:Z:PAO89685_pass__2264ba8c_afee3a87_14.pod5       sm:f:-741.842   sd:f:0.00795814 sv:Z:pa dx:i:0  RG:Z:afee3a87585a5c58b78955ac2f01d681f6359a75_dna_r10.4.1_e8.2_400bps_sup@v5.0.0        NM:i:74 ms:i:345        AS:i:312   nn:i:0   de:f:0.140299   tp:A:S  cm:i:3  s1:i:45 MD:Z:11^GATGGAT2^A14C7A3^GTA9G3^ATTC15A4^ATTGATGA9^T0G16A4^TGA8G2G7A0A5^GAATG7G0A4G2G12A5A25^ATAC7A14G0G5T0C15A4C1T0G7G6A8A6A8G10A7C2G4A3A4G3       rl:i:1959
```
2. Only tested on ONT data, but in theory compatible with other sequencing platforms, such as PacBio sequencing.


## Validation script
This tool is validated by a single-thread Python script in [here](./min_qv.py) using `pysam`. The validation script calculates the minimum read QV of a given BAM file and outputs the value to Standard output.
```
# build environment
conda create -n valid-env python=3.7;
conda activate valid-env;
pip install pysam;

# validate BamQVFilter
python min_qv.py test.bam;
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


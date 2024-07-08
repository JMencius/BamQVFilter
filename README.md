# BamQVFilter
## Installation
Download a ready-to-use binary from the release
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

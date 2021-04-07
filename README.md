[![install with conda](https://anaconda.org/micknudsen/snpahoy/badges/version.svg)](https://anaconda.org/micknudsen/snpahoy) ![CI](https://github.com/micknudsen/snpahoy/workflows/CI/badge.svg?branch=master) [![Coverage Status](https://coveralls.io/repos/github/micknudsen/snpahoy/badge.svg?branch=master)](https://coveralls.io/github/micknudsen/snpahoy?branch=master)

# SNP Ahoy!

Just a little tool for checking ID SNPs. It works in both germline and somatic modes as described in the sections below. By default, only sites with at least `30X` coverage are considered, and sites with major allele frequency greater than or equal to `95%` are considered homyzygote.

```
$ snpahoy --help
Usage: snpahoy [OPTIONS] COMMAND [ARGS]...

Options:
  --minimum_coverage INTEGER      Only consider SNP positions with a least this
                                  coverage  [default: 30]

  --minimum_base_quality INTEGER  Only count bases with at least this quality
                                  [default: 1]

  --homozygosity_threshold FLOAT  Consider a SNP position homozygote if
                                  frequency of most common allele is this or
                                  higher  [default: 0.95]

  --help                          Show this message and exit.

Commands:
  germline
  somatic
```

## Germline Mode

To run in germline mode, simply provide a BAM file using the `--bam_file` option.

```
$ snpahoy germline --help
Usage: snpahoy germline [OPTIONS]

Options:
  --bed_file PATH          BED file with SNP postions  [required]
  --bam_file PATH          BAM file (must be indexed)  [required]
  --output_json_file PATH  JSON output file  [required]
  --help                   Show this message and exit.
```

The output JSON file contains input information, genotypes at all SNP positions, and a summary. In case a SNP is not genotyped (as for the `chrY` ones in the example below), the empty string is reported as genotype.

```
{
    "input": {
        "files": {
            "bed-file": "snps.bed",
            "bam-file": "germline.bam"
        },
        "settings": {
            "minimum-coverage": 30,
            "homozygosity-threshold": 0.95
        }
    },
    "output": {
        "genotypes": {
            "chr1:4789323": "CC",
            "chr1:4895801": "CC",
            "chr1:7374482": "TT",
            ...
            "chrY:20768865": "",
            "chrY:23164803": ""
        },
        "summary": {
            "snps-total": 1041,
            "snps-genotyped": 1040,
            "heterozygotes-fraction": 0.4404,
            "mean-maf-homozygote-sites": 0.0006,
            "mean-off-genotype-frequency": 0.0006
        }
    }
}
```

## Somatic Mode

To run in somatic mode, provide tumor and germline BAM files using the `--tumor_bam_file` and `--germline_bam_file` options.

```
$ snpahoy somatic --help
Usage: snpahoy somatic [OPTIONS]

Options:
  --bed_file PATH           BED file with SNP postions  [required]
  --tumor_bam_file PATH     Tumor BAM file (must be indexed)  [required]
  --germline_bam_file PATH  Germline BAM file (must be indexed)  [required]
  --output_json_file PATH   JSON output file  [required]
  --help                    Show this message and exit.
```

Output is similar to that in germline mode. Only sites which are genotyping in both tumor and germline are used, and the homozygote sites used in mean MAF calculations are the homozygote sites in the germline sample.

```
{
    "input": {
        "files": {
            "bed-file": "snps.bed",
            "tumor-bam-file": "tumor.bam",
            "germline-bam_file": "germline.bam"
        },
        "settings": {
            "minimum-coverage": 30,
            "homozygosity-threshold": 0.95
        },
        "output": {
            "tumor-genotypes": { ... },
            "germline-genotypes": { ... }
        },
        "summary": {
            "snps-total": 1041,
            "snps-genotyped": 1040,
            "heterozygotes-fraction-tumor": 0.4337,
            "heterozygotes-fraction-germline": 0.4404,
            "mean-maf-homozygote-sites-tumor": 0.0007,
            "mean-maf-homozygote-sites-germline": 0.0006,
            "mean-off-genotype-frequency-tumor": 0.0009,
            "mean-off-genotype-frequency-germline": 0.0006
        }
    }
```

This tool is developed with with the [MSK IMPACT](https://doi.org/10.1016/j.jmoldx.2014.12.006) panel in mind. Suggesed cut-offs for identifying sample swap or contamination are `0.55` for heterozygotes fractions and `0.01` for mean MAFs.

## Installation

The recommended way to install `snpahoy` is by using conda:

```
$ conda install -c micknudsen snpahoy
```
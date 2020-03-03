[![install with conda](https://img.shields.io/badge/install%20with-conda-brightgreen.svg?style=flat)](https://conda.anaconda.org/micknudsen) [![Build Status](https://travis-ci.org/micknudsen/snpahoy.svg?branch=master)](https://travis-ci.org/micknudsen/snpahoy) [![Coverage Status](https://coveralls.io/repos/github/micknudsen/snpahoy/badge.svg?branch=master)](https://coveralls.io/github/micknudsen/snpahoy?branch=master)

# SNP Ahoy!

Just a little tool for checking ID SNPs. It works in both germline and somatic modes as described in the sections below. By default, only sites with at least `30X` coverage are considered, and sites with major allele frequency greater than or equal to `95%` are considered homozygous.

```
$ snpahoy --help
Usage: snpahoy [OPTIONS] COMMAND [ARGS]...

Options:
  --minimum_coverage INTEGER      Only consider SNP positions with a lest this
                                  coverage in both tumor and normal  [default:
                                  30]
  --homozygosity_threshold FLOAT  Consider a SNP position homozygote if
                                  frequency of most common allele is this or
                                  higher  [default: 0.95]
  --help                          Show this message and exit.

Commands:
  germline
  somatic
```

## Somatic Mode

See example output below. The homozygote sites used in mean MAF calculations are the homozygote sites in the normal sample. Suggested cut-offs for the [MSK IMPACT](https://doi.org/10.1016/j.jmoldx.2014.12.006) panel are `0.55` for heterozygotes fractions and `0.01` for mean MAFs.

```
{
    "input": {
        "files": {
            "bed-file": "my_id_snps.bed",
            "normal-bam_file": "my_normal.bam",
            "tumor-bam-file": "my_tumor.bam"
        },
        "settings": {
            "minimum-coverage": 30,
            "homozygosity-threshold": 0.95
        }
    },
    "output": {
        "summary": {
            "snps-total": 1042,
            "snps-genotyped": 1005,
            "heterozygotes-fraction-normal": 0.4239,
            "heterozygotes-fraction-tumor": 0.4201,
            "mean-maf-homozygote-sites-normal": 0.0041,
            "mean-maf-homozygote-sites-tumor": 0.0035
        }
    }
}
```

## Germline Mode

Here be dragons!

## Installation

The recommended way to install `snpahoy` is by using conda:

```
$ conda install -c micknudsen snpahoy
```
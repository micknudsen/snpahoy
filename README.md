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

To run in germline mode, simply provide a BAM/CRAM file using the `--bam_file` option.

```
$ snpahoy germline --help
Usage: snpahoy germline [OPTIONS]

Options:
  --bed_file PATH              BED file with SNP postions  [required]
  --bam_file PATH              BAM file (must be indexed)  [required]
  --reference_fasta_file PATH  Reference FASTA file for CRAM files
  --output_json_file PATH      JSON output file  [required]
  --help                       Show this message and exit.
```

The output JSON file contains input information, genotypes at all SNP positions, and a summary. In case a SNP is not genotyped (as for the `chrY` ones in the example below), the empty string is reported as genotype.

```
{
    "input": {
        "settings": {
            "minimum-coverage": 30,
            "minimum-base-quality": 1,
            "homozygosity-threshold": 0.95
        },
        "files": {
            "bed-file": "snps.bed",
            "bam-file": "germline.bam"
        }
    },
    "output": {
        "details": {
            "chr1:13866176": {
                "depth": 982,
                "counts": {
                    "A": 3,
                    "C": 484,
                    "G": 0,
                    "T": 495
                },
                "alleles": {
                    "major": {
                        "base": "T",
                        "frequency": 0.5041
                    },
                    "minor": {
                        "base": "C",
                        "frequency": 0.4929
                    }
                },
                "genotype": "CT",
                "off_genotype_frequency": 0.0031
            },
            "chr1:22261384": {
                "depth": 323,
                "counts": {
                    "A": 0,
                    "C": 0,
                    "G": 323,
                    "T": 0
                },
                "alleles": {
                    "major": {
                        "base": "G",
                        "frequency": 1.0
                    },
                    "minor": {
                        "base": "",
                        "frequency": 0.0
                    }
                },
                "genotype": "GG",
                "off_genotype_frequency": 0.0
            },

            (...)

            "chrY:9401049": {
                "depth": 0,
                "counts": {
                    "A": 0,
                    "C": 0,
                    "G": 0,
                    "T": 0
                },
                "alleles": {
                    "major": {
                        "base": "",
                        "frequency": 0.0
                    },
                    "minor": {
                        "base": "",
                        "frequency": 0.0
                    }
                },
                "genotype": "",
                "off_genotype_frequency": 0.0
            }
        },
        "summary": {
            "snps": {
                "total": 1041,
                "genotyped": 1016
            },
            "heterozygotes-fraction": 0.4744,
            "mean-maf-homozygote-sites": 0.0022,
            "mean-off-genotype-frequency": 0.0023
        }
    }
}
```

## Somatic Mode

To run in somatic mode, provide tumor and germline BAM/CRAM files using the `--tumor_bam_file` and `--germline_bam_file` options.

```
$ snpahoy somatic --help
Usage: snpahoy somatic [OPTIONS]

Options:
  --bed_file PATH              BED file with SNP postions  [required]
  --tumor_bam_file PATH        Tumor BAM file (must be indexed)  [required]
  --germline_bam_file PATH     Germline BAM file (must be indexed)  [required]
  --reference_fasta_file PATH  Reference FASTA file for CRAM files
  --output_json_file PATH      JSON output file  [required]
  --help                       Show this message and exit.
```

Output is similar to that in germline mode. Only sites which are genotyping in both tumor and germline are used, and the homozygote sites used in mean MAF calculations are the homozygote sites in the germline sample.

```
{
    "input": {
        "settings": {
            "minimum-coverage": 30,
            "minimum-base-quality": 1,
            "homozygosity-threshold": 0.95
        },
        "files": {
            "bed-file": "snps.bed",
            "tumor-bam-file": "tumor.bam",
            "germline-bam-file": "germline.bam"
        },
        "output": {
            "details": {
                "tumor": { ... },
                "germline": { ... }
            }
        },
        "summary": {
            "snps": {
                "total": 1041,
                "genotyped": 1000
            },
            "tumor": {
                "heterozygotes-fraction": 0.474,
                "mean-maf-homozygote-sites": 0.0019,
                "mean-off-genotype-frequency": 0.0019
            },
            "germline": {
                "heterozygotes-fraction": 0.474,
                "mean-maf-homozygote-sites": 0.0022,
                "mean-off-genotype-frequency": 0.0019
            }
        }
    }
}
```

This tool is developed with the [MSK IMPACT](https://doi.org/10.1016/j.jmoldx.2014.12.006) panel in mind. Suggested cut-offs for identifying sample swap or contamination are `0.55` for heterozygotes fractions and `0.01` for mean MAFs.

## Installation

The recommended way to install `snpahoy` is by using conda:

```
$ conda install -c micknudsen snpahoy
```

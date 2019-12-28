[![install with conda](https://img.shields.io/badge/install%20with-conda-brightgreen.svg?style=flat)](https://conda.anaconda.org/micknudsen) [![Build Status](https://travis-ci.org/micknudsen/snpahoy.svg?branch=master)](https://travis-ci.org/micknudsen/snpahoy) [![Coverage Status](https://coveralls.io/repos/github/micknudsen/snpahoy/badge.svg?branch=master)](https://coveralls.io/github/micknudsen/snpahoy?branch=master)

# SNP Ahoy!

Calculates various measures for identifying sample swaps and contamination in paired tumor-normal sequencing.

```
% snpahoy --help
usage: snpahoy [-h] --bed_file BED_FILE --tumor_bam_file TUMOR_BAM_FILE
               --normal_bam_file NORMAL_BAM_FILE --output_json_file
               OUTPUT_JSON_FILE [--minimum_coverage MINIMUM_COVERAGE]
               [--homozygosity_threshold HOMOZYGOSITY_THRESHOLD]

optional arguments:
  -h, --help            show this help message and exit
  --bed_file BED_FILE   BED file with SNP postions
  --tumor_bam_file TUMOR_BAM_FILE
                        Tumor BAM file. Must be indexed.
  --normal_bam_file NORMAL_BAM_FILE
                        Normal BAM file. Must be indexed.
  --output_json_file OUTPUT_JSON_FILE
                        JSON output file
  --minimum_coverage MINIMUM_COVERAGE
                        Only consider SNP positions with a lest this coverage
                        in both tumor and normal (default: 30).
  --homozygosity_threshold HOMOZYGOSITY_THRESHOLD
                        Consider a SNP position homozygote if frequency of
                        most common allele is this or higher (default: 0.95).
```

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

The simplest way to install `snpahoy` is by using conda:

```
$ conda install -c micknudsen snpahoy
```
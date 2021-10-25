import click
import json
import os

from collections import defaultdict
from typing import Dict, List

from pysam import AlignmentFile

from snpahoy.core import Genotyper, SNP
from snpahoy.parsers import get_snps
from snpahoy.parsers import parse_bed_file
from snpahoy.utilities import count_heterozygotes
from snpahoy.utilities import mean_minor_allele_frequency
from snpahoy.utilities import mean_off_genotype_frequency


def get_counts(alignment: AlignmentFile, chromosome: str, position: int, minimum_base_quality: int) -> Dict[str, int]:
    coverage = alignment.count_coverage(contig=chromosome, start=position, stop=position + 1, quality_threshold=minimum_base_quality)
    return {'A': coverage[0][0], 'C': coverage[1][0], 'G': coverage[2][0], 'T': coverage[3][0]}


def get_details(snps: List[SNP]):

    return {
        snp.__str__(): {
            'depth': snp.depth,
            'counts': {base: snp.count(base) for base in list('ACGT')},
            'alleles': {
                'major': {
                    'base': snp.major_allele if snp.count(snp.major_allele) > 0 else '',
                    'frequency': float('%.4f' % snp.major_allele_frequency())
                },
                'minor': {
                    'base': snp.minor_allele if snp.count(snp.minor_allele) > 0 else '',
                    'frequency': float('%.4f' % snp.minor_allele_frequency())
                }
            },
            'genotype': snp.genotype if snp.genotype else '',
            'off_genotype_frequency': float('%.4f' % snp.off_genotype_frequency()) if snp.genotype else 0.0
        }
        for snp in snps}


@click.group()
@click.option('--minimum_coverage', default=30, show_default=True, help='Only consider SNP positions with a least this coverage')
@click.option('--minimum_base_quality', default=1, show_default=True, help='Only count bases with at least this quality')
@click.option('--homozygosity_threshold', default=0.95, show_default=True, help='Consider a SNP position homozygote if frequency of most common allele is this or higher')
@click.pass_context
def client(ctx, minimum_coverage, homozygosity_threshold, minimum_base_quality):

    ctx.obj['genotyper'] = Genotyper(minimum_coverage=minimum_coverage,
                                     homozygosity_threshold=homozygosity_threshold)

    results = {}

    results['input'] = defaultdict(dict)
    results['output'] = defaultdict(dict)

    results['input']['settings'] = {'minimum-coverage': minimum_coverage,
                                    'minimum-base-quality': minimum_base_quality,
                                    'homozygosity-threshold': homozygosity_threshold}

    ctx.obj['results'] = results


@client.command()
@click.option('--bed_file', type=click.Path(), required=True, help='BED file with SNP postions')
@click.option('--tumor_bam_file', type=click.Path(), required=True, help='Tumor BAM file (must be indexed)')
@click.option('--germline_bam_file', type=click.Path(), required=True, help='Germline BAM file (must be indexed)')
@click.option('--output_json_file', type=click.Path(), required=True, help='JSON output file')
@click.pass_context
def somatic(ctx, bed_file, tumor_bam_file, germline_bam_file, output_json_file):

    results = ctx.obj['results']

    results['input']['files'] = {'bed-file': os.path.basename(bed_file)}
    results['input']['files']['tumor-bam-file'] = os.path.basename(tumor_bam_file)
    results['input']['files']['germline-bam-file'] = os.path.basename(germline_bam_file)

    with open(bed_file, 'rt') as f:
        snp_coordinates = parse_bed_file(f.read().splitlines())

    germline_snps = get_snps(coordinates=snp_coordinates,
                             genotyper=ctx.obj['genotyper'],
                             get_counts=lambda chromosome, position: get_counts(alignment=AlignmentFile(germline_bam_file),
                                                                                chromosome=chromosome,
                                                                                position=position,
                                                                                minimum_base_quality=results['input']['settings']['minimum-base-quality']))

    tumor_snps = get_snps(coordinates=snp_coordinates,
                          genotyper=ctx.obj['genotyper'],
                          get_counts=lambda chromosome, position: get_counts(alignment=AlignmentFile(tumor_bam_file),
                                                                             chromosome=chromosome,
                                                                             position=position,
                                                                             minimum_base_quality=results['input']['settings']['minimum-base-quality']))

    results['output']['details'] = {'tumor': get_details(snps=tumor_snps),
                                    'germline': get_details(snps=germline_snps)}

    # Only consider SNPs which are genotyped in both germline and tumor sample.
    genotyped_snp_pairs = []
    for germline_snp, tumor_snp in zip(germline_snps, tumor_snps):
        if germline_snp.genotype and tumor_snp.genotype:
            genotyped_snp_pairs.append({'germline': germline_snp, 'tumor': tumor_snp})

    number_of_heterozygotes_germline = count_heterozygotes(snps=[pair['germline'] for pair in genotyped_snp_pairs])
    number_of_heterozygotes_tumor = count_heterozygotes(snps=[pair['tumor'] for pair in genotyped_snp_pairs])

    # Homozygote positions are those at which the germline sample is homzygote.
    germline_snps_at_homozygote_positions = [pair['germline'] for pair in genotyped_snp_pairs if pair['germline'].is_homozygote()]
    tumor_snps_at_homozygote_positions = [pair['tumor'] for pair in genotyped_snp_pairs if pair['germline'].is_homozygote()]

    results['output']['summary'] = {'snps': {'total': len(snp_coordinates), 'genotyped': len(genotyped_snp_pairs)}}

    if genotyped_snp_pairs:
        results['output']['summary']['tumor'] = {
            'heterozygotes-fraction': float('%.4f' % (number_of_heterozygotes_tumor / len(genotyped_snp_pairs))),
            'mean-maf-homozygote-sites': float('%.4f' % mean_minor_allele_frequency(snps=tumor_snps_at_homozygote_positions)),
            'mean-off-genotype-frequency': float('%.4f' % mean_off_genotype_frequency(snps=[pair['tumor'] for pair in genotyped_snp_pairs]))
        }
        results['output']['summary']['germline'] = {
            'heterozygotes-fraction': float('%.4f' % (number_of_heterozygotes_germline / len(genotyped_snp_pairs))),
            'mean-maf-homozygote-sites': float('%.4f' % mean_minor_allele_frequency(snps=germline_snps_at_homozygote_positions)),
            'mean-off-genotype-frequency': float('%.4f' % mean_off_genotype_frequency(snps=[pair['germline'] for pair in genotyped_snp_pairs]))
        }

    with open(output_json_file, 'w') as json_file_handle:
        json.dump(results, json_file_handle, indent=4)


@client.command()
@click.option('--bed_file', type=click.Path(), required=True, help='BED file with SNP postions')
@click.option('--bam_file', type=click.Path(), required=True, help='BAM file (must be indexed)')
@click.option('--output_json_file', type=click.Path(), required=True, help='JSON output file')
@click.pass_context
def germline(ctx, bed_file, bam_file, output_json_file):

    results = ctx.obj['results']

    results['input']['files'] = {'bed-file': os.path.basename(bed_file)}
    results['input']['files']['bam-file'] = os.path.basename(bam_file)

    with open(bed_file, 'rt') as f:
        snp_coordinates = parse_bed_file(f.read().splitlines())

    snps = get_snps(coordinates=snp_coordinates,
                    genotyper=ctx.obj['genotyper'],
                    get_counts=lambda chromosome, position: get_counts(alignment=AlignmentFile(bam_file),
                                                                       chromosome=chromosome,
                                                                       position=position,
                                                                       minimum_base_quality=results['input']['settings']['minimum-base-quality']))

    genotyped_snps = [snp for snp in snps if snp.genotype]

    results['output']['details'] = get_details(snps=snps)
    results['output']['summary'] = {'snps': {'total': len(snps), 'genotyped': len(genotyped_snps)}}

    number_of_heterozygotes = count_heterozygotes(snps=genotyped_snps)
    homozygote_snps = [snp for snp in snps if snp.is_homozygote()]

    if genotyped_snps:
        results['output']['summary']['heterozygotes-fraction'] = float('%.4f' % (number_of_heterozygotes / len(genotyped_snps)))
        results['output']['summary']['mean-maf-homozygote-sites'] = float('%.4f' % mean_minor_allele_frequency(snps=homozygote_snps))
        results['output']['summary']['mean-off-genotype-frequency'] = float('%.4f' % mean_off_genotype_frequency(snps=homozygote_snps))

    with open(output_json_file, 'w') as json_file_handle:
        json.dump(results, json_file_handle, indent=4)


def run():
    client(obj={})

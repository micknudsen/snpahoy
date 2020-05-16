import click
import json
import os

from collections import defaultdict
from typing import Dict

from pysam import AlignmentFile

from snpahoy.core import Genotyper
from snpahoy.parsers import get_snps
from snpahoy.parsers import parse_bed_file
from snpahoy.utilities import count_heterozygotes
from snpahoy.utilities import mean_minor_allele_frequency
from snpahoy.utilities import mean_off_genotype_frequency


def get_counts(alignment: AlignmentFile, chromosome: str, position: int) -> Dict[str, int]:
    coverage = alignment.count_coverage(contig=chromosome, start=position, stop=position + 1)
    return {'A': coverage[0][0], 'C': coverage[1][0], 'G': coverage[2][0], 'T': coverage[3][0]}


@click.group()
@click.option('--minimum_coverage', default=30, show_default=True, help='Only consider SNP positions with a lest this coverage')
@click.option('--homozygosity_threshold', default=0.95, show_default=True, help='Consider a SNP position homozygote if frequency of most common allele is this or higher')
@click.pass_context
def client(ctx, minimum_coverage, homozygosity_threshold):

    ctx.obj['genotyper'] = Genotyper(minimum_coverage=minimum_coverage,
                                     homozygosity_threshold=homozygosity_threshold)

    results = {}

    results['input'] = defaultdict(dict)
    results['output'] = defaultdict(dict)

    results['input']['settings'] = {'minimum-coverage': minimum_coverage,
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

    results['input']['files']['tumor-bam-file'] = os.path.basename(tumor_bam_file)
    results['input']['files']['germline-bam_file'] = os.path.basename(germline_bam_file)

    with open(bed_file, 'rt') as f:
        snp_coordinates = parse_bed_file(f.read().splitlines())

    germline_snps = get_snps(coordinates=snp_coordinates,
                             genotyper=ctx.obj['genotyper'],
                             get_counts=lambda chromosome, position: get_counts(alignment=AlignmentFile(germline_bam_file), chromosome=chromosome, position=position))

    tumor_snps = get_snps(coordinates=snp_coordinates,
                          genotyper=ctx.obj['genotyper'],
                          get_counts=lambda chromosome, position: get_counts(alignment=AlignmentFile(tumor_bam_file), chromosome=chromosome, position=position))

    tumor_genotypes = {}
    for snp in tumor_snps:
        tumor_genotypes[snp.__str__()] = snp.genotype if snp.genotype else ''
    results['output']['tumor-genotypes'] = tumor_genotypes

    germline_genotypes = {}
    for snp in germline_snps:
        germline_genotypes[snp.__str__()] = snp.genotype if snp.genotype else ''
    results['output']['germline-genotypes'] = germline_genotypes

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

    results['output']['summary'] = {'snps-total': len(snp_coordinates),
                                    'snps-genotyped': len(genotyped_snp_pairs)}

    if genotyped_snp_pairs:
        results['output']['summary']['heterozygotes-fraction-tumor'] = float('%.4f' % (number_of_heterozygotes_tumor / len(genotyped_snp_pairs)))
        results['output']['summary']['heterozygotes-fraction-germline'] = float('%.4f' % (number_of_heterozygotes_germline / len(genotyped_snp_pairs)))
        results['output']['summary']['mean-maf-homozygote-sites-tumor'] = float('%.4f' % mean_minor_allele_frequency(snps=tumor_snps_at_homozygote_positions))
        results['output']['summary']['mean-maf-homozygote-sites-germline'] = float('%.4f' % mean_minor_allele_frequency(snps=germline_snps_at_homozygote_positions))
        results['output']['summary']['mean-off-genotype-frequency-tumor'] = float('%.4f' % mean_off_genotype_frequency(snps=[pair['tumor'] for pair in genotyped_snp_pairs]))
        results['output']['summary']['mean-off-genotype-frequency-germline'] = float('%.4f' % mean_off_genotype_frequency(snps=[pair['germline'] for pair in genotyped_snp_pairs]))

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
                    get_counts=lambda chromosome, position: get_counts(alignment=AlignmentFile(bam_file), chromosome=chromosome, position=position))

    genotypes = {}
    genotyped_snps = []

    for snp in snps:
        if snp.genotype:
            genotypes[snp.__str__()] = snp.genotype
            genotyped_snps.append(snp)
        else:
            genotypes[snp.__str__()] = ''

    results['output']['genotypes'] = genotypes

    results['output']['summary'] = {'snps-total': len(snps),
                                    'snps-genotyped': len(genotyped_snps)}

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

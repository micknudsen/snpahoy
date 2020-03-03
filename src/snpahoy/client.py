import click
import json
import os

from collections import defaultdict
# from statistics import mean
from typing import Dict


from pysam import AlignmentFile

from snpahoy.core import Genotyper
from snpahoy.parsers import get_snps
from snpahoy.parsers import parse_bed_file


def get_counts(alignment: AlignmentFile, chromosome: str, position: int) -> Dict[str, int]:
    coverage = alignment.count_coverage(contig=chromosome, start=position, stop=position + 1)
    return {'A': coverage[0][0], 'C': coverage[1][0], 'G': coverage[2][0], 'T': coverage[3][0]}


@click.group()
@click.option('--bed_file', type=click.Path(), required=True, help='BED file with SNP postions')
@click.option('--output_json_file', type=click.Path(), required=True, help='JSON output file')
@click.option('--minimum_coverage', default=30, show_default=True, help='Only consider SNP positions with a lest this coverage in both tumor and normal')
@click.option('--homozygosity_threshold', default=0.95, show_default=True, help='Consider a SNP position homozygote if frequency of most common allele is this or higher')
@click.pass_context
def client(ctx, bed_file, output_json_file, minimum_coverage, homozygosity_threshold):

    ctx.obj['minimum_coverage'] = minimum_coverage
    ctx.obj['homozygosity_threshold'] = homozygosity_threshold
    ctx.obj['output_json_file'] = output_json_file

    genotyper = Genotyper(minimum_coverage=ctx.obj['minimum_coverage'],
                          homozygosity_threshold=ctx.obj['homozygosity_threshold'])

    with open(bed_file, 'rt') as f:
        snp_coordinates = parse_bed_file(f.read().splitlines())

    ctx.obj['snp_coordinates'] = snp_coordinates
    ctx.obj['genotyper'] = genotyper

    results = {}

    results['input'] = defaultdict(dict)
    results['output'] = defaultdict(dict)

    results['input']['files'] = {'bed-file': os.path.basename(bed_file)}
    results['input']['settings'] = {'minimum-coverage': ctx.obj['minimum_coverage'],
                                    'homozygosity-threshold': ctx.obj['homozygosity_threshold']}

    ctx.obj['results'] = results


@client.command()
@click.option('--normal_bam_file', type=click.Path(), required=True, help='Normal BAM file. Must be indexed.')
@click.option('--tumor_bam_file', type=click.Path(), required=True, help='Tumor BAM file. Must be indexed.')
@click.pass_context
def somatic(ctx, tumor_bam_file, normal_bam_file):

    results = ctx.obj['results']

    results['input']['files']['normal-bam_file'] = os.path.basename(normal_bam_file)
    results['input']['files']['tumor-bam-file'] = os.path.basename(tumor_bam_file)

    normal_snps = get_snps(coordinates=ctx.obj['snp_coordinates'],
                           genotyper=ctx.obj['genotyper'],
                           get_counts=lambda chromosome, position: get_counts(alignment=AlignmentFile(normal_bam_file), chromosome=chromosome, position=position))

    tumor_snps = get_snps(coordinates=ctx.obj['snp_coordinates'],
                          genotyper=ctx.obj['genotyper'],
                          get_counts=lambda chromosome, position: get_counts(alignment=AlignmentFile(tumor_bam_file), chromosome=chromosome, position=position))

    normal_genotypes = {}
    for snp in normal_snps:
        normal_genotypes[snp.__str__()] = snp.genotype if snp.genotype else ''
    results['output']['normal-genotypes'] = normal_genotypes

    tumor_genotypes = {}
    for snp in tumor_snps:
        tumor_genotypes[snp.__str__()] = snp.genotype if snp.genotype else ''
    results['output']['tumor-genotypes'] = tumor_genotypes

    # Only consider SNPs which are genotyped in both normal and tumor sample.
    genotyped_snp_pairs = []
    for normal_snp, tumor_snp in zip(normal_snps, tumor_snps):
        if normal_snp.genotype and tumor_snp.genotype:
            genotyped_snp_pairs.append({'normal': normal_snp, 'tumor': tumor_snp})

    # def count_heterozygotes(sample: str) -> int:
    #     """Exactly as advertized. Counts the number of heterozygote sites."""
    #     return len([pair for pair in genotyped_snp_pairs if pair[sample].is_heterozygote()])

    # def mean_minor_allele_frequency_at_homozygote_sites(sample: str) -> float:
    #     """Computes the mean minor allele frequency at sites which are homozygote in the NORMAL sample."""
    #     return mean([pair[sample].minor_allele_frequency() for pair in genotyped_snp_pairs if pair['normal'].is_homozygote()])

    results['output']['summary'] = {'snps-total': len(ctx.obj['snp_coordinates']),
                                    'snps-genotyped': len(genotyped_snp_pairs)}

    # if genotyped_snp_pairs:
    #     results['output']['summary']['heterozygotes-fraction-normal'] = float('%.4f' % (count_heterozygotes(sample='normal') / len(genotyped_snp_pairs)))
    #     results['output']['summary']['heterozygotes-fraction-tumor'] = float('%.4f' % (count_heterozygotes(sample='tumor') / len(genotyped_snp_pairs)))
    #     results['output']['summary']['mean-maf-homozygote-sites-normal'] = float('%.4f' % mean_minor_allele_frequency_at_homozygote_sites('normal'))
    #     results['output']['summary']['mean-maf-homozygote-sites-tumor'] = float('%.4f' % mean_minor_allele_frequency_at_homozygote_sites('tumor'))

    with open(ctx.obj['output_json_file'], 'w') as json_file_handle:
        json.dump(results, json_file_handle, indent=4)


@client.command()
@click.option('--bam_file', type=click.Path(), required=True, help='BAM file. Must be indexed.')
@click.pass_context
def germline(ctx, bam_file):

    results = ctx.obj['results']

    results['input']['files']['bam-file'] = os.path.basename(bam_file)

    snps = get_snps(coordinates=ctx.obj['snp_coordinates'],
                    genotyper=ctx.obj['genotyper'],
                    get_counts=lambda chromosome, position: get_counts(alignment=AlignmentFile(bam_file), chromosome=chromosome, position=position))

    genotypes = {}
    for snp in snps:
        genotypes[snp.__str__()] = snp.genotype if snp.genotype else ''
    results['output']['genotypes'] = genotypes

    results['output']['summary'] = {'snps-total': len(snps),
                                    'snps-genotyped': len([snp for snp in snps if snp.genotype])}

    with open(ctx.obj['output_json_file'], 'w') as json_file_handle:
        json.dump(results, json_file_handle, indent=4)


def run():
    client(obj={})

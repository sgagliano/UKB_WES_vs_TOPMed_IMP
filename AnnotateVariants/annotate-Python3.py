##add annotations for each variant outputted by `correlation.py` script
##annotations from https://useast.ensembl.org/info/genome/variation/prediction/predicted_data.html

import argparse
import gzip
import pysam

argparser = argparse.ArgumentParser(description = 'Adds annotations from VCF/BCF.')
argparser.add_argument('-i', '--in', metavar = 'file', dest = 'in_file', required = True, help = 'Input file with variants.')
argparser.add_argument('-c', '--consequences', metavar = 'file', dest = 'in_CSQ_VCF', required = True, help = 'Input VCF/BCF file with CSQ INFO field.')
argparser.add_argument('-o', '--out', metavar = 'file', dest = 'out_file', required = True, help = 'Output file.')

so_terms = [
    'transcript_ablation',
    'splice_acceptor_variant',
    'splice_donor_variant',
    'stop_gained',
    'frameshift_variant',
    'stop_lost',
    'start_lost',
    'transcript_amplification',
    'inframe_insertion',
    'inframe_deletion',
    'missense_variant',
    'protein_altering_variant',
    'splice_region_variant',
    'incomplete_terminal_codon_variant',
    'start_retained_variant',
    'stop_retained_variant',
    'synonymous_variant',
    'coding_sequence_variant',
    'mature_miRNA_variant',
    '5_prime_UTR_variant',
    '3_prime_UTR_variant',
    'non_coding_transcript_exon_variant',
    'intron_variant',
    'NMD_transcript_variant',
    'non_coding_transcript_variant',
    'upstream_gene_variant',
    'downstream_gene_variant',
    'TFBS_ablation',
    'TFBS_amplification',
    'TF_binding_site_variant',
    'regulatory_region_ablation',
    'regulatory_region_amplification',
    'feature_elongation',
    'regulatory_region_variant',
    'feature_truncation',
    'intergenic_variant',
    'LoF'
]

if __name__ == '__main__':
    args = argparser.parse_args()

    with gzip.open(args.in_file, 'rt') as ifile, pysam.VariantFile(args.in_CSQ_VCF) as ivcf, gzip.open(args.out_file, 'wt') as ofile:
        header = ifile.readline().rstrip().split()
        if any(x not in ['CHROM', 'POS', 'REF', 'ALT', 'N_GT', 'IMP_AF', 'GT_AF', 'DOSE_AF', 'IMP_R2', 'GT_vs_GT', 'GT_vs_DS'] for x in header):
            raise Exception('Wrong header!')


        vcf_meta = ivcf.header.info.get('CSQ', None)
        if vcf_meta is None:
            raise Exception('Missing CSQ INFO field!')

        csq_header = vcf_meta.description.split(':', 1)[-1].strip().split('|')
        has_chr_prefix = list(ivcf.header.contigs)[0].startswith('chr')

        ofile.write('{}\t{}\n'.format('\t'.join(header), '\t'.join(so_terms)))

        for line in ifile:
            fields = dict(zip(header, line.rstrip().split()))
            if has_chr_prefix and not fields['CHROM'].startswith('chr'):
                chrom = 'chr' + fields['CHROM']
                variant_name = f"chr{fields['CHROM']}_{fields['POS']}_{fields['REF']}_{fields['ALT']}"
            elif not has_chr_prefix and fields['CHROM'].startswith('chr'):
                chrom = fields['CHROM'][3:]
                variant_name = f"{fields['CHROM'][3:]}_{fields['POS']}_{fields['REF']}_{fields['ALT']}"
            else:
                chrom = fields['CHROM']
                variant_name = f"{fields['CHROM']}_{fields['POS']}_{fields['REF']}_{fields['ALT']}"

            pos = int(fields['POS'])

            consequences = set()
            for record in ivcf.fetch(chrom, max(0, pos - 1), pos + 1):
                if record.pos == pos and record.ref == fields['REF'] and record.alts[0] == fields['ALT']: # assuming annoation VCF has only bi-allelic entry per line
                    for transcript_consequence in record.info['CSQ']:
                        transcript_consequence = dict(zip(csq_header, transcript_consequence.split('|')))
                        if transcript_consequence['BIOTYPE'] != 'protein_coding':
                            continue
                        for x in transcript_consequence['Consequence'].split(','):
                            consequences.update(x.split('&'))
                        if transcript_consequence['LoF'] == 'HC':
                            consequences.add('LoF')
                    break

            ofile.write('{}\t{}\n'.format( '\t'.join(fields[h] for h in header), '\t'.join(str(int(t in consequences)) for t in so_terms)))

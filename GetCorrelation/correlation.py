import argparse
import pysam
import numpy
import gzip

argparser = argparse.ArgumentParser(description = 'Computes correlation between imputed dosages and "exact" genotypes.')
argparser.add_argument('-i', '--imputed', metavar = 'file', dest = 'in_imputed_VCF', type = str, required = True, help = 'Input VCF/BCF with imputed genotypes.')
argparser.add_argument('-g', '--genotyped', metavar = 'file', dest = 'in_genotyped_VCF', type = str, required = True, help = 'Input VCF/BCF with "real" genotypes.')
argparser.add_argument('-s', '--samples', metavar = 'file', dest = 'in_samples_file', type = str, required = True, help = 'Input samples list.')
argparser.add_argument('-c', '--chromosome', metavar = 'name', dest = 'chromosome', type = str, required = True, help = 'Chromosome name.')
argparser.add_argument('-b', '--begin', metavar = 'number', dest = 'begin', type = int, required = True, help = 'Start position (bp).')
argparser.add_argument('-e', '--end', metavar = 'number', dest = 'end', type = int, required = True, help = 'End position (bp).')
argparser.add_argument('-o', '--output', metavar = 'file', dest = 'out_file', type = str, required = True, help = 'Output file name (will be gzipp\'ed).')


def cache_records(ifile, chrom, start, stop):
    cache = {}
    print('Load into cache region {}:{}-{}'.format(chrom, start, stop))
    for record in ifile.fetch(chrom, start - 1, stop):
        if record.pos < start: # avoid indels which start before the region
            continue
        if len(record.alts) != 1:
            raise Exception('Multi-allelic entries are not supported.')
        variant_name = '{}_{}_{}'.format(record.pos, record.ref, record.alts[0])
        cache[variant_name] = record
    return cache


if __name__ == '__main__':
    args = argparser.parse_args()

    # read samples into array
    samples = []
    with open(args.in_samples_file, 'r') as ifile:
        for line in ifile:
            sample = line.rstrip()
            if sample:
                samples.append(sample)

    print('Loaded {} sample(s).'.format(len(samples)))

    with pysam.VariantFile(args.in_imputed_VCF, 'r') as imputed_file, pysam.VariantFile(args.in_genotyped_VCF, 'r') as genotyped_file, gzip.open(args.out_file, 'w') as ofile:
        print('{} sample(s) in {}'.format(len(imputed_file.header.samples), args.in_imputed_VCF))
        print('{} sample(s) in {}'.format(genotyped_file.header.samples, args.in_genotyped_VCF))

        imputed_file.subset_samples(samples)
        genotyped_file.subset_samples(samples)

        if len(imputed_file.header.samples) != len(genotyped_file.header.samples):
            raise Exception('Different number of samples in input VCFs/BCFs.')

        print('Subsetted {} sample(s).'.format(len(imputed_file.header.samples)))

        imputed_chrom = args.chromosome
        has_prefix = list((imputed_file.header.contigs))[0].startswith('chr')
        if has_prefix and not imputed_chrom.startswith('chr'):
            imputed_chrom = 'chr' + imputed_chrom
        elif not has_prefix and imputed_chrom.startswith('chr'):
            imputed_chrom = imputed_chrom[3:]

        genotyped_chrom = args.chromosome
        has_prefix = list((genotyped_file.header.contigs))[0].startswith('chr')
        if has_prefix and not genotyped_chrom.startswith('chr'):
            genotyped_chrom = 'chr' + genotyped_chrom
        elif not has_prefix and genotyped_chrom.startswith('chr'):
            genotyped_chrom = genotyped_chrom[3:]

        use_ds = 'DS' in imputed_file.header.formats
        if not use_ds:
            print('No DS format meta information found. DS will not be used.')

        ofile.write('CHROM\tPOS\tREF\tALT\tN_GT\tIMP_AF\tGT_AF\tDOSE_AF\tIMP_R2\tGT_vs_GT\tGT_vs_DS\n'.encode())

        genotyped_cache = {}
        genotyped_cache_stop = 0
        gt_vs_gt = 'nan'
        gt_af = 'nan'
        gt_vs_ds = 'nan'
        dose_af = 'nan'
        for imputed_record in imputed_file.fetch(imputed_chrom, args.begin - 1, args.end):
            if imputed_record.pos < args.begin: # avoid overlapping indels which start before the region
                continue
            variant_name = '{}_{}_{}'.format(imputed_record.pos, imputed_record.ref, imputed_record.alts[0]) # imputed entries are always bi-allelic
            if genotyped_cache_stop < imputed_record.pos:
                genotyped_cache = cache_records(genotyped_file, genotyped_chrom, imputed_record.pos, imputed_record.pos + 100000)
                genotyped_cache_stop = imputed_record.pos + 100000

            genotyped_record = genotyped_cache.get(variant_name, None)
            if genotyped_record is None:
                continue

            real_gt_array = []
            imputed_gt_array = []
            imputed_ds_array = []
            for sample, genotyped_fmt in genotyped_record.samples.iteritems():
                if None in genotyped_fmt['GT']:
                    continue
                else:
                    imputed_fmt = imputed_record.samples[sample]
                    real_gt_array.append(sum(genotyped_fmt['GT']))
                    imputed_gt_array.append(sum(imputed_fmt['GT']))
                    if use_ds:
                        imputed_ds_array.append(imputed_fmt['DS'])

            imp_af = imputed_record.info.get('AF', '.')
            imp_r2 = imputed_record.info.get('R2', '.')
            n_gt = len(real_gt_array)
            if n_gt > 0:
                gt_af = sum(real_gt_array) / (2 * len(real_gt_array))
                if use_ds:
                    dose_af = sum(imputed_ds_array) / (2 * len(imputed_ds_array))
            else:
                gt_af = 'nan'
                if use_ds:
                    dose_af = 'nan'

            gt_vs_gt = numpy.corrcoef(real_gt_array, imputed_gt_array)[0, 1]
            if use_ds:
                gt_vs_ds = numpy.corrcoef(real_gt_array, imputed_ds_array)[0, 1]
            ofile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(args.chromosome, imputed_record.pos, imputed_record.ref, imputed_record.alts[0], n_gt, imp_af, gt_af, dose_af, imp_r2, gt_vs_gt, gt_vs_ds).encode())

    print('Done')

import argparse

''' argument parsing '''
parser = argparse.ArgumentParser()
parser.add_argument('-i',
                    help='input file name',
                    default='/scratch0/battle-fs1/GTEx_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf')
parser.add_argument('-f1',
                    help='field number (starting from 1) to cut from',
                    default=1)
parser.add_argument('-f2',
                    help='field number (starting from 1) to cut to',
                    default=3)
parser.add_argument('-sep',
                    help='separator',
                    default='\t')
parser.add_argument('-batch',
                    help='maximum number of lines to process together',
                    default=100000)
parser.add_argument('-comment',
                    help='comment char',
                    default='#')
parser.add_argument('-o',
                    help='output file name',
                    default='results/pycut.txt')

args = parser.parse_args()

fn = args.i
field_start= args.f1
field_end= args.f2
sep = args.sep
batch = args.batch
comment_char = args.comment
outfn = args.o


# fn = '/scratch0/battle-fs1/GTEx_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf'
# field_start=1  # 1-based
# field_end=3    # 1-based
# sep = '\t'
# outfn = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/tmp/vcf_info.txt'
# batch = 100000

idx1 = field_start-1  # python index idx1:idx2
idx2 = field_end      # python index idx1:idx2

infh = open(fn, 'r')

def write_lines(lines, fn, mode='a'):
  outfh = open(fn, mode)
outfh.write('\n'.join(lines))
if(len(lines)>0):
  outfh.write('\n')
outfh.close()

write_lines([], outfn, 'w')

parsed_lines = []
for i, line in enumerate(infh):
  if line.startswith(comment_char):
  continue
line = line.rstrip()
parsed_lines.append(sep.join(line.split(sep)[idx1:idx2]))
if i % batch == 0:
  write_lines(parsed_lines, outfn, mode='a')
parsed_lines = []
print ''.join(["parsed upto line ",  str(i)])
#if (i+1) % 50000 == 0:
#    break

infh.close()

if len(parsed_lines) > 0:
  write_lines(parsed_lines, outfn)


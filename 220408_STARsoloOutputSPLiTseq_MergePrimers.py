import pandas as pd
import argparse
import os
import sys

def get_args():
	parser = argparse.ArgumentParser()

	#parser.add_argument('-s', dest='samfile',
	#	help='SAM file with Split-seq '+\
	#		 'barcode+UMI information in the BAM tags as created by STARsolo')
	parser.add_argument('--suffix', dest='suff', default=None,
		help='Suffix to add to cell barcodes. Useful if merging separate Split-seq experiments '+
		'that might have overlapping barcodes otherwise.')
	#parser.add_argument('-o', dest='oprefix',
	#	help='Output file path/prefix')
	parser.add_argument('-bcf', dest='bcfile',
		help='Path/Name of barcode file')
	parser.add_argument('-sinfo', dest='splinfo',help='File containing barcode sample association')

	args = parser.parse_args()
	return args

def get_bc1_matches(bcfile):
	bc_file = bcfile #Default: 'bc_8nt_v1.csv'
	bc_df = pd.read_csv(bc_file, index_col=0, names=['bc'])
	bc_df['well'] = [i for i in range(0, 48)]+[i for i in range(0, 48)]
	bc_df['primer_type'] = ['dt' for i in range(0, 48)]+['randhex' for i in range(0, 48)]
	bc_df = bc_df.pivot(index='well', columns='primer_type', values='bc')
	bc_df = bc_df.rename_axis(None, axis=1).reset_index()
	bc_df.rename({'dt': 'bc1_dt', 'randhex': 'bc1_randhex'}, axis=1, inplace=True)
	return bc_df

def get_bc_sample_match(sinfofile):
	si_file = sinfofile
	si_df = pd.read_csv(si_file, index_col=0, names=['bc','sample'])
	return si_df

def main():
	args = get_args()
	ifile = sys.stdin
	#oprefix = args.oprefix
	suff = args.suff
	bcfile = args.bcfile
	sinfofile = args.splinfo
	
	bc_df = get_bc1_matches(bcfile)
	#fname = '{}_merged_primers.sam'.format(oprefix)
	print(bc_df.head())
	si_df = get_bc_sample_match(sinfofile)
	print(si_df.head())
	#samples = si_df['sample'].unique().tolist()
	
		
	
	#ofile = open(fname,'w')
	header = open('header.txt','w')
	
	for line in ifile:
		if line.startswith('@'):
			header.write(line)
		else:
			line = line.strip().split('\t')
			read_bc = line[22]
			read_bc = read_bc.split(':')
			r_bc = read_bc[-1]
			if r_bc!='-':
				bc = r_bc.split('_')
				bc3 = bc[0]
				bc2 = bc[1]
				bc1 = bc[2]
				if bc1 in bc_df.bc1_randhex.tolist():
					bc1 = bc_df.loc[bc_df.bc1_randhex==bc1, 'bc1_dt'].values[0]
				bc = bc3+'_'+bc2+'_'+bc1
				cr_tag = 'CR:Z:{}'.format(bc)
				cb_tag = 'CB:Z:{}'.format(bc)
				line[15] = cr_tag
				line[22] = cb_tag
			else:
				r_bc = line[15].split(':')[-1]
				bc1 = r_bc.split('_')[2]
			line = '\t'.join(line)+'\n'
			if bc1 in si_df.bc.tolist():
				sample = si_df.loc[si_df.bc==bc1, 'sample'].values[0]
				of = '{}_merged_primers.sam'.format(sample)
				ofile = open(of,'a')
				ofile.write(line)
				ofile.close()
					
			
			
	header.close()
	
	
if __name__ == '__main__': main()


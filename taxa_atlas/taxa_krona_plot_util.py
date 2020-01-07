import os
import sys
import argparse

def main():
	parser = argparse.ArgumentParser(description="Format consensus_based_on_*.txt file to visualize with Krona chart")
	parser.add_argument("-i","--input_file", help="Input file",required=True)
	parser.add_argument("-o","--op_file", help="Output file",required=True)
	args = parser.parse_args()

	counts = {}
	with open(args.input_file) as f:
		for line in f:
			val = line.strip().split('\t')
			if val[1] in counts:
				counts[val[1]] += 1
			else:
				counts[val[1]] = 1
	fw = open(args.op_file, 'w')
	for key in counts:
		val = '\t'.join(key.split(';'))
		fw.write(str(counts[key]) +'\t' + val + '\n')


if __name__ == '__main__':
	main()
import os
import sys
import numpy as np
taxonomy = {} 
def main():
	#read taxonomy
	with open(str(sys.argv[1])) as f:
		for line in f:
			l= line.strip()
			lsplit =  l.split('\t')
			taxonomy[lsplit[0]] = lsplit[1]
			
	#Read outlier file

	with open(str(sys.argv[2])) as f1:
		for line in f1:
			l = line.strip() 
			lsplit = l.split('\t')
			if len(lsplit) <= 1:
				continue
			else:
				outliers =  lsplit[1].split(';')
				names =[]
				for i in outliers:
					database_name = taxonomy[i].split(';')
					while len(database_name) < 7:
						database_name.append('')
						print database_name
					# print database_name
					names.append(database_name)

				names1 = np.array(names)
				# print names1, names1.shape[1]
				final_name = []
				if len(names1.shape) == 1:
					continue
				# print names1
				# print len(names1.shape), lsplit[0]

				for i in xrange(0, names1.shape[1]):
					if len(np.unique(names1[:,i])) == 1:
						final_name.append(np.unique(names1[:,i])[0])	
					else:
						break
				print lsplit[0] + '\t' + ';'.join(final_name)
				# break
				# taxon_assignment = os.path.commonprefix(names)
				# print lsplit[0] + '\t' + taxon_assignment
				

if __name__ == "__main__":
	main()
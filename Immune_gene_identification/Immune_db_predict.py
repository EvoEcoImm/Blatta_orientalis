#! /usr/bin/python3


import argparse
import re

def parser():
	args=argparse.ArgumentParser()
	args.add_argument('-i','--input',help='The name of the Hmmsearch output_file.')
	args.add_argument('-x','--blastx',help='The name of the blastx output_file')
	args.add_argument('-p','--blastp',help='The name of the blastp output_file')
	args.add_argument('-o','--output',help='The name of the results output file')
	args=args.parse_args()
	return args


hmmfile = open(parser().input,'rt')
blastxfile= open(parser().blastx,'rt')
blastpfile= open(parser().blastp, 'rt')
outfile = open(parser().output,'wt')

'''edit hmm output file "tab"'''
def hmm_refine_dict(file):
	entries= [re.split("  *", entry) for entry in file.read().split('\n') if "#" not in entry and entry]
	'''spec=[(entry[0],entry[2],float(entry[4]), float(entry[7])) for entry in (entry for entry in entries) if float(entry[4])<=0.00001 and float(entry[7])<=0.001] file from online'''
	spec=[(entry[0].split('::')[0],entry[2],float(entry[4]),float(entry[7]),entry[18],entry[19],entry[0].split('::')[3],entry[0].split('::')[1]) for entry in (entry for entry in entries) if float(entry[4])<=0.00001 and float(entry[7]) <=0.001]  
	'''file from trinity'''
	spec.sort(key=lambda spec:(spec[2],spec[3]),reverse=True)
	dict={entry[0]:entry[1:] for entry in spec}
	return dict

'''open orthofamily_gene name file'''
with open("Immunedbgene2family",'rt') as file:
	gene2family_dict={entry[0]:entry[1] for entry in (entry.strip('\n').split('\t') for entry in file.readlines())}


''' edit blastp output file'''
def blastp_refine_dict(file, g2f_dict):
	blastplist=[entry.strip('\n').split('\t') for entry in file.readlines()]
	blastp_re=[(entry[0].split('::')[3],entry[1],float(entry[10])) for entry in blastplist]
	blastp_re.sort(key=lambda blastp_re:blastp_re[2],reverse=True)
	dict={entry[0]:(g2f_dict[entry[1]],entry[2]) for entry in blastp_re if entry[2]<=0.001}
	return dict

'''edit blastx output file'''
def blastx_refine_dict(file,g2f_dict):
	blastxlist=[entry.strip('\n').split('\t') for entry in file.readlines()]
	blastx_re=[(entry[0],entry[1],float(entry[-2])) for entry in blastxlist if float(entry[-2])<=0.001]
	blastx_re.sort(key=lambda blastx_re:blastx_re[2],reverse=True)
	dict={entry[0]:(g2f_dict[entry[1]],entry[2]) for entry in blastx_re}
	return dict

'''combine hmm and blastp output'''

def combine_blastxp(hmm_refine_dict, blastp_refine_dict, blastx_refine_dict,combine_dict):
	for entry in hmm_refine_dict.keys():
		if hmm_refine_dict[entry][-2] in blastp_refine_dict.keys() and hmm_refine_dict[entry][-1] in blastx_refine_dict.keys():
			if hmm_refine_dict[entry][0]==blastp_refine_dict[hmm_refine_dict[entry][-2]][0] and hmm_refine_dict[entry][0]==blastx_refine_dict[hmm_refine_dict[entry][-1]][0]:
				combine_dict[entry]=hmm_refine_dict[entry]+(blastp_refine_dict[hmm_refine_dict[entry][-2]][-1],blastx_refine_dict[hmm_refine_dict[entry][-1]][-1])
	return combine_dict


'''write out'''
def output_write_out(file, dict):
	for gene in dict.keys():
		file.write(gene+'\t'+('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % dict[gene])+'\n')
	file.close()


hmm_dict=hmm_refine_dict(hmmfile)
blastp_dict=blastp_refine_dict(blastpfile, gene2family_dict)
blastx_dict=blastx_refine_dict(blastxfile, gene2family_dict)
final_dict={}
combine_blastxp(hmm_dict, blastp_dict, blastx_dict,final_dict)
output_write_out(outfile,final_dict)




'''
for entry in entries:
	slentry=[entry[0],entry[2],float(entry[4]),float(entry[7])]
	if slentry[2]<=0.00001 and slentry[3]<=0.001:
		spec1.append(slentry)
		if slentry[0] in sldict1.keys() and slentry[3] < sldict1[slentry[0]]:
			sldict1[slentry[0]]=slentry[1:]
		sldict1[slentry[0]]=slentry[1:]

for name in sldict.keys():
	speci=[sitem for sitem in spec if name in sitem]
	if len(speci)>1:
		for i in range(len(speci)-1):
			if speci[i][2] < sldict[name][1]:
				sldict[name]=speci[i][1:]




totalname=[entry[0] for entry in spec]
namecount={entryname:totalname.count(entryname) for entryname in totalname}
namecount={}
for name in sldict.keys():
	if totalname.count(name) > 1:
		namecount[name]=totalname.count(name)

sorted(a, key=lambda a:a[0])[:10]
sorted(a, key=lambda a:a[2])[:10]
dict={}
for entry in a:
	k,v = entry[0], entry[1:]
	dict[k]=v
speci=[s for s in a if "BGER006296-PA" in s] '''

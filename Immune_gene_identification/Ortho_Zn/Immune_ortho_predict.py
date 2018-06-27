with open("Znev_immune_gene",'rt') as file:
	zn_dict={entry[2]:entry[:2] for entry in (entry.strip('\n').split('\t') for entry in file.readlines())}


with open("Orthogroups.csv") as file:
	ortho=[entry.strip('\r\n').replace('Bo_','').split('\t') for entry in file.readlines()]


c=[]
for i in zn_dict.keys():
	for j in ortho:
		if i in j[1] and j not in c:
			c.append(j)

bo=','.join([entry[2] for entry in c]).replace(" ","").split(',')


bodict={}
for i in bo:
	for j in ortho:
		if i in j[2]:
			bodict[i]=j[1].split(', ')


for i in bodict.keys():
	for j in zn_dict.keys():
		if j in bodict[i]:
			bodict[i].insert(bodict[i].index(j)+1,zn_dict[j][0])



with open("../Longest_orfs.pep",'rt') as file:
	name=[[entry1[0],entry1[1],entry1[3],entry1[4]] for entry1 in (entry.strip('>').replace('::',' ').split() for entry in file.readlines() if ">" in entry)]



for gene,iso,trans,ntype in name:
	if trans in bodict.keys():
		bodict[trans].extend([iso,ntype,gene])


finallist=[entry[::-1] for entry in bodict.values()]
finallist.sort(key=lambda finallist:(finallist[0],finallist[3]))
finallist=['\t'.join(entry) for entry in finallist]


with open("ortho_Immune",'wt') as file:
	for i in finallist:
		file.write("%s\n" %i )


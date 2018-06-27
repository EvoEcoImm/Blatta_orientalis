"""immune gene prediction derived from orthogroup"""

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

""" annotation of immune gene from database in trinotate"""

with open("bo_Immunedb",'rt') as file:
	a=[entry.strip('\n').split('\t') for entry in file.readlines()]
	a1=[entry[6] for entry in a]
	Immunedb_dict={entry[6]:entry for entry in a}


with open("../bo_trinotate_annotation_report.xls",'rt') as file:
	ctrinotate=[[entry[1],entry[2],entry[6],entry[7],entry[4]] for entry in (entry.strip('\n').split('\t') for entry in file.readlines()[1:]) if (entry[6]!="." or entry[7]!=".") and entry[4].split('::')[3] in a1]



trinotate_edi=[]
for trans,blastx,blastp,pfam,transid in ctrinotate:
	d=[trans]
	idedi=[transid.split('::')[3]]
	if blastp !=".":
		blastp1=[blastp.split('^')[0],blastp.split('^')[4].replace('E:',''),blastp.split('^')[5].replace('RecName: Full=','').replace(';',''),blastp.split('^')[6].split(';')[1].replace(' ','')]
	if blastp ==".":
		blastp1=['.']*4
	if pfam !=".":
		pfam1=[entry.split('^') for entry in pfam.replace('E:','').split('`')]
		pfam1.sort(key=lambda pfam1:pfam1[-1])
		pfams=[pfam1[0][0],pfam1[0][2],pfam1[0][4]]
	if pfam ==".":
		pfams=['.']*3
	d.extend(idedi+blastp1+pfams+[blastx])
	trinotate_edi.append(d)


for i in trinotate_edi:
	if i[1] in Immunedb_dict.keys():	
		Immunedb_dict[i[1]].extend(i[2:])
		key=Immunedb_dict[i[1]][7]
		Immunedb_dict[key]=Immunedb_dict.pop(i[1],key)

"""combine all immune gene to manually curate"""

with open("ortho_Immune",'rt') as file:
	ortho_Immune_dict={entry[2]:entry for entry in (entry1.strip('\n').split('\t') for entry1 in file.readlines())}

for i in ortho_Immune_dict.keys():
	if i in Immunedb_dict.keys():
		Immunedb_dict[i].extend(ortho_Immune_dict[i])
	else:
		Immunedb_dict[i]=ortho_Immune_dict[i]

allimmunegene=Immunedb_dict.values()
allimmunegene.sort(key=lambda allimmunegene:(allimmunegene[1],allimmunegene[2]))


allimmune=['\t'.join(entry) for entry in allimmunegene]
with open('all_immune','wt') as file:
	for i in allimmune:
		
			file.write("%s\n" %i)









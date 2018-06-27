import re
from collections import Counter
from ete3 import NCBITaxa
ncbi=NCBITaxa()



with open("bo.match",'rt') as file:
	bonr= [[entry[0],entry[11],entry[12]] for entry in (entry.strip('\n').split('\t') for entry in file.readlines()) if entry[11]!=""]





bonredi=[[re.sub("_i\d*","", entry[0]),entry[1],float(entry[2])] for entry in bonr]
bonredi.sort(key=lambda bonredi:(bonredi[0],bonredi[2]),reverse=True)
genedict={entry[0]:(entry[1].split(';')[0],entry[2]) for entry in bonredi}

gene_Metazoa=[]
gene_Other=[]
for a,b in genedict.items():
	name=ncbi.get_lineage(b[0])
	c=ncbi.get_taxid_translator(name)
	d=[c[iname] for iname in name]
	if "Metazoa" in d:
		gene_Metazoa.append([a,d[-1],b[1]])
	else:
		gene_Other.append([a,d[-1],b[1]])


gene_all=gene_Metazoa+gene_Other
gene_tax=[entry[1].split()[0] for entry in gene_all]


with open("nr_gene_taxname_evalue",'wt') as file:
	for a,b,i in gene_all:
		file.write("%s\t%s\t%.2e\n" % (a,b,i))





with open("nr_gene_contamination",'wt') as file:
	for a,b,i in gene_Other:
		file.write("%s\t%s\t%.2e\n" % (a,b,i))




genuscount=[[a,b] for a,b in Counter(gene_tax).items()]
genuscount.sort(key=lambda genuscount:genuscount[1],reverse=True)
with open("nr_genus_statistic",'wt') as file:
	for a,b in genuscount:
		file.write("%s\t%d\n" % (a,b))



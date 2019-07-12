#!/bin/env python
"Identifies directories containing STRINGTIE and collate the FPKM, TPM and coverage values in respective csv files"
import sys
sourcepath='/'.join(sys.argv[0].split('/')[0:-1])
exec(open(sourcepath+'/../HEADER/python-header.py').read())
exec(open(sourcepath+'/../src/bookkeeping.py').read())

systnames=[f.split('/')[1] for f in glob.glob('./*/STRINGTIE')]

dirs=['./stringtie-count-tables','./stringtie-count-pkl']
subprocess.call(['mkdir','-pv',dirs[0],dirs[1]])

tabledict,maxlen,maxsyst={},0,''

print ('\n Reading stringtie generated files for: ')
for syst in systnames:
	print ('\t--> '+syst)
	table = pd.read_table('./'+syst+'/STRINGTIE/count.tab')
	tabledict[syst]=table.drop_duplicates('Gene ID',keep='first')
	if len(tabledict)>maxlen:
		maxlen=len(tabledict)
		maxsyst=syst

table = tabledict[maxsyst]
targetdict = make_dict(table["Gene ID"])
columns=table.columns[0:6]
coverage,fpkm,tpm=table[columns].copy(deep=True),table[columns].copy(deep=True),table[columns].copy(deep=True)
	
for syst in systnames:
	table=tabledict[syst]
	coverage.loc[:,syst] = populate_table(coverage["Gene ID"],table["Gene ID"],table.Coverage,targetdict)
	fpkm.loc[:,syst] = populate_table(fpkm["Gene ID"],table["Gene ID"],table.FPKM ,targetdict)
	tpm.loc[:,syst] = populate_table(tpm["Gene ID"],table["Gene ID"],table.TPM,targetdict)
	
print ('\n')
#filter out PAR_Y and mitochondrial DNA
fpkm = filter_entries(fpkm,'Gene ID','PAR_Y',keep=False)
fpkm = filter_entries(fpkm,'Reference','chrM',keep=False)
tpm = filter_entries(tpm,'Gene ID','PAR_Y',keep=False)
tpm = filter_entries(tpm,'Reference','chrM',keep=False)
coverage=filter_entries(coverage,'Gene ID','PAR_Y',keep=False)
coverage=filter_entries(coverage,'Reference','chrM',keep=False)


#dump all rawdata to a csv file and also to a pickle 
fpkm.to_csv(dirs[0]+'/gene-fpkm.csv')
fpkm.to_pickle(dirs[1]+'/gene-fpkm.pkl')
tpm.to_csv(dirs[0]+'/gene-tpm.csv')
tpm.to_pickle(dirs[1]+'/gene-tpm.pkl')
coverage.to_csv(dirs[0]+'/gene-coverage.csv')
coverage.to_pickle(dirs[1]+'/gene-coverage.pkl')

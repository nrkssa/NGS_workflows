#!/bin/env python
"Identifies directories containing CUFFNORM and collate the gene_FPKM gene_count in respective csv files"

import sys
from joblib import Parallel,delayed
import multiprocessing as mp
import argparse,argcomplete
sourcepath='/'.join(sys.argv[0].split('/')[0:-1])
exec(open(sourcepath+'/../src/python-header.py').read())
exec(open(sourcepath+'/../src/bookkeeping.py').read())

cuffnormdict={  'fpkm':{'genes':{'attr':'genes.attr_table','table':'genes.fpkm_table'}, 
						'isoforms':{'attr':'isoforms.attr_table','table':'isoforms.fpkm_table'}},
				'count':{'genes':{'attr':'genes.attr_table','table':'genes.count_table'}, 
						'isoforms':{'attr':'isoforms.attr_table','table':'isoforms.count_table'}}
			}

def process_count_data(dinfo,cuffdict,syslocation,systnames,outdir='./count-tables',exclude=[],round=10,procno=0):
	mode,feature = dinfo
	refdict=cuffdict[mode][feature]
	print ('Processor:',procno," : ",mode, feature)
	outfile=mode+'-'+feature
	
	for nn, (sysloc,sysname) in enumerate(zip(syslocation,systnames)):
		if nn==0:
			fpkmtable = pd.read_csv(sysloc+'/'+refdict['attr'],usecols=(0,3,4,6),sep='\t')
			targetdict = make_dict(fpkmtable.tracking_id)
			
		fpkmval= pd.read_csv(sysloc+'/'+refdict['table'],sep='\t')
		fpkmtable[sysname] = populate_table(fpkmtable.tracking_id,fpkmval.tracking_id,
												fpkmval[fpkmval.columns[1]],targetdict)
		
	fpkmtable = filter_entries(fpkmtable,'tracking_id','PAR_Y',keep=False)
	fpkmtable = filter_entries(fpkmtable,'locus','chrM',keep=False)
	if len(exclude)>0: fpkmtable = filter_entries(fpkmtable,'tracking_id',exclude,keep=False)

	fpkmtable[systnames] = fpkmtable[systnames].round(round)	
	visibleflag = [False,False,True,False]+[True for _ in systnames]

	#dump all rawdata to a csv file and also to a pickle 
	fpkmtable.to_csv(outdir+'/'+outfile+'.csv')
	fpkmtable.to_json(outdir+'/'+outfile+'.json',orient='table')
	write_datatable_json(fpkmtable,outdir+'/'+outfile,visible=visibleflag)
	
	print ('\t \t--> Completed: ',mode, feature)
	return fpkmtable

if __name__== "__main__":
	parser=argparse.ArgumentParser(description="Assemble the results of CUFFNORM into raw counts and fpkm")
	parser.add_argument('-e','--exclude',dest='exclude',nargs='+',
						help="List of gene_ids to exclude: -e/--exclude list1 list2",default=[])
	parser.add_argument('-np','--ncore',dest='ncore',help="-np/--ncore <number of processors>: default=all",
						default=mp.cpu_count(),type=int)
	parser.add_argument('-r','--round',dest='round',help="-r/--round <number of digits to round>: default=None",
						default=5,type=int)
	parser.add_argument('-o','--outdir',dest='outdir',help="-o/--outdir: default=./count_tables",
						default='./count-tables')
	parser.add_argument('-i','--input',dest='input',help="-i/--input: default=./",default='./')
	argcomplete.autocomplete(parser)
	args=parser.parse_args()

	args.input = args.input.rstrip('/')


	syslocation = sorted(glob.glob(args.input+'/*/CUFFNORM'))
	systnames=[f.split('/')[-2] for f in syslocation]

	print (syslocation)
	
	excludegenes=[]
	for elist in args.exclude: excludegenes += open(elist,'r').read().split()
	if not os.path.isdir(args.outdir): subprocess.call(['mkdir','-pv',args.outdir])
	dinfos=[[mode,feature]  for mode in ['fpkm','count'] for feature in ['genes','isoforms']]
	num_cores = min(len(dinfos),args.ncore)

	jobs = []
	for _nn, dinfo in enumerate(dinfos):
		_args=(dinfo,cuffnormdict,syslocation,systnames)
		_kwargs = dict(outdir=args.outdir,exclude=excludegenes,round=args.round,procno=_nn)
		jobs.append(mp.Process(target=process_count_data,args=_args,kwargs=_kwargs))
	for _job in jobs: _job.start()
	for _job in jobs: _job.join()
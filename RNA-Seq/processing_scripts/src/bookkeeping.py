def get_alloverlaps(dlen):
	'''
	get all possible overlaps between dlen datasets
	'''
	pairs = [[i] for i in range(dlen)]
	tpairs = [[i,j] for i in range(dlen-1) for j in range(i+1,dlen)]
	pairs += tpairs

	interval=np.arange(3,dlen+1)

	for inter in interval:
		ttpairs = []
		for pair in tpairs:
			for i in range(dlen):
				temp = pair + [i]
				if len(set(temp)) == len(temp):
					olapflag = True
					if len(ttpairs) > 0:
						for tp1 in ttpairs:
							if len(set(temp+tp1))==inter: olapflag=False
					if olapflag: ttpairs.append(temp)
		pairs += ttpairs
		tpairs = ttpairs
	return pairs

def dump_command(command,outdir='./',filename='command.sh'):
	'''
	Dump the executed command to a file
	'''
	headstr1='#!/bin/bash'
	headstr2='module load anaconda/3.6'
	
	try:
		curstr = open(outdir+'/'+filename,'r').read()
		fp = open(outdir+'/'+filename,'a')
		if not headstr1 in curstr: fp.write(headstr1 + '\n')
		if not headstr2 in curstr: fp.write(headstr2 + '\n')
		if not command in curstr: fp.write(command + '\n')
		fp.close()
	except:
		fp = open(outdir+'/'+filename,'w')
		fp.write(headstr1 + '\n')
		fp.write(headstr2 + '\n')
		fp.write(command + '\n')
		fp.close()

def prepare_dir(dirname,overwrite=False):
	if (not os.path.isdir(dirname)) or overwrite:
		subprocess.call(['mkdir','-pv',dirname])

def convert_to_absolute_path(objname):
	'''
	convert relative paths to absolute paths
	'''

	import os
	if type(objname)==str: objname = os.path.abspath(objname)
	if type(objname)==list:
		for _nn,_xx in enumerate(objname): objname[_nn] = os.path.abspath(_xx)
	return objname


## DESeq 
def compute_intersection(data):
	'''
	 Extract the rows that intersect across all columns of the dataframe data	 
	'''
	
	inter=data[0]
	for d1 in data[1:]:
		inter=inter.intersection(d1)
	return inter


def compute_cross_intersection(data):
	'''
	compute the number of rows that overlap between each column in the dataframe
	'''
	
	cross_inter=[]
	for i in range(len(data)-1):
		for j in range(i+1,len(data)):
			cross_inter.append(len(data[i].intersection(data[j])))

	if len(data)>2:
		inter=data[0]
		for d1 in data[1:]:
			inter=inter.intersection(d1)
		cross_inter.append(len(inter))

	return cross_inter


def filterframe(data,colval):
	'''
	perform a logical operation on a column in a dataframe
	pass all filters as [colname, operator, filterval]
	'''
	temp=data.copy(deep=True)
	for c1,op,c2 in colval:
		temp=temp[op(temp[c1],c2)]
	return temp

##-----------------------------------------------------------------------------

def make_dict(data):
	'''
    create a dictionary of an array/list such that the entries are numbered numerically
    Input: array/list/dataframe column of the form [e1,e2,...]
    Output: dictionary {e1:1,e2:2,.......}
    '''
	
	ddict={}
	for i,j in enumerate(data):
		ddict[j] = i
	return (ddict)


def filter_entries(frame,column,string,keep=False,show=False):
	'''
	Function to filter entries in a dataframe. The filtering string can either be a str or a list.
	Depending on the type of entry different types of filtering will be performed. 
	The keep flag if True will return the rows that match the search string, if False will return the complement
	'''

	if keep:
		if show: print ("keeping only entries with ",string, 'in field ', column)
		if type(string)==str: frame = frame[frame[column].str.contains(string)]
		if type(string)==list: frame = frame[frame[column].isin(string)]
	else:
		if show: print ("removing all entries with ",string, 'in field ', column)
		if type(string)==str: frame = frame[~frame[column].str.contains(string)]
		if type(string)==list: frame = frame[~frame[column].isin(string)]
	return frame


def extract_values(ddict,name,value):
	'''
	matches the entries in name to a dictionary and populates an array with value
	'''
	
	nn = 0
	finalarr=np.zeros(len(ddict)) 
	for n,v in zip(name,value):
		try:
			finalarr[ddict[n]] = v
		except:
			nn += 1
			pass
	if nn>0:
		print (str(nn)+'/'+str(len(name)), 'genes did not find a dict entry')
	return (finalarr)

def populate_table(targetname,sourcename,sourceval,targetdict):
	returnarr = np.zeros(len(targetname))
	if (targetname.equals(sourcename)):
		returnarr = sourceval.values
	else:
		returnarr = extract_values(targetdict,sourcename,sourceval)
							
	return returnarr

def write_datatable_json(dataframe,filename,cssid='datatable',elementid="tableinput",fixpos=None,script='local',title=None,visible=[]):

	'''
	Convert a dataframe into a html file. Here the data is stored as a JSON file which is rendered
	in the html file using a JAVASCRIPT
	'''

	visibleoff = ['0']                                              # for the index
	for nn,visibleflag in enumerate(visible):
		if not visibleflag: visibleoff.append(str(nn+1))            # set the visibility of the columns based on the boolean strings passed
	
	if not fixpos: fixpos = visibleoff[-1]
	
	jsonfile = filename+'.json'
	fp=open(filename+'.html','w')
	
	fp.write('<html>' + '\n')
	fp.write('<head>' + '\n')
	if script.lower()=='local':
		fp.write('<link rel="stylesheet" href="https://cbi-asang-009.pmacs.upenn.edu/NGS/templates/datatables/DataTables/css/dataTables.jqueryui.css">' + '\n')
		fp.write('<link rel="stylesheet" href="https://cbi-asang-009.pmacs.upenn.edu/NGS/templates/datatables/DataTables/css/jquery.dataTables.min.css">' + '\n')
		fp.write('<link rel="stylesheet" href="https://cbi-asang-009.pmacs.upenn.edu/NGS/templates/datatables/Buttons/css/buttons.dataTables.min.css">' + '\n')
		fp.write('<link rel="stylesheet" href="https://cbi-asang-009.pmacs.upenn.edu/NGS/templates/jquery/jquery-ui.css">' + '\n')
		fp.write('<link rel="stylesheet" href="https://cbi-asang-009.pmacs.upenn.edu/NGS/templates/datatables/display-tables.css">' + '\n')
		fp.write('<script src="https://cbi-asang-009.pmacs.upenn.edu/NGS/templates/jquery/jquery.min.js"></script>' + '\n')
		fp.write('<script src="https://cbi-asang-009.pmacs.upenn.edu/NGS/templates/datatables/datatables.min.js"></script>' + '\n')
		fp.write('<script src="https://cbi-asang-009.pmacs.upenn.edu/NGS/templates/datatables/DataTables/js/dataTables.bootstrap4.min.js"></script>' + '\n')
		fp.write('<script src="https://cbi-asang-009.pmacs.upenn.edu/NGS/templates/datatables/Scroller/js/dataTables.scroller.min.js"></script>' + '\n')
		fp.write('<script src="https://cbi-asang-009.pmacs.upenn.edu/NGS/templates/datatables/Buttons/js/dataTables.buttons.min.js"></script>' + '\n')
		fp.write('<script src="https://cbi-asang-009.pmacs.upenn.edu/NGS/templates/datatables/Buttons/js/buttons.flash.min.js"></script>' + '\n')
		fp.write('<script src="https://cbi-asang-009.pmacs.upenn.edu/NGS/templates/datatables/JSZip/jszip.min.js"></script>' + '\n')
		fp.write('<script src="https://cbi-asang-009.pmacs.upenn.edu/NGS/templates/datatables/pdfmake/pdfmake.min.js"></script>' + '\n')
		fp.write('<script src="https://cbi-asang-009.pmacs.upenn.edu/NGS/templates/datatables/pdfmake/vfs_fonts.js"></script>' + '\n')
		fp.write('<script src="https://cbi-asang-009.pmacs.upenn.edu/NGS/templates/datatables/Buttons/js/buttons.html5.min.js"></script>' + '\n')
		fp.write('<script src="https://cbi-asang-009.pmacs.upenn.edu/NGS/templates/datatables/Buttons/js/buttons.print.min.js"></script>' + '\n')
		fp.write('<script src="https://cbi-asang-009.pmacs.upenn.edu/NGS/templates/datatables/Buttons/js/buttons.colVis.min.js"></script>' + '\n')
		fp.write('<script src="https://cbi-asang-009.pmacs.upenn.edu/NGS/templates/datatables/FixedColumns/js/dataTables.fixedColumns.min.js"></script>' + '\n')
		fp.write('<script src="https://cbi-asang-009.pmacs.upenn.edu/NGS/templates/datatables/FixedHeader/js/dataTables.fixedHeader.min.js"></script>' + '\n')
		fp.write('<script src="https://cbi-asang-009.pmacs.upenn.edu/NGS/templates/datatables/DataTables/js/dataTables.jqueryui.js"></script>' + '\n')
	else:
		fp.write('<script src="https://code.jquery.com/jquery-3.3.1.js"></script> ' + '\n')
		fp.write('<script src="https://cdn.datatables.net/1.10.16/js/jquery.dataTables.min.js"></script>' + '\n')
		fp.write('<link rel="stylesheet" href="https://cdn.datatables.net/1.10.16/css/jquery.dataTables.min.css">' + '\n')
		fp.write('<script src="https://cdn.datatables.net/buttons/1.5.1/js/dataTables.buttons.min.js"></script> ' + '\n')
		fp.write('<script src="https://cdn.datatables.net/buttons/1.5.1/js/buttons.flash.min.js"></script> ' + '\n')
		fp.write('<script src="https://cdnjs.cloudflare.com/ajax/libs/jszip/3.1.3/jszip.min.js"></script> ' + '\n')
		fp.write('<script src="https://cdnjs.cloudflare.com/ajax/libs/pdfmake/0.1.32/pdfmake.min.js"></script> ' + '\n')
		fp.write('<script src="https://cdnjs.cloudflare.com/ajax/libs/pdfmake/0.1.32/vfs_fonts.js"></script> ' + '\n')
		fp.write('<script src="https://cdn.datatables.net/buttons/1.5.1/js/buttons.html5.min.js"></script> ' + '\n')
		fp.write('<script src="https://cdn.datatables.net/buttons/1.5.1/js/buttons.print.min.js"></script> ' + '\n')
		fp.write('<script src="https://cdn.datatables.net/buttons/1.5.1/js/buttons.colVis.min.js"></script> ' + '\n')
		fp.write('<link rel="stylesheet" href="http://cbi-asang-009.pmacs.upenn.edu/NGS/templates/datatables/DataTables-1.10.16/css/jquery.dataTables.min.css">' + '\n')
		fp.write('<link rel="stylesheet" href="https://cdn.datatables.net/buttons/1.5.1/css/buttons.dataTables.min.css">' + '\n')
		fp.write('<link rel="stylesheet" href="http://cbi-asang-009.pmacs.upenn.edu/NGS/templates/display-tables.css">' + '\n')
		fp.write('<link rel="stylesheet" href="https://code.jquery.com/ui/1.12.1/themes/base/jquery-ui.css">' + '\n')
		fp.write('<link rel="stylesheet" href="https://code.jquery.com/ui/1.12.1/themes/base/jqueryui.min.css">' + '\n')
		fp.write('<script src="https://cdn.datatables.net/fixedcolumns/3.2.4/js/dataTables.fixedColumns.min.js"></src>' + '\n')
		fp.write('<script src="https://cdn.datatables.net/1.10.16/js/dataTables.jqueryui.min.js"></src>' + '\n')
		
		
	fp.write('<meta charset=utf-8 />' + '\n')
	fp.write('</head>' + '\n')
	fp.write('<div class="row">' + '\n')
	fp.write('\t<div class="container">' + '\n')
	if title != None: fp.write('<h3>'+title +'</h3><br>\n')
	fp.write('\t<table id="'+cssid+'" style="width:100%;text-align:center" class=" display">' + '\n')    #table 
	fp.write('\t\t<thead>' + '\n')                                                                       #header
	fp.write('\t\t<tr style="text-align: center;">' + '\n')
	fp.write('\t\t<th>index</th>' + '\n')
	
	for ccc in dataframe.columns:
		fp.write('\t\t<th> ' +ccc +' </th>' + '\n')
	fp.write('\t\t</tr>' + '\n')
	fp.write('\t\t</thead>' + '\n')
	
	fp.write('\t\t<tfoot>' + '\n')                                                                       #header
	fp.write('\t\t<tr style="text-align: center;">' + '\n')
	fp.write('\t\t<th>index</th>' + '\n')
	
	for ccc in dataframe.columns:
		fp.write('\t\t<th> ' +ccc +' </th>' + '\n')
	fp.write('\t\t</tr>' + '\n')
	fp.write('\t\t</tfoot>' + '\n')
	fp.write('\t</table>' + '\n')
	fp.write('\t</div>' + '\n')
	fp.write('</div>' + '\n')
	
	#add a filter function to filter out based on gene names
	fp.write('<script>'+'\n')
	fp.write('\t $(document).ready(function() {'+'\n')

	fp.write('\t\t // Setup - add a text input to each footer cell' + '\n')
	fp.write('\t\t $("#'+cssid+' tfoot th").each( function () {  ' + '\n')
	fp.write('\t\t\t var title = $("#'+cssid+' tfoot th").eq( $(this).index() ).text(); ' + '\n')
	fp.write('\t\t\t $(this).html( \'<input type="text" placeholder="Search \'+title+\'" />\' );' + '\n')
	fp.write('\t\t } );' + '\n')


	fp.write('\t var columns=[{"data" : "index"},')
	fp.write(','.join(['{"data" : "'+ccc+'"}' for ccc in dataframe.columns])+'];'+'\n')
	fp.write('\t var table= $("#'+cssid+'").DataTable( {' + '\n')
	fp.write('\t\t dom: \'Bfrtlip\',' + '\n')
	fp.write('\t\t lengthMenu: [[10, 20, 25, 50, 75, -1], [10, 20, 25, 50, 75, "All"]],'+'\n')
	fp.write('\t\t colReorder: true,'+'\n')
	fp.write('\t\t columnDefs: [{"targets":[ '+', '.join(visibleoff) +' ], "visible": false}],'+'\n')
	fp.write('\t\t search: {"regex": true},'+'\n')
	fp.write('\t\t buttons: [{extend:"colvis",text:"show/hide columns"},{extend:"copyHtml5", exportOptions:{columns: ":visible"}},{extend:"csvHtml5",title:"Data_Export",exportOptions:{columns: ":visible"}},{extend:"excelHtml5",title:"Data_Export",exportOptions:{columns: ":visible"}},{extend:"pdfHtml5",title:"Data_Export",exportOptions:{columns: ":visible"}} ],' + '\n')
	fp.write('\t\t processing: true,' + '\n')
	fp.write('\t\t bserverSide: false,' + '\n')
	fp.write('\t\t ajax: { url:"' + jsonfile.split('/')[-1]+'", dataType:"json"},'+'\n')
	fp.write('\t\t columns: columns,' + '\n')	
	fp.write('\t\t autoWidth: true,' + '\n')	
	fp.write('\t\t deferRender: true,' + '\n')	
	fp.write('\t\t fixedColumns: {' + '\n')	
	fp.write('\t\t\t leftColumns:"'+str(fixpos)+'"' + '\n')	
	fp.write('\t\t },' + '\n')	
	fp.write('\t\t fixedHeader: {' + '\n')	
	fp.write('\t\t\t header: true,' + '\n')	
	fp.write('\t\t\t footer: false' + '\n')	
	fp.write('\t\t },' + '\n')	
	fp.write('\t\t scrollX: true,' + '\n')	
	fp.write('\t\t scrollY: "60vh",' + '\n')	
	fp.write('\t\t lengthChange: true,' + '\n')	
	fp.write('\t\t stateSave: true' + '\n')
	fp.write('\t } );' + '\n')
	
	
	fp.write('\t\t // make the footer visible and override the scrollbar' + '\n')
	fp.write("\t\t $( table.table().container() ).on( 'keyup change', 'tfoot input', function () {" + '\n')
	fp.write('\t\t\t table ' + '\n')
	fp.write('\t\t\t\t .column( $(this).parent().index()+\':visible\') ' + '\n')
	fp.write('\t\t\t\t  .search( this.value ) ' + '\n')
	fp.write('\t\t\t\t .draw(); ' + '\n')
	fp.write('\t\t }); ' + '\n')


	fp.write('\t } );' + '\n')
	
	fp.write('</script>'+'\n')
	fp.write('</html>')
	fp.close()
	dataframe.to_json(jsonfile,orient='table')


def write_datatable_html(dataframe,filename,cssid='datatable',elementid="tableinput",genepos=0,script='local'):
	'''
	convert a dataframe to a datatable
	'''
	fp=open(filename,'w')
	fp.write('<html>' + '\n')
	fp.write('<head>' + '\n')
	
	if script.lower()=='local':
		fp.write('<script src="https://cbi-asang-009.pmacs.upenn.edu/NGS/templates/jquery/jquery-3.3.1.min.js"></script>' + '\n')
		fp.write('<script src="https://cbi-asang-009.pmacs.upenn.edu/NGS/templates/datatables/DataTables/js/jquery.dataTables.min.js"></script>' + '\n')
		fp.write('<script src="https://cbi-asang-009.pmacs.upenn.edu/NGS/templates/datatables/Scroller/js/dataTables.scroller.min.js"></script>' + '\n')
		fp.write('<link rel="stylesheet" href="https://cbi-asang-009.pmacs.upenn.edu/NGS/templates/datatables/DataTables/css/jquery.dataTables.min.css">' + '\n')
		fp.write('<link rel="stylesheet" href="https://cbi-asang-009.pmacs.upenn.edu/NGS/templates/display-tables.css">' + '\n')
		fp.write('<!-- ' + '\n')
		fp.write('<script src="https://cdn.datatables.net/1.10.16/js/jquery.dataTables.min.js"></script>' + '\n')
		fp.write('<link rel="stylesheet" href="https://cdn.datatables.net/1.10.16/css/jquery.dataTables.min.css">' + '\n')
		fp.write('<script src="https://code.jquery.com/jquery-3.3.1.js"></script> ' + '\n')
		fp.write('-->' + '\n')
	else:
		fp.write('<script src="https://cdn.datatables.net/1.10.16/js/jquery.dataTables.min.js"></script>' + '\n')
		fp.write('<link rel="stylesheet" href="https://cdn.datatables.net/1.10.16/css/jquery.dataTables.min.css">' + '\n')
		fp.write('<script src="https://code.jquery.com/jquery-3.3.1.js"></script> ' + '\n')
		
	fp.write('</head>' + '\n')
	
	fp.write('<div class="row">' + '\n')
	fp.write('<div class="container">' + '\n')
	table = re.sub('border="1"','id="'+str(cssid)+'"style="width:100%"',dataframe.to_html(classes='display'))
	fp.write(table +'\n')
	fp.write('</div>' + '\n')
	fp.write('</div>' + '\n')
	
	#add a filter function to filter out based on gene names
	fp.write('<script>'+'\n')
	fp.write('$(document).ready(function() {'+'\n')
	fp.write("$('#"+cssid+"').DataTable();"+'\n')
	fp.write('  } ); '+'\n')
	fp.write('</script>'+'\n')
	fp.write('</html>')
	fp.close()


def write_formatted_table_custom(dataframe,filename,cssid='datatable',elementid="tableinput",genepos=0):
	
	'''
	Pure HTML based table and search: works for small tables but is not suitable when loading gene matrices that contain 10s of thousands of rows
	'''
	
	fp=open(filename,'w')
	fp.write('<html>' + '\n')
	fp.write('<head>' + '\n')
	css=open('/home/ramn/workdir/NGS-pipeline/scripts/templates/display-tables.css').read()
	fp.write(css+'\n')
	fp.write('</head>' + '\n')

	fp.write('<div class="row">' + '\n')
	fp.write('<div class="search-container">' + '\n')
	fp.write('<input type="text" id="'+elementid+'" onkeyup="SearchGene()" placeholder="Search for gene name/pattern">' + '\n')
	fp.write('</div>' + '\n')
	fp.write('</div>' + '\n')
	
	
	fp.write('<div class="row">' + '\n')
	fp.write('<div class="container">' + '\n')
	table = re.sub('border=','id="'+str(cssid)+'" border=',dataframe.to_html(classes='table table-striped table-hover'))
	fp.write(table +'\n')
	fp.write('</div>' + '\n')
	fp.write('</div>' + '\n')
	
	
	#add a filter function to filter out based on gene names
	fp.write('<script>'+'\n')
	fp.write('function SearchGene() {'+'\n')
	fp.write('var input, filter, table, tr, td, i;'+'\n')
	fp.write('input = document.getElementById("'+elementid+'");'+'\n')
	fp.write('filter = input.value.toUpperCase();'+'\n')
	fp.write('table = document.getElementById("'+cssid+'");'+'\n')
	fp.write('tr = table.getElementsByTagName("tr");'+'\n')
	fp.write(''+'\n')
	fp.write('if (filter.length>1) {'+'\n')
	fp.write('for (i = 0; i < tr.length; i++) {'+'\n')
	fp.write('  td = tr[i].getElementsByTagName("td")['+str(genepos)+'];'+'\n')
	fp.write('  if (td) {'+'\n')
	fp.write('    if (td.innerHTML.toUpperCase().indexOf(filter) > -1) {'+'\n')
	fp.write('      tr[i].style.display = "";'+'\n')
	fp.write('    } else {'+'\n')
	fp.write('      tr[i].style.display = "none";'+'\n')
	fp.write('    }'+'\n')
	fp.write('  } '+'\n')
	fp.write('  } '+'\n')
	fp.write('  } '+'\n')
	fp.write('else {'+'\n')
	fp.write('for (i = 0; i < tr.length; i++) {'+'\n')
	fp.write('      tr[i].style.display = "";'+'\n')
	fp.write('  } '+'\n')
	fp.write('}'+'\n')
	fp.write('}'+'\n')
	fp.write('</script>'+'\n')
	fp.write('</html>')
	fp.close()

def run_mp_jobs(jobs,queue=None):
	for _job in jobs: _job.start()
	if queue: _results=[queue.get() for _job in jobs]
	results=[[] for _ in jobs]
	for _xx in _results:
		for _xxx,_yyy in _xx.items():
			results[_xxx] = _yyy
	return results
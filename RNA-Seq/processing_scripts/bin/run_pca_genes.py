def compute_PCA(dframe,columns,coverage=0.975,ncomponents=None):
    '''
   This function performs a Principal Component Analysis of the gene expression data.Dimensionality reduction is performed 
   along the gene axis. For a gene matrix with M samples and N genes the input is expected to be a NXM matrix.
   PCA analysis will return a NXm matrix where m is determined by 
      -- the value of the coverage variance (coverage=0.975 implies use all PC modes that contains 97.5% of the signal) or
      -- the number of components prescribed (ncomponents=2 compute PC only for the top two modes)
    '''
    if ncomponents==None:
        pca=PCA(coverage)
    else:
        pca=PCA(n_components=ncomponents)
    
    pca.fit(dframe[columns].T)
    variance=pca.explained_variance_ratio_*100
    pca_transformed = pca.transform(dframe[columns].T)
    pca_frame = pd.DataFrame(pca_transformed)
    pca_frame.columns = ['PC'+str(i+1) for i in range(np.shape(pca_transformed)[1])]
    pca_frame = pca_frame.set_index(columns)
    return [pca_frame,variance]

def plot_pca_static(pca,sampgroup,sampindex,filename='PCA-analysis'):
    "pass pca as a dataframe whose length matches columns. Sampledict will provide a unique name to each entry in the column"
    mew,ms=1,8
    fig=plt.figure(figsize=(6.0,4.0))
    plt.style.use('seaborn')
    grid=gs.GridSpec(100,100)
    ax=fig.add_subplot(grid[0:95,5:70])
    for xx,yy in sampgroup.items():
        marker,color=markers[sampindex[xx]%nmarker],pcacolor[sampindex[xx]%npcacolor]
        ax.plot(pca.loc[yy,'PC1'],pca.loc[yy,'PC2'],color=color,marker=marker,mew=mew,ms=ms,label=xx,ls='None',mfc='None')
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    ax.set_title('PCA clustering of samples')
    ax.legend(loc=(1.01,0.0),markerscale=0.5,ncol=1,columnspacing=0.1)
    fig.savefig(filename+'.svg')
    fig.savefig(filename+'.pdf')
    del fig


def plot_pca_interactive(pca,variance,sampgroup,sampindex,filename='PCA-analysis',width=1200,height=900):
    ''' 
    plot an interactive html of the PCA analysis. This function generates three plots
        1. A 3D plot showing PC1, PC2 and PC3
        2. A 2D plot showing PC1 and PC2
        3. A 2D plot showing explained variance vs PC mode
    '''
    symbols =['circle',"circle-open","square","square-open","diamond","diamond-open","cross","x"]
    nsymbols = len(symbols)
    #'square','x','triangle-left','diamond','triangle-up','triangle-right',
    #  'triangle-down','star-diamond','cross','octagon','star','pentagon','hexagon']

    layout = go.Layout(
        width=width,height=height,
        xaxis=dict(title='PC1',domain=[0.0,0.1],anchor='x'),
        yaxis=dict(title='PC2',domain=[0.0,1.0]),
        scene=dict( 
            domain=dict(x=[0.4, 1.0],y=[0.6, 1.0]),
            xaxis=dict(showbackground=False,backgroundcolor='rgb(230, 230,230)',title='PC1'),
            yaxis=dict(showbackground=False,backgroundcolor='rgb(230, 230,230)',title='PC2'),
            zaxis=dict(showbackground=False,backgroundcolor='rgb(230, 230,230)',title='PC3'),
            camera=dict(up=dict(x=500000,y=500000,z=0),),
            aspectratio = dict( x=2, y=2, z=1.0 ),
            aspectmode = 'manual'
        ),
        xaxis2=dict(title='PC1',domain=[0.0,0.45]),
        yaxis2=dict(title='PC2',domain=[0.0,0.45],anchor='x2'),
        xaxis3=dict(title='Principal Component',domain=[0.55,1.0]),
        yaxis3=dict(title='Explained Variance',domain=[0.0,0.45],anchor='x3'),

        legend=dict(x=0.0,y=1.0, traceorder='normal',font=dict(family='sans-serif', size=12, color='#000' ),
        bgcolor='#E2E2E2', bordercolor='#FFFFFF',borderwidth=2),

    )

    cvariance = np.cumsum(variance)
    trace1=[go.Scatter(x=pca.loc[yy,'PC1'],y= pca.loc[yy,'PC2'],mode='markers',xaxis='x2',yaxis='y2',name=xx,
                    marker=dict(symbol=symbols[sampindex[xx]%nsymbols],size=16,color=pcacolor[sampindex[xx]%npcacolor]), text=pca.loc[yy].index) for  xx,yy in sampgroup.items()]
    if len(variance)>=3:
        trace2 = [go.Scatter3d(x=pca.loc[yy,'PC1'],y= pca.loc[yy,'PC2'],z=pca.loc[yy,'PC3'],mode='markers',name=xx,
                    marker=dict(symbol=symbols[sampindex[xx]%nsymbols],size=10,color=pcacolor[sampindex[xx]%npcacolor]), text=pca.loc[yy].index,showlegend=False) for  xx,yy in sampgroup.items()]
    else:
        trace2 = [go.Scatter3d(x=pca.loc[yy,'PC1'],y= pca.loc[yy,'PC2'],z=[0]*len(variance),mode='markers',name=xx,
                    marker=dict(symbol=symbols[sampindex[xx]%nsymbols],size=10,color=pcacolor[sampindex[xx]%npcacolor]), text=pca.loc[yy].index,showlegend=False) for  xx,yy in sampgroup.items()]

    trace3 = [go.Scatter(x=['PC'+str(i+1) for i in range(len(cvariance))],y=cvariance,mode='markers+lines',showlegend=False,xaxis='x3',yaxis='y3',
                    text=[int(cv) for cv in cvariance],name='Explained variance',marker=dict(size=12),line=dict(width=2))]

    fig = go.Figure(data=trace1+trace2+trace3, layout=layout)
    save_plotly_html(fig,filename,width=width,height=height)


#!/bin/env python
import sys,os
sourcepath=os.path.realpath('/'.join(sys.argv[0].split('/')[0:-1]))
exec(open(sourcepath+'/../HEADER/python-header.py').read())
exec(open(sourcepath+'/../src/plotter.py').read())
parser = argparse.ArgumentParser(description="compute and plot PCA of the FPKM matrix")
parser.add_argument('-m','--matrix',dest='gmt',help="-m/--matrix gmtfile",default='./count-tables/fpkm-genes.csv')
parser.add_argument('-o','--odir',dest='outdir',help="-o/--odir <output dir>: write results to this folder",default='./PCA')
parser.add_argument('--fpkm_lowcut',dest='fpkm_lowcut',help="--fpkm_lowcut <cutoff value>",default=0.0,type=float)
parser.add_argument('--fpkm_upcut',dest='fpkm_upcut',help="--fpkm_upcut <cutoff value>",default=10000.0,type=float)
args = parser.parse_args()

if not os.path.isdir(args.outdir): subprocess.call(['mkdir','-pv',args.outdir])
fpkm=pd.read_csv(args.gmt)
fpkm = fpkm[fpkm.locus.str.contains('chr')]          # consider only the known chromosomes
#filcols=fpkm.columns[0:5].append(fpkm.columns[5:])
#fpkm=fpkm[filcols]
columns=fpkm.columns[5:]
fpkm = fpkm[(fpkm[columns]>=args.fpkm_lowcut).all(axis=1)]
fpkm = fpkm[(fpkm[columns]<=args.fpkm_upcut).all(axis=1)]



delim=[c for c in columns[0] if not c.isalnum()][-1]
samples=([delim.join(f.split(delim)[0:-1]) for f in columns ])
sampledict=OD([(s1,s2) for s1,s2 in zip(columns,samples)])
sampgroup={xx:[] for xx in list(set(samples))}
sampindex = {xx:ii for ii,xx in enumerate(list(set(samples)))}
for xx,yy in zip(columns,samples):
    sampgroup[yy].append(xx)
pca,variance = compute_PCA(fpkm,columns,ncomponents=len(columns))

plot_pca_static(pca,sampgroup,sampindex,filename=args.outdir+'/PCA-genes')
plot_pca_interactive(pca,variance,sampgroup,sampindex,filename=args.outdir+'/PCA-genes')

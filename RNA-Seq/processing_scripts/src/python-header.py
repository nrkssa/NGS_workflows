#!/bin/env python

import os,sys,glob,subprocess,re,copy
import pickle as pkl
import pandas as pd
import pandas.plotting as ppl
import numpy as np

from pathlib import Path
import argparse,argcomplete

from bs4 import BeautifulSoup
import re,operator,time,regex
import json

import collections
from collections import OrderedDict as OD
from functools import reduce

#parallelization
from joblib import Parallel,delayed
import multiprocessing

#plotting
import matplotlib as mpl
from matplotlib import gridspec as gs
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, LinearLocator, AutoLocator, MaxNLocator
import seaborn as sb
from matplotlib_venn import venn2
import matplotlib_venn as venn
from matplotlib import cm
mpl.use('tkagg')

import mpld3
import cmocean
plt.rcParams['svg.fonttype'] = 'none'


#plotly
from plotly import offline as py
import plotly.tools as tls
import plotly.graph_objs as go
import plotly.figure_factory as ff
from plotly import tools
import cufflinks as cf

PLOTLY_CMAPS=['Blackbody','Bluered','Blues','Earth','Electric','Greens','Greys','Hot',
            'Jet','Picnic','Portland','Rainbow','RdBu','Reds','Viridis','YlGnBu','YlOrRd']
condition_colors={'VEH':'#bdbdbd','VEHICLE':'#bdbdbd','DHT':'#e41a1c','THZ':'#4daf4a','MDV':'#4daf4a','DHT+ENZ':'#4daf4a','JQ1':'#377eb8','BCL':'#984ea3'}

#bokeh
from bokeh.palettes import Set1_9 as colors
from bokeh.palettes import Category20_20 as pcacolor
from bokeh.models import LinearColorMapper,ColorMapper,ColumnDataSource,HoverTool
from bokeh.layouts import column,row,gridplot,GridSpec
from bokeh.plotting import figure,output_file, show,save

#colorlover
#import colorlover as cl

#scipy and sklearn
from scipy.cluster.hierarchy import linkage,dendrogram
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
from scipy.optimize import minimize
import scipy.stats as stats

#tqdm
import tqdm as tqdm
from tqdm import tqdm_notebook as progress
from tqdm import tqdm as progressn

#jit
from numba.decorators import jit

markers = ["o","s","v","d","p","h","^","<",">","1","2","3","4","8","*","H","+","D","|","_"]
nmarker,npcacolor,ncolor = len(markers),len(pcacolor),len(colors)

# logo color scheme for plot motif sequences
logocolorscheme= { 'basepairing':{'G': 'blue','C': 'blue','T': 'darkorange','A': 'darkorange','U': 'darkorange'},
			          'meme':{'G': 'orange', 'A': 'red', 'C': 'blue', 'T': 'darkgreen'},
				'classic': {'G': 'orange','T': 'red','U': 'red','C': 'blue','A': 'darkgreen'}}

CMOCEAN_CMAPS={'haline':cmocean.cm.haline,
			   'turbid':cmocean.cm.turbid,
			   'speed':cmocean.cm.speed,
			   'tempo':cmocean.cm.tempo,
			   'gray':cmocean.cm.gray,
			   'phase':cmocean.cm.phase,
			   'balance':cmocean.cm.balance,
			   'delta':cmocean.cm.delta,
			   'curl':cmocean.cm.curl}



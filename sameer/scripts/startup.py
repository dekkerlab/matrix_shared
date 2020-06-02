import sys
import os
import importlib as imp
import time
import glob
import json
import h5py
from IPython.core.display import display, HTML
display(HTML("<style>.container {width:80% !important; }</style/"))

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.gridspec import GridSpec
mpl.rc('animation', html='html5')

import numpy as np
import pandas as pd
import scipy 
import scipy.sparse as sp
from scipy import signal, ndimage
import multiprocess as mp
from collections import defaultdict, namedtuple

import cooler
import cooltools
import bioframe
import mirnylib.plotting

#path = '/net/levsha/share/sameer/analysis_code/'
def add_to_path(main_path):
    paths = []
    excluded = []
    for item in os.walk(main_path):
        exclude = False
        path = item[0]
        files = item[2]

        if 'setup.py' in files:
            excluded.append(path)

        if np.any([exc in path for exc in excluded]):
            exclude = True

        if not exclude and path.find('/_') == -1 and path.find('/.') == -1:
            paths.append(path)

    for item in paths:
        sys.path.insert(-1, item)

add_to_path('/net/levsha/share/sameer/github/mirnylab-experimental/sameer/')


def get_custom_cmap(num=None, alpha=False):
    total_list = ["#FF4A46", "#008941", "#006FA6", "#A30059", "#7A4900", "#0000A6", "#000000", "#FFFF00", 
     "#1CE6FF", "#FF34FF", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87", "#5A0007", 
     "#809693", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80", "#61615A", "#BA0900", 
     "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100", "#000035", "#7B4F4B", 
     "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F", "#372101", "#FFB500", "#A079BF", 
     "#CC0744", "#C0B9B2", "#001E09", "#00489C", "#6F0062", "#0CBD66", "#456D75", "#B77B68",
     "#7A87A1", "#788D66", "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", 
     "#0086ED", "#886F4C", "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", 
     "#FF913F", "#938A81", "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", 
     "#C8A1A1", "#1E6E00", "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", 
     "#D790FF", "#9B9700", "#549E79", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", 
     "#922329", "#5B4534", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C"]
    
    if num is None:
        sublist = total_list
    else:
        sublist = total_list[0:num]
#    else:
#	sublist = total_list
    
    c = []
    for hex in sublist:
        c.append(mcolors.to_rgb(hex))
        
    if alpha:
        color_list = []
        for color in c:
            r,g,b = color
            new_color = (r,g,b,1)
            color_list.append(new_color)
            new_color = (r,g,b,0.75)
            color_list.append(new_color)
    else:
        color_list = c
        
    return color_list


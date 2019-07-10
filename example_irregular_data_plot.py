import json
import os
import argparse

import numpy as np
from numpy import linspace, meshgrid

import matplotlib.pylab as plt
import pandas as pd
import tqdm
from inference.gp_tools import GpRegressor
from matplotlib import cm, ticker
from matplotlib.lines import Line2D
from pandas.io.json import json_normalize
from plotly.graph_objs import (Contour, Heatmap, Layout, Scatter, Scatter3d,
                               Surface)
from plotly.offline import download_plotlyjs, plot
from scipy.interpolate import griddata, interp2d


def grid(x, y, z, resX=100, resY=100): 
    "Convert 3 column data to matplotlib grid" 
    xi = linspace(min(x), max(x), resX) 
    yi = linspace(min(y), max(y), resY) 
    Z = griddata((x, y), z, (xi[None,:], yi[:,None]),method='linear') 
#     Z = griddata(x, y, z, xi, yi) #using matplotlib griddata
    X, Y = meshgrid(xi, yi) 
    return X, Y, Z 

def make_2d_surface_trace(gp_mu_folded,x_gp,y_gp,z_axis_title):
        trace = Contour(
            z=gp_mu_folded,
            x=x_gp,
            y=y_gp,
            colorscale='Viridis',
            colorbar={'title':z_axis_title,
                      'titleside':'right'},
            opacity=0.9,
            line= dict(width=1,smoothing=0.85),
            contours=dict(
                # showlines=False,
                # showlabels=False,
                coloring= "heatmap",
                start=min(gp_mu),
                end=max(gp_mu),
                size=0.05,
                        labelfont=dict(
                            size=15,
                        ),
                )
        )
        return trace

def make_2d_scatter_trace(x,y,z,text_values, name=''):
        trace = Scatter(x=x, 
                        y=y,
                        mode = 'markers',
                        name = name,
                        hoverinfo='text' ,
                        text=text_values,
                        marker = {'size':6,
                                       'color':'red',
                                       'symbol':'cross'} 
        )
        return trace


def generate_2d_layout(title,x_axis_name,y_axis_name,x,y):
    layout = Layout(xaxis={'title':x_axis_name,
                           'showline':True,
                           'range':[min(x),max(x)],
                           'mirror':"ticks",
                           'ticks':"outside",
                           "linewidth":3,
                           "tickwidth": 3,
                           'titlefont':dict(size=40),
                            'tickfont':dict(size=40)                       
                           },
                    yaxis={'title':y_axis_name,
                           'showline':True,
                           'range':[min(y),max(y)],
                           'mirror':"ticks",
                           'ticks':"outside",
                           "linewidth":3,
                           "tickwidth": 3,
                           'titlefont':dict(
                                size=40,
                            ),
                            'tickfont':dict(
                                size=40,
                            )
                          },
                    margin=dict(l=150, r=10, b=120, t=120),
                    hovermode = 'closest',
                    title=title,
                    titlefont=dict(size=40)   ,                 
                    legend= {'y': 0.7,#.1,
                             'x': 0.1,
                             #'orientation': "h",
                             'font':dict(size=40),
                            },
                    )
    return layout


def make_2d_surface_contours(x, y, z, title='', x_axis_name='', y_axis_name='', z_axis_name=''):
        fig = plt.figure()
        ax = plt.axes()
        X, Y, Z = grid(x, y, z, resX=187, resY=187) 

        ax.scatter(x,y,c='red',s=2,zorder=5)
        contours = plt.contour(X, Y, Z, colors='black', vmin=0.5, vmax=max(z))
        plt.contourf(X, Y, Z, 100)
        

        ax.set_xlim([min(x),max(x)])
        ax.set_ylim([min(y),max(y)])                                

        plt.title(title)
        plt.xlabel(x_axis_name)
        plt.ylabel(y_axis_name)
        cbar = plt.colorbar()
        cbar.ax.set_ylabel(z_axis_name)    
        plt.tight_layout()
        plt.show()
        plt.savefig('plot.png')



parser = argparse.ArgumentParser(description='csg conversion tool.')

parser.add_argument('-i','--input', 
                    help = 'Filename to read in',
                    required=True)

parser.add_argument('-x','--x_axis', 
                    help = 'Material to plot of x axis',
                    required=True)                    

args = parser.parse_args()


filename = args.input
with open(os.path.join(filename)) as f:
    results = json.load(f)

results_df = json_normalize(data=results)
df_filtered_by_mat = results_df

# df_filtered_by_mat = results_df[results_df['fw_material']=='eurofer'&&results_df['armour_material']=='tungsten']

x = list(df_filtered_by_mat['fw_thickness'])
y = list(df_filtered_by_mat['armour_thickness'])
z = list(df_filtered_by_mat['leakage_neutron_current.value'])
z_e = list(df_filtered_by_mat['leakage_neutron_current.std_dev'])
labels = [str(i)+'+-'+str(j) for i,j in zip(z,z_e)]

if len(x) < 40:

    coords = list(zip(x,y))

    GP = GpRegressor(coords, z, y_err=z_e)

    x_gp = linspace(start=min(x), stop=max(x), num=len(x)*0.5)
    y_gp = linspace(start=min(y), stop=max(y), num=len(y)*0.5)

    coords_gp = [ (i,j) for i in x_gp for j in y_gp ]

    gp_mu, gp_sigma = GP(coords_gp)
    gp_mu_folded = np.reshape(gp_mu,(len(x_gp),len(y_gp))).T
            
    traces= []
    traces.append(make_2d_surface_trace(gp_mu_folded,x_gp,y_gp,'Neutron Leakage Current (Particle/Source Particle)'))
    traces.append(make_2d_scatter_trace(x,y,z,labels))

    layout = generate_2d_layout('my plot','x axis title','y axis title',x,y)
    print(traces)
    plot({'data':traces,
        'layout':layout},
        filename='example_irregular_plot.html')    

make_2d_surface_contours(x,y,z,'my plot',args.x_axis,'arm','leak')

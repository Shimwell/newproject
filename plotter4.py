
from plotly.offline import download_plotlyjs, plot
from plotly.graph_objs import Scatter, Layout
import json
import pandas as pd
from pandas.io.json import json_normalize

with open('sm2v3ng_results_muir.json') as f:
    results = json.load(f)

results_df = json_normalize(data=results)

def make_plot(x_axis_name, y_axis_name, colour_code_name):

    x_axis_values = results_df[x_axis_name]
    y_axis_values = results_df[y_axis_name]
    colour_code_values = results_df[colour_code_name]

    # print(colour_code_values)

    traces = (Scatter(x=x_axis_values,
                      y=y_axis_values,
                      mode = 'markers',
                      marker=dict(color=colour_code_values,
                                  colorscale='Viridis',
                                  colorbar=dict(
                                      title='Neutron Leakage Current (Particle/Source Particle)'),
                                  showscale=True),
                      ))

    
    layout = {'title':'Neutron Leakage',
              'xaxis':{'title':'First Wall thickness (cm)'},
              'yaxis':{'title':'First Wall Armour thickness (cm)'},
    }

    plot({'data':[traces],
          'layout':layout},
          filename='data_output.html')

make_plot(x_axis_name='fw_thickness',
          y_axis_name='armour_thickness',
          colour_code_name='leakage_neutron_current.value')
#simple first wall model
#this model will provide justification for the project we are going to perform
#it will demonstrate that parameters of the first wall have an impact on the number of neutrons
#and the spectrum of neutrons that pass into the breeder blanket, affecting its function

#the simple model will consist of;
#a spherical shell representing the first wall
#a central, monoenergetic, isotropic neutron source
#tallies 'before' and 'after' the first wall shell to measure neutron current and spectrum
#this will then be extended further to make it more complex


#this first model will investigate how first wall thickness of a single material affects the
#neutron current that passes into the breeder blanket

#need to construct geometry
#need to specify first wall material - eurofer
#need to specify source
#need to construct tallies

#however, first, we must find a justification for the size of the shell used and the thicknesses investigated
#AT THE MOMENT, WE WILL JUST USE A SPHERE OF 1M, WITH A THICKNESS OF 10CM SO THAT WE CAN JUST GET THE CODE WORKING
#I.E WE NEED TO REMEMBER TO CHANGE THIS 




import openmc
import matplotlib.pyplot as plt
import os
import numpy as np
from tqdm import tqdm

#define materials

mats = openmc.Materials()

eurofer = openmc.Material(name='eurofer')
eurofer.set_density('g/cm3', 7.75)
eurofer.add_element('Fe', 89.067, percent_type='wo')
eurofer.add_element('C', 0.11, percent_type='wo')
eurofer.add_element('Mn', 0.4, percent_type='wo')
eurofer.add_element('Cr', 9.0, percent_type='wo')
eurofer.add_element('Ta', 0.12, percent_type='wo')
eurofer.add_element('W', 1.1, percent_type='wo')
eurofer.add_element('N', 0.003, percent_type='wo')
eurofer.add_element('V', 0.2, percent_type='wo')
mats.append(eurofer)

#define a function that basically forms the geometry and runs the simulation according to that geometry

def make_geometry_tallies(batches,nps,inner_radius,thickness):
    
    first_wall_inner_surface = openmc.Sphere(r=inner_radius)
    first_wall_outer_surface = openmc.Sphere(r=inner_radius + thickness, boundary_type='vacuum')

    first_wall = +first_wall_inner_surface & -first_wall_outer_surface
    first_wall = openmc.Cell(region=first_wall)
    first_wall.fill = eurofer
    
    inner_vac_cell = -first_wall_inner_surface
    inner_vac_cell = openmc.Cell(region=inner_vac_cell)

    universe = openmc.Universe(cells=[first_wall,inner_vac_cell])
    geom = openmc.Geometry(universe)

    geom.export_to_xml('geometry')

    vox_plot = openmc.Plot()
    vox_plot.type = 'voxel'
    vox_plot.width = (15, 15, 15)
    vox_plot.pixels = (200, 200, 200)
    vox_plot.filename = 'plot_3d'
    vox_plot.color_by = 'material'
    vox_plot.colors = {eurofer: 'blue'}
    plots = openmc.Plots([vox_plot])
    plots.export_to_xml()

    openmc.plot_geometry()

    os.system('openmc-voxel-to-vtk plot_3d.h5 -o plot_3d.vti')
    os.system('paraview plot_3d.vti')

    sett = openmc.Settings()
    sett.batches = batches
    sett.inactive = 0
    sett.particles = nps
    sett.run_mode = 'fixed source'

    source = openmc.Source()
    source.space = openmc.stats.Point((0,0,0))
    source.angle = openmc.stats.Isotropic()
    source.energy = openmc.stats.Discrete([14.08e6], [1])
    #source.energy = openmc.stats.Muir(e0=14080000.0, m_rat=5.0, kt=20000.0)
    sett.source = source

    sett.export_to_xml('settings.xml')

    #tallies
    particle_filter = openmc.ParticleFilter([1])
    surface_filter_front = openmc.SurfaceFilter(first_wall_inner_surface)
    surface_filter_rear = openmc.SurfaceFilter(first_wall_outer_surface)
    bins = openmc.mgxs.GROUP_STRUCTURES['VITAMIN-J-175']
    #think will need to change this
    energy_filter = openmc.EnergyFilter(bins)

    tallies = openmc.Tallies()

    tally = openmc.Tally(name='incident_neutron_current')
    tally.filters=[surface_filter_front,particle_filter]
    tally.scores = ['current']
    tallies.append(tally)

    tally = openmc.Tally(name='leakage_neutron_current')
    tally.filters = [surface_filter_rear,particle_filter]
    tally.scores = ['current']
    tallies.append(tally)

    tally = openmc.Tally(name='incident_neutron_spectrum')
    tally.filters = [surface_filter_rear,particle_filter,energy_filter]
    tally.scores = ['flux']
    tallies.append(tally)

    tally = openmc.Tally(name='leakage_neutron_spectrum')
    tally.filters = [surface_filter_front,particle_filter,energy_filter]
    tally.scores = ['flux']
    tallies.append(tally)

    model = openmc.model.Model(geom, mats, sett, tallies)
    model.run()

    sp = openmc.StatePoint('statepoint.'+str(batches)+'.h5')

    #we want to retrieve our tallies, but we now want to save them in a .json file
    #therefore, we setup the json file to recieve the tally data
    #for now, we will simply get the json file to get the neutron current
    #first, we specify the 'general' parameters about the setup that we want the .json file to recieve

    json_output = {'inner_radius':inner_radius,
                   'thickness':thickness}

                   #i.e. these are the general parameters about the setup that we want the json file to recieve

    #however, we also want the json file to retrieve the data from the tallies

    #first, we want to retrieve the neutron current at the inner and outer surfaces
    tallies_to_retrieve = ['incident_neutron_current', 'leakage_neutron_current']
    for tally_name in tallies_to_retrieve:
        tally = sp.get_tally(name=tally_name)

        df = tally.get_pandas_dataframe()
        #defining something that stands for dataframe, need to investigate this
        #its basically something that we use to obtain the mean value and the std deviation value of the tally
        tally_result = df['mean'].sum()
        tally_std_dev = df['std. dev.'].sum()

        json_output[tally_name] = {'value': tally_result,
                                   'std_dev': tally_std_dev}


    #next we wnat to retrieve the neutron spectra data at the inner and outer surfaces of the shell

    spectra_tallies_to_retrieve = ['incident_neutron_spectrum','leakage_neutron_spectrum']
    for spectra_name in spectra_tallies_to_retrieve:
        spectra_tally = sp.get_tally(name=spectra_name)
        spectra_tally_result = [entry[0][0] for entry in spectra_tally.mean]
        spectra_tally_std_dev = [entry[0][0] for entry in spectra_tally.std_dev]
        #print(spectra_tally_result)

    return json_output
    
    #note at the moment, this data is not being 'saved' anywhere
    #i.e. is NOT PART OF THE .JSON FILE
    
results = []
num_simulations = 5

for i in tqdm(range(0,num_simulations)):
    thickness=np.linspace(1,5,6)
    #here 
    for t in thickness:
        result = make_geometry_tallies(batches=2,
                                       inner_radius=1,
                                       thickness=t,
                                       nps=1)
        results.append(result)

#having run the simulations, we now need to output the .json file

output_filename = 'simulation_results.json'
with open(output_filename, mode='w', encoding='utf-8') as f:
    json.dump(results, f)
























#print(first_wall_inner_surface)
#thickness=np.linspace(0,1,6)
#first_wall_inner_surface = openmc.Sphere(r=1)
#first_wall_outer_surface = openmc.Sphere (first_wall_inner_surface + thickness)

#instead of defining the geometry this way, we should change this so we can define the thickness


#define cells

#first_wall = +first_wall_inner_surface & -first_wall_outer_surface
#first_wall = openmc.Cell(region=first_wall)
#first_wall.fill = eurofer

#universe = openmc.Universe(cells=[first_wall])
#geom = openmc.Geometry(universe)
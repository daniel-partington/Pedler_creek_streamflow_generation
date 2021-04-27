# -*- coding: utf-8 -*-
"""
Created on Sun Feb 16 17:45:14 2020

@author: guti0021
"""

# coding: utf-8

# #  Analysis of HGS output files for Pedler Creek Model


# In[1]:

#------------------------------------------------------------------------------
#----------------------------Import packages-----------------------------------
#------------------------------------------------------------------------------

#    First let's get the relevant packages required for the analysis
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import datetime
from matplotlib import style
from matplotlib.ticker import MaxNLocator

# In[2]: Setup plot style
plt.rc('axes', linewidth=0.4)
plt.rcParams.update({'font.size': 11})
plt.rcParams.update({'legend.labelspacing':0.1})
plt.rc('font', **{'sans-serif' : 'Arial','family' : 'sans-serif'})
plt.rcParams['lines.linewidth'] = 1
style.use('default')

# In[4]:

p_j = os.path.join

# Script path
Dir = os.getcwd()

# Model path

soils_dir = r'R:\SE-EphemeralRivers\Pedler_model\run_models\2_CoarseMesh_scenarios\Test_soils2'

models_dirs = {}

for i in range(1,5):
    print(i)
    key = 'Scenario{}'.format(str(i))
    models_dirs[key] = p_j(soils_dir, key)


models_dirs['incised_10m'] =  r'R:\SE-EphemeralRivers\Pedler_model\run_models\2_CoarseMesh_scenarios\Test_river_incision2\Scenario_10m'

# To save final figures
figures2 = r'c:\Users\guti0021\OneDrive - Flinders\2nd_paper\HESS\Files\figures'


# In[7]: Create functions to easily load the files
def get_header(fname, header_line=1):
    '''Function to get the header from the output files:
       First 'split' gets rid of the '=' and then replace the double
       quotations to get a list with str format
       Function takes:
           fname = file name
           header_line = line location of the header
    '''
    with open(fname, 'r') as f:
        header = f.readlines()[header_line].split('=')[1].replace('"','').replace('\n', '').split(',')
    return header

def hgs_wat_bal_style_to_df(fname, ftype='wat_bal'):
    '''Function loads the output file into a pd DF specifying the header line
       and skiping appropiate rows depending on the file. Function takes:
           fname = file name
           ftype = file type. Options are 'wat_bal', 'hydrograph' or 'obs_well'
    '''
    if ftype == 'hydrograph':
        header_line = 0
        skiprows = 2
    elif ftype == 'wat_bal':
        header_line = 1
        skiprows = 3
    else:
        header_line = 1
        skiprows = 2
    # End if
    header = get_header(fname, header_line=header_line)
    df = pd.read_csv(fname, delim_whitespace=True, index_col=0,
                     skiprows=skiprows, names=header)
    df.dropna(inplace=True, axis=1, how='all')
    df.dropna(inplace=True, axis=0)
    return df

# In[] Load the WATER BALANCE files

wat_bal_dict = {}
hmc_dict = {}
sand = {}
clay = {}
loam = {}
for key in models_dirs:

    model_dir= models_dirs[key]

    prefix = "Pedler"
    wat_bal_file = prefix + "o.water_balance.dat"
    hmc_hydro = prefix +'o.hmc_Catchment_outlet_hydrograph.dat'

    #----- Scenario Characteristics
    if 'Test_soils' in model_dir:
        print('Soil properties scenarios')

        if 'Scenario4' in model_dir:
            sand[key] = 0.314
            loam[key] = 0.108
            clay[key] = 0.0009
            print(key)
        elif 'Scenario3' in model_dir:
            sand[key]= 1.06
            loam[key]= 0.0624
            clay[key]= 0.0009
            print('Scenario 3')
        elif 'Scenario2' in model_dir:
            sand[key]= 0.314
            loam[key]= 0.0624
            clay[key]= 0.0009
            print('Scenario 2')
        elif 'Scenario1' in model_dir:
            sand[key]= 1.06
            loam[key]= 0.108
            clay[key]= 0.0009
            print('Scenario 1')

    else:
        key = 'Scenario8'
        print('River incision Scenario')
        sand[key]= 1.06
        loam[key]= 0.108
        clay[key]= 0.0009

    #-----Read the WB file and set up a folder to store plots
    os.chdir(model_dir)

    wat_bal_dict[key] = hgs_wat_bal_style_to_df(wat_bal_file, ftype='wat_bal')

    hmc_dict[key] = hgs_wat_bal_style_to_df(hmc_hydro, ftype='hydrograph')

    figures = os.path.join(model_dir, 'figures')
    if not os.path.exists(figures):
        os.mkdir(figures)

#---- Analyse WB from day 5 to avoid large errors at the begining of the run
for df in wat_bal_dict:
    wat_bal_dict[df] = wat_bal_dict[df].loc[5:]

# In[]
wat_bal_volumes = {}
for df in wat_bal_dict:

    wat_bal_dict[df]['dt'] = [wat_bal_dict[df].index.tolist()[0]] + \
                                 list(np.diff(wat_bal_dict[df].index))
    wb_ref = wat_bal_dict[df].copy()
    wb_ref = wb_ref.multiply(wat_bal_dict[df]['dt'].tolist(), axis=0)
    wat_bal_volumes[df] = wb_ref

summary_dict = {}
for scenario in wat_bal_dict:
    print (scenario)

    df = wat_bal_dict[scenario]
    summary = {}
    for key in df:
        summary[key] = df[key].sum()
        print(key)

    summary_dict[scenario] = summary

#===========================================================Summary in volumes
summary_dict_vol = {}
for scenario in wat_bal_dict:
    print (scenario)

    df = wat_bal_volumes[scenario]
    summary = {}
    for key in df:
        summary[key] = df[key].sum()

    summary_dict_vol[scenario] = summary

# In[9]: First, we can look at the overall water balance
'''
Key to the water balance:
 - Time
 - FluidTransfer 1     _(Sink)_
 - rainfall            _(Source)_
 - ET                  _(Sink)_
 - crit depth          _(Sink)_
 - NET1 Sources/Sinks  _(Sum of the sources and sinks)_
 - PM
 - Overland
 - NET2 Accumulation   _(Sum of the rates of change of storage in PM + OLF)
 - ERROR (NET1-NET2)
 - Error rel
 - Error percent
 - Infilt              _(Infiltration into the subsurface from the surface)_
 - Exfilt              _(Exfiltration from the subsurface to the surface)_
 - Tot_ET
 - Delta_Stor-Int


The water balance shows all the components in the water balance when
considering the whole modelling domain. That is, it is a representation of:
            dS/dt = Q_in - Q_out
   where:
   ds/dt is the sum of the terms PM and Overland

   Q_in is the inflow to the system which is essentially
        the precipitation in this model

   Q_out is the outflow from the system from the critical depth boundary
         (CritDepth), evapotranspiration (Total_ET) and the fluid transfer
         (FluidTransfer_1)
'''

# In[ ]: Specify colors, labels, and  lines style/width to use in all plots

colours = {'crit_depth':'#660066', 'Tot_ET':'#74d600',
           'FluidTransfer_1':'#ffa700','PM':'#66b2b2',
           'Overland':'#0057e7','ERROR (NET1-NET2)':'#ae1d1d',
           'rainfall':'#0a1c5c', 'Infilt': '#a3c1ad', 'Exfilt': '#1ebbd7',
           'Canopy_evap':'#77ab59', 'Surf_evap':'#36802d',
           'PM_evap':'#234d20', 'PM_trans':'#c9df8a',
           'NET2 Accumulation':'cyan'}

labels = {'crit_depth':'Critical depth', 'Tot_ET':'Total ET',
           'FluidTransfer_1':'Fluid transfer','PM':'Porous media',
           'Overland':'Overlandflow','ERROR (NET1-NET2)':'Error',
           'rainfall':'Rainfall', 'Infilt': 'Infilt', 'Exfilt': 'Exfilt',
           'Canopy_evap':'Canopy evap', 'Surf_evap':'Surface evap',
           'PM_evap':'PM evap', 'PM_trans':'PM trans',
           'NET2 Accumulation':'Change in Storage' }

linesw = {'crit_depth':1.1, 'Tot_ET':1.1, 'FluidTransfer_1':1.1, 'PM':1.5,
           'Overland':1.1, 'ERROR (NET1-NET2)':1.1, 'rainfall':0.5,
           'Infilt':1.1, 'Exfilt':1.1,'Canopy_evap':1.1, 'Surf_evap':1.1,
           'PM_evap':1.1, 'PM_trans':1.1, 'NET2 Accumulation':1.1 }

linesty = {'crit_depth':'-', 'Tot_ET':'-', 'FluidTransfer_1':'-', 'PM':'-',
           'Overland':'-', 'ERROR (NET1-NET2)':':', 'rainfall':':',
           'Infilt':'-', 'Exfilt':'-', 'Canopy_evap':'-', 'Surf_evap':'-',
           'PM_evap':'-', 'PM_trans':'-',  'NET2 Accumulation':'-'}

scolours = {'Scenario8':'#0084ff', 'Scenario1':'#44bec7',
            'Scenario2':'#ffc300', 'Scenario4':'#fa3c4c',
            'Scenario3':'#d696bb'}

closep = plt.close
a = 'all'

# In[ ]: Plot the CUMULATIVE VALUES for the principal components of the WB

#==============================================================================
#==============================CUMULATIVE VALUES===============================
#==============================================================================
components = ['rainfall','Tot_ET','FluidTransfer_1', 'crit_depth',
              'PM', 'Overland', 'NET2 Accumulation']

scenarios_ordered = ['Scenario1', 'Scenario2', 'Scenario3',
                     'Scenario4', 'Scenario8']

fig = plt.figure(figsize=(8, 9))

for i , component in enumerate(components):
    print (i, component)

    ax = fig.add_subplot(4, 2, i + 1)

    for df in scenarios_ordered:
        print (df)

        colour = scolours[df]
        labelx = df

        wat_bal_dict[df]['dt'] = [wat_bal_dict[df].index.tolist()[0]] + \
                                 list(np.diff(wat_bal_dict[df].index))
        wb_ref = wat_bal_dict[df].copy()
        wb_ref = wb_ref.multiply(wat_bal_dict[df]['dt'].tolist(), axis=0)
        wb_ref[component].cumsum().plot(ax=ax, color=colour,
                                                      lw=1, ls='-',
                                                      label=labelx,
                                                      fontsize=10)

        if i == 6:
            ax.legend(loc="center left", frameon=0, fancybox=0,
                      bbox_to_anchor=(1.15, .45), fontsize=11,
                      labelspacing=1.2)

    ax.set_xlabel('')
    ax.text(0.40, 0.85, labels[component], transform=ax.transAxes, fontsize=11,
            bbox=dict(facecolor='white', edgecolor='black', alpha=0.5,
            boxstyle='round,pad=0.15', linewidth=0.5))

# Adjust subplots spacing
left  = 0.15
right = 0.98
bottom = 0.1
top = 0.97
wspace = 0.35
hspace = 0.2
plt.subplots_adjust(left=left, bottom=bottom, right=right,
                    top=top, wspace=wspace, hspace=hspace)

# X and Y Axis labels
fig.text(0.3, 0.03, 'Time (d)',  fontsize=12)
fig.text(0.0050, 0.65, 'Cumulative volumes (m$^3$)',  fontsize=12, rotation=90)

# Save figure
fig.savefig(os.path.join(figures2, 'Cum_values_ALL.png'),
            bbox_inches = "tight", dpi=600)
# In[ ]: Get the values for the table
scenarios_ordered = ['Scenario8', 'Scenario1', 'Scenario2',
                     'Scenario3', 'Scenario4']

for i, dfx in enumerate (scenarios_ordered):
    print (i, dfx)
    df = wat_bal_dict[dfx]

    # Get the values as porcentage for the pie chart
    total = (abs(summary_dict_vol[dfx]['PM']) +
             abs(summary_dict_vol[dfx]['Overland']) +
             abs(summary_dict_vol[dfx]['Tot_ET']) +
             abs(summary_dict_vol[dfx]['FluidTransfer_1']) +
             abs(summary_dict_vol[dfx]['crit_depth']) +
             abs(summary_dict_vol[dfx]['ERROR (NET1-NET2)']))

    for item in components:
        val = ((abs(summary_dict_vol[dfx][item])) / total) * 100
        print(item, round(val, 2))

# In[ ]: Plot the Water balance with a pie chart showing the percentage that
#        each component represents

#==============================================================================
#================================WATER BALANCE=================================
#==============================================================================
components = ['PM', 'Overland', 'Tot_ET','FluidTransfer_1',
              'crit_depth', 'ERROR (NET1-NET2)']

linesw = {'crit_depth':0.8, 'Tot_ET':0.8, 'FluidTransfer_1':1.2, 'PM':1.1,
           'Overland':0.8, 'ERROR (NET1-NET2)':0.8, 'rainfall':0.3,
           'Infilt':0.8, 'Exfilt':0.8,'Canopy_evap':0.8, 'Surf_evap':0.8,
           'PM_evap':0.8, 'PM_trans':0.8, 'NET2 Accumulation':0.8}

scenarios_ordered = ['Scenario8', 'Scenario1']

figx = 17/2.54  # figure width
figy = 15/2.54  # figure height
x0 = 0.17  # left edge (0-1)
y0 = 0.08   # bottom edge (0-1)
xs = 0.016 # horizontal space between panels (0-1)
ys = 0.025  # vertical space between panels (0-1)
fig = plt.figure(figsize=(figx, figy))

# ---------------------------------------------------Top plot
xw = 0.65 # width of panel (0-1)
yh = 0.15   # height of panel(0-1)
ax0 = plt.axes([x0, 0.80, xw, yh])


df = wat_bal_dict["Scenario8"]
df['rainfall'].plot(ax=ax0, color='#325b84', fontsize=8,
                    label='Rainfall', linewidth=0.7)

# Set axis limits
ax0.invert_yaxis()
ax0.set_ylim((ax0.get_ylim()[0])+200000, 0)
ax0.set_yticks([0, 500000])

# Set axis labels
ax0.set_xticklabels('')
ax0.set_xlabel('')

ax0.legend(loc=4, fancybox=0, frameon=0, fontsize=8)

ax0.tick_params(direction='out', length=3, width=0.5)

for i, dfx in enumerate (scenarios_ordered):
    print (i, dfx)
    # -----------------------------------------------Left Plots
    xw1 = 0.65 # width of panel (0-1)
    yh1 = 0.335   # height of panel(0-1)
    ax = plt.axes([x0, y0+(ys*i)+(yh1*i), xw1, yh1])

    df = wat_bal_dict[dfx]

    min_val = df['Tot_ET'].min()
    max_val = 1000000
    diff = max_val - min_val

    for component in components:
        colour = colours[component]
        labelx = labels[component]
        df[component].plot(ax=ax, color=colour,
                                   linewidth=linesw[component],
                                   ls=linesty[component],
                                   label=labelx, fontsize=8)

    ax.set_ylim(min_val * 10 , max_val)
    ax.tick_params(direction='out', length=3, width=0.5)

    if dfx == 'Scenario1':
        ax.set_ylabel('Flux (m$^3$/d)', fontsize=10)
        ax.get_yaxis().set_label_coords(-0.2,-0.1)


    if dfx == 'Scenario8':
        ax.set_xlabel('Time (d)', fontsize=10)
    else:
        ax.set_xlabel('')
        ax.set_xticklabels('')

    # ---------------------------------------------------------------Right Plot
    # Get the total flux amounts for each of the components of the WB
    #summary = {}
    #for key in df.keys():
    #    summary[key] = df[key].sum()
    #    print(key)

    idx_start = 0
    idx_end = 1183

    # Get the values as porcentage for the pie chart
    total = (abs(summary_dict_vol[dfx]['PM']) +
             abs(summary_dict_vol[dfx]['Overland']) +
             abs(summary_dict_vol[dfx]['Tot_ET']) +
             abs(summary_dict_vol[dfx]['FluidTransfer_1']) +
             abs(summary_dict_vol[dfx]['crit_depth']) +
             abs(summary_dict_vol[dfx]['ERROR (NET1-NET2)']))

    # Get the elements needed to call the pie chart
    values = []
    labels_p = []
    colours_p = []

    for item in components:
        val = ((abs(summary_dict_vol[dfx][item])) / total) * 100
        print(item, round(val, 2))

        values.append(val)
        #label = labels[item] + ' ({})'.format("%.0f%%" % val)
        labels_p.append(labels[item])
        colours_p.append(colours[item])

    # -----------------------------------------------Right Plots
    # Call the pie subplot
    xw2 = 0.2  # width of panel (0-1)
    yh2 = 0.2  # height of panel(0-1)
    ax = plt.axes([0.85, y0+(ys*i)+(yh1*i)-0.01, xw2, yh2])

    wedges, texts  = ax.pie(values, colors=colours_p, startangle=90,
                            wedgeprops={'width':0.5, 'edgecolor':'white'})
    ax.axis('equal')


    if dfx == 'Scenario8':
        ax.text(0.05, 1.48,  '{} (River Incised)'.format(dfx),
                transform=ax.transAxes, fontsize=8)

    else:
        ax.text(0.23, 1.48,  '{}'.format(dfx),
                transform=ax.transAxes, fontsize=8)

    ax.text(0.18, 1.36,  'Soil     Ksat [m/d]', color= "#696969",
            transform=ax.transAxes, fontsize=7)
    ax.text(0.18, 1.24, 'Sand       {}'.format(sand[dfx]), color= "#696969",
            transform=ax.transAxes, fontsize=7)
    ax.text(0.18, 1.12,  'Loam      {}'.format(loam[dfx]), color= "#696969",
            transform=ax.transAxes, fontsize=7)
    ax.text(0.18, 1.0, 'Clay        {}'.format(clay[dfx]), color= "#696969",
            transform=ax.transAxes, fontsize=7)

ax.legend(wedges, labels_p, loc="center left", frameon=0, fancybox=0,
          bbox_to_anchor=(0.03, 2.15), fontsize=8, labelspacing=0.3)

fig.savefig(os.path.join(figures2, 'Water_balance_{}.png'.format(dfx)),
                         bbox_inches= "tight", dpi=600)
plt.close()

# In[ ]:
#==============================================================================
#================================ET COMPONENTS=================================
#==============================================================================
components = ["PM_evap", 'Surf_evap','Canopy_evap', 'PM_trans']

for i, dfx in enumerate (scenarios_ordered):
    print (i, dfx)

    df = wat_bal_dict[dfx]

    # Get the total flux amounts for each of the components of the WB
    summary = {}
    for key in df.keys():
        summary[key] = df[key].sum()


    total = summary['Canopy_evap'] + summary['Surf_evap'] + \
            summary['PM_evap'] + summary['PM_trans']

    # If a component == 0 then don't plot
    components_plot= {key: summary.get(key) for key in components}
    components_plot= {k:v for k,v in components_plot.items() if v != 0}

    # Sort keys by smaller to larger values so values dont get cover over
    components_plot= [k for k, v in sorted(components_plot.items(),
    key=lambda item: item[1])]

    for item in components_plot:
        val = (summary[item] / total * 100)
        print (item, val)

# In[ ]:
#==============================================================================
#================================ET COMPONENTS=================================
#==============================================================================
components = ["PM_evap", 'Surf_evap','Canopy_evap', 'PM_trans']

scenarios_ordered = ['Scenario8', 'Scenario1']

figx = 17/2.54  # figure width
figy = 13/2.54  # figure height
x0 = 0.17  # left edge (0-1)
y0 = 0.08   # bottom edge (0-1)
xs = 0.016 # horizontal space between panels (0-1)
ys = 0.03  # vertical space between panels (0-1)
fig = plt.figure(figsize=(figx, figy))

for i, dfx in enumerate (scenarios_ordered):
    print (i, dfx)

    df = wat_bal_dict[dfx]

    # Get the total flux amounts for each of the components of the WB
    summary = {}
    for key in df.keys():
        summary[key] = df[key].sum()
        print(key)

    total = summary['Canopy_evap'] + summary['Surf_evap'] + \
            summary['PM_evap'] + summary['PM_trans']

    # If a component == 0 then don't plot
    components_plot= {key: summary.get(key) for key in components}
    components_plot= {k:v for k,v in components_plot.items() if v != 0}

    # Sort keys by smaller to larger values so values dont get cover over
    components_plot= [k for k, v in sorted(components_plot.items(),
    key=lambda item: item[1])]

    # ---------------------------------------------------------------Left Plot
    xw1 = 0.65 # width of panel (0-1)
    yh1 = 0.4   # height of panel(0-1)
    ax = plt.axes([x0, y0+(ys*i)+(yh1*i), xw1, yh1])

    alpha_p = [1, 0.7, 0.6]

    for ix, component in enumerate(components_plot):
        colour = colours[component]
        labelx = labels[component]
        df[component].plot(ax=ax, color=colour,
                                   linewidth=0.5,
                                   ls=linesty[component],
                                   label=labelx, fontsize=8, alpha=alpha_p[ix])

    if  dfx == "Scenario1":
        ax.set_ylabel('Flux (m$^3$/d)', fontsize=10)
        ax.get_yaxis().set_label_coords(-0.15,-0.15)
    #ax.set_ylim(min_val * 10 , max_val)

    if dfx == 'Scenario8':
        ax.set_xlabel('Time (d)', fontsize=10)
    else:
        ax.set_xlabel('')
        ax.set_xticklabels('')

    ax.tick_params(direction='out', length=3, width=0.5)
    # ---------------------------------------------------------------Right Plot
    # Call the pie subplot
    xw2 = 0.2  # width of panel (0-1)
    yh2 = 0.2  # height of panel(0-1)
    ax2 = plt.axes([0.83, y0+(ys*i)+(yh1*i)+0.07, xw2, yh2])

    # Get the elements needed to call the pie chart
    values = []
    labels_p = []
    colours_p = []

    for item in components_plot:
        val = (summary[item] / total * 100)
        values.append(val)
        labels_p.append(labels[item])
        colours_p.append(colours[item])


    wedges, texts, autotext  = ax2.pie(values, colors=colours_p,
                                  startangle=90, autopct= '%.0f%%',
                                  wedgeprops={'width':0.5, 'edgecolor':'white'},
                                  pctdistance=0.75  )
    ax2.axis('equal')

    plt.setp(autotext, **{'color':'white','weight':'bold', 'fontsize':6})

    for text in autotext:
        print (text)
        if text.get_text() == '-0%':
            text.set_visible(0)

    if dfx == 'Scenario8':
        ax2.legend(wedges, labels_p, loc="center left", frameon=0, fancybox=0,
                   labelspacing=0.08, bbox_to_anchor=(0.0, -0.23), fontsize=8)

        ax.text(1.03, 0.95,  '{} (River incised)'.format(dfx),
                transform=ax.transAxes, fontsize=8)

    else:
        ax.text(1.09, 0.95,  '{}'.format(dfx),
                transform=ax.transAxes, fontsize=8)

    ax.text(1.07, 0.88,  'Soil     Ksat [m/d]', color= "#696969",
            transform=ax.transAxes, fontsize=7)
    ax.text(1.07, 0.82, 'Sand       {}'.format(sand[dfx]), color= "#696969",
            transform=ax.transAxes, fontsize=7)
    ax.text(1.07, 0.76,  'Loam      {}'.format(loam[dfx]), color= "#696969",
            transform=ax.transAxes, fontsize=7)
    ax.text(1.07, 0.70, 'Clay        {}'.format(clay[dfx]), color= "#696969",
            transform=ax.transAxes, fontsize=7)

fig.savefig(os.path.join(figures2, 'ET_components{}.png'.format(dfx)),
                             bbox_inches= "tight", dpi=600)

plt.close()

# In[]
df = wat_bal_dict['Scenario1'].loc[71.5: 77.5].copy()

fig = plt.figure(figsize=(14, 8))
ax = fig.add_subplot(1, 1, 1)

df.dt.plot(ax=ax, drawstyle='steps', color='#b61c1c')

ax2=ax.twinx()

df.rainfall.plot(ax=ax2)

ax.set_xlim(74.5, 77)
ax.set_ylim(2,6)

ax2.invert_yaxis()
ax2.set_ylim((60000, 0))
ax2.set_yticks([0, 30000, 50000])

ax.set_ylabel('Cummulative dt', fontsize=12, labelpad=20)

ax2.set_ylabel('Precipitation inputs (m$^3$/d)', fontsize=12, labelpad=20)

ax.set_xlabel('Time (d)', fontsize=12, labelpad=20)


fig.savefig(os.path.join(figures2, 'Convergence_w_rainfall_inputs.png'),
            bbox_inches = "tight", dpi=600)

# In[ ]: Plot the CUMULATIVE VALUES for the ET components

#==============================================================================
#============================CUMULATIVE ET VALUES==============================
#==============================================================================
components = ['Tot_ET', "PM_evap", 'Surf_evap', 'PM_trans']

scenarios_ordered = ['Scenario1', 'Scenario2', 'Scenario3',
                     'Scenario4', 'Scenario8']


figx = 14 # figure width in cm
figy = 10  # figure height in cm
fig = plt.figure(figsize=(figx, figy))

for i, df in enumerate (scenarios_ordered):
    print (i, df)

    # ---------------------------------------------------------------Left Plot
    ax = fig.add_subplot(5, 1, i + 1)

    for ix , component in enumerate(components):
        print (component)

        colour = colours[component]
        labelx = labels[component]

        wat_bal_dict[df]['dt'] = [wat_bal_dict[df].index.tolist()[0]] + \
                                 list(np.diff(wat_bal_dict[df].index))
        wb_ref = wat_bal_dict[df].copy()
        wb_ref = wb_ref.multiply(wat_bal_dict[df]['dt'].tolist(), axis=0)
        wb_ref[component].cumsum().plot(color=colour,lw=1, ls='-',label=labelx,
                                                      fontsize=8)

        if i == 3:
            ax.legend(loc="center left", frameon=0, fancybox=0,
                          bbox_to_anchor=(1., .45), fontsize=12, labelspacing=1.2)


    ax.set_xlabel('')
    ax.text(0.45, 0.85, df, transform=ax.transAxes,  fontsize=11,
            bbox=dict(facecolor='white', edgecolor='black', alpha=0.5,
            boxstyle='round,pad=0.15', linewidth=0.5))

# Adjust subplots spacing
left  = 0.05
right = 0.85
bottom = 0.1
top = 0.97
wspace = 0.2
hspace = 0.2
plt.subplots_adjust(left=left, bottom=bottom, right=right,
                    top=top, wspace=wspace, hspace=hspace)

# X and Y Axis labels
fig.text(0.5, 0.03, 'Time (d)',  fontsize=12)
fig.text(0.015, 0.65, 'Cumulative volumes (m$^3$)',  fontsize=12, rotation=90)

# Save figure
fig.savefig(os.path.join(figures2, 'Cum_ET_comp{}.png'),
            bbox_inches = "tight", dpi=600)
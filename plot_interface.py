import glob
import os

import matplotlib.pyplot as plt
import numpy as np

"""
paraview filter Plot over line
"""

from parameter import *
result_type = '.csv'
result_folder = 'results/'
#os.chdir(result_folder)
filename_root = '*'


def add_subplot_axes(ax,rect,axisbg='w'):
    #rect: [xmin, ymin, width_ratio, height ratio]
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)    
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    subax = fig.add_axes([x,y,width,height],axisbg=axisbg)
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax

fig = plt.figure()
#ax1 = fig.add_subplot(111)
ax1 = plt.gca()

"""
os.chdir(result_folder)
files = glob.glob( filename_root + result_type)
xlabel =
for f in files:
    names = f[:-len(result_type)].split('__')
    param_value = float(names[-1])
    param_name = names[-2].replace('_', ' ')
    label = param_name + ' = ' + str(param_value)
    Tvalues = np.loadtxt(f)
    plot(Tvalues, label = label)
"""
length_scale = friction_heat_length
param_name = "Chip distance from the tool tip / friction heat length"
f = result_folder + 'T_case2_profile.csv'
data = np.loadtxt(f, delimiter = ',', skiprows=1)
T = data[:, 1]
x = data[:, 0] / length_scale

ax1.plot(x, T, 'r', label = "uniform")
xlim = plt.gca().get_xlim()


if nonuniform_friction_heat_configuration == 1:
    f = result_folder + "T_case2_linear_with_cutter_material.csv"
    Tlim = [0, 1600] #plt.gca().get_ylim()
    extra_label = "constant thermal properties"
else:
    f = result_folder + 'T_case2_2d_profile.csv'
    Tlim = [0, 1400] #plt.gca().get_ylim()
    extra_label = "2D geomtry"

data = np.loadtxt(f, delimiter = ',', skiprows=1)
T = data[:, 1]
x = data[:, 0] / length_scale
ax1.plot(x, T, 'g--', label = extra_label)

#######################################

#assert case_id = 2, "check parameter.py for setting"
assert nonuniform_friction_heat, "check parameter.py for setting"

# get zone line
x_start = friction_heat_start / length_scale
x_uniform_end = uniform_friction_heat_length / length_scale
x_end = actual_friction_heat_length / length_scale
xh = np.array([0, x_uniform_end, x_end]) + x_start
qh = np.array([uniform_heat_ratio, uniform_heat_ratio, 0])
max_q_ratio = uniform_heat_ratio * 1.1


if nonuniform_friction_heat_configuration == 1:
    f2= result_folder + "T_case2_nonuniform_profile.csv"
else:
    f2= result_folder + "T_case2_nonuniform2_profile.csv"

data2 = np.loadtxt(f2, delimiter = ',', skiprows=1)
T2 = data2[:, 1]
x2 = data2[:, 0] / length_scale
ax1.plot(x2, T2, '-.b', label = "nonuniform")

rect = [0.4, 0.05, 0.3, 0.3]
ax2 = add_subplot_axes(ax1,rect)

ax2.set_ylim([0,max_q_ratio])
ax2.plot(xh, qh)
ax2.set_title('non-uniform distribution \n of friction heating')
#ax2.set_xlabel('')
ax2.set_ylabel('ratio of heat load')

ax1.vlines(xh[0], *Tlim,  linestyles = 'dashed', label = "start of friction heating")
ax1.vlines(xh[1], *Tlim,  linestyles = 'dashed', label = "end of uniform heating")
ax1.vlines(xh[-1], *Tlim,  linestyles = 'dashed', label = "end of nonuniform heating")
ax1.set_xlim([0, 5.5])

ax1.set_xlabel(param_name)
ax1.set_ylabel(r'Interface temperature ($\degree$C)')
ax1.set_title("Tool-chip interface temperature profile")
leg = ax1.legend()
if leg:
    leg.draggable()

saved_filename = result_folder + "interface_temperature_profile_{}.png".format(nonuniform_friction_heat_configuration)
plt.savefig(saved_filename)
plt.show()


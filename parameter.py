from __future__ import print_function, division
import sys
from math import cos, sin, atan, pi
from parameter_case import *  #default parameter could be overwritten by sys.argv
from material import *

validation_datafile = "validation.csv"
is_batch_mode = True
using_debug = False
using_salome = True
is_preprocessing = True # for non-parallel, it can be done all togother

# salome can not parse the parameter this way!
if not using_salome:
    if  (len(sys.argv)>=3):
        print("=============argv================")
        print(sys.argv)  # mpirun will be ignored,   only arg after python are captured

        is_batch_mode = True
        param_name = sys.argv[1]
        param_value = float(sys.argv[2])
        if param_name == "case_id":
            case_id = param_value
            case_id_supplied = True
        elif param_name == "cutting_speed":
            cutting_speed = param_value
            print("!!! cutting_speed !!! is provided by command line")
        elif param_name == "C_0":
            C_0 = param_value
        elif param_name == "delta":
            delta = param_value
        else:
            print("parameter {} is not recognized".format(param_name))
            sys.exit()
        print("==== recoginised parameter name {} and value {} =====".format(param_name, param_value))

    if  (len(sys.argv)>3):
        result_filename =  sys.argv[3]

def defined(x):
    return x in locals() or x in globals()


def get_shear_angle(feed_thickness, chip_thickness, cutter_angle_v):
    #chip_thickness/(cos(a_phi)*cos(a_v) + sin(a_v)*sin(a_phi)) = feed_thickness/sin(a_phi)
    a = sin(cutter_angle_v*pi/180)
    x = cos(a) / (chip_thickness/feed_thickness - sin(a))
    shear_angle = atan(x) /pi * 180
    return shear_angle

if not defined('case_id_supplied'):
    case_id_supplied = True; case_id = 2 # default to case 2, which has the best matching for exp and analytical solution

###############################################
result_folder = './output/'
datafile_root = 'result'
validation_tol = 0.02  # should be less than 2.5%
extracting_data = False
using_3D = False
using_2D = not using_3D
using_workpiece_extra_extusion = False and using_3D
#if using_2D:
#    cutter_thickness = 1  # overriding 

exporting_foam = False

#using_mapping_bc = 1;   # otherwise fillet
using_mapping_bc = using_salome and True
if using_mapping_bc:
    using_double_shear_heat_layer = True
else:
    using_double_shear_heat_layer = False

using_chip_cutter_mapping_bc = not has_friction_transition_zone # keep it False, it does not work quite as expected
using_chip_cutter_transition = has_friction_transition_zone  # another name define in material.py

has_convective_velocity = True
considering_radiation = False # neglectible impact on interface temperature

element_degree = 2  # will second order help easing velocity unsmoothing,  yes!
#using_nonlinear_thermal_properties = False  # controlled in material.py
using_stab = False # stabilization causes less mass heat flux, while no stab ; totally diff contour
IP_alpha = 2


####################################################
using_experimental_data = True
result_filename =  'results/results_testing.csv'
#result_filename =  'results/results_speed_linear_exp.csv'
#result_filename =  'results/results_speed_nonlinear_exp.csv'
#result_filename =  'results/results_speed_nonlinear_tube_exp.csv'
#result_filename =  'results/results_delta_nonlinear_exp.csv'
#result_filename =  'results/results_C_0_nonlinear_exp.csv'
#result_filename =  'results/results_2D_nonlinear_exp.csv'
#result_filename =  'results/results_validation_nonlinear_exp.csv'
#result_filename =  'results/results_validation_nonlinear.csv'

#if "meshing_only" in sys.argv or "preprocessing" in sys.argv:
#    is_preprocessing = True

#replace `double#doubleZZdouble##` line to parameter by sed in batch mode
##ZZ##

###############################################
if defined('case_id_supplied') and case_id_supplied:
    import csv
    input = open(validation_datafile)
    reader = csv.DictReader(input)
    dict_list = []
    for line in reader:
        dict_list.append(line)
    input.close()
    rd = dict_list [int(case_id)-1]
    for k in rd:
        rd[k] = float(rd[k])
    #print(rd)
    """
    import pandas as pd  # salome's python2 has no pandas
    d = pd.read_csv(validation_datafile)
    print(d.index)
    rd = d.iloc[case_id-1] # loc[] return series, loc() return iterator
    print(rd, type(rd), rd.index)
    """
    assert rd['case_id'] == case_id
    #print(rd['feed_thickness'], type(rd['feed_thickness']))
    feed_thickness = rd['feed_thickness'] * 0.001
    k_AB = rd['k_AB'] * 1e6
    #cutting_speed = rd['cutting_speed'] / 60.0
    #chip_friction_distance = rd['h'] * 0.001
    #= rd['']
    cutter_angle_v = rd['cutter_angle_v']  #deg
    chip_thickness = rd['chip_thickness']  * 0.001

    #shear_angle = get_shear_angle(feed_thickness, chip_thickness, cutter_angle_v)
    if using_experimental_data:
        F_cutting = rd['F_c_exp']
        F_thrush = rd['F_t_exp']
    else:
        F_cutting = rd['F_cutting']
        F_thrush = rd['F_thrush']

    shear_angle = get_shear_angle(feed_thickness, chip_thickness, cutter_angle_v)
    #print("calculated shear angle in deg = ", shear_angle)
    l_AB = feed_thickness/sin(shear_angle*pi/180)
    if not defined('chip_friction_distance'):
        chip_friction_distance = rd['h'] * 0.001
    if not defined('cutting_speed'):
        cutting_speed = rd['cutting_speed'] / 60.0
    if not defined('C_0'):
        C_0 = rd['C_0']
    shear_heat_thickness = l_AB / C_0
    if not defined('delta'):
        delta = rd['delta']
    friction_heat_thickness = chip_thickness * delta
    print("friction_heat_thickness, shear_heat_thickness", friction_heat_thickness, l_AB, shear_heat_thickness)
else:
    #print( feed_thickness, chip_thickness, cutter_angle_v, shear_angle, cutting_speed, F_cutting, F_thrush)
    shear_angle = get_shear_angle(feed_thickness, chip_thickness, cutter_angle_v)
    print("calculated shear angle in deg = ", shear_angle)
    l_AB = feed_thickness/sin(shear_angle*pi/180)


###############################################
shear_heat_volume = l_AB * shear_heat_thickness * cutter_thickness # m3, ignore angle
if using_double_shear_heat_layer:
    friction_heat_start = shear_heat_thickness/2.0/sin((90-cutter_angle_v+shear_angle)*pi/180)  #double layers
else:
    friction_heat_start = shear_heat_thickness/sin((90-cutter_angle_v+shear_angle)*pi/180)
friction_heat_length = chip_friction_distance - friction_heat_start
pressing_heat_volume = pressing_heat_width * pressing_heat_thickness * cutter_thickness
chip_volume = chip_thickness * cutter_thickness * chip_length # rough value  for a straight chip


# uniform at the beginning, with higher density, then drop linearly
nonuniform_friction_heat = False  # currently assume nonuniform_friction_thickness == False
nonuniform_friction_heat_configuration = 1
if nonuniform_friction_heat:
    if nonuniform_friction_heat_configuration == 1:
        uniform_length_ratio = 0.5  # total heating length increased, but the heat ratio kept
        uniform_heat_ratio = 1.0
        actual_friction_heat_length = friction_heat_length * ( uniform_length_ratio + (1.0 - uniform_length_ratio) * 2)
    else:
        uniform_length_ratio = 0.3333333333  # total heating length increased, but the heat ratio kept
        uniform_heat_ratio = 1.5
        actual_friction_heat_length = friction_heat_length
    actual_friction_heat_distance = friction_heat_start + actual_friction_heat_length
    uniform_friction_heat_length= friction_heat_length * (uniform_length_ratio)  # assuming linear drop to zero
else:
    actual_friction_heat_length = friction_heat_length


nonuniform_friction_thickness = True
if nonuniform_friction_thickness:
    friction_end_thickness = 5e-6
    actual_friction_heat_thickness = (friction_heat_thickness + friction_end_thickness) * 0.5
    friction_heat_volume = actual_friction_heat_length * actual_friction_heat_thickness * cutter_thickness
else:
    friction_heat_volume = actual_friction_heat_length * friction_heat_thickness * cutter_thickness

nominal_friction_heat_volume = friction_heat_length * friction_heat_thickness * cutter_thickness
print("calculated volume for shear, friction and chip:", shear_heat_volume, friction_heat_volume, chip_volume)

"""
uniform_length_ratio = 0.3
uniform_friction_heat_distance = friction_heat_start + friction_heat_length * uniform_length_ratio
uniform_heat_ratio = 2.0/(1 + uniform_length_ratio)
"""

###############################################
angle_phi = shear_angle;  # alias name
holder_angle_v = cutter_angle_v+5;

#chip_sliding_distance = chip_friction_distance  # velocity expression
tool_chip_interface_length = chip_friction_distance # tool chip interface length
assert tool_chip_interface_length*1.5<chip_length
chip_speed = cutting_speed * (feed_thickness/chip_thickness)
#if workpiece is a disc
frequency = cutting_speed / disc_R / 2.0 / pi
direction = -1
omega = cutting_speed / (disc_R+0.5*feed_thickness)

#################################################




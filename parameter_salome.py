from __future__ import print_function, division
from parameter import *
from math import cos, sin, tan, pi, fabs

#fillet_r = feed_thickness * 0.25
#a_phi = shear_angle * pi/180

#change geometry tolerance to 1e-9
#this file will not reload until salome program restart

salome_mesh_output = '/tmp/Mesh_1.med'
mapping_bc_gap = 5e-6 # should be small after debug

friction_mapping_gap_volume = mapping_bc_gap * cos(cutter_angle_v*pi/180) * actual_friction_heat_length * cutter_thickness
print("friction_mapping_gap_volume = ", friction_mapping_gap_volume)

# 0 is rectangle section, 1 disc, 2 work is very thin, 3 for mapping BC
using_fillet_shear_zone = not using_mapping_bc

extrusion_thickness = cutter_thickness #
workpiece_extra_extusion_thickness = cutter_thickness * 1
workpiece_min_z = - workpiece_extra_extusion_thickness

pressing_heat_x = pressing_heat_width

# sometime, quad is not divided into tri
using_quad_mesh = False   # failed to convert all cell into tri
#and not using_nonlinear_thermal_properties  # quad at chip as possible, but not convergent for nonliner model

if using_3D:
    line_segment_nb_default = 20
    line_segment_nb_small = 2  #
else:
    line_segment_nb_default = 30  # case 6,8  can not converge
    line_segment_nb_small = 1  # 3 make low quality cells which may lead to divergence
line_segment_nb_vertical = 2

# workpiece_type_id defined in parameter_case.py
# if workpiece is not a full disc (type 1)
work_left_x = feed_thickness * -5
work_right_x = feed_thickness * 5  # > 5 , does not affect Tmax
if workpiece_type_id == 2:  # no bottom region, only chip
    work_bottom_y = 0
elif workpiece_type_id == 0:  # rectangle workpiece as in this paper
    work_bottom_y = feed_thickness * -5  # insentitive to this value, fixed as -5
    # orthagonal cutting, tube,  adding extra tube thickness, will affect heat parittion
else:
    pass

#chip_shear_width = math.sqrt(chip_start_x*chip_start_x + chip_start_y*chip_start_y); # == l_AB

p_chip_turning_x = - l_AB * cos(shear_angle *pi/180)
p_chip_turning_y = feed_thickness

#
chip_start_x, chip_start_y = p_chip_turning_x, p_chip_turning_y

#if using_fillet_shear_zone:  #used in velocity_expression
fillet_angle = 90 - cutter_angle_v
fillet_r = shear_heat_thickness / (cos(angle_phi*pi/180)*2*sin(fillet_angle/2.0*pi/180))   # shearing_layer_thickness

cutting_fillet_hole = False
_ratio = (pi/4 + cutter_angle_v*pi/360) - sin(fillet_angle/2.0*pi/180) * cos(fillet_angle/2.0*pi/180)
hole_volume =  cutter_thickness * fillet_r * fillet_r * (1 - sin(cutter_angle_v*pi/360) / 2.0 - _ratio)

fillet_shift_x = - fillet_r*tan(fillet_angle/2*pi/180)
fillet_shift_y = fillet_r
fillet_p2_shift_x = -fillet_shift_x * sin(cutter_angle_v*pi/180)
fillet_p2_shift_y = -fillet_shift_x * cos(cutter_angle_v*pi/180)

turning_center_x = p_chip_turning_x + fillet_shift_x
turning_center_y = p_chip_turning_y + fillet_shift_y

_ll = shear_heat_thickness/cos((shear_angle - cutter_angle_v)*pi/180)
if using_double_shear_heat_layer:
    _ll = _ll/2
    work_shear_thickness_shift_x = shear_heat_thickness/2.0/sin(angle_phi*pi/180) * -1
shear_thickness_shift_x = _ll * sin(cutter_angle_v*pi/180)
shear_thickness_shift_y = _ll * cos(cutter_angle_v*pi/180)
#shear_thickness_distance = shear_heat_thickness / sin((90-cutter_angle_v-angle_phi)*pi/180)
#
if using_3D:
    z_mid = cutter_thickness*0.5
else:
    z_mid = 0
# probe point 
point_shear_surface_center_coordinates = [shear_thickness_shift_x-0.5*l_AB*cos((shear_angle)*pi/180), 
                            shear_thickness_shift_y + fabs(0.5*l_AB*sin((shear_angle)*pi/180)), z_mid]

# heat zone is inside chip, it is 2.5% higher?
_l = friction_heat_thickness / cos((shear_angle - cutter_angle_v)*pi/180)
friction_b_shift_x = _l * cos(( angle_phi)*pi/180) * -1  # bottom
friction_b_shift_y = _l * sin(( angle_phi)*pi/180)
if not nonuniform_friction_thickness:  # assume rect heating zone, 
    friction_t_shift_x = friction_b_shift_x
    friction_t_shift_y = friction_b_shift_y
else:  # top end has a smaller thickness, which is several um
    _l = friction_end_thickness / cos((shear_angle - cutter_angle_v)*pi/180)
    friction_t_shift_x = _l * cos(( angle_phi)*pi/180) * -1  # top
    friction_t_shift_y = _l * sin(( angle_phi)*pi/180)

if nonuniform_friction_heat:
    chip_ccc_x = (actual_friction_heat_distance) * sin(cutter_angle_v*pi/180)
    chip_ccc_y = (actual_friction_heat_distance) * cos(cutter_angle_v*pi/180)
else:
    chip_ccc_x = (chip_friction_distance) * sin(cutter_angle_v*pi/180)
    chip_ccc_y = (chip_friction_distance) * cos(cutter_angle_v*pi/180)

#print(chip_ccc_x, chip_ccc_y)
chip_top_x = chip_length * sin(cutter_angle_v*pi/180)
chip_top_y = chip_length * cos(cutter_angle_v*pi/180)
print(chip_top_x, chip_top_y)
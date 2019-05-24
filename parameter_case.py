############### could be overwritten #############
"""
case_id = 2
feed_thickness = 0.00015;

chip_thickness = 0.00031;  # also feed speed  meter per rov, it is a derived value
shear_angle = 25.1659836334;

chip_friction_distance = 0.0006;  # lp <  sliding contact length
chip_sliding_distance = 0.0006;  #  h calc from equation (16)

#cutting force measurement data?
cutter_angle_v = 5;  # degree
cutting_speed =200/60.0;

F_cutting = 433;
F_thrush = 171;  # ivester case 2

delta = 0.055  # case 2,   delta = S1 / chip_thickness
C_0 = 4.2
k_AB = 5.84e8  # shearing yielding stress, depends on temperature

shear_heat_thickness = 4e-5;  # shearing_layer_thickness
friction_heat_thickness = 1e-5;
"""

##############################################
# 0 is rectangle section,  5X5 feed_thickness,  using_extra_extrusion_thickness
# 1 disc, full disc
# 2 work is very thin, only chip thickness
# 3 extra work thickness, to simulate 
workpiece_type_id = 0; 
is_straight_chip = 1;

#later read from a table csv file, default to case 2
cutter_thickness = 0.0016;  #kasper data, b=0.004

##################################

pressing_heat_thickness = 2e-6;
pressing_heat_width = 1e-5;
has_pressing_zone = 0;


##################################

# shearing calc from cutting force, * beta (work -> heat conversion coeff)
# friction_force is known, but CoF is not known
#F_cutting = 596;
#F_thrush= 324;  # jasper case 4

# http://www.carbidedepot.com/CTAPR123B-P78366.aspx,  8-10 mm
cutter_H = 0.005;  # half inch
cutter_W =0.005;  # half inch

cutter_angle_h = 5;  #degree, not used in this program yet

holder_H = 0.012;  # 0.75 inch, 
holder_W =0.012;
holder_thickness = cutter_thickness;
holder_top_y = holder_H;

# even rect work is part of disc, used as length scale
disc_R = 1e-1;
disc_thickness = cutter_thickness;  #same with cutter now

#this is should be set
#chip_length = 0.001;  # it should be defined later
# only if chip is an arc
chip_end_angle = 90;  # relative to chip arc center
chip_radius = disc_R;  #avg of innner and outer radius


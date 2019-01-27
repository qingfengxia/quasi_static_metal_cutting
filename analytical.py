from __future__ import print_function, division
import math
from math import cos, sin, pi
from parameter import *

# Lalwani 2009
# ============== other control parameter===============


#different way to cal
def get_friction_heat_source():
    angle_F = math.atan2(F_thrush, F_cutting)*180/math.pi # atan2(y, x)
    angle_lambda = angle_F + cutter_angle_v
    angle_theta = angle_phi + angle_F
    F_reaction = (F_cutting*F_cutting + F_thrush*F_thrush)**0.5
    _tmp = (angle_theta - angle_phi + cutter_angle_v) * math.pi/180
    F_secondary = F_reaction * math.sin(_tmp)  # sign !
    print(angle_theta, angle_F)
    F_AB_1 = F_reaction * math.cos(angle_theta*math.pi/180)  # shearing force, not correct
    print("F_AB_1, F_secondary", F_AB_1, F_secondary)
    friction_heat = (chip_speed * F_secondary)  #  * taylor_quinn_coeff
    return friction_heat

def get_shear_heat_source(angle_phi, k_AB):
    l_AB = feed_thickness/math.sin(shear_angle*math.pi/180)
    F_AB = k_AB  * cutter_thickness * l_AB
    #V_AB equation Fig 3 Lalvani 2009
    V_AB = cutting_speed * math.cos(cutter_angle_v*math.pi/180) / math.cos((angle_phi - cutter_angle_v)*math.pi/180)
    shear_heat = (F_AB * V_AB) * taylor_quinn_coeff
    return shear_heat

def get_heat_source():
    shear_heat = get_shear_heat_source(angle_phi, k_AB)
    friction_heat = get_friction_heat_source()
    print('heat generated from primary zone and secondary zone', shear_heat, friction_heat)
    return shear_heat, friction_heat

mass_flow_rate = metal_density * cutting_speed * feed_thickness * cutter_thickness
#thermal_number = metal_density * metal_cp_0 * cutting_speed * feed_thickness / metal_k_0

def get_heat_partition(thermal_number):
    cond = thermal_number*math.tan(angle_phi*math.pi/180.0)
    if cond > 10:
        heat_partition_work = 0.3 - 0.15*math.log10(cond)
    else:
        heat_partition_work = 0.5 - 0.35*math.log10(cond)
    #print("heat_partition_work = ", heat_partition_work)
    return heat_partition_work

T_tol = 0.1
relaxation_coeff = 0.3
def get_T_AB(shear_heat):
    #which reference temperature should be used to calc material property?
    T_AB = T_reference
    T_AB_new = T_AB + 10
    while math.fabs(T_AB_new - T_AB) > T_tol:
        # temperature loop
        T_AB = T_AB + (T_AB_new - T_AB ) * relaxation_coeff  # relaxation_coeff
        Cp = metal_cp_f(T_AB)
        K = metal_k_f(T_AB)

        thermal_number = metal_density * Cp * cutting_speed * feed_thickness / K
        heat_partition_work = get_heat_partition(thermal_number)
        delta_T_shear_zone = shear_heat * (1-heat_partition_work) / (mass_flow_rate * Cp)

        T_AB_new = T_reference + delta_T_shear_zone
    return T_AB

def get_T_int(T_AB, friction_heat):
    #thermal_number = metal_density * metal_cp_0 * cutting_speed * feed_thickness / metal_k_0

    T_chip = T_AB
    T_chip_new = T_chip + 10
    while math.fabs(T_chip_new - T_chip) > T_tol:
        T_chip = T_chip + (T_chip_new - T_chip ) * relaxation_coeff
        Cp = metal_cp_f(T_chip)
        K = metal_k_f(T_chip)
        thermal_number = metal_density * Cp * cutting_speed * feed_thickness / K

        delta_T_friction_zone = friction_heat / (mass_flow_rate * Cp)
        T_chip_new = T_AB + delta_T_friction_zone
        #print("T_chip loop:", T_AB, T_chip, delta_T_friction_zone,  thermal_number)
        ksi = 0.9
        power10 = 0.06 - 0.195 *  delta * math.sqrt(thermal_number * chip_thickness/feed_thickness) \
                        + 0.5 * math.log10(thermal_number * chip_thickness/tool_chip_interface_length)
        #print("math.pow(10, power10)", math.pow(10, power10))

    # end of T_chip loop
    T_int= T_AB + ksi * math.pow(10, power10) * delta_T_friction_zone   # avg interface temperature
    return T_int

#print("delta_T_shear_zone, delta_T_friction_zone", delta_T_shear_zone, delta_T_friction_zone)

# there is another paper to predict this temperature, List. G et al, 2012
"""
l_c = tool_chip_interface_length
w = cutter_thickness
A_bar = 2/math.pi * (math.log(2*w/l_c) + w/l_c/3 + 0.5)
R_s = 1.0/(1.0 + 0.754*((k_ratio)/A_bar)/math.sqrt( metal_diffusivity_0 / (chip_speed*l_c)) )  #partition coeff
delta_T_friction_zone_1 = ksi*(1-R_s) * shear_heat/(l_c*w)  *l_c * A_bar / (2*metal_k_0)
T_friction_interface_1 = T_ambient + delta_T_shear_zone + delta_T_friction_zone_1
print("A_bar, R_s, delta_T_friction_zone_1 (alt) = ", A_bar, R_s, delta_T_friction_zone_1 )
"""

def get_analytical_T():
    print('============= analytical prediction =============')
    Re = omega * disc_R * disc_R / air_kinematic_viscosity
    #print(metal_k_0/(metal_cp_0*metal_density))
    Pe = omega * disc_R * disc_R / (metal_k_0/(metal_cp_0*metal_density))
    print("At the cutting position: Re = {}, Pe = {}".format(Re, Pe))
    #print("T_aB/delta_T_friction_zone = ", delta_T_friction_zone)
    shear_heat, friction_heat = get_heat_source()
    T_shear_interface = get_T_AB(shear_heat)
    T_friction_interface = get_T_int(T_shear_interface, friction_heat)
    print("T_shear_interface = ", T_shear_interface)
    print("T_friction_interface = ", T_friction_interface)
    return T_shear_interface, T_friction_interface

get_analytical_T()
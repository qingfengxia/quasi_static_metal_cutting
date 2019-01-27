# Lalvani et al. 2009
import math
from math import cos, sin, tan, pi, fabs
import numpy as np

#from parameter import *
#from material import *

############### could be overwritten #############
#case_id = 2
feed_thickness = 0.00015;
cutter_thickness = 0.0016;  #kasper data, b=0.004
#chip_thickness = 0.00031;  # also feed speed  meter per rov, it is a derived value
cutter_angle_v = 5;  # degree
cutting_speed =200/60.0;

#######################################
#material 1045 steel
T_reference = 25
taylor_quinn_coeff = 0.9  # 1045 steel
#conductivity = 420+ 0.504*T, cp = 52.21 -0.0281*T,  density = 8000 constant
metal_density = 8000.0
k_ratio = 1.0  #tool to work thermal diffusivity
metal_cp_0 = 420.0
metal_k_0 = 52.0
metal_capacity_0 = metal_cp_0*metal_density
metal_diffusivity_0 = 52.21/metal_capacity_0

metal_cp_f = lambda T: (420 + 0.504*T)
metal_k_f = lambda T: (52.21 -0.0281*T)
metal_capacity_f = lambda T: (420*8000 + 0.504*8000*T)
metal_diffusivity_f = lambda T: (52.21/8000/420 - 0.0281/8000/420*T)

def strain_hardening_coeff(epsilon):
    #1045 steel
    A = 5.5531e8
    B = 6.008e8
    n = 0.234
    return n*B * epsilon**n / (A + B * epsilon**n)

def normal_stress(epsilon, epsilon_dot, T):
    #1045 steel
    A = 5.5531e8
    B = 6.008e8
    n = 0.234
    C = 0.0134
    m = 1
    T_m = 1460  # Celsius degree
    epsilon_dot_0 = 1# reference strain, Material behaviour in conditions similar to metal cutting: flow stress in the primary shear zone
    T_factor = (1 - ((T - T_reference)/(T_m - T_reference))**m)
    return (A + B * epsilon**n) * (1 + C*math.log(epsilon_dot/epsilon_dot_0)) * T_factor

mass_flow_rate = metal_density * cutting_speed * feed_thickness * cutter_thickness  # cosnt for a specific cutting case
#thermal_number = metal_density * metal_cp_0 * cutting_speed * feed_thickness / metal_k_0

delta_0, C_0_0, angle_phi_0 = 0.055, 4.2, 25.8 # case 2
_debug = False
T_tol = 0.2
nc = 20
T_w = T_reference
# this loop will calc chip_thickness and angle_phi

def calc(delta, C_0, angle_phi, final = False):
    l_AB = feed_thickness/sin(angle_phi*pi/180)
    V_AB = cutting_speed * cos(cutter_angle_v*pi/180) / cos((angle_phi - cutter_angle_v)*pi/180)

    gamma_AB = 0.5 * cos(cutter_angle_v*pi/180) / cos((angle_phi - cutter_angle_v)*pi/180)/sin(angle_phi*pi/180)

    epsilon_dot_AB = C_0 * V_AB / l_AB / (3**0.5) # avg strain-rate, C_0 is strain rate constant
    epsilon_AB = gamma_AB/(3**0.5)

    T_AB = T_w
    T_AB_new = T_AB + 10
    relaxation_coeff = 0.3
    while fabs(T_AB_new - T_AB) > T_tol:
        # temperature loop
        T_AB = T_AB + (T_AB_new - T_AB ) * relaxation_coeff  # relaxation_coeff
        Cp = metal_cp_f(T_AB)
        K = metal_k_f(T_AB)

        sigma_normal = normal_stress(epsilon_AB, epsilon_dot_AB, T_AB)
        k_AB = sigma_normal/(3.0**0.5)
        F_AB = k_AB * cutter_thickness * l_AB

        shear_heat = (F_AB * V_AB) * taylor_quinn_coeff
        thermal_number = metal_density * Cp * cutting_speed * feed_thickness / K
        cond = thermal_number*math.tan(angle_phi*math.pi/180.0)
        #rint("T_AB loop:", T_AB, thermal_number, cond)
        if cond > 10:
            heat_partition_work = 0.3 - 0.15*math.log10(cond)
        else:
            heat_partition_work = 0.5 - 0.35*math.log10(cond)
        delta_T_shear_zone = shear_heat * (1-heat_partition_work) / (mass_flow_rate * Cp)

        T_AB_new = T_w + delta_T_shear_zone
    #loop to find T_AB

    chip_thickness = feed_thickness * cos((angle_phi-cutter_angle_v)*math.pi/180) / sin(angle_phi*math.pi/180)
    chip_speed = cutting_speed * (feed_thickness/chip_thickness)
    n_eq = strain_hardening_coeff(epsilon_AB)
    tan_theta = 1 + 2 * (pi/4 - angle_phi/180*pi)  - C_0*n_eq
    angle_theta = math.atan(tan_theta) * 180 / pi

    angle_F = angle_theta - angle_phi
    angle_lambda = angle_theta - angle_phi + cutter_angle_v

    F_reaction = F_AB / cos(angle_theta*pi/180)
    F_secondary = F_reaction * math.sin(angle_lambda*pi/180)   # frictional force on tool_chip
    friction_heat = (F_secondary * chip_speed) * taylor_quinn_coeff  # different from that paper
    F_normal = F_reaction * math.cos(angle_lambda*pi/180)
    F_cutting = F_reaction * cos(angle_F * pi/180)
    F_thrush = F_reaction * sin(angle_F * pi/180)

    tool_chip_interface_length = feed_thickness * sin(angle_theta*math.pi/180) \
        / (cos(angle_lambda*math.pi/180) * sin(angle_phi*math.pi/180)) \
        * (1 + C_0 * n_eq/3/(1 + 2*(pi/4 - angle_phi*math.pi/180) - C_0 * n_eq)) # h Eq (18)
    #print(chip_thickness, tool_chip_interface_length)
    #print(angle_theta, F_secondary)
    assert tool_chip_interface_length > 0 and F_reaction>0

    T_chip = T_AB
    T_chip_new = T_chip + 10
    while fabs(T_chip_new - T_chip) > T_tol:
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

    tau_int = F_secondary / (tool_chip_interface_length * cutter_thickness)
    gamma_M = tool_chip_interface_length/(delta * chip_thickness)# total total shear strain,  List 2012 has diff equation min(Cof*normal_stress, )
    gamma_int = 2*gamma_AB + 0.5 * gamma_M  #max shear strain at the frictional interface
    epsilon_int = gamma_int / (3.0**0.5)
    epsilon_dot_int = chip_speed / (delta * chip_thickness) / 3.0**0.5
    k_int = normal_stress(epsilon_int, epsilon_dot_int, T_int) / (3.0**0.5)

    sigma_normal = F_normal / (tool_chip_interface_length*cutter_thickness)
    sigma_normal_1 = k_AB * (1 + pi/2 - 2*cutter_angle_v*pi/180 - 2*C_0*n_eq)  # Eq 16
    sigma_diff = fabs(sigma_normal_1 - sigma_normal)

    k_diff = fabs(k_int - tau_int)

    if final:
        h = tool_chip_interface_length
        return [delta, C_0, chip_thickness, epsilon_AB, epsilon_dot_AB, epsilon_int, epsilon_dot_int, T_AB, T_int, k_AB, n_eq, angle_theta, sigma_normal_1, h, F_cutting, F_thrush, angle_phi]
    else:
        return fabs(F_cutting), sigma_diff, k_diff

def run(delta_, C_0_, angle_phi_, relative_step = 0.005):
    results = []
    delta_tol = fabs(delta_) * relative_step
    C_0_tol = fabs(C_0_) * relative_step
    angle_phi_tol = fabs(angle_phi_) * relative_step

    delta_range = np.arange(delta_-nc*delta_tol, delta_+nc*delta_tol, delta_tol)
    C_0_range = np.arange(C_0_-nc*C_0_tol, C_0_+nc*C_0_tol, C_0_tol)
    angle_phi_range = np.arange(angle_phi_-nc*angle_phi_tol, angle_phi_+nc*angle_phi_tol, angle_phi_tol)

    delta_final = delta_
    F_cutting_min = 1e10 # very big number
    for delta in delta_range: # delta = friction_heat_thickness / chip_thickness
        C_0_final = C_0_
        sigma_diff_min = 1e10
        for C_0 in C_0_range:              #C_0 = l_AB/shear_heat_thickness
            angle_phi_final = angle_phi_
            k_diff_min = 1e10
            for angle_phi in angle_phi_range:
                if _debug: print("=== looping: ", delta, C_0, angle_phi, " === ")
                #
                F_cutting, sigma_diff, k_diff = calc(delta, C_0, angle_phi)
                if k_diff < k_diff_min:
                    k_diff_min = k_diff
                    angle_phi_final = angle_phi
                results.append([delta, C_0, angle_phi, F_cutting, sigma_diff, k_diff])
                #print("=== looping: ", delta, C_0, angle_phi, " final: === ", delta_final, C_0_final, angle_phi_final)
            #end of angle_phi loop, minimum k_diff
            F_cutting, sigma_diff, k_diff = calc(delta, C_0, angle_phi_final)
            if sigma_diff < sigma_diff_min:
                sigma_diff_min = sigma_diff
                C_0_final = C_0
        #end of C_0 loop
        F_cutting, sigma_diff, k_diff = calc(delta, C_0_final, angle_phi_final)
        if F_cutting < F_cutting_min:
            F_cutting_min = F_cutting
            delta_final = delta
        #finalise delta vlue when F_cutting is minimum
    #end of delta loop
    #ret = sorted(results, key=lambda d: d[3]*d[4]*d[5])
    #print(ret[-1])
    print("min(F_cutting_min, sigma_diff_min, k_diff_min)", F_cutting_min, sigma_diff_min, k_diff_min)
    return delta_final, C_0_final, angle_phi_final

if False:
    delta_final, C_0_final, angle_phi_final = run(delta_0, C_0_0, angle_phi_0)
    results = calc(delta_final, C_0_final, angle_phi_final, final=True)
    print(results)
else:
    import scipy.optimize as optimize
    def f(params):
        delta, C_0, angle_phi = params
        res = calc(delta, C_0, angle_phi, final = False)
        return res[0] * res[1] * res[2]
    initial_guess = [0.052, 4.2, 25.8]
    bounds = [(0.035, 0.75), (3,6), (15, 45)]
    result = optimize.minimize(f, initial_guess, method = 'TNC', bounds = bounds)
    if result.success:
        fitted_params = result.x
        print(fitted_params)
    else:
        raise ValueError(result.message)

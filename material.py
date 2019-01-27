using_nonlinear_thermal_properties = True
using_nonlinear_loop = False # not using_nonlinear_thermal_properties  # controlled in material.py
has_friction_transition_zone = True
has_pressing_zone = False
###################################
# work piece tube or disc, define in gmsh *.geo file
default_subdomain_id = 0

work_subdomain_id = 1
cutter_subdomain_id = 2
holder_subdomain_id = 3
chip_subdomain_id = 4  # line + arc,  once this value changes, velocity cpp code needs updated!
shear_subdomain_id = 5  # primary zone
friction_subdomain_id = 6  # secondary zone
subdomain_number = 7  # holder  is in used, inc zero as default

# cutter_chip_transition, different thermal properties to model thermal resistance
if has_friction_transition_zone:
    friction_transition_subdomain_id = subdomain_number
    subdomain_number  += 1
if has_pressing_zone:
    pressing_subdomain_id = subdomain_number
    subdomain_number  += 1

work_htc_boundary_id = work_subdomain_id
cutter_htc_boundary_id = cutter_subdomain_id
chip_htc_boundary_id = chip_subdomain_id
holder_htc_boundary_id = holder_subdomain_id

chip_end_boundary_id = subdomain_number + 1
holder_top_boundary_id = subdomain_number + 2
mapping_boundary_id = subdomain_number + 3
mapped_shear_boundary_id = subdomain_number + 4
mapped_friction_boundary_id = subdomain_number + 5
work_bottom_boundary_id = subdomain_number + 6  #give it a big number, bigger than gmsh surface and volume id

################# thermla boundary #####################
air_kinematic_viscosity = 1e-5

emissivity = 0.75
T_ambient = 0
T_holder_top = T_ambient + 25
T_work_inlet = T_ambient + 25   # ref paper:
T_reference = T_work_inlet
T_reference_material = 25 # Celsius degree

using_high_HTC = False
if not using_high_HTC:
    htc_work = 100  # FIXME function of speed
    htc_chip = 100
    htc_cutter = 50
    htc_holder = 50
else:
    htc_work = 1000  # FIXME function of speed
    htc_chip = 1000
    htc_cutter = 500
    htc_holder = 50

##############################################
#material 1045 steel
taylor_quinn_coeff = 0.9  # 1045 steel

metal_cp_f = lambda T: (420 + 0.504*(T-T_reference_material))
metal_k_f = lambda T: (52.21 -0.0281*(T-T_reference_material))
metal_capacity_f = lambda T: (420*8000 + 0.504*8000*(T-T_reference_material))
#metal_diffusivity_f = lambda T: (52.21/8000/420 - 0.0281/8000/420*(T-T_reference_material))  not correct

#conductivity = 420+ 0.504*T, cp = 52.21 -0.0281*T,  density = 8000 constant
metal_density = 8000.0
k_ratio = 1.0  #tool to work thermal diffusivity
#T_reference = 800  # diff cutter material can change the max temperature by 4.5%
metal_cp_0 = metal_cp_f(T_reference)
metal_k_0 = metal_k_f(T_reference)
metal_capacity_0 = metal_cp_0*metal_density
metal_diffusivity_0 = metal_k_0/metal_capacity_0
T_reference = T_work_inlet

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
    T_m = 1460
    epsilon_dot_0 = 1 # reference strain
    T_w = T_reference_material
    T_factor = (1 - ((T - T_w)/(T_m - T_w))**m)
    return (A + B * epsilon**n) * (1 + C*log(epsilon_dot/epsilon_dot_0)) * T_factor

if using_nonlinear_thermal_properties:
    metal_cp = metal_cp_f
    metal_k = metal_k_f
    metal_capacity = metal_capacity_f
    #metal_diffusivity = metal_diffusivity_f
else:
    metal_cp = metal_cp_0
    metal_k = metal_k_0
    metal_capacity = metal_capacity_0
    metal_diffusivity = metal_diffusivity_0

# material 1045 steel
material_work = {'thermal_conductivity':  metal_k, 'density': metal_density, 'specific_heat_capacity': metal_cp, 'capacity': metal_capacity,
'dk_dT': -0.0281,  'dc_dT': 0.504*8000,
'thermal_expansion_coefficient': 11.2e-6}



material_friction_transition = material_work.copy()
material_friction_transition['thermal_conductivity'] = metal_k_0/2.0
# not checked
material_cutter = {'thermal_conductivity': 85 , 'density': 15200, 'specific_heat_capacity': 280, 'capacity': 15200*280,
    'thermal_expansion_coefficient': 4.5e-6}
#https://www.makeitfrom.com/material-properties/Tungsten-Carbide-WC

material_list = [material_work] * subdomain_number
##cutter has higher k, but contact resitance is not considered, treat cutter as workpiece material to consider that
material_list[cutter_subdomain_id] = material_cutter
material_list[friction_transition_subdomain_id] = material_friction_transition
#if has_pressing_zone:
#    material_list.append((pressing_subdomain_id, material_work))
material_k_list = [t['thermal_conductivity'] for t in material_list]
material_cp_list = [t['specific_heat_capacity'] for t in material_list]
material_capacity_list = [t['capacity'] for t in material_list]
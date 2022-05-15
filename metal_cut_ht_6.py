"""
solver developing notes
version 6.0
updated on Jan 4 2019

updated on May 14 2022, updated code for dolfin v2019.1


## version history

+ HTC equation from paper
+ velocity field for all in cpp velocity_expression
+ material as a list
+ heat partition: surface integration
+ third zone heat source. add in V3
+ mesh refinement -> done in salome_parameter.py
+ friction heating distribution controlled in parameter.py
+ tune velocity field, make it tangential with boundary  -> mapping boundary pair in V6
+ v6 remove all gmsh foam related code, use salome for all cases
+ frictional heat zone is marked out to use different thermal properties (linear model)
 since May 2019, paper revision
+ 3D holder shape, but mapping boundary bc pair as periodic_boundary failed
+ nonlinear looping is working now, but assuming same dc_dT for all zones
+ WC material properties updated

## observation or just information

* meshio-convert MED lib version too low to support
* scalarTransportFoam:  property, U, heat shource
* beta, ratio of heat from mechanical work
* radiation: still much small compared with HTC
* fenics does not support mesh mixed of tetra and prism
* Fenics v2017.1 fails to converge for nonlinear model, but v2017.2 is fine

## parallel computation is possible but not necessary

* DG0 for material expression and heat source, not working with MPI
* parallel not working, interpolate(velocity expression)  -> save serial result into HDF5File or XMDF

## bug fixed
- salome mesh exported volume id is changed in gmsh -> directly edit *.msh ascii file/ use salome mesh tool
- salome_mesh can not control sudomain, boundary export sequence (id)
- salome_mesh 3D, nan error, for this 2 outlets convection problem (solved, bad mesh quality)


- holder:  thermal boundary to lathe fixing (temperature unknown)  -> this is negletiable/important, since the heat loss ratio is very small compared
- workpiece bottom: temperature unknown -> this is negletiable/important,

###  no stabilization_method is used in this solver
- Internal panelty does not work with nonlinear solver, https://fenicsproject.org/qa/3974/expression-in-interior-facets/
- capacity is fixed as constant, otherwise, ulferror for internal panelty items  (run without internal panelty)

## exp validation cases

* geometry of exp study: Ivester 2000
* cutter material, Kennametal grade K68     WC/Co
* cutter holder, Kennametal CTAPR 123B
* work geometry, tube, Tubes-152.4 mm overhang,  101.6 mm  diameter,  1.6 mm  wall thickness

Lalwani 2009 table6 Aisi 1045 steel, done
Karpart 2006 table2 and table3, but nonlinear material property is not found

-## issues for nonlinear material

- judge from parameter study, the biggest error (deviation from experimental interface temperature)
 is caused by the assumption of rectangular friction zone

- heat generation != volume integral, estimated error 1~1.5%
  reason: metal_k/dT does not affect, but metal_capacity/dT make a difference, refine mesh can reduce this error a bit

- cutter material thermal properties uses metal work 's, error 0.5% estimated from linear model
  : thermal contact can decrease heat transfer into cutter,
  while regarding cutter material as workpiece can also decrease that,
  the Tungsten Carbide's thermal property as functions are unknown,
  so just do nothing on this.

"""

from __future__ import print_function, division
import math, sys, os.path
#from subprocess import run_command
import numpy as np

#=============================
try:
    import dolfin
    #from mshr import *
    from dolfin import *
    print("found fenics version: ", dolfin.__version__)
except ImportError as e:
    print(e)
    print('Fenics is not installed, exit')
    exit(-1)

try:
    ver = dolfin.dolfin_version().split('.')
except:
    ver = dolfin.__version__.split('.')

if int(ver[0]) <= 2017 and int(ver[1])<2:
    UserExpression = Expression
    using_VTK = True
else:
    using_VTK = False

# dolfin v2018  Rename ERROR, CRITICAL etc. to LogLevel.ERROR, LogLevel.CRITICAL.
# dolfin v2018  Rename mpi_comm_world() to MPI.comm_world.
if int(ver[0]) >= 2018:
    set_log_level(LogLevel.WARNING) # depend on dolfin version, this line may failed
    if dolfin.MPI.comm_world.size > 1:
        using_MPI = True
    else:
        using_MPI = False
else:
    set_log_level(WARNING) # depend on dolfin version, this line may failed
    if dolfin.MPI.size(dolfin.mpi_comm_world()) > 1:
        using_MPI = True
    else:
        using_MPI = False

# in case if you have not install Fenics and FenicsSolver to python path
sys.path.append('/media/sf_OneDrive/gitrepo/FenicsSolver')
from FenicsSolver import ScalarTransportSolver

#set_log_level(3)
#set_log_level(TRACE)
#solver.parameters['monitor_convergence'] = verbose
#solver.parameters['report'] = verbose
parameters['std_out_all_processes'] = False
#parameters['ghost_mode'] = 'shared_facet'
#===============================

def defined(x):
    return x in locals() or x in globals()


from dolfin_utilities import mark_boundary_for_subdomain, convert_to_hdf5_mesh_file
'''
#shearing zone velocity step has big impact on chipend flux out.
def mark_boundary_for_subdomain(mesh, subdomains):
    boundaries = MeshFunction("size_t", mesh, mesh.geometry().dim()-1)
    boundaries.set_all(0)
    D = mesh.topology().dim()
    mesh.init(D-1,D) # Build connectivity between facets and cells
    #cellmap = boundaries.entity_map(2)
    for f in facets(mesh):
        for c in cells(f):
            #print(f.index(), c.index())
            boundaries[f.index()] = subdomains[c.index()]

    return boundaries

def convert_to_hdf5_mesh_file(filename):
    assert  filename[-4:] == ".xml"
    filename_base = filename[:-4]
    mesh = Mesh(filename)
    dim = mesh.geometry().dim()
    hdf = HDF5File(mesh.mpi_comm(), filename_base + ".h5", "w")
    hdf.write(mesh, "/mesh")
    subdomain_file = filename_base + "_physical_region.xml"
    if os.path.exists(subdomain_file):
        subdomains = MeshFunction("size_t", mesh, subdomain_file)
        #plot(subdomains)
        hdf.write(subdomains, "/subdomains")
    bmeshfile =filename_base + "_facet_region.xml"
    if os.path.exists(bmeshfile):
        boundaries = MeshFunction("size_t", mesh, bmeshfile)
    else:
        boundaries = mark_boundary_for_subdomain(mesh, subdomains)
    hdf.write(boundaries, "/boundaries")
    print("XML mesh files have been converted into HDF5 file")
'''

################geometry parameter for gmsh################
using_gmsh = False  # gmsh only works for  metal_cut_ht_5
using_salome = not using_gmsh  # does not work since it is rect work shape


if using_salome:
    from parameter_salome import *
else:
    raise RuntimeError('gmsh meshing has dropped/deleted in this version, use salome to mesh')
    from parameter_gmsh import *  # current not completed
    sys.exit(-1)

from velocity_expression import velocity_code
from mesh_utilities import generate_salome_mesh, convert_salome_mesh_to_dolfin, convert_salome_mesh_to_foam
from analytical import get_heat_source, get_analytical_T, get_heat_partition
############### derived parameters,  math function has diff name ! #################

if using_salome:  # now mesh output folder are unified
    meshfolder = 'salome_mesh'
    if using_2D:
        meshfolder = 'salome_mesh'
    if using_mapping_bc:
        meshfolder = 'salome_mesh'
else:
    meshfolder = 'gmsh_mesh_v4'  # gmsh4
mesh_file_root = meshfolder + "/metal_cut"

if is_preprocessing:
    if using_MPI:
        raise NotImplementedError()
    #also save some field into hdf5, which can not been run in parallel
    if using_salome:
        #generate_salome_mesh('salome_mesh_metal_cut_1.py')
        generate_salome_mesh('salome_mesh_metal_cut.py')
        convert_salome_mesh_to_dolfin(mesh_file_root + '.xml')
        if exporting_foam:
            if using_mapping_bc:
                foam_data_folder = '/media/sf_OneDrive/gitrepo/icoThermoFoam/metalCutMapping/'
            else:
                foam_data_folder = '/media/sf_OneDrive/gitrepo/icoThermoFoam/metalCut/'
            #export and conver the mesh from gmshToFoam
            convert_salome_mesh_to_foam(foam_data_folder)
            sys.exit()  # no need to carry on, since fenics will not support hex mesh or mixed cell in xml file

    convert_to_hdf5_mesh_file(mesh_file_root + '.xml')
#MPI: load from hdf5 file
mesh_file = meshfolder + "/metal_cut.h5"
print("mesh file name", mesh_file)

########### empirical formula ###########



#######################################
# Sub domain for Periodic boundary condition
# facet near origin is not selected, why?
tol = 2e-7
radian_tol = 1e-3  # this is required for minus cutting angle
_a_phi = pi - shear_angle * pi/180
class TurningInterface(SubDomain):
    # chip size shear plane
    def inside(self, x, on_boundary):
        if near(x[0], 0) and near(x[1], 0) and on_boundary:
            #print(x)  # if both x, y are zero,  atan2() returns zero
            return True
        else:
            angle = math.atan2(x[1], x[0])
            return bool( math.fabs(angle - _a_phi) < radian_tol   \
                                and x[1] >= -tol and x[0] <= tol and on_boundary)

_a_cutter_v = pi/2.0 - cutter_angle_v * pi/180
class FrictionalInterface(SubDomain):
    def inside(self, x, on_boundary):
        if near(x[0]-mapping_bc_gap, 0) and near(x[1], 0) and on_boundary:
            #print(x)  # if both x, y are zero,  atan2() returns zero
            return True
        else:
            angle = math.atan2(x[1], x[0]-mapping_bc_gap)
            return bool( math.fabs(angle - _a_cutter_v) < radian_tol and x[1] <= chip_ccc_y+tol \
                                and x[1] >= -tol and on_boundary)

if using_chip_cutter_mapping_bc:  # frictional interface
    class PeriodicBoundary(SubDomain):
        # Left boundary is "target domain" G
        def inside(self, x, on_boundary):
            #
            if (near(x[0]+mapping_bc_gap, 0) or near(x[0], 0)) and near(x[1], 0) and on_boundary:
                #print(x)
                return True
            if using_3D and x[2]<0:
                return False
            angle = math.atan2(x[1], x[0]+mapping_bc_gap)
            angle2 = math.atan2(x[1], x[0])
            return bool(((math.fabs(angle - _a_phi) < radian_tol and x[1] >= -tol and x[0]+mapping_bc_gap <= tol) \
                            or (math.fabs(angle2 - _a_cutter_v) < radian_tol and x[1] <= chip_ccc_y+tol and x[1] >= -tol)) \
                        and on_boundary)

        # Map x in right boundary (H) to y in left boundary (G)
        def map(self, x, y):
            y[0] = x[0] - mapping_bc_gap
            y[1] = x[1]
            if using_3D:
                y[2] = x[2]
else:  # shear interface
    class PeriodicBoundary(SubDomain):
        # Left boundary is "target domain" G
        def inside(self, x, on_boundary):
            #
            if near(x[0]+mapping_bc_gap, 0) and near(x[1], 0) and on_boundary:
                #print(x)
                return True
            angle = math.atan2(x[1], x[0]+mapping_bc_gap)
            return bool( math.fabs(angle - _a_phi) < radian_tol and x[1] >= -tol and x[0]+mapping_bc_gap <= tol \
                        and on_boundary)

        # Map x in right boundary (H) to y in left boundary (G)
        def map(self, x, y):
            y[0] = x[0] - mapping_bc_gap
            y[1] = x[1]
            if using_3D:
                y[2] = x[2]

if using_mapping_bc:
    pbc = PeriodicBoundary()
    print('=== using mapping boundary ===')
else:
    pbc = None

###########################################


def solve_ht():
    #  error for inf temperature field if all BC are HTC
    bcs = {
            "work": { 'boundary_id': work_htc_boundary_id, 'values': {
                                'temperature': {'variable': "temperature", 'type': 'HTC', 'value': htc_work, 'ambient': T_ambient } } },
            "cutter":  {'boundary_id': cutter_htc_boundary_id, 'values': {
                                #'temperature': {'variable': "temperature", 'type': 'Dirichlet', 'value': T_ambient + 200} } },
                                'temperature': {'variable': "temperature", 'type': 'HTC', 'value': htc_cutter, 'ambient': T_ambient } } },
            "holder":  {'boundary_id': holder_htc_boundary_id, 'values': {
                                'temperature': {'variable': "temperature", 'type': 'HTC', 'value': htc_holder, 'ambient': T_ambient} } },
                                #'temperature': {'variable': "temperature", 'type': 'Dirichlet', 'value': T_ambient + 20} } },
            "chip":  { 'boundary_id': chip_htc_boundary_id, 'values': {
                                'temperature': {'variable': "temperature", 'type': 'HTC', 'value': htc_chip, 'ambient': T_ambient } } },
            #other are zero heat flux
    }

    settings = {'solver_name': 'ScalarEquationSolver',
                    'mesh': mesh_file, 'function_space': None,  'fe_degree': element_degree, 'periodic_boundary': pbc,
                    'boundary_conditions': bcs,
                    'body_source': None,
                    'initial_values': {'temperature': T_ambient},
                    'material': material_work,
                    'solver_settings': {
                        'transient_settings': {'transient': False, 'starting_time': 0, 'time_step': 0.01, 'ending_time': 0.03},
                        'reference_values': {'temperature': T_ambient},
                        'solver_parameters': {"relative_tolerance": 1e-8,  # mapping to solver.parameters of Fenics
                                                        'absolute_tolerance': 1E-9,
                                                        "maximum_iterations": 500,
                                                        "relaxation_parameter": 0.3,  # case 6 , case 8 failed to converge
                                                        "monitor_convergence": True,  # print to console
                                                        },
                        },
                    # solver specific settings
                    'scalar_name': 'temperature',
                    "convective_velocity": None,
                    #'radiation_settings': {'ambient_temperature': T_ambient-20, 'emissivity': 0.9}
                    }
    if using_stab:  # stab is not necessary if mapping boundary is used
        settings['advection_settings'] = {'stabilization_method': 'IP', 'alpha': IP_alpha, 'Pe':  Pe/1000.0}
    if considering_radiation:
        settings['radiation_settings']= {'ambient_temperature': T_ambient, 'emissivity': emissivity}

    import time  # dolfin has time function
    ts_before_preprocessing = time.time()
    solver = ScalarTransportSolver.ScalarTransportSolver(settings)

    Q = solver.function_space
    vector_degree = element_degree+1
    V = VectorFunctionSpace(solver.mesh, 'CG', vector_degree)
    # Expression has been deprecated, user UserExpression instead
    v_e = UserExpression(cppcode=velocity_code, degree=vector_degree)  # degree, must match function space
    v_e.subdomain_id = solver.subdomains
    velocity = Function(V)
    velocity.interpolate(v_e)  # does not work for MPI
    if has_convective_velocity:
        #solver.convective_velocity = Constant((1, 1, 0))
        solver.convective_velocity = velocity

    #DG for conductivity and heat source, save to hdf5?
    DG0 = FunctionSpace(solver.mesh, 'DG', 0)  # can not mixed higher order CG and DG?
    unity = Function(DG0)
    unity.vector()[:] = 1 #interpolate(Expression('1', element_degree), DG0)


    # Numerical simulations of advection-dominated scalar mixing with applications to spinal CSF flow and drug transport
    if not using_nonlinear_thermal_properties:
        k_f = get_material_function(DG0, solver.subdomains, material_k_list)
        solver.material['thermal_conductivity'] = k_f
        capacity_f = get_material_function(DG0, solver.subdomains, material_capacity_list)
        solver.material['capacity'] = capacity_f

    #################################################
    shear_heat, friction_heat = get_heat_source()
    shear_heat_density = shear_heat / shear_heat_volume  # W/m3
    friction_heat_density = friction_heat / nominal_friction_heat_volume
    print("heating volume and intensity; ", shear_heat_volume, friction_heat_volume, shear_heat_density, friction_heat_density)
    if True:  # test passed, non-uniform heat source can be implemented later
        class HeatSourceExpression(Expression):
            def eval_cell(self, values, x, ufc_cell):
                cindex = ufc_cell.index
                did = solver.subdomains[cindex]
                if did == shear_subdomain_id:
                    values[0] = shear_heat_density
                elif did == friction_subdomain_id:
                    if nonuniform_friction_heat:
                        distance = math.sqrt(x[0]*x[0] + x[1]*x[1])
                        _start = (friction_heat_start + uniform_friction_heat_length)
                        if distance <= _start:
                            values[0] = friction_heat_density * uniform_heat_ratio
                        else:
                            _ratio_l = (distance -_start)/(actual_friction_heat_distance - _start)
                            _ratio_h = uniform_heat_ratio * (1 - _ratio_l)
                            values[0] = friction_heat_density * _ratio_h
                    if nonuniform_friction_thickness:
                        distance = math.sqrt(x[0]*x[0] + x[1]*x[1])
                        _start = friction_heat_start
                        _ratio_l = (distance -_start)/actual_friction_heat_length
                        _ratio_h = 1.0 / ( 1.0 - _ratio_l * ( 1.0 - friction_end_thickness/friction_heat_thickness))
                        values[0] = friction_heat_density * _ratio_h
                    else:
                        values[0] = friction_heat_density
                else:
                    values[0] = 0
        # FIXME: cause error of UFLException: Discontinuous type Coefficient must be restricted.
        v_s = HeatSourceExpression(degree = element_degree)  # degree, must match function space
        heat_source = Function(DG0)
        heat_source.interpolate(v_s)  # does not work for MPI
        solver.body_source = heat_source
    else:
        heat_source_list = [0] * subdomain_number
        heat_source_list[shear_subdomain_id] = shear_heat_density
        heat_source_list[friction_subdomain_id] = friction_heat_density
        heat_source = get_material_function(DG0, solver.subdomains, heat_source_list)  # does not work for MPI
        solver.body_source =  heat_source
    '''
    solver.boundary_facets.set_all(0)
    for facet in facets(solver.mesh):
        for cell in cells(facet):
            #print(f.index(), c.index())
            if facet.exterior() and solver.subdomains[cell.index()] < 5:
                solver.boundary_facets[facet.index()] = solver.subdomains[cell.index()]
    '''

    if exporting_foam and False:  # this is not supported in public repo
        #preset boundary value
        #DG vector space
        #v_e = Expression(cppcode=velocity_code, degree=vector_degree)  # degree, must match function space?
        #v_e.subdomain_id = solver.subdomains
        #print(velocity_code)

        DGV = VectorFunctionSpace(solver.mesh, 'DG', 0)
        velocity = Function(DGV)
        velocity.interpolate(v_e)  # does not work for MPI
        print(' chip volocity = ({}, {}, {})'.format(chip_sliding_vel_x, chip_sliding_vel_y, 0))
        print(' work volocity = ({}, {}, {})'.format(cutting_speed, 0, 0))

        sys.path.append('/media/sf_OneDrive/gitrepo/VTKToFoam')
        from VTKToFoam import set_internal_field, write_field
        foam_data_folder = '0/'
        #foam_data_folder = '0.origin/'

        #print('size of DG0 vector numpy array', velocity.vector().array().shape)  # correct 1D array
        #print('size of DG0 scalar numpy array', heat_source.vector().array().shape)
        # need to regerate mesh and transform mesh even if boundary, subdomain id changed.
        foam_velocity_data_file = foam_data_folder + foam_data_folder + 'U'
        set_internal_field(foam_velocity_data_file, velocity.vector().array(), 'vector',  foam_data_folder + '0.origin/U')

        # icoThermoFoam TEqn  accepts  power/capacity as source
        #set_internal_field(foam_data_folder + 'S', heat_source.vector().array()/metal_capacity_0, 'scalar', foam_base_folder + 'S')  # if not existing, write_field()
        print("set heat source density in system/")
        sys.exit()

    #
    holder_top_boundary = AutoSubDomain(lambda x, on_boundary: \
                                        near(x[1], holder_top_y) and on_boundary)
    holder_top_boundary.mark(solver.boundary_facets, holder_top_boundary_id)
    solver.boundary_conditions['holder_top'] = { 'boundary_id': holder_top_boundary_id, 'values': {
                    'temperature': {'variable': "temperature", 'type': 'Dirichlet', 'value': T_holder_top } } }
    #
    if is_straight_chip:
        chip_end_boundary = AutoSubDomain(lambda x, on_boundary: \
                    near(x[1], chip_top_y) and on_boundary)
    else:
        chip_end_boundary = AutoSubDomain(lambda x, on_boundary: \
                    near(x[0], chip_center_x) and (x[1]>0) and on_boundary)

    chip_end_boundary.mark(solver.boundary_facets, chip_end_boundary_id)
    if not has_convective_velocity:
        solver.boundary_conditions['chipend_htc'] = { 'boundary_id': chip_end_boundary_id, 'values': {
                    'temperature': {'variable': "temperature", 'type': 'HTC', 'value':  htc_chip, 'ambient':T_ambient } } }

    if has_convective_velocity:
        if using_mapping_bc: # test if bc is found
            # those are used in post-processing, PeriodicBoundary should be done in function space creation
            mapping_periodic_boundary = PeriodicBoundary()
            mapping_periodic_boundary.mark(solver.boundary_facets, mapping_boundary_id, check_midpoint = True)
            #
            shear_interface = TurningInterface()  # chip side
            shear_interface.mark(solver.boundary_facets, mapped_shear_boundary_id, check_midpoint = True)
            #
            if using_chip_cutter_mapping_bc:
                friction_interface = FrictionalInterface()  # cutter side
                friction_interface.mark(solver.boundary_facets, mapped_friction_boundary_id, check_midpoint = True)

        if workpiece_type_id == 0 or workpiece_type_id == 2:  # rect work shape
            work_bottom_boundary = AutoSubDomain(lambda x, on_boundary: near(x[1], work_bottom_y) and on_boundary)
            work_bottom_boundary.mark(solver.boundary_facets, work_bottom_boundary_id)
            #
            work_left_boundary_id = work_bottom_boundary_id+1
            work_left_boundary = AutoSubDomain(lambda x, on_boundary: near(x[0], work_left_x) and x[1]<feed_thickness+DOLFIN_EPS and on_boundary)
            work_left_boundary.mark(solver.boundary_facets, work_left_boundary_id)
            #
            work_right_boundary = AutoSubDomain(lambda x, on_boundary: near(x[0], work_right_x) and x[1]<feed_thickness+DOLFIN_EPS and on_boundary)
            work_right_boundary.mark(solver.boundary_facets, work_bottom_boundary_id+2)
            #right are natural boundary, zero gradient, or it does not mater
            solver.boundary_conditions['work_inlet'] = { 'boundary_id': work_left_boundary_id, 'values': {
                                        'temperature': {'variable': "temperature", 'type': 'Dirichlet', 'value': T_work_inlet } } }
            #solver.boundary_conditions['work_bottom'] = { 'boundary_id': work_bottom_boundary_id, 'values': {
            #                            'temperature': {'variable': "temperature", 'type': 'Dirichlet', 'value': T_ambient } } }
        elif workpiece_type_id == 1:
            disc_center_y = - disc_R
            work_clamp_boundary_id = 8
            work_clamp_boundary = AutoSubDomain(lambda x, on_boundary: \
                between(x[0], (-disc_R*0.05, disc_R*0.05)) \
                and between(x[1], (disc_center_y-disc_R*0.1, disc_center_y+disc_R*0.1)) and on_boundary)
            work_clamp_boundary.mark(solver.boundary_facets, work_clamp_boundary_id)
            solver.boundary_conditions['work_clamp'] = { 'boundary_id': work_clamp_boundary_id, 'values': {
                                        'temperature': {'variable': "temperature", 'type': 'Dirichlet', 'value': T_ambient } } }
        else:
            raise NotImplementedError('this work piece shape id is not supported')


    def update_material_property(DG0, subdomains, material_values, T_DG0, cutter_f, metal_f):
        # it takes very long time to complete if set individually
        did = np.asarray(subdomains.array(), dtype=np.int32)
        u = Function(DG0)
        u.vector()[:]  = metal_f(T_DG0.vector().get_local())  # return numpy.array
        #print(type(u_new), u_new.size, u.vector().array().size)

        u.vector()[did == cutter_subdomain_id]  = cutter_f(T_DG0.vector().get_local())[did == cutter_subdomain_id]
        #u.vector()[did == cutter_subdomain_id] = cutter_value
        """
        for i, v in enumerate(did):
            if (v == cutter_subdomain_id):
                u.vector()[i] = cutter_value
            else:
                u.vector()[i] = metal_f(T_DG0.vector().array()[i])
        """
        return u

    ts_after_preprocessing = time.time()
    # my own nonliner loop, updating material property in each iteration
    if using_nonlinear_loop:  # 3D looping is possible
        assert not using_MPI  # setting material property may not working in parallel
        loop_i = 0
        loop_N = 10
        T_old = Function(solver.function_space)  #default to zero, Yes
        while loop_i < loop_N:
            T = solver.solve()
            T_mean_diff = np.mean(T_old.vector()[:] - T.vector()[:])
            print("=== nonlinear loop ", loop_i, " ======= ", T_mean_diff)
            if  math.fabs(T_mean_diff) < 0.1:
                break
            T_DG0 = project(T, DG0)
            capacity_field = update_material_property(DG0, solver.subdomains, material_capacity_list, T_DG0, cutter_capacity_f, metal_capacity_f)
            solver.material['capacity'] = capacity_field
            k_field = update_material_property(DG0, solver.subdomains, material_k_list, T_DG0, cutter_k_f, metal_k_f)
            solver.material['thermal_conductivity'] = k_field
            #solver.solve() has considered dc/dT        ignore dk/dT
            #to trigger this condition: if 'dc_dT' in self.material:

            T_old.vector()[:] = T.vector()[:]
            loop_i += 1
            #ofile = File(result_folder + "k.pvd")  # not for parallel
            #ofile << solver.material['thermal_conductivity']
    else:
        T = solver.solve()

    ts_before_postprocessing = time.time()
    if using_MPI:
        mf = HDF5File(mpi_comm_world(),  result_folder + "metal_cutting_result" +'.h5', 'w')
        #mf = XDMFFile(mpi_comm_world(), result_name+'.xdmf')
        mf.write(T, "T")
        mf.write(V.mesh(), "mesh")
        #post processing is not possible in parallel?
        sys.exit()
    else:
        #interpolate the higher order solution onto a finer linear space
        T.rename("Temperature (C)", "temperature contour")
        ofile = File(result_folder + "T.pvd")  # not for parallel
        ofile << T

    ############################### post-processing#####################

    is_simulation_valid = True
    error_message = ''
    dx= Measure("dx", subdomain_data=solver.subdomains)
    print("================= summary =================")
    print("shear_heat, friction_heat by definition", shear_heat, friction_heat)
    volume_shear_zone = assemble(unity*dx(shear_subdomain_id))
    volume_friction_zone = assemble(unity*dx(friction_subdomain_id))
    print("volume_shear_zone, volume_friction_zone by integral:", volume_shear_zone, volume_friction_zone)

    heat_source = solver.body_source
    heat1 = assemble(heat_source*dx(shear_subdomain_id))
    heat2 = assemble(heat_source*dx(friction_subdomain_id))
    heat_chip = assemble(heat_source*dx(chip_subdomain_id))
    print("shear_heat, friction_heat, heating by chip (should be zero) by integration", heat1, heat2, heat_chip)
    total_heat = assemble(heat_source*dx)
    if using_3D:
        heat_ratio = total_heat/(shear_heat + friction_heat)
    else:
        #total_heat = total_heat*cutter_thickness
        heat_ratio = total_heat*cutter_thickness/(shear_heat + friction_heat)
    print("total_heat, ratio of heat integral to the defined", total_heat,  heat_ratio)
    if math.fabs(heat_ratio - 1) > validation_tol:
        error_message += "heat generation != heat source integration"
        is_simulation_valid = False

    print("========== heat loss from HTC =================")
    ds= Measure("ds", subdomain_data=solver.boundary_facets)
    # radiation is not added into cooling
    cooling_frictinal_zone = assemble(htc_work*(T-T_ambient)*ds(friction_subdomain_id))
    cooling_shear_zone = assemble(htc_work*(T-T_ambient)*ds(shear_subdomain_id))
    print("cooling from friction and shear zone surface",  cooling_frictinal_zone, cooling_shear_zone)
    cooling_work = assemble(htc_work*(T-T_ambient)*ds(work_htc_boundary_id))
    cooling_chip = assemble(htc_chip*(T-T_ambient)*ds(chip_htc_boundary_id))
    cooling_tool = assemble(htc_cutter*(T-T_ambient)*ds(cutter_htc_boundary_id) +
                                          htc_holder*(T-T_ambient)*ds(holder_htc_boundary_id))
    cooling_total_HTC = cooling_work + cooling_chip + cooling_tool + cooling_frictinal_zone + cooling_shear_zone
    cooling_HTC_all = assemble(htc_work*(T-T_ambient)*ds)  # not all surface are HTC
    print("convective cooling_work, cooling_chip, cooling_tool and holder, cooling_HTC_all",\
                cooling_work, cooling_chip, cooling_tool, cooling_HTC_all)
    print("ratio of HTC loss to generation", cooling_total_HTC/total_heat)
    total_heat_loss = cooling_total_HTC

    if considering_radiation:
        Stefan_constant = 5.670367e-8  # W/m-2/K-4
        m_ = emissivity * Stefan_constant
        radiation_outflux = - m_*(T_ambient**4 - pow(T, 4))
        R_work = assemble(radiation_outflux*ds(work_htc_boundary_id))
        R_chip = assemble(radiation_outflux*ds(chip_htc_boundary_id))
        R_tool = assemble(radiation_outflux*ds(cutter_htc_boundary_id) +
                                              radiation_outflux*ds(holder_htc_boundary_id))
        R_total = R_work + R_chip + R_tool
        R_all_surfaces = assemble(radiation_outflux*ds)
        total_heat_loss += R_total
        print("radiative: cooling_work, cooling_chip, cooling_tool and holder", R_work, R_chip, R_tool, R_all_surfaces)
        if math.fabs(R_total / R_all_surfaces- 1) > 0.3:
            error_message += "radiation sum != radiation integration"
            is_simulation_valid = False
        print("ratio of heat generation radiative cooling",  R_total/total_heat, R_all_surfaces/total_heat)
    else:
        R_total = 0


    T_DG0 = interpolate(T, DG0)
    did = np.asarray(solver.subdomains.array(), dtype=np.int32)
    T_DG0.vector()[did != shear_subdomain_id] = 0.0
    Tmax_shear_zone = np.max(T_DG0.vector().array())
    # defined in salome_parameter.py, failed in 2D for v2017.2
    try:
        p_probe = Point(*point_shear_surface_center_coordinates)
        T_probe = T(p_probe)
    except:
        T_probe = 0

    Tmean_shear_zone = assemble(T*dx(shear_subdomain_id))/volume_shear_zone
    Tmean_friction_zone = assemble(T*dx(friction_subdomain_id))/volume_friction_zone
    Tmean_chip_zone = assemble(T*dx(chip_subdomain_id))/chip_volume
    Tmax_friction_zone = np.max(T.vector().array())

    T_analytical_shear_zone,  T_analytical_friction_zone = get_analytical_T()
    print('================ temperature prediction by FEA ==================')
    print('Tmean_shear_zone = ', Tmean_shear_zone)
    print('Tmax_shear_zone = ', Tmax_shear_zone, "T_probe = ", T_probe)
    print('Tmean_chip_zone = ', Tmean_chip_zone)
    print('Tmean_friction_zone = ', Tmean_friction_zone)
    print('Tmax_friction_zone = ', Tmax_friction_zone)


    if has_convective_velocity:
        #chip material flow taking away from system, how to get T_chip_end, (chip_cross_area*chip_velocity) is known
        #capacity = metal_density*metal_cp
        chip_speed = cutting_speed * (feed_thickness/chip_thickness)
        #QC = capacity * chip_speed* (cutter_thickness * chip_thickness)
        #print("outflow heat = %f * delta_T"%QC)
        #T_chip_end_center = T[p_chip_end_xyz]
        surface_normal = FacetNormal(solver.mesh)

        #area_shear_zone = assemble(unity*ds(shear_subdomain_id))
        #area_friction_zone = assemble(unity*ds(friction_subdomain_id))
        #print('area_shear_zone = ', area_shear_zone, "should ==", 2*l_AB * shear_heat_thickness + shear_heat_thickness*cutter_thickness)
        #print('area_friction_zone = ', area_friction_zone, "should ==", 2*actual_friction_heat_length * friction_heat_thickness)

        if using_nonlinear_loop:
            _capacity = solver.material['capacity']
            _conductivity = solver.material['thermal_conductivity']
        elif using_nonlinear_thermal_properties:  #CG, cutter has same material with metal,  no difference
            #T_DG0 = interpolate(T, DG0)
            #_capacity = update_material_property(DG0, solver.subdomains, material_capacity_list, T_DG0, material_cutter['capacity'], metal_capacity_f)
            #_conductivity = update_material_property(DG0, solver.subdomains, material_k_list, T_DG0, material_cutter['capacity'], metal_k_f)
            _capacity = metal_capacity_f(T)
            _conductivity = metal_k_f(T)
        else:
            _capacity = Constant(metal_density*metal_cp_0)
            _conductivity = metal_k_0
        _density = metal_density

        system_heat_capacity = assemble((T-Constant(T_reference)) * _capacity*_density*dx)
        print("characteristic time = ", system_heat_capacity/total_heat)  # too big not quite usful

        if using_mapping_bc:
            #validation
            print('============ mapping bc validation ================')
            if using_3D:
                _area_theory = l_AB*cutter_thickness
            else:
                _area_theory = l_AB
            area_mapping = assemble(unity*ds(mapping_boundary_id))
            print('area mapping bounadry and theroetical (should be equal): ', area_mapping, _area_theory)
            area_mapped_shear = assemble(unity*ds(mapped_shear_boundary_id))

            if using_workpiece_extra_extusion:
                if math.fabs(area_mapping / (_area_theory) - 1) > validation_tol:
                    error_message += "\n mapping interface area is not equal\n"
                    is_simulation_valid = False

            M_mapping = assemble(dot(surface_normal, velocity)*ds(mapping_boundary_id))  # not equal, why?
            print('mass flux from mapping bounadries vs in theory (should be equal): ', M_mapping, cutting_speed*cutter_thickness*feed_thickness)
            Q_mapping_bottom = assemble(dot(surface_normal, velocity) * (T-T_reference)*_capacity*ds(mapping_boundary_id))
            Q_mapping_top = assemble(dot(surface_normal, velocity) * (T-T_reference)*_capacity*ds(mapped_shear_boundary_id))
            print('heat flux from mapping bounadries (should be equal): ', Q_mapping_top, Q_mapping_bottom)
            if using_workpiece_extra_extusion:
                if math.fabs(math.fabs(Q_mapping_top / Q_mapping_bottom) - 1) > validation_tol:
                    error_message += "\n heat flux on mapping interfaces area is not equal\n"
                    is_simulation_valid = False

            Q_mapping_work = assemble(_conductivity*dot(grad(T), surface_normal) * ds(mapping_boundary_id)) * -1
            Q_mapping_chip = assemble(_conductivity*dot(grad(T), surface_normal) * ds(mapped_shear_boundary_id)) * -1
            print("conductive heat transfer, Q_mapping_work, Q_mapping_chip", Q_mapping_work, Q_mapping_chip)
            if using_chip_cutter_mapping_bc:
                area_mapped_friction = assemble(unity*ds(mapped_friction_boundary_id))
                print('area for a chip-cutter fricitonal interface  and theroetical value (should be equal): ', \
                        area_mapped_friction, actual_friction_heat_length*cutter_thickness)  # 2D has diff area
                print('total area for mapping (should be equal to), and sum of cutter and work', \
                        area_mapping, area_mapped_friction+area_mapped_shear)

        print('============ heat with mass flow ================')
        # heat loss due to mass flow, should be zero
        _tmp = dot(surface_normal, velocity) * (T-Constant(T_reference)) * _capacity
        Q_friction = assemble(_tmp*ds(friction_subdomain_id))
        print("heat flow with mass out of system from friction zone is:", Q_friction)
        Q_shearing = assemble(_tmp*ds(shear_subdomain_id))
        print("heat flow with mass out of system from sheairing zone is:", Q_shearing)

        Q_chip = assemble(_tmp*ds(chip_htc_boundary_id))
        print("heat flow with mass out of system from chip htc boundary is:", Q_chip)
        Q_work = assemble(_tmp*ds(work_htc_boundary_id))
        print("heat flow with mass out of system from work htc boundary (should be zero)is:", Q_work)
        Q_cutter = assemble(_tmp*ds(cutter_htc_boundary_id))
        print("heat flow with mass out of system from cutter htc is:", Q_cutter)
        Q_holder = assemble(_tmp*ds(holder_htc_boundary_id))
        print("heat flow with mass out of system from holder htc is:", Q_holder)
        Q_chipend = assemble(_tmp*ds(chip_end_boundary_id))
        print("heat flow with mass from  chipend", Q_chipend)
        if has_friction_transition_zone:  # no using mapping bc for friction contact
            Q_friction_transition_subdomain = assemble(_tmp*ds(friction_transition_subdomain_id))
            print("heat flow with mass from friction_transition_subdomain", Q_friction_transition_subdomain)
        Q_ds0 = assemble(_tmp*ds(0))
        print("heat flow with mass from  ds(0), should be zero", Q_ds0)
        total_heat_loss += Q_chipend  # Q_material_out_sum,  excluding the workpiece outlet

        """
        print('============ mass flow ================')
        M_shearing = assemble(dot(surface_normal, velocity)*ds(shear_subdomain_id))
        print("mass flow out of system from shearing zone boundary  is:", M_shearing)
        M_friction = assemble(dot(surface_normal, velocity)*ds(friction_subdomain_id))
        print("mass flow out of system from friction boundary is:", M_friction)
        M_work = assemble(dot(surface_normal, velocity)*ds(work_htc_boundary_id))
        print("mass flow out of system from work  is:", M_work)
        M_chip = assemble(dot(surface_normal, velocity)*ds(chip_htc_boundary_id))
        print("mass flow out of system from chip HTC boundary is:", M_chip)
        M_chipend = assemble(dot(surface_normal, velocity)*ds(chip_end_boundary_id))
        print("mass flow from chip end in theory VS integral at chipend", chip_speed * (cutter_thickness * chip_thickness), M_chipend)
        """

        Q_material_out = assemble(_tmp*ds)
        if workpiece_type_id == 0:  #rectangle
            #Q_bottom = cutting_speed * capacity *ds(work_bottom_boundary_id)
            Q_inflow = assemble( (T-Constant(T_reference))* cutting_speed * _capacity * ds(work_bottom_boundary_id + 1))
            #* cutter_thickness * (chip_start_y - work_bottom_y)
            Q_outflow = assemble( (T-Constant(T_reference))* cutting_speed * _capacity * ds(work_bottom_boundary_id + 2))
            print('inflow from left boundary and outflow from right bounadry is', Q_inflow, Q_outflow)

            Q_bottom = assemble(_conductivity*dot(grad(T), surface_normal)*ds(work_bottom_boundary_id)) * -1 #flow out

            if using_workpiece_extra_extusion:
                Q_flowing_out_work = (Q_material_out  - Q_chipend) + Q_bottom
            else:
                Q_flowing_out_work = (Q_outflow - Q_inflow) + Q_bottom
            print('heat flow out from work oulet and bottom is:', Q_flowing_out_work)
            total_heat_loss += Q_flowing_out_work
        else:  # disc clamping
            Q_clamp = assemble(_conductivity*dot(grad(T), surface_normal)*ds(work_clamp_boundary_id)) * -1 #flow out
            print('heat flux from clamp bounadry is', Q_clamp)
            Q_flowing_out_work = Q_clamp + cooling_work
            total_heat_loss += Q_flowing_out_work

        Q_material_out_sum = Q_chip + Q_work + Q_cutter + Q_holder + Q_chipend + Q_outflow + Q_inflow
        print("heat flow with mass out of system from all boundary byintegral (ds) and sum:",  Q_material_out, Q_material_out_sum)
        if math.fabs(Q_material_out_sum/Q_material_out - 1.0) > validation_tol*0.2:
            error_message += "\n heat flow with mass out of system from all boundary byintegral (ds) and sum is not equal\n"
            is_simulation_valid = False

        print('=============== conduction heat loss ===========')
        Q_conduction_all = assemble(_conductivity*dot(grad(T), surface_normal) * ds) * -1  # surface_normal pointing out
        Qk_chipend = assemble(_conductivity*dot(grad(T), surface_normal)*ds(chip_end_boundary_id)) * -1# flow out
        print('conductive heat flux from bottom and chipend bounadries', Q_bottom, Qk_chipend)   #should be set zero-flux
        Qk_holder_top = assemble(_conductivity*dot(grad(T), surface_normal) * ds(holder_top_boundary_id)) * -1 #flow out
        Q_conduction_work = assemble(_conductivity*dot(grad(T), surface_normal) * ds(work_subdomain_id)) * -1
        Q_conduction_chip = assemble(_conductivity*dot(grad(T), surface_normal) * ds(chip_subdomain_id)) * -1
        Q_conduction_tool = assemble(_conductivity*dot(grad(T), surface_normal) * ds(holder_subdomain_id)) * -1 +\
                                        assemble(_conductivity*dot(grad(T), surface_normal) * ds(cutter_subdomain_id)) * -1
        Qk_bottom = assemble(_conductivity*dot(grad(T), surface_normal) * ds(work_bottom_boundary_id)) * -1
        Qk_inlet= assemble(_conductivity*dot(grad(T), surface_normal) * ds(work_bottom_boundary_id + 1)) * -1
        Qk_outlet= assemble(_conductivity*dot(grad(T), surface_normal) * ds(work_bottom_boundary_id + 2)) * -1
        print('conductive heat at holder top, Qk_bottom, Qk_inlet, Qk_outlet, conduction all is',\
                    Qk_holder_top, Qk_bottom, Qk_inlet, Qk_outlet, Q_conduction_all)
        print('Q_conduction_tool, Q_conduction_work, Q_conduction_chip',\
                  Q_conduction_tool, Q_conduction_work, Q_conduction_chip)
        total_heat_loss += Qk_holder_top

        print('=============== heat loss ratios ===========')
        Q_baseline = total_heat_loss
        print("convection, radiation, Q_flowing_out_work, Q_bottom, Q_chipend, Q_holder_top")
        _ratios = np.array([cooling_total_HTC, R_total, Q_flowing_out_work, Q_bottom, Q_chipend, Qk_holder_top])/Q_baseline
        print(_ratios)

        Q_total_passing_friction_interface = (Qk_holder_top + cooling_tool)
        partition_coeff_shear = 1.0 -  Q_flowing_out_work/ shear_heat
        partition_coeff_friction = 1.0 -  Q_total_passing_friction_interface/ friction_heat
        print("partition_coeff_shear, partition_coeff_friction = ", partition_coeff_shear, partition_coeff_friction)

    print("ratio of total heat loss {}, heat generation {}, ratio {}"\
            .format(total_heat_loss, total_heat, total_heat_loss/total_heat))
    if math.fabs(total_heat_loss/total_heat - 1.0) > validation_tol:
        error_message += "\n heat generaton and heat loss are not equal\n"
        is_simulation_valid = False

    ts_after_post_processing = time.time()
    print('=============== time consumption ===========')
    print('preprocessing time', ts_after_preprocessing - ts_before_preprocessing)
    print('assembling and LA solving time', ts_before_postprocessing - ts_after_preprocessing)

    if using_debug:
        #plot(solver.boundary_facets)
        bfile = File(result_folder + "boundary.pvd")  # not working for parallel
        bfile << solver.boundary_facets
    if using_debug:
        ofile = File(result_folder + "HeatSource.pvd")  # not for parallel
        ofile << solver.body_source   #it is possible to save DG0 data

    if has_convective_velocity and using_debug:
        ofile = File(result_folder + "U.pvd")  # not for parallel, using hdf5 for parallel
        ofile << velocity

    if not is_batch_mode:
        if using_VTK:
            plot(solver.boundary_facets, title = "boundary id")
            plot(solver.subdomains, title = "subdomain id")
            plot(heat_source, title = "heat source density (W/m3)")
            #plot(solver.material['specific_heat_capacity'], title = "specific_heat_capacity")  # can not plot DG
            plot(velocity, title = "pseudo convective velocity")
            interactive()
            solver.plot()  # plt.legend() does not work: warnings.warn("No labelled objects found. "
        else:
            import matplotlib.pyplot as plt
            solver.plot()
            plt.show()
            #`paraview --data=T.pvd`

    if is_simulation_valid == False:
        print("======== the simulation is INVALID ! ==================")
        print(error_message)
        #if not using_workpiece_extra_extusion:
        #    sys.exit()
    else:
        print("======== the simulation is completed sucessfully ==================")


    thermal_number = metal_density * metal_cp_f(T_analytical_shear_zone) * cutting_speed * feed_thickness / metal_k_f(T_analytical_shear_zone)
    partition_coeff_eq = get_heat_partition(thermal_number)
    #also save thermal_number
    from datetime import datetime
    time_stamp  = datetime.utcnow().isoformat()
    parameters = [time_stamp, case_id, workpiece_type_id,
        feed_thickness, chip_thickness, cutter_angle_v, shear_angle, cutting_speed, F_cutting, F_thrush,
        shear_heat_thickness, friction_heat_thickness, tool_chip_interface_length,
        T_ambient,  T_analytical_shear_zone, T_analytical_friction_zone,
        Tmean_shear_zone, Tmax_shear_zone, Tmean_friction_zone, Tmax_friction_zone,
        thermal_number, partition_coeff_eq, partition_coeff_shear, partition_coeff_friction, is_simulation_valid]

    output = ",".join([str(v) for v in parameters])
    with open(result_filename, 'a') as wf:
        wf.write(output)
        wf.write('\n')

    if extracting_data:
        #sys.path.append('/opt/fenicstools/')  # this package must be installed
        from fenicstools import Probes
        dist = np.linspace(0, tool_chip_interface_length, 20)
        probing_pts = np.array([(math.sin(cutter_angle_v*pi/180)*l, math.cos(cutter_angle_v*pi/180)*l, 0) for l in dist])
        probes = Probes(probing_pts.flatten(), solver.function_space)
        probes(T)  # evaluate f at all probing points
        extracted_datafile = result_folder + datafile_root + '__' + param_name + '__' + str(param_value) + '.csv'
        np.savetxt(probes.array(), extracted_datafile)
        print(probes.array())


#################### utility functions #######################


def get_material_function(DG0, sub_domains, material_values):
    "DeprecationWarning: GenericVector.array() is being deprecated, use GenericVector.get_local()"
    # material_values is a list, index is given by sub_domains (CellFunction)
    # Create material function and assign values (vectorized using numpy)
    mat_func = Function(DG0)
    tmp = np.asarray(sub_domains.array(), dtype=np.int32)
    #assert len(material_values) == np.max(tmp) + 1
    # it assume DG0 always map to CellFunction
    mat_func.vector()[:] = np.choose(tmp, material_values)
    return mat_func


def get_material_property(DG0, subdomains, material_values):
    #material_values = (10., 100.)
    helper = np.asarray(subdomains.array(), dtype=np.int32)
    tcond = subdomains.copy()
    tcond.vector()[:] = np.choose(helper, material_values)  #should get a deep copy
    #
    V = DG0 # FunctionSpace(mesh, "DG", 0)
    u = Function(V)

    dm = V.dofmap()
    for cell in cells(DG0.mesh()):
       helper[dm.cell_dofs(cell.index())] = tcond[cell]

    u.vector()[:] = np.choose(helper, material_values)
    return u

# ratio of HTC to radiation on unit surface


solve_ht()

# -*- coding: utf-8 -*-

###
### This file is generated automatically by SALOME v8.5.0 with dump python functionality
###
import os
import sys
import salome

salome.salome_init()
theStudy = salome.myStudy

import salome_notebook
notebook = salome_notebook.NoteBook(theStudy)

# to find parameter_salome in the current folder
try:
    os.chdir('/media/sf_OneDrive/gitrepo/quasi_static_metal_cut')  # change to working dir
except:
    print('First of all, os.chdir to the folder where salome_parameter.py is located')
    sys.exit()
    #os.chdir('/media/OneDrive/Fenics/metal_cut/')
from math import sin, cos, tan, pi

from parameter_salome import *  # restart salome if any change in this parameter file

#/opt/SALOME-8.5.0-UB16.04-SRC/runSalome -t -b /myScript
## the sequence of gmsh options is important
#gmsh4 -format msh2   -o metal_cut.msh -save Mesh_1.med
#dolfin-convert metal_cut.msh metal_cut.xml
##########################################
#derived parameter



###
### GEOM component
###


import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS


geompy = geomBuilder.New(theStudy)

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
OZ_minus = geompy.MakeVectorDXDYDZ(0, 0, -1)

#============================
p_origin = geompy.MakeVertex(0, 0, 0)
x0 = 0
y0 = 0
p_chip_turning = geompy.MakeVertex(p_chip_turning_x, p_chip_turning_y, 0)

#Geometry preferences default tol -9
#geom_tolerance = 1e-11
#shape1_lt = geompy.LimitTolerance(shape1, geom_tolerance)
#shape2_lt = geompy.LimitTolerance(shape2, geom_tolerance)

p_work_left_top = geompy.MakeVertex(work_left_x, feed_thickness, 0)
p_work_pressing_y0 = geompy.MakeVertex(pressing_heat_x, 0, 0)

p_cutter_bottom = p_work_pressing_y0
#p_cutter_bottom = geompy.MakeVertex(0, 1e-6, 0)
p_cutter_H = geompy.MakeVertex(cutter_H*math.tan(holder_angle_v*math.pi/180), cutter_H, 0)
p_cutter_W = geompy.MakeVertex(cutter_W, cutter_W*math.tan(cutter_angle_h*math.pi/180), 0)
#p_cutter_HW = geompy.MakeVertex(cutter_W, cutter_H, 0)

p_ccc = geompy.MakeVertex(chip_ccc_x, chip_ccc_y, 0)
p_chip_left_mid = geompy.MakeVertex(p_chip_turning_x + chip_ccc_x, p_chip_turning_y + chip_ccc_y, 0)
p_top_o = geompy.MakeVertex(chip_top_x, chip_top_y, 0)
p_top_m = geompy.MakeVertex(chip_top_x+friction_b_shift_x, chip_top_y, 0) # intermediate point to enable quad mesh, firction heat thickness, 
p_top_t = geompy.MakeVertex(chip_top_x - chip_thickness/cos(cutter_angle_v*pi/180), chip_top_y, 0)

Line_chip_top = geompy.MakePolyline([p_top_t,  p_top_m, p_top_o])
Line_chip_right = geompy.MakeLineTwoPnt(p_ccc, p_top_o)
p_friction_top = geompy.MakeVertex(chip_ccc_x + friction_t_shift_x, chip_ccc_y + friction_t_shift_y, 0)
Line_chip_middle = geompy.MakePolyline([p_chip_left_mid, p_friction_top])

split_tool_ = [Line_chip_middle]
if using_fillet_shear_zone:  # fillet
    p_fillet_origin_center = geompy.MakeVertex(fillet_shift_x, fillet_shift_y, 0)
    p_fillet_origin_1 = geompy.MakeVertex(fillet_shift_x, 0, 0)
    p_fillet_origin_2 = geompy.MakeVertex(fillet_p2_shift_x, fillet_p2_shift_y, 0)
    fillet_origin = geompy.MakeArcCenter(p_fillet_origin_center, p_fillet_origin_1, p_fillet_origin_2, False)
    #tolerance ?
    p_fillet_t_center = geompy.MakeVertex(p_chip_turning_x + fillet_shift_x, p_chip_turning_y + fillet_shift_y, 0)
    p_fillet_t_1 = geompy.MakeVertex(p_chip_turning_x + fillet_shift_x, p_chip_turning_y, 0)
    p_fillet_t_2 = geompy.MakeVertex(p_chip_turning_x + fillet_p2_shift_x, p_chip_turning_y + fillet_p2_shift_y, 0)
    fillet_turning = geompy.MakeArcCenter(p_fillet_t_center, p_fillet_t_1, p_fillet_t_2, False)
    #doc(geompy.MakeArcCenter)

    #uisng subshape of arc, vertex, otherwise the arc and line are not closed
    fillet_origin_vertex_3 = geompy.GetSubShape(fillet_origin, [3])
    fillet_turning_vertex_3 = geompy.GetSubShape(fillet_turning, [3])

    p_chip_shear_end_o = fillet_origin_vertex_3
    p_chip_shear_end_t = fillet_turning_vertex_3

    p_friction_bottom = geompy.MakeVertex(fillet_p2_shift_x + friction_b_shift_x, fillet_p2_shift_y + friction_b_shift_y, 0)
    p_cutter_mapping_bottom = geompy.MakeVertex(fillet_p2_shift_x + mapping_bc_gap, fillet_p2_shift_y, 0)
    #p_friction_top is defined at top level at the beginning
    Line_chip_left = geompy.MakePolyline([p_chip_shear_end_t, p_chip_left_mid, p_top_t])
    Line_chip_friction_interface = geompy.MakeLineTwoPnt(p_chip_shear_end_o, p_ccc)

    Line_to_1 = geompy.MakeLineTwoPnt(p_fillet_t_1, p_fillet_origin_1)
    Line_to_2 = geompy.MakePolyline([p_chip_shear_end_t, p_friction_bottom, p_chip_shear_end_o])
    Line_work_top_left = geompy.MakeLineTwoPnt(p_work_left_top, p_fillet_t_1)

    Line_work_y0 = geompy.MakeLineTwoPnt(p_origin, p_fillet_origin_1)  #
    Line_cutter_work_joint = geompy.MakeLineTwoPnt(p_chip_shear_end_o, p_origin)
    #===================================
    if cutting_fillet_hole:
        cutter_work_joint_outline = [Line_cutter_work_joint, fillet_origin, Line_work_y0] #cut a hole in the outer loop, it works!  keep it as a shear heating zone
    else:
        cutter_work_joint_outline = []
        split_tool_ += [fillet_origin, Line_work_y0]  #Line_cutter_work_joint, 
    chip_outline = [ fillet_turning, Line_chip_left, Line_chip_top, Line_chip_right]
else:
    Line_turning = geompy.MakeLineTwoPnt(p_chip_turning, p_origin)  #mapping bc on chip side

    p_friction_bottom = geompy.MakeVertex(shear_thickness_shift_x + friction_b_shift_x, shear_thickness_shift_y + friction_b_shift_y, 0)
    p_cutter_mapping_bottom = geompy.MakeVertex(shear_thickness_shift_x + +mapping_bc_gap, shear_thickness_shift_y, 0)

    p_chip_shear_end_o =  geompy.MakeVertex(shear_thickness_shift_x, shear_thickness_shift_y, 0)
    p_chip_shear_end_t = geompy.MakeVertex(p_chip_turning_x + shear_thickness_shift_x, p_chip_turning_y + shear_thickness_shift_y, 0)
    #Line_shear_o = geompy.MakeLineTwoPnt(p_chip_shear_end_o, p_origin)
    Line_shear_interface = geompy.MakePolyline([p_chip_shear_end_t, p_friction_bottom, p_chip_shear_end_o, p_origin])

    Line_chip_left = geompy.MakePolyline([p_chip_turning, p_chip_shear_end_t, p_chip_left_mid, p_top_t])
    Line_chip_friction_interface = geompy.MakeLineTwoPnt(p_chip_shear_end_o, p_ccc)
    # mapping bc for work size, double shear layers
    p_wm_top = geompy.MakeVertex(p_chip_turning_x - mapping_bc_gap, p_chip_turning_y, 0)
    p_wm_bottom = geompy.MakeVertex(-mapping_bc_gap, y0 , 0)
    Line_wm_bc = geompy.MakeLineTwoPnt(p_wm_top, p_wm_bottom)  # workpiece side

    p_work_sh_o = geompy.MakeVertex(work_shear_thickness_shift_x - mapping_bc_gap, 0, 0)
    _work_wm_shift = p_chip_turning_x + work_shear_thickness_shift_x - mapping_bc_gap
    p_work_sh_t = geompy.MakeVertex(_work_wm_shift, feed_thickness, 0)
    Line_work_sh = geompy.MakePolyline([p_wm_bottom, p_work_sh_o, p_work_sh_t]) # split work and shear zone

    Line_work_top_left = geompy.MakeLineTwoPnt(p_work_left_top, p_wm_top)
    Line_wm_bottom_gap = geompy.MakeLineTwoPnt(p_wm_bottom, p_origin)
    cutter_work_joint_outline = []
    chip_outline = [Line_turning, Line_chip_left, Line_chip_top, Line_chip_right]

if using_chip_cutter_mapping_bc or using_chip_cutter_transition:
    p_cutter_mapping_top = geompy.MakeVertex(chip_ccc_x + mapping_bc_gap, chip_ccc_y, 0)
    #p_cutter_mapping_bottom = geompy.MakeVertex(shear_thickness_shift_x+mapping_bc_gap, shear_thickness_shift_y , 0)
    #cutter_mapping_bottom = geompy.MakeVertex(x0+mapping_bc_gap, y0 , 0)

    Line_cutter_mapping_bc = geompy.MakeLineTwoPnt(p_cutter_mapping_bottom, p_cutter_mapping_top)  # cutter side
    #Line_chip_friction_interface = geompy.MakeLineTwoPnt(p_ccc, p_origin)
    Line_chip_cutter_mapping_bottom_gap = geompy.MakeLineTwoPnt(p_chip_shear_end_o, p_cutter_mapping_bottom)
    Line_chip_cutter_mapping_top_gap = geompy.MakeLineTwoPnt(p_ccc, p_cutter_mapping_top)
    geompy.addToStudy(Line_chip_cutter_mapping_bottom_gap, 'Line_chip_cutter_mapping_bottom_gap' )
    #geompy.addToStudy(Line_chip_cutter_mapping_top_gap, 'Line_chip_cutter_mapping_top_gap' )
    geompy.addToStudy(Line_cutter_mapping_bc, 'Line_cutter_mapping_bc' )
    
    chip_outline += [Line_chip_cutter_mapping_top_gap]
    split_tool_ += [Line_chip_cutter_mapping_bottom_gap, Line_cutter_mapping_bc]
    Line_cutter_left = geompy.MakeLineTwoPnt(p_cutter_mapping_top, p_cutter_H)
else:
    Line_cutter_left = geompy.MakeLineTwoPnt(p_ccc, p_cutter_H)

#make the work piece profile
if workpiece_type_id == 2:  # work piece is only feed_thickness thick, with fillet structure
    assert using_fillet_shear_zone
    p_work_left_bottom = geompy.MakeVertex(work_left_x, 0, 0)
    Line_work_left = geompy.MakeLineTwoPnt(p_work_left_top, p_work_left_bottom)
    #Line_work_bottom = geompy.MakeLineTwoPnt(p_work_left_bottom, p_origin)
    Line_work_bottom = geompy.MakeLineTwoPnt(p_work_left_bottom, p_fillet_origin_1)  # fillet

    work_outline = [Line_work_bottom, Line_work_left, Line_work_top_left]

elif workpiece_type_id == 0:  # with a rect work piece, with or without extra workpiece in Z axis

    p_work_left_y0 = geompy.MakeVertex(work_left_x, y0 , 0)
    #p_work_right_y0 = geompy.MakeVertex(work_right_x, y0 , 0)
    p_work_left_bottom = geompy.MakeVertex(work_left_x, work_bottom_y , 0)
    p_work_right_top = geompy.MakeVertex(work_right_x, y0, 0)
    p_work_pressing_right = geompy.MakeVertex(pressing_heat_x + pressing_heat_thickness, y0, 0)  # shift to make sure near the pressing point, mesh size on the right is not too big
    p_work_right_bottom = geompy.MakeVertex(work_right_x, work_bottom_y, 0)
    p_work_bottom_pressing = geompy.MakeVertex(pressing_heat_x, work_bottom_y, 0)
    p_work_bottom_t = geompy.MakeVertex(0, work_bottom_y, 0)  #
    p_work_bottom_o = geompy.MakeVertex(pressing_heat_x, work_bottom_y, 0)

    Line_work_left = geompy.MakePolyline([p_work_left_top, p_work_left_y0, p_work_left_bottom])
    Line_work_top_right = geompy.MakePolyline([p_work_right_top, p_work_pressing_right, p_work_pressing_y0])
    Line_work_right = geompy.MakePolyline([p_work_right_bottom, p_work_right_top])

    Line_work_cutter_interface = geompy.MakeLineTwoPnt(p_cutter_bottom, p_origin)
    split_tool_ += [Line_work_cutter_interface]
    if using_fillet_shear_zone:
        #pass
        Line_work_top_left = geompy.MakeLineTwoPnt(p_work_left_top, p_fillet_t_1)
        Line_work_bottom = geompy.MakePolyline([p_work_left_bottom, p_work_bottom_t, p_work_bottom_o, p_work_right_bottom])
        work_outline = [Line_work_top_right, Line_work_right, Line_work_bottom, Line_work_left, Line_work_top_left]
        split_tool_ += [Line_to_1, Line_to_2]  # Line_work_y0 may be not used, fillet_origin,
    else: # mapping
        Line_work_top_left = geompy.MakePolyline([p_work_left_top, p_work_sh_t, p_wm_top])
        Line_work_bottom = geompy.MakePolyline([p_work_left_bottom, p_work_bottom_t, p_work_bottom_o, p_work_bottom_pressing, p_work_right_bottom])
        work_outline = [Line_work_top_right, Line_work_right, Line_work_bottom, Line_work_left, Line_work_top_left, Line_wm_bc, Line_wm_bottom_gap]
        cutter_work_joint_outline = []
        #
        split_tool_ += [Line_shear_interface]  #shear_end and chip interface
        if using_double_shear_heat_layer:
            split_tool_ += [Line_work_sh]
    ################### ########################
    if using_workpiece_extra_extusion:
        p_work_extra_right_top = geompy.MakeVertex(work_right_x, feed_thickness, 0)
        Line_extra_workpiece = geompy.MakePolyline([p_wm_bottom, p_work_right_top, p_work_extra_right_top, p_wm_top])
        Face_extra_work = geompy.MakeFaceWires([Line_extra_workpiece,  Line_wm_bc], 1)
        # may be totally seperate solid

    Line_work_left_y0_interface = geompy.MakeLineTwoPnt(p_work_left_y0, p_origin)
    split_tool_ += [Line_work_left_y0_interface]  # split workpiece into 2 zones for better meshing
else:
    pass


#######################################################
# 3 lines to split friction heating layer
p_pressing_shift_left = geompy.MakeVertex(x0, y0-pressing_heat_thickness, 0)
p_pressing_shift_right = geompy.MakeVertex(pressing_heat_x, y0-pressing_heat_thickness, 0)
Line_pressing_polyline = geompy.MakePolyline([p_origin, p_pressing_shift_left, p_pressing_shift_right, p_work_pressing_y0])
#if has_pressing_zone:
split_tool_ += [Line_pressing_polyline] 

#Line_friction_top = geompy.MakeLineTwoPnt(p_ccc, p_friction_top)
#Line_friction_right = geompy.MakeLineTwoPnt(p_friction_top, p_friction_bottom)
#Line_friction_bottom = geompy.MakeLineTwoPnt(fillet_origin_vertex_3, p_friction_bottom)
Line_friction_polyline = geompy.MakePolyline([p_ccc, p_friction_top, p_friction_bottom, p_chip_shear_end_o])
geompy.addToStudy(Line_friction_polyline, 'Line_friction_polyline')

Line_cutter_bottom = geompy.MakeLineTwoPnt(p_cutter_bottom, p_cutter_W)
Line_cutter_HW = geompy.MakeLineTwoPnt(p_cutter_H, p_cutter_W)
geompy.addToStudy( Line_cutter_left, 'Line_cutter_left')
geompy.addToStudy( Line_cutter_HW, 'Line_cutter_HW')
#geompy.addToStudy( , '')

if using_3D_holder:
    holder_z = 0 # cutter_thickness
else:
    holder_z = 0
holder_top_y = holder_H
p_holder_H = geompy.MakeVertex(holder_H*math.tan(holder_angle_v*math.pi/180), holder_H, holder_z)
p_holder_W = geompy.MakeVertex(holder_W, holder_W*math.tan(cutter_angle_h*math.pi/180), holder_z)
p_holder_HW = geompy.MakeVertex(holder_W, holder_H, holder_z)
p_holder_bottom = geompy.MakeVertex(holder_bottom_x, holder_bottom_y, holder_z)

Line_holder_top = geompy.MakeLineTwoPnt(p_holder_H, p_holder_HW)
Line_holder_right = geompy.MakeLineTwoPnt(p_holder_W, p_holder_HW)
Line_cutter_holder_left = geompy.MakeLineTwoPnt(p_holder_bottom, p_cutter_H)
Line_cutter_holder_bottom = geompy.MakeLineTwoPnt(p_holder_bottom, p_cutter_W)
"""
if using_3D_holder:
    Line_holder_bottom = geompy.MakeLineTwoPnt(p_holder_W, p_holder_bottom)
    Line_holder_left = geompy.MakeLineTwoPnt(p_holder_H, p_holder_bottom)
"""
Line_holder_bottom = geompy.MakeLineTwoPnt(p_holder_W, p_cutter_W)
Line_holder_left = geompy.MakeLineTwoPnt(p_holder_H, p_cutter_H)

## #############################################
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )

#geompy.addToStudy( p_work_left_top, 'p_work_left_top' )
#geompy.addToStudy( p_work_right_top, 'p_work_right_top' )
#geompy.addToStudy( p_work_left_bottom, 'p_work_left_bottom' )
#geompy.addToStudy( p_work_right_bottom, 'p_work_right_bottom' )
geompy.addToStudy( p_origin, 'p_origin' )

if using_fillet_shear_zone:
    #geompy.addToStudy( p_chip_turning, 'p_chip_turning' )
    geompy.addToStudy( p_fillet_t_center, 'p_fillet_t_center' )
    geompy.addToStudy( p_fillet_t_1, 'p_fillet_t_1' )
    geompy.addToStudy( p_fillet_t_2, 'p_fillet_t_2' )
    geompy.addToStudy( fillet_turning, 'fillet_turning' )

    geompy.addToStudy( p_fillet_origin_center, 'p_fillet_origin_center' )
    geompy.addToStudy( p_fillet_origin_1, 'p_fillet_origin_1' )
    geompy.addToStudy( p_fillet_origin_2, 'p_fillet_origin_2' )
    geompy.addToStudy( fillet_origin, 'fillet_origin' )

    geompy.addToStudyInFather( fillet_origin, fillet_origin_vertex_3, 'fillet_origin:vertex_3' )
    geompy.addToStudyInFather( fillet_turning, fillet_turning_vertex_3, 'fillet_turning:vertex_3' )
    #
    geompy.addToStudy( Line_work_y0 , 'Line_work_y0 ' )
    geompy.addToStudy( Line_to_1, 'Line_to_1' )
    geompy.addToStudy( Line_to_2, 'Line_to_2' )

geompy.addToStudy( Line_chip_top, 'Line_chip_top')
geompy.addToStudy( Line_chip_friction_interface, 'Line_chip_friction_interface')
geompy.addToStudy( Line_chip_right, 'Line_chip_right')
geompy.addToStudy( Line_chip_left, 'Line_chip_left')

#geompy.addToStudy( Line_friction_top, 'Line_friction_top')
#geompy.addToStudy( Line_friction_right, 'Line_friction_right')
geompy.addToStudy( Line_friction_polyline, 'Line_friction_polyline')
##################################

if using_3D_holder:

    Face_1 = geompy.MakeFaceWires( cutter_work_joint_outline + work_outline + chip_outline
        + [Line_cutter_left, Line_cutter_HW, Line_cutter_bottom], 1)  #

    geompy.addToStudy( Face_1, 'Face_1' )
    split_tool = split_tool_ + [Line_chip_friction_interface, Line_friction_polyline]
    split_tool += [Line_cutter_holder_left, Line_cutter_holder_bottom]
    Partition_1 = geompy.MakePartition([Face_1], split_tool, [], [], geompy.ShapeType["FACE"], 0, [], 0)

    #make holder face, not shared with cutter
    Face_holder = geompy.MakeFaceWires([Line_holder_left, Line_holder_top, Line_holder_right,  Line_holder_bottom, Line_cutter_HW], 1)
    geompy.addToStudy( Face_holder, 'Face_holder' )

else:
    Face_1 = geompy.MakeFaceWires( cutter_work_joint_outline + work_outline + chip_outline
        + [Line_cutter_left, Line_holder_left, Line_holder_top, Line_holder_right,  Line_holder_bottom, Line_cutter_bottom], 1)  #

    geompy.addToStudy( Face_1, 'Face_1' )
    split_tool = split_tool_ + [Line_cutter_HW, Line_chip_friction_interface, Line_friction_polyline]
    Partition_1 = geompy.MakePartition([Face_1], split_tool, [], [], geompy.ShapeType["FACE"], 0, [], 0)


if using_fillet_shear_zone:  #for some bug, need to split again!
    Partition_1 = geompy.MakePartition([Partition_1], split_tool, [], [], geompy.ShapeType["FACE"], 0, [], 0)
geompy.addToStudy( Partition_1, 'Partition_1' )

#capture two workpiece face, and extrusion to minus z direction

def near(pa, pb):
    if type(pa) == type(1.0) and type(pb) == type(1.0):
        return math.fabs(pa-pb) < 1e-6
    else:
        A = geompy.PointCoordinates(pa)
        B = geompy.PointCoordinates(pb)
        dist2 = sum([(A[i]-B[i])*(A[i]-B[i]) for i in range(3)])
        return dist2 < 1e-5

def edge_length(edge):
    [pa, pb] = geompy.ExtractShapes(edge, geompy.ShapeType["VERTEX"], True)
    A = geompy.PointCoordinates(pa)
    B = geompy.PointCoordinates(pb)
    dist2 = sum([(A[i]-B[i])*(A[i]-B[i]) for i in range(3)])
    return dist2**0.5

def get_shape_list_by_center_condition(geom, cond, shape_type = "EDGE"):
    edgelist = geompy.ExtractShapes(geom, geompy.ShapeType[shape_type], True)
    elist = []
    for i,e in enumerate(edgelist):
        cg = geompy.PointCoordinates(geompy.MakeCDG(e))  #center of gravity, only for solid? return type?
        if cond(cg):  #z-coordinate
            elist.append(e)
            #geompy.addToStudyInFather( geom, e, "Edge_{}".format(i))
            # geompy.getCoordinate(vertex)  # return a 3elm tuple?
    return elist

def add_group_by_center_condition(geom, cond, shape_type, group_name):
    sl = get_shape_list_by_center_condition(geom, cond, shape_type = shape_type)
    shape_group = geompy.CreateGroup(geom, geompy.ShapeType[shape_type])
    geompy.UnionList(shape_group, sl)
    geompy.addToStudyInFather(geom, shape_group, group_name)
    return shape_group

def get_edge_list_by_length_condition(geom, cond):
    edgelist = geompy.ExtractShapes(geom, geompy.ShapeType["EDGE"], True)
    elist = []
    for i,e in enumerate(edgelist):
        cg = geompy.PointCoordinates(geompy.MakeCDG(e))  #center of gravity, only for solid? return type?
        l = edge_length(e)
        if cond(cg, l):  #z-coordinate
            elist.append(e)
            #geompy.addToStudyInFather( geom, e, "Edge_{}".format(i))
            # geompy.getCoordinate(vertex)  # return a 3elm tuple?
    return elist

#if not using_3D:
#explosion to get subshape

#facelist = geompy.ExtractShapes(Partition_1, geompy.ShapeType["FACE"], True)
#for i,f in enumerate(facelist):
#    geompy.addToStudyInFather( Partition_1, f, "Face_{}".format(i))

"""
props = geompy.BasicProperties(box)
print "\nBox 100x30x100 Basic Properties:"
print " Wires length: ", props[0]
print " Surface area: ", props[1]
print " Volume      : ", props[2]
"""

_angle_cutter = math.pi/2.0 - cutter_angle_v*math.pi/180.0
_angle_phi = math.pi - angle_phi*math.pi/180.0
_angle_shear = math.pi - angle_phi*math.pi/180.0
if using_3D:
    ###################################
    Extrusion_1 = geompy.MakePrismVecH(Partition_1, OZ, extrusion_thickness)
    geompy.addToStudy( Extrusion_1, 'Extrusion_1' )
    body_list = [Extrusion_1]
    if using_3D_holder:
        cond_faces_holder = lambda p: (p[1] > cutter_H*0.3 and p[0] > cutter_W*0.3 and math.fabs(p[2]) < 1e-5)
        Face_holder_cutter = get_shape_list_by_center_condition(Extrusion_1, cond_faces_holder, "FACE")
        faces_holder = [Face_holder] + Face_holder_cutter
        if len(faces_holder) < 2:
            print('cutter and holder overlapping face is not identified, exit')
            sys.exit(-1)
        #faces_holder_merged = geompy.MakeFuseList(faces_holder, True, True)  # holder and cutter are not connected
        faces_holder_merged = geompy.MakePartition(faces_holder, [], [], [], geompy.ShapeType["FACE"], 0, [], 0)
        # two holder zones are not connected
        Extrusion_holder = geompy.MakePrismVecH(faces_holder_merged, OZ_minus, holder_thickness)
        geompy.addToStudy( Extrusion_holder, 'Extrusion_holder')
        body_list.append(Extrusion_holder)
    if using_workpiece_extra_extusion:
        #print('type of Extrusion_1', type(Extrusion_1)) #GEOM_Object compound
        #faces_workpiece = []
        #cond_faces_workpiece = lambda p: (p[1] < 0 or (p[0] < work_left_x * 0.3 and p[1]<feed_thickness*0.7)) and math.fabs(p[2]) < 1e-5
        cond_faces_workpiece = lambda p: (p[1] < 0 or  math.atan2(p[1], p[0]) > _angle_shear) and math.fabs(p[2]) < 1e-5
        faces_workpiece = get_shape_list_by_center_condition(Extrusion_1, cond_faces_workpiece, "FACE")
        #faces_workpiece.append(Face_extra_work)  # not conected to other part, need a fusion
        Extrusion_3 = geompy.MakePrismVecH(Face_extra_work, OZ_minus, math.fabs(workpiece_extra_extusion_thickness))
        Extrusion_2 = geompy.MakePrismVecH(geompy.MakeCompound(faces_workpiece), OZ_minus, math.fabs(workpiece_extra_extusion_thickness))
        geompy.addToStudy( Extrusion_3, 'Extrusion_3' )
        geompy.addToStudy( Extrusion_2, 'Extrusion_2' )
        #Fuse_1 = geompy.MakeFuseList([Extrusion_3, Extrusion_2], True, True)
        #geompy.addToStudy( Fuse_1, 'Fuse_1' )
        body_list.append(Extrusion_2)
    if len(body_list) > 1:
        Compound_1 = geompy.MakePartition(body_list, Face_holder_cutter, [], [], geompy.ShapeType["SOLID"], 0, [], 0)
        #Compound_1 = geompy.MakeCompound(body_list)  # MakeCompound: not connnected
    else:
        Compound_1 = Extrusion_1
    geompy.addToStudy(Compound_1, 'Compound_1')
    domain = Compound_1
    #is there a way to detect whether a point is in solid?
    vlist = geompy.ExtractShapes(domain, geompy.ShapeType["SOLID"], True)
    solid_type = "SOLID"
    boundary_type = "FACE"
else:
    domain = Partition_1
    #is there a way to detect whether a point is in solid?
    vlist = geompy.ExtractShapes(domain, geompy.ShapeType["FACE"], True)
    solid_type = "FACE"
    boundary_type = "EDGE"


shear_heat_zones = []
work_zones = []
cutter_zones = []
chip_zones = []
holder_zones = []
cutter_chip_transition_zones = []
_debug = False

for i,v in enumerate(vlist):
    cm = geompy.MakeCDG(v)
    if cm is None:
        raise RuntimeError("can not makeCDG for the solid")
    else:
        coords = geompy.PointCoordinates(cm)
        if _debug: print("\nCentre of gravity for the solid_{} is".format(i), coords)
    props = geompy.BasicProperties(v)
    if using_3D:
        _solid_volume = props[2]
    else:
        _solid_volume = props[1]*extrusion_thickness

    volume_near = lambda v1, v2:  math.fabs(v2/v1 - 1) < 0.05
    _angle = math.atan2(coords[1], coords[0])
    if _debug: print("Volume of solid_{} is:".format(i), _solid_volume)

    if volume_near(_solid_volume, friction_heat_volume) and _angle > _angle_cutter and math.fabs(_angle - _angle_cutter) < 0.5:
        Solid_friction = v
        print("identified as Solid_friction")
    elif using_workpiece_extra_extusion and using_3D and coords[2] < 0:
        work_zones.append(v)
        print("identified as Solid_work extra extension zone")
    elif volume_near(_solid_volume, shear_heat_volume)  and math.fabs(_angle - _angle_shear) < 0.1:
        shear_heat_zones.append(v)
        print("identified as Solid_shear single zone")
    elif using_fillet_shear_zone and volume_near(_solid_volume, hole_volume):
        cutter_zones.append(v)
        print("identified as hole and added to cutter zones")
    elif volume_near(_solid_volume, pressing_heat_volume) and coords[1] < 0:
        if has_pressing_zone:
            Solid_pressing = v
        else:
            work_zones.append(v)
        print("identified as pressing zone")
    elif volume_near(_solid_volume, friction_mapping_gap_volume) and _angle < _angle_cutter:
        if using_chip_cutter_transition:
            cutter_chip_transition_zones.append(v)  # diff thermal properties
            print("identified as transition layer between cutter and chip")
        elif using_chip_cutter_mapping_bc:
            print("identified as mapping gap between cutter and chip, ignore it")
        else:
            pass

    elif volume_near(_solid_volume, shear_heat_volume*0.5):
        shear_heat_zones.append(v)  #using_double_shear_heat_layer
        print("identified as Solid_shear double-zone")
    elif _angle > _angle_cutter and _angle < _angle_phi  and coords[1] > feed_thickness:
        chip_zones.append(v)
        print("identified as Solid_chip")
    else:
        if coords[1] < 0 or ( _angle > _angle_shear): #
            work_zones.append(v)
            print("identified as Solid_work")
        elif ( _angle < _angle_cutter) and coords[1] < cutter_H*0.6:
            if using_3D:
                if coords[2] > 0:
                    cutter_zones.append(v)
                    print("identified as Solid_cutter")
                else:
                    holder_zones.append(v)
                    print("identified as Solid_holder overlapping with cutter")
            else:
                cutter_zones.append(v)
                print("identified as Solid_cutter")
        elif coords[1] > cutter_H * 0.6:  # is that works for 
            holder_zones.append(v)
            print("identified as Solid_holder")
        else:
            print("Volume of solid_{} is:".format(i), _solid_volume)
            print("\nCentre of gravity for the solid_{} is".format(i), coords)
            raise RuntimeError("solid is not identified")

Solid_work = geompy.CreateGroup(domain, geompy.ShapeType[solid_type])
geompy.UnionList( Solid_work, work_zones)
Solid_shear = geompy.CreateGroup(domain, geompy.ShapeType[solid_type])
geompy.UnionList( Solid_shear, shear_heat_zones)
Solid_cutter= geompy.CreateGroup(domain, geompy.ShapeType[solid_type])
geompy.UnionList( Solid_cutter, cutter_zones)
Solid_holder= geompy.CreateGroup(domain, geompy.ShapeType[solid_type])
geompy.UnionList( Solid_holder, holder_zones)
Solid_chip= geompy.CreateGroup(domain, geompy.ShapeType[solid_type])
geompy.UnionList(Solid_chip, chip_zones)
Solid_cutter_chip_transition= geompy.CreateGroup(domain, geompy.ShapeType[solid_type])
geompy.UnionList(Solid_cutter_chip_transition, cutter_chip_transition_zones)

geompy.addToStudyInFather( domain, Solid_shear, 'Solid_shear' )
geompy.addToStudyInFather( domain, Solid_chip, 'Solid_chip' )
geompy.addToStudyInFather( domain, Solid_friction, 'Solid_friction' )
geompy.addToStudyInFather( domain, Solid_work, 'Solid_work' )
geompy.addToStudyInFather( domain, Solid_cutter, 'Solid_cutter' )
geompy.addToStudyInFather( domain, Solid_holder, 'Solid_holder' )
geompy.addToStudyInFather( domain, Solid_cutter_chip_transition, 'Solid_cutter_chip_transition' )

# independent of 2D or 3D, in same case,  heat thickness is not counted as short edge
velist = get_edge_list_by_length_condition(domain, lambda p, l: math.fabs(l) < small_edge_length_threshold  and math.fabs(p[2]) < 1e-8 )
ShortEdgeGroup = geompy.CreateGroup(domain, geompy.ShapeType["EDGE"])
geompy.UnionList(ShortEdgeGroup, velist)
geompy.addToStudyInFather(domain, ShortEdgeGroup, 'ShortEdgeGroup' )

cond_triangle = lambda p: math.fabs(p[2]) < 1e-6 and (coords[1]<0 or math.atan2(p[1], p[0]) < _angle_cutter)
fl_triangle = get_shape_list_by_center_condition(domain, cond_triangle, shape_type = 'FACE')  # for mesh control non-quad
FaceGroupTriangleMesh = geompy.CreateGroup(domain, geompy.ShapeType["FACE"])
geompy.UnionList(FaceGroupTriangleMesh, fl_triangle)
geompy.addToStudyInFather( domain, FaceGroupTriangleMesh, 'FaceGroupTriangleMesh' )

#2D also have z coordinates? friction heat zone may be triangle shape
cond_quad = lambda p: math.fabs(p[2]) < 1e-6 and math.atan2(p[1], p[0]) < math.pi * 0.75 and math.atan2(p[1], p[0]) > _angle_cutter
FaceGroupQuadMesh = add_group_by_center_condition(domain, cond_quad, "FACE", 'FaceGroupQuadMesh')

if using_3D:
    velist = get_shape_list_by_center_condition(domain, lambda p: math.fabs(p[2] - 0.5*extrusion_thickness) < 1e-5)
    VerticalEdgeGroup = geompy.CreateGroup(domain, geompy.ShapeType["EDGE"])
    if using_3D_holder:
        velist += get_shape_list_by_center_condition(domain, lambda p: math.fabs(p[2] - holder_middle_z) < 1e-5)
    geompy.UnionList(VerticalEdgeGroup, velist)
    geompy.addToStudyInFather(domain, VerticalEdgeGroup, 'VerticalEdgeGroup' )
    if using_workpiece_extra_extusion:
        cond_wve = lambda p: near(p[2], workpiece_min_z*0.5)
        edges_wve = add_group_by_center_condition(domain, cond_wve, "EDGE", 'extra_workpiece_vertical_edges')

if using_3D:
    on_boundary = lambda p: (near(p[2], workpiece_min_z*0.5) or near(p[2], extrusion_thickness*0.5))
else:
    on_boundary = lambda p: True

cond_inlet = lambda p:  on_boundary(p) and near(p[0], work_left_x) and p[1]< feed_thickness
cond_outlet = lambda p:  on_boundary(p) and near(p[0], work_right_x) and p[1]< feed_thickness
cond_chipend = lambda p: on_boundary(p) and near(p[1], chip_top_y)
cond_holder_top = lambda p: on_boundary(p) and near(p[1], holder_top_y)
cond_work_bottom = lambda p: on_boundary(p) and near(p[1], work_bottom_y)
if using_mapping_bc:
    cond_mapping_work = lambda p: on_boundary(p) and math.fabs(math.atan2(p[1], p[0]+mapping_bc_gap) - _angle_shear) < 0.001
    mapping_work = add_group_by_center_condition(domain, cond_mapping_work, "FACE", 'mapping_work')
    cond_mapping_chip = lambda p: on_boundary(p) and math.fabs(math.atan2(p[1], p[0]) - _angle_shear) < 0.001
    mapping_chip = add_group_by_center_condition(domain, cond_mapping_chip, "FACE", 'mapping_chip')

inlet = add_group_by_center_condition(domain, cond_inlet, boundary_type, 'inlet')
outlet = add_group_by_center_condition(domain, cond_outlet, boundary_type, 'outlet')
chipend = add_group_by_center_condition(domain, cond_chipend, boundary_type, 'chipend')
holder_top = add_group_by_center_condition(domain, cond_holder_top, boundary_type, 'holder_top')
work_bottom = add_group_by_center_condition(domain, cond_work_bottom, boundary_type, 'work_bottom')

#################################################################################

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(True)

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New(theStudy)


if using_3D:
    Mesh_1 = smesh.Mesh(domain)
    Regular_1D = Mesh_1.Segment()
    Number_of_Segments_default = Regular_1D.NumberOfSegments(line_segment_nb_default)#
    smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
    MEFISTO_2D = Mesh_1.Triangle(algo=smeshBuilder.MEFISTO)
    smesh.SetName(MEFISTO_2D.GetAlgorithm(), 'MEFISTO_2D')
    Prism_3D = Mesh_1.Prism()
    smesh.SetName(Prism_3D.GetAlgorithm(), 'Prism_3D')

    #submesh the face partition
    if using_3D_holder:
        # still not working
        flist = get_shape_list_by_center_condition(Extrusion_1, lambda p: math.fabs(p[2]) < 1e-5)
        Z0FaceGroup = geompy.CreateGroup(Extrusion_1, geompy.ShapeType["FACE"])
        geompy.UnionList(Z0FaceGroup, flist)
        geompy.addToStudyInFather(Extrusion_1, Z0FaceGroup, 'Z0FaceGroup' )
        face_domain = Z0FaceGroup
    else:
        face_domain = Partition_1
        # this is important to make shear mapping bc pair
        Regular_1D_1 = Mesh_1.Segment(geom=face_domain)
        status = Mesh_1.AddHypothesis(Number_of_Segments_default, face_domain)
        MEFISTO_2D_1 = Mesh_1.Triangle(algo=smeshBuilder.MEFISTO,geom=face_domain)
        Sub_mesh_1 = Regular_1D_1.GetSubMesh()

    Regular_1D_2 = Mesh_1.Segment(geom=ShortEdgeGroup)
    Number_of_Segments_single = Regular_1D_2.NumberOfSegments(line_segment_nb_small)
    Sub_mesh_2 = Regular_1D_2.GetSubMesh()

    Regular_1D_3 = Mesh_1.Segment(geom=VerticalEdgeGroup)
    Number_of_Segments_vertical = Regular_1D_3.NumberOfSegments(line_segment_nb_vertical)
    Sub_mesh_3 = Regular_1D_3.GetSubMesh()
    smesh.SetName(Number_of_Segments_vertical, 'Number of Segments_vertical')

    if using_workpiece_extra_extusion:
        Regular_1D_3_1 = Mesh_1.Segment(geom=edges_wve)
        Number_of_Segments_vertical = Regular_1D_3_1.NumberOfSegments(line_segment_nb_vertical)
        Sub_mesh_3_1 = Regular_1D_3_1.GetSubMesh()
        smesh.SetName(Sub_mesh_3_1, 'Sub-mesh_3_1')

    #MEFISTO_2D_2 = Mesh_1.Triangle(algo=smeshBuilder.MEFISTO,geom=FaceGroupTriangleMesh)
    #Sub_mesh_4 = MEFISTO_2D_2.GetSubMesh()

    ## Set names of Mesh objects
    smesh.SetName(Number_of_Segments_single, 'Number_of_Segments_small')
    smesh.SetName(Number_of_Segments_default, 'Number_of_Segments_default')

    smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
    smesh.SetName(Sub_mesh_1, 'Sub-mesh_1')
    smesh.SetName(Sub_mesh_3, 'Sub-mesh_3')
    smesh.SetName(Sub_mesh_2, 'Sub-mesh_2')

    if using_quad_mesh:
        Quadrangle_2D = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE,geom=FaceGroupQuadMesh)
        Sub_mesh_4 = Quadrangle_2D.GetSubMesh()
        smesh.SetName(Sub_mesh_4, 'Sub-mesh_4')

    mesh_solid_type= SMESH.VOLUME
    mesh_boundary_type = SMESH.FACE
else:
    Mesh_1 = smesh.Mesh(domain)
    Regular_1D = Mesh_1.Segment()
    Number_of_Segments_default = Regular_1D.NumberOfSegments(line_segment_nb_default)
    MEFISTO_2D = Mesh_1.Triangle(algo=smeshBuilder.MEFISTO)

    Regular_1D_2 = Mesh_1.Segment(geom=ShortEdgeGroup)
    Number_of_Segments_single = Regular_1D_2.NumberOfSegments(line_segment_nb_small)
    Sub_mesh_2 = Regular_1D_2.GetSubMesh()

    smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
    smesh.SetName(Sub_mesh_2, 'Sub-mesh_2')

    # this working for meshing, but solver does not converge
    if using_quad_mesh:
        Quadrangle_2D = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE,geom=FaceGroupQuadMesh)
        Sub_mesh_4 = Quadrangle_2D.GetSubMesh()
        smesh.SetName(Sub_mesh_4, 'Sub-mesh_4')

    mesh_solid_type= SMESH.FACE
    mesh_boundary_type = SMESH.EDGE


#isDone = Mesh_1.SetMeshOrder( [ [ Sub_mesh_2, Sub_mesh_4, Sub_mesh_1, Sub_mesh_3 ] ])  ## cause error?
isDone = Mesh_1.Compute()



## group
Solid_work_1 = Mesh_1.GroupOnGeom(Solid_work,'Solid_work',mesh_solid_type)
Solid_cutter_1 = Mesh_1.GroupOnGeom(Solid_cutter,'Solid_cutter',mesh_solid_type)
Solid_holder_1 = Mesh_1.GroupOnGeom(Solid_holder,'Solid_holder',mesh_solid_type)
Solid_chip_1 = Mesh_1.GroupOnGeom(Solid_chip,'Solid_chip',mesh_solid_type)
Solid_shear_1 = Mesh_1.GroupOnGeom(Solid_shear,'Solid_shear',mesh_solid_type)
Solid_friction_1 = Mesh_1.GroupOnGeom(Solid_friction,'Solid_friction',mesh_solid_type)
Solid_cutter_chip_transition_1 = Mesh_1.GroupOnGeom(Solid_cutter_chip_transition,'Solid_cutter_chip_transition',mesh_solid_type)

## here _1 means group not topo objects
smesh.SetName(Solid_work_1, 'Solid_work')
smesh.SetName(Solid_cutter_1, 'Solid_cutter')
smesh.SetName(Solid_holder_1, 'Solid_holder')
smesh.SetName(Solid_chip_1, 'Solid_chip')
smesh.SetName(Solid_shear_1, 'Solid_shear')
smesh.SetName(Solid_friction_1, 'Solid_friction')
smesh.SetName(Solid_cutter_chip_transition_1, 'Solid_cutter_chip_transition')

##
if using_3D:
    if exporting_foam:
        
        inlet_1 = Mesh_1.GroupOnGeom(inlet,'inlet',SMESH.FACE)
        outlet_1 = Mesh_1.GroupOnGeom(outlet,'outlet',SMESH.FACE)
        chipEnd_1 = Mesh_1.GroupOnGeom(chipend,'chipEnd',SMESH.FACE)
        workBottom_1 = Mesh_1.GroupOnGeom(work_bottom,'workBottom',SMESH.FACE)
        holderTop_1 = Mesh_1.GroupOnGeom(holder_top,'holderTop',SMESH.FACE)
        #chipEnd_1 = Mesh_1.GroupOnGeom(HTC_moving,'HTC_moving',SMESH.FACE)

        smesh.SetName(inlet_1, 'inlet')
        smesh.SetName(outlet_1, 'outlet')
        smesh.SetName(chipEnd_1, 'chipEnd')
        smesh.SetName(workBottom_1, 'workBottom')
        smesh.SetName(holderTop_1, 'holderTop')
        #smesh.SetName(HTC_moving_1, 'HTC_moving')

        if using_mapping_bc:
            mapping_work_1 = Mesh_1.GroupOnGeom(mapping_work,'mappingWork',SMESH.FACE)
            smesh.SetName(mapping_work_1, 'mapping_work')
            mapping_chip_1 = Mesh_1.GroupOnGeom(mapping_chip,'mappingChip',SMESH.FACE)
            smesh.SetName(mapping_chip_1, 'mapping_chip')
    else:  #split, dolfin only support tetra
        Mesh_1.SplitVolumesIntoTetra( Mesh_1, 1)
else:
    Mesh_1.QuadTo4Tri(Mesh_1)  # how to split 2D?
    pass

isDone = Mesh_1.Compute()

#[ Solid_chip_1, Solid_shear_1, Solid_friction_1, Solid_work_1, Solid_cutter_1, inlet_1, chipEnd_1, Solid_holder_1] = Mesh_1.GetGroups()


#Mesh_1.ExportMED( '/media/sf_OneDrive/Fenics/metal_cut/salome_mesh/Mesh_dolfin.med', 0, SMESH.MED_V2_2, 1, None ,1)
Mesh_1.ExportMED(salome_mesh_output, 0, SMESH.MED_V2_2, 1, None ,1)
#can not replace the media file, jsut can not replace it in shared_folder
'''
try:
    if exporting_foam:
        #file should not exist!
        Mesh_1.ExportMED( '/tmp/Mesh_foam.med', 0, SMESH.MED_V2_2, 1, None ,1)
    else:
        Mesh_1.ExportMED( '/media/sf_OneDrive/Fenics/metal_cut/salome_mesh/Mesh_dolfin.med', 0, SMESH.MED_V2_2, 1, None ,1)
    #pass
except:
    print 'ExportToMED() failed. Invalid file name?'
'''

print('==== completed mesh saving ====')

#meshio-convert Mesh_1.med test.xdmf  # support quad mesh
if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(True)
else:
    from runSalome import *
    myArgs={}
    myArgs["killall"] = True
    kill_salome(myArgs)
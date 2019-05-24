## Introduction

This is code for the paper to be published:
*Quasi-static 3D and 2D FEM nonlinear thermal model for high speed metal cutting*

Journal paper: preprint link to be add later?

## License (Academic Creative)

It is free to use my code for non-commercial purposes, but please refer to my *paper* to ack this work.
For commercial usage, contact me to relicense it, I can do contract coding/consultancy for industries.

## key features

see more details in my paper <link to be added later>

1. nonlinear material properties and realistic geometry
2. quasi-static reduced order model, it is very fast
3. special way to deal with unphysical heat flux at the material-air boundary
4. heat partition for work-chip, tool-chip can be calcuated directly from this multi-body FEM model
5. this is an example of automated engineering simulation, see my ppt on Fenics18 conference: <>

```bash
#/opt/SALOME-8.5.0-UB16.04-SRC/runSalome -t -b /myScript
## the sequence of gmsh options is important
#gmsh4 -format msh2   -o metal_cut.msh -save Mesh_1.med
#dolfin-convert metal_cut.msh metal_cut.xml
```

## Limitation in this nonlinear model

There is no known issue for constant thermal properties (linear model), but some issues are spotted for nonlinear material.
For the thermal modelling of metal cutting, large error comes from key geometries (the primary and secondary zones dimension prediction) and material properties model. I list the issues below, hopefully you can shed some light on thise minor issues:

1. 1% error in heat convervation (heat generation != heat loss form all surfaces) for nonlinear material model
    It does not mean 1% error in maximum temperature, it is just the upper limit. This error does not drop with speed, mesh fineness?

2. Cutter material (can be higher than workpiece metal) is treated as metal material. There is a trend of reduced interface temperature and increased heat into cutter, 0.5% error is evaluated in linear model
    Thermal contact resistance is ignored in this model (for very high normal stress), which will slightly increase interface temperature
    user-defined nonlinear loop could improve this situation.

3. It is possible to make the solver run in parallel, but I have not made the effort for this.
   However, it is not necessary as I can complete the computation on a single core within several minutes for a quasi-static 3D model and several seconds for a quasi-static 2D model.

4. I try to make the code more readable, even python code is very readable compared with Fortran.
    but it is not a software project, I have no time to do documentation.

5. I tried to implement the algorith of `Oxley_JC_model.py` in one paper, but I have some issues in convergence for my specific implementation, for this multivariale optimisation problem.
    That is useful for mechcanical-thermal coupling, as my future work.

## Test

This code (no GUI code is involved) is develoepd on ubuntu 16.04, but it should be anable to run on any POSIX platform.

As this is not a commecial product/complete software, this solver can not run on your PC out of box.

First of all, install all dependencies listed below, then fix any problem related with filename and path name, good luck!

Secondly, run the solver
    edit `parameter.py` to switch on and off features.
    edit `material.py` to switch on and off nonlinear material models and replace material.
    `python metal_cut_ht_6.py`

Finally, using paraview to view the result file `T.pvd`

### Dependencies

+ python2 or pyhton3 with Fenics installed
+ Fenics  install from PPA, it will install all other essential python packages like numpy and matplotlib
+ FenicsSolver   My Multiphysics solver based on Fenics
             <> I have not uploaded to pypi, it is still in heavy develoepment stage.
+ salome platform
+ gmsh4

### Extra setup

+ gmsh4 and salome path is hardcoded in my script,may need to replace it
+ `FenicsSolver` should also be included in your  python`sys.path` need to fix it.
+ create folder_structure:  `mkdir output && mkdir salome_mesh && mkdir results`

## Ackowledgement

This is an independent reseach during my employment in the University of Oxford, while it is not a derivative work from my current reseach project, partially supported by industrial partners.

Thanks are given to my supervisor, Prof Gillespie, who is also one the authors of the paper, for his financial support and career guidance at the University of Oxford.

I appreciate support from all my family members, esp my wife, Ms J Wang, to free me from house work to complete this work.
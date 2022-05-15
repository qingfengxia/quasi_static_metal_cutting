# Reduced-order thermal modelling for multi-body system (metal cutting)

## Introduction

This is code for the journal paper published:
Qingfeng Xia, David R.H. Gillespie, Quasi-static  finite element modelling of thermal distribution and heat partitioning  for the multi-component system of high speed metal cutting, 
*Journal of Materials Processing Technology*, 2019,  https://doi.org/10.1016/j.jmatprotec.2019.116389

## License (Academic Creative)

It is free to use my code for non-commercial purposes, but please cite my *paper* to ack this work.
For commercial usage, contact me to relicense it, I can do contract coding/consultancy for industries.

## key features

see more details in my paper <https://doi.org/10.1016/j.jmatprotec.2019.116389>
if you can not access from the publisher, here is the pre-print
<https://www.researchgate.net/publication/332018765_Quasi-static_finite_element_modelling_of_thermal_distribution_and_heat_partitioning_for_the_multi-component_system_of_high_speed_metal_cutting_accepted>

1. nonlinear material properties and realistic geometry
2. quasi-static reduced order model, it is very fast, quasi-static 2D simulation can be completed in 3 seconds.
3. special way to deal with unphysical heat flux at the material-air boundary, see details in my paper
4. heat partition for work-chip, tool-chip can be calculated directly from this multi-body FEM model
5. this is an example of automated engineering simulation, see my ppt on Fenics18 conference and my `CAE_pipeline` repo
<https://github.com/qingfengxia/CAE_pipeline>


## Install software required

### Python3 + numpy scipy matplotlib + vtk

The python3 program (I used in our published paper quasi-static thermal modelling) is based on Fenics, which is Linux only. I run it on ubuntu 18.04 in virtualbox, while you can have ubuntu 20.04 in WSL2 + Xwindows on windows 10.  If you have ubuntu 18.04/20.04 somewhere, you can run it without change anything hopefully.

### mehser: gmsh Salome

Salome could be needed,  it is a total open source workflow but with automated parameter study.

https://www.salome-platform.org/?page_id=15

gmsh4 can be download from gmsh website if on ubutnu 18.04. Here is the link: https://gmsh.info/bin/Linux/gmsh-4.10.2-Linux64.tgz

On Ubuntu 20.04, [gmsh4 is the default version](https://ubuntu.pkgs.org/20.04/ubuntu-universe-arm64/gmsh_4.4.1+ds1-2build1_arm64.deb.html) `sudo apt-get install gmsh`

### Fenics and FenicsSolver

Recently Fenics upgraded to Fenics-X, the next generation, which will break the program. If stay with Ubuntu 20.04,  it is still the original fenics. compatible with FenicsSolver (a high-level wrapper of fenics library) I wrote and published on my github.

`sudo apt-get install fenics`
do not forget `dolfin-bin`

This link above should be useful to run the code, except install fenics 2019.1 from official repo, not from PPA.

Install `FenicsSolver`: a high-level wrapper of fenics library

### Extra setup (for some hardcoded path in scripts)

+ gmsh4 and salome path is hardcoded in my script, so you need to put them on PATH
+ `FenicsSolver` should also be included in your python`sys.path` need to fix it.
+ create folder_structure:  `mkdir output && mkdir salome_mesh && mkdir results`

## Workflow pipeline

### general workflow

CAE_pipeline repo https://github.com/qingfengxia/CAE_pipeline

### use case for metal cut in this repo

To run the code, some basic knowledge of python3 is needed, as your new env would be different such as file path.

0. install all the software required, as stated in previous section

1. set parameters,  usually, there is a `parameter.py` 
    edit `parameter.py` to switch on and off features.
    edit `material.py` to switch on and off nonlinear material models and replace material.

2. generate mesh using  gmsh or Salome (generate geometry and mesh)

   please edit the var `gmsh_app` and `salome_app` path in`mesh_utitilies.py`

   Here, we run salome in batch mode instead of GUI mode, myScript will be played to generate mesh.`/opt/SALOME-8.5.0-UB16.04-SRC/runSalome -t -b myScriptToGenerateMesh.py`  

   On Ubuntu 18.04, the download tar name will be  `SALOME-x.y.0-UB18.04-SRC`,  salome 8.x is recommended, since I have not tested 9.x.   The python `myScriptToGenerateMesh.py` in metal cut case, is `salome_mesh_metal_cut.py`

     `gmsh` is also needed to convert salome mesh into the med format, then `dolfin-convert` convert the mesh into vtk format for Fenics

 ```sh
## the sequence of gmsh options is important,  gmsh 4.x used, 
gmsh -format msh2   -o metal_cut.msh -save /tmp/Mesh_1.med
dolfin-convert metal_cut.msh metal_cut.xml
# actually 2 files are generated:   metal_cut_physical_region.xml  metal_cut.xml
 ```

3. solve by  FenicsSolver (to solve and save vtk file)
   there could be some batch shell script to drive parameter study

   `python3 metal_cut_ht_6.py`

4. post processing python matplotlib to plot result (postprocessing) and also paraview 
   using paraview to view the result file `T.pvd`

## Test platforms

This code (no GUI code is involved) is developed on ubuntu 16.04/18.04, but it should run on any POSIX platform.

As this is not a commercial product/complete software, this solver may not run on your PC out of box.

+ pyhton3.6 with Fenics installed, python2 would not be support any more
+ Fenics  v2017.2-v2019.1 tested, it will install all other essential python packages like numpy and matplotlib
+ FenicsSolver   My Multiphysics solver based on Fenics 
+ salome platform 8.5
+ gmsh4


### Test without salome installed (mesh can be generated but program not working)
using dolfin 2017.x is the choice.

mesh has been generated by salome 9.7 and uploaded to github repo

`~/SALOME-9.7.0/salome -t -b salome_mesh_metal_cut.py `

Without installing salome, it is possible to test the program (installation) by skip the preprocessing stage (mesh generation by salome), while the solving may fail if the mesh is not for the parameter.py

```sh
# assuming using the default setup in parameter.py, is_preproceesing = False
python3 metal_cut_ht_6.py
```
https://fenicsproject.org/olddocs/dolfin/2019.1.0/python/

### Parameter study by generate parameter.py from parameter_template.py

here is the script for speed parameter study

```sh
ORIGIN="parameter_template.py"
SCRIPT="parameter.py"
VAR="cutting_speed"
#first of all, comment out the var assign if after anchor line (##ZZ##) in parameter_template
#also check result file name

#double quote to expand variable into string
for VALUE in 0.1 0.2 0.4 0.6 0.8 1.0 1.5 2.0 2.5 3.0 3.5 4.0; do
#for VALUE in 2.5 3.0 3.5 4.0; do
    newline="$VAR=$VALUE"
    echo "processing parameter:$newline"
    sed "s/##ZZ##/$newline/" $ORIGIN > $SCRIPT
    #is_batch_mode = False
    python3 metal_cut_ht_6.py
done
```

## Limitation in this nonlinear model

There is no known issue for constant thermal properties (linear model), but some issues are spotted for nonlinear material.
For the thermal modelling of metal cutting, large error comes from key geometries (the primary and secondary zones dimension prediction) and material properties model. I list the issues below, hopefully you can shed some light on thise minor issues:

1. 1% error in heat conservation (heat generation != heat loss form all surfaces) for nonlinear material model
   It does not mean 1% error in maximum temperature, it is just the upper limit. This error does not drop with speed, mesh fineness?

2. Cutter material (can be higher than workpiece metal) is treated as metal material. There is a trend of reduced interface temperature and increased heat into cutter, 0.5% error is evaluated in linear model
   Thermal contact resistance is ignored in this model (for very high normal stress), which will slightly increase interface temperature
   user-defined nonlinear loop could improve this situation.

3. It is possible to make the solver run in parallel, but I have not made the effort for this.
   However, it is not necessary as I can complete the computation on a single core within several minutes for a quasi-static 3D model and several seconds for a quasi-static 2D model.

4. I try to make the code more readable, even python code is very readable compared with Fortran.
   but it is not a software project, I have no time to do documentation.

5. I tried to implement the algorithm of `Oxley_JC_model.py` in one paper, but I have some issues in convergence for my specific implementation, for this multivariale optimisation problem.
   That is useful for mechcanical-thermal coupling, as my future work.

## Acknowledgement

This is an independent reseach during my employment in the University of Oxford, while it is not a derivative work from my current reseach project, partially supported by industrial partners.

Thanks are given to my supervisor, Prof Gillespie, who is also one the authors of the paper, for his financial support and career guidance at the University of Oxford.

I appreciate support from all my family members, esp my wife, Ms J Wang, to free me from house work to complete this work.
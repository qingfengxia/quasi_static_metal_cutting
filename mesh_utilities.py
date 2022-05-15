#utf8
from __future__ import print_function, division
import subprocess
import os
import os.path
import sys
import time

default_tmp_mesh_filename = '/tmp/Mesh_1.med'
# salome 9.x may be too new for this script, some version released in 2018 is better
# salome_app = '/media/sf_OneDrive/Salome-9.8.0/salome'
salome_app = '/home/qingfeng/SALOME-9.7.0/salome'
gmsh_app = '/media/sf_OneDrive/gmsh-4.10.2-Linux64/bin/gmsh'

def run_command(comandlist):
    has_error = False
    try:
        print(comandlist)
        p = subprocess.Popen(comandlist, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, error = p.communicate()
        print(output)  # stdout is still cut at some point but the warnings are in stderr and thus printed :-)
        print(error)  # stdout is still cut at some point but the warnings are in stderr and thus printed :-)
    except:
        print('Error executing: {}\n'.format(comandlist))
        has_error = True
    return has_error

def defined(x):
    return x in locals() or x in globals()


def generate_salome_mesh(mesh_gen_script, smesh_file = default_tmp_mesh_filename):
    # this function is quite specific to my pc setup
    if not os.path.exists(salome_app):
        print('Error: salome executable is not found in the specified path\n', salome_app)
        sys.exit(-1)
    salome_cmd = "{} -t -b {}".format(salome_app, mesh_gen_script)
    os.system(salome_cmd)  # run_command will bring up GUI, but os.system does not!
    # there is salome shutdown command inside the script
    #"%PYTHONBIN%" "%KERNEL_ROOT_DIR%\bin\salome\runSalome.py" -t -u myScript.py
    #%PYTHONBIN% "%KERNEL_ROOT_DIR%\bin\salome\killSalome.py"
    # salome_cmd = ["/opt/SALOME-8.5.0-UB16.04-SRC/salome", "-t", script]
    # run_command(salome_cmd)
    check_mtime(smesh_file)

def convert_salome_mesh_to_dolfin(output_dolfin_mesh, smesh_file = default_tmp_mesh_filename):
    gmsh_filename = output_dolfin_mesh[:-4] + ".msh"
    cmdline = "{} -format msh2 -o {} -save {}".format(gmsh_app, gmsh_filename, smesh_file)
    run_command(cmdline)  # check the return value 'has_error' does not always work so check the output timestamp
    check_mtime(gmsh_filename)
    run_command("dolfin-convert {} {}".format(gmsh_filename, output_dolfin_mesh))
    check_mtime(output_dolfin_mesh)

def convert_salome_mesh_to_foam(output_foam_case_folder, smesh_file = default_tmp_mesh_filename):
    gmsh_filename = smesh_file[:-4] + ".msh"
    cmdline = "{} -format msh2 -o {} -save {}".format(gmsh_app, gmsh_filename, smesh_file)
    run_command(cmdline)
    check_mtime(gmsh_filename)
    run_command("gmshToFoam -case {} {}".format(output_foam_case_folder, gmsh_filename))
    check_mtime(output_foam_case_folder + '/constant/polyMesh/boundary')

def from_python_file_to_parameter_string(python_file):
    #from StringIO import StringIO
    #from io import stringIO
    parameter_string = []
    with open(python_file, "r") as inf:
        lines = inf.readlines()
        for i,l in enumerate(lines):
            parameter_string.append(l.replace('#' , '//'))
    return ''.join(parameter_string)

def generate_gmsh_mesh(mesh_file_root, mesh_parameter_string):
    # if mesh_parameter_string is not provide, just use provided geo file
    if mesh_parameter_string:
        template_file = mesh_file_root + "_template.geo"
        with open(template_file) as inf:
            with open(mesh_file_root + ".geo", "w") as outf:
                outf.write(mesh_parameter_string)
                outf.write(inf.read())
    # use gmsh4 instead gmsh3 here
    gmshcmd = ['{} - -match -tol 1e-12 - {}.geo'.format(gmsh_app, mesh_file_root)]
    if(run_command(gmshcmd)):
        sys.exit()
    check_mtime(mesh_file_root + ".msh")

def convert_salome_mesh_to_dolfin(mesh_file_root):
    convertcmd= ['dolfin-convert {}.msh {}.xml'.format(mesh_file_root, mesh_file_root)]
    if (run_command(convertcmd) ):
        sys.exit()
    check_mtime(mesh_file_root + ".xml")

def check_mtime(filename):
    modified_time = os.path.getmtime(filename)
    current_time = time.time() # in second since epoch
    second_delta = current_time - modified_time
    if second_delta > 200:
        print('file `{}` modified time is more than {} seconds, mesh conversion failed?'.format(filename, second_delta))
        print("file last modified at %s" % time.ctime(modified_time))
        sys.exit()

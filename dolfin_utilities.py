from __future__ import print_function, division
import os.path

try:
    #from mshr import *
    from dolfin import *
except:
    print('Fenics is not installed, exit')
    exit()


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
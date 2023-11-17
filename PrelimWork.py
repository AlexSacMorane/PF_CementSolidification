#-------------------------------------------------------------------------------
# Librairies
#-------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import vtk
from vtk.util.numpy_support import vtk_to_numpy

#-------------------------------------------------------------------------------
# User
#-------------------------------------------------------------------------------

namefile = 'vtk/PF_Cement_Solidification_other_000_0.vtu'

#-------------------------------------------------------------------------------
# Work
#-------------------------------------------------------------------------------

# load a vtk file as input
reader = vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName(namefile)
reader.Update()

#Grab a scalar from the vtk file
my_vtk_array = reader.GetOutput().GetPointData().GetArray("phi")

#Get the coordinates of the nodes and the scalar values
nodes_nummpy_array = vtk_to_numpy(nodes_vtk_array)
my_numpy_array = vtk_to_numpy(my_vtk_array )

x,y,z = nodes_nummpy_array[:,0],\
        nodes_nummpy_array[:,1],\
        nodes_nummpy_array[:,2]

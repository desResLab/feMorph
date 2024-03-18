import numpy as np
import matplotlib.pyplot as plt
from vtk import vtkXMLUnstructuredGridReader
from vtk import vtkXMLPolyDataReader
from vtk.util.numpy_support import vtk_to_numpy

def exctract_time(name):
  name = name.replace('velocity_','')
  try:
    res = float(name)
  except:
    print('Invalid floating point number for velocity time')
    exit(-1)
  return res

# =================
# READ FROM VTU/VTP
# =================

# Read nodes and elements from VTU
def read_nd_el(vtu_file):
  reader = vtkXMLUnstructuredGridReader()
  reader.SetFileName(vtu_file)
  reader.Update()
  nd = vtk_to_numpy(reader.GetOutput().GetPoints().GetData())
  # WARNING: Assumes the elements are all TET4
  # el = vtk_to_numpy(reader.GetOutput().GetCells().GetData()).reshape((-1,5))[:,1:,]
  # WARNING: Assumes the elements are all HEXA8
  el = vtk_to_numpy(reader.GetOutput().GetCells().GetData()).reshape((-1,9))[:,1:,]
  return nd,el

# Read bct file
def read_bct_vtp(bct_vtp_file):
  reader = vtkXMLPolyDataReader()
  reader.SetFileName(bct_vtp_file)
  reader.Update()
  bct = reader.GetOutput()
  # Get global node ID
  inlet_nodes = vtk_to_numpy(bct.GetPointData().GetScalars('GlobalNodeID'))
  # Get Nodal Coordinates
  inlet_coords = vtk_to_numpy(bct.GetPoints().GetData())

  vel_time = [0.0,10000.0]
  vels     = []
  inlet_vels = np.zeros_like(inlet_coords)
  for loopA in range(len(inlet_vels)):
    y = inlet_coords[loopA,1]
    inlet_vels[loopA,0] = 1.0
    inlet_vels[loopA,1] = 0.0
    inlet_vels[loopA,2] = 0.0
  vels.append(inlet_vels)
  vels.append(inlet_vels)
  
  return inlet_nodes,vel_time,vels

def read_wall(ext_file):
  reader = vtkXMLPolyDataReader()
  reader.SetFileName(ext_file)
  reader.Update()
  # Get global node ID
  nd_id = vtk_to_numpy(reader.GetOutput().GetPointData().GetScalars('GlobalNodeID'))
  return nd_id

def read_outlet(outlet_file):
  reader = vtkXMLPolyDataReader()
  reader.SetFileName(outlet_file)
  reader.Update()
  # Get global node ID
  nd_id = vtk_to_numpy(reader.GetOutput().GetPointData().GetScalars('GlobalNodeID'))
  return nd_id

# ===============
# WRITE FUNCTIONS
# ===============

def write_problem(out_file_name):
  f = open(out_file_name,"w")
  # Problem Type
  f.write('PROBLEM EXPLICITVMS\n')
  # Nodal DOFs
  f.write('NODEDOF 3\n')  
  # Timestep
  f.write('TIMESTEP 5.0e-3 5000 10\n')
  f.close()

def write_nd(nodes,out_file_name):
  f = open(out_file_name,"a")
  for loopA in range(len(nodes)):
    # NODE 1 0.000000e+00 0.000000e+00 0.000000e+00
    f.write('NODE %d %.8e %.8e %.8e\n' % (loopA+1,nodes[loopA,0],nodes[loopA,1],nodes[loopA,2]))
  f.close()

def write_el(elements,out_file_name):
 f = open(out_file_name,"a")
 for loopA in range(len(elements)):
   # ELEMENT QUAD4 1 1 1 2 103 102
   # ELEMENT TET4 %d 1 %d %d %d %d
   f.write('ELEMENT HEXA8 %d 1 %d %d %d %d %d %d %d %d\n' % (loopA+1,elements[loopA,0]+1,
                                                                     elements[loopA,1]+1,
                                                                     elements[loopA,2]+1,
                                                                     elements[loopA,3]+1,
                                                                     elements[loopA,4]+1,
                                                                     elements[loopA,5]+1,
                                                                     elements[loopA,6]+1,
                                                                     elements[loopA,7]+1))
 f.close()

def write_dir(wall_nodes,out_file_name):
 f = open(out_file_name,"a")
 for loopA in range(len(wall_nodes)):
   # NODEDIRBC 7979 0.0 0.0 0.0 0
   f.write('NODEDIRBC %d 0.0 0.0 0.0\n' % (wall_nodes[loopA]))
 f.close()

# Prescribed velocity
def write_inlet(inlet_nodes,vel_time,vels,out_file_name):
  f = open(out_file_name,"a")
  # Loop over velocity nodes
  for loopA in range(len(inlet_nodes)):
    curr_node = inlet_nodes[loopA]
    # Loop over time steps
    for loopB in range(len(vel_time)):
      curr_time = vel_time[loopB]
      curr_vel = vels[loopB][loopA]
      # NODEDIRBC 7979 0.0 0.0 0.0 1 
      f.write('NODEVEL %d %.8e %.8e %.8e %.8e\n' % (curr_node,curr_vel[0],curr_vel[1],curr_vel[2],curr_time))
  f.close()

def write_outlet(outlet_nodes,out_file_name):
 f = open(out_file_name,"a")
 for loopA in range(len(outlet_nodes)):
   f.write('NODEPRES %d 0.0\n' % (outlet_nodes[loopA]))
 f.close()

# # Write properties to file
def write_prop(out_file_name):
 f = open(out_file_name,"a")
 # double rho = model->vmsProps[0];
 # double mu  = model->vmsProps[1];
 # double lam = model->vmsProps[2];
 # double cI  = model->vmsProps[3];
 # Reynolds number 100: lid velocity U = 1 m·s-1, L = 1 m, ρ = 10 kg·m-3, µ = 0.1 Pa·s
 # Reynolds number 1,000: lid velocity U = 1 m·s-1, L = 1 m, ρ = 10 kg·m-3, µ = 0.01 Pa·s
 # Reynolds number 10,000: lid velocity U = 1 m·s-1, L = 1 m, ρ = 10 kg·m-3, µ = 0.001 Pa·s
 # OR
 # Reynolds number 100: lid velocity U = 10 m·s-1, L = 1 m, ρ = 1 kg·m-3, µ = 0.1 Pa·s
 # Reynolds number 1,000: lid velocity U = 10 m·s-1, L = 1 m, ρ = 1 kg·m-3, µ = 0.01 Pa·s
 # Reynolds number 10,000: lid velocity U = 10 m·s-1, L = 1 m, ρ = 1 kg·m-3, µ = 0.001 Pa·s
 f.write('VMSPROPS 10.0 0.01 0.05 4.0 0.0\n')
 f.close()
 
# =========
# MAIN CODE
# =========
if __name__ == "__main__":

  # Set file name
  vtu_file         = 'fluid.vtu'
  inlet_bct_file   = 'wall_top.vtp'
  wall_file_left   = 'wall_left.vtp'
  wall_file_right  = 'wall_right.vtp'
  wall_file_bottom = 'wall_bottom.vtp'

  # No outlet boundary conditions
  out_file_name = 'model.txt'

  # Read Qtys
  nodes,elements            = read_nd_el(vtu_file)
  inlet_nodes,vel_time,vels = read_bct_vtp(inlet_bct_file)
  wall_nodes_left           = read_wall(wall_file_left)
  wall_nodes_right          = read_wall(wall_file_right)
  wall_nodes_bottom         = read_wall(wall_file_bottom)

  # Write to model file
  write_problem(out_file_name)
  write_prop(out_file_name)
  write_nd(nodes,out_file_name)
  write_el(elements,out_file_name)
  write_dir(wall_nodes_left,out_file_name)
  write_dir(wall_nodes_right,out_file_name)
  write_dir(wall_nodes_bottom,out_file_name)

  # Write inlet velocities
  write_inlet(inlet_nodes,vel_time,vels,out_file_name)



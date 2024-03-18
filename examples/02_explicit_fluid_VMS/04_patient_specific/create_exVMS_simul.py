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
  el = vtk_to_numpy(reader.GetOutput().GetCells().GetData()).reshape((-1,5))[:,1:,]
  return nd,el

# Read bct file
def read_bct_vtp(bct_vtp_file):
  reader = vtkXMLPolyDataReader()
  reader.SetFileName(bct_vtp_file)
  reader.Update()
  bct = reader.GetOutput()
  # Get global node ID
  inlet_nodes = vtk_to_numpy(bct.GetPointData().GetScalars('GlobalNodeID'))
  # Find velocities in PointData
  vel_time = []
  vel_idx  = []
  vels     = []
  for loopA in range(bct.GetPointData().GetNumberOfArrays()):
    res_name = bct.GetPointData().GetArray(loopA).GetName()
    if('velocity' in res_name):
      curr_time = exctract_time(res_name)
      vel_time.append(curr_time)
      vel_idx.append(loopA)
      vels.append(vtk_to_numpy(bct.GetPointData().GetArray(loopA)))
      # Set initial velocities to zero for ramp
      if(curr_time == 0.0):
        vels[-1][:] = 0.0
  return inlet_nodes,vel_idx,vel_time,vels

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
  f.write('NODEDOF 4\n')
  # Timestep
  f.write('TIMESTEP 0.0001 10 1\n')
  f.write('DUALSTEP 1.0e-3 0\n')
  f.close()

def write_nd(nodes,out_file_name):
  f = open(out_file_name,"a")
  for loopA in range(len(nodes)):
    # NODE 1 0.000000e+00 0.000000e+00 0.000000e+00
    f.write('NODE %d %8.3e %8.3e %8.3e\n' % (loopA+1,nodes[loopA,0],nodes[loopA,1],nodes[loopA,2]))
  f.close()

def write_el(elements,out_file_name):
 f = open(out_file_name,"a")
 for loopA in range(len(elements)):
   # ELEMENT QUAD4 1 1 1 2 103 102
   f.write('ELEMENT TET4 %d 1 %d %d %d %d\n' % (loopA+1,elements[loopA,0]+1,elements[loopA,1]+1,elements[loopA,2]+1,elements[loopA,3]+1))
 f.close()

def write_dir(wall_nodes,out_file_name):
 f = open(out_file_name,"a")
 for loopA in range(len(wall_nodes)):
   # NODEDIRBC 7979 0.0 0.0 0.0 0
   f.write('NODEDIRBC %d 0.0 0.0 0.0\n' % (wall_nodes[loopA]))
 f.close()

# Prescribed velocity
def write_inlet(inlet_nodes,vel_idx,vel_time,vels,out_file_name):
  f = open(out_file_name,"a")
  # Loop over velocity nodes
  for loopA in range(len(inlet_nodes)):
    curr_node = inlet_nodes[loopA]
    # Loop over time steps
    for loopB in range(len(vel_time)):
      curr_time = vel_time[loopB]
      curr_vel = vels[loopB][loopA]
      # NODEDIRBC 7979 0.0 0.0 0.0 1 
      f.write('NODEVEL %d %8.3e %8.3e %8.3e %8.3e\n' % (curr_node,curr_vel[0],curr_vel[1],curr_vel[2],curr_time))
  f.close()

def write_outlet(outlet_nodes,out_file_name):
 f = open(out_file_name,"a")
 for loopA in range(len(outlet_nodes)):
   f.write('NODEPRES %d 0.0\n' % (outlet_nodes[loopA]))
 f.close()

# # Write properties to file
def write_prop(out_file_name):
 f = open(out_file_name,"a")
 # density, viscosity, epsilon^-1, c1, c2
 f.write('VMSPROPS 1.06 0.04 400.0 4.0 2.0\n')
 f.close()
 
# =========
# MAIN CODE
# =========
if __name__ == "__main__":

  # Set file name
  vtu_file      = 'mesh/aorta.vtu'
  bct_vtp_file  = 'mesh/bct.vtp'
  wall_file     = 'mesh/walls_combined.vtp'
  outlet_file   = 'mesh/outlets_combined.vtp'
  out_file_name = 'model.txt'

  # Read Qtys
  print('Reading nodes and cells...')
  nodes,elements                    = read_nd_el(vtu_file)
  print('Reading bct...')
  inlet_nodes,vel_idx,vel_time,vels = read_bct_vtp(bct_vtp_file)
  print('Reading walls...')
  wall_nodes                        = read_wall(wall_file)
  print('Reading outlets...')  
  outlet_nodes                      = read_outlet(outlet_file)

  # # Write to model file
  print('Writing problem...')
  write_problem(out_file_name)
  print('Writing nodes...')
  write_nd(nodes,out_file_name)
  print('Writing elements...')
  write_el(elements,out_file_name)
  print('Writing walls...')
  write_dir(wall_nodes,out_file_name)

  # Write inlet velocities
  print('Writing inlets...')
  write_inlet(inlet_nodes,vel_idx,vel_time,vels,out_file_name)

  # Write outlet pressure conditions
  print('Writing outlets...')
  write_outlet(outlet_nodes,out_file_name)

  # Write propo
  print('Writing properties...')
  write_prop(out_file_name)

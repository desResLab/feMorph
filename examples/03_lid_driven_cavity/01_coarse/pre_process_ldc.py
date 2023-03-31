import numpy as np
import matplotlib.pyplot as plt
from vtk import vtkXMLUnstructuredGridReader
from vtk import vtkXMLPolyDataReader
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import vtk

# Read nodes and elements from VTU
def read_nd_el(vtu_file):
  reader = vtkXMLUnstructuredGridReader()
  reader.SetFileName(vtu_file)
  reader.Update()
  nd = vtk_to_numpy(reader.GetOutput().GetPoints().GetData())
  # WARNING: Assumes the elements are all TET4
  el = vtk_to_numpy(reader.GetOutput().GetCells().GetData()).reshape((-1,5))[:,1:,]
  return nd,el

# Read and find mapping betwee
def proc_gID(il_file,model_nd,model_el):
  reader = vtkXMLPolyDataReader()
  reader.SetFileName(il_file)
  reader.Update()
  # Extract node coordinates
  il_nd = vtk_to_numpy(reader.GetOutput().GetPoints().GetData())
  # WARNING: Assumes the elements are all TRIANGLES
  il_el= vtk_to_numpy(reader.GetOutput().GetPolys().GetData()).reshape((-1,4))[:,1:,]
  # Init Global ID
  il_gID = np.zeros(len(il_nd))

  # Loop through the inlet nodes and find the correspondence  
  for loopA in range(len(il_nd)):
    il_coords = il_nd[loopA]
    for loopB in range(len(model_nd)):
      model_coords = model_nd[loopB]
      dist = np.sqrt(np.sum((model_coords-il_coords)**2))
      if(dist < 1e-8):
        # Node number is zero based in feMorph!!!
        il_gID[loopA] = loopB+1

  return il_nd,il_el,il_gID


def process_faces(nd,el,gID):

  points = vtk.vtkPoints()
  cells  = vtk.vtkCellArray()

  # insert points
  for loopA in range(len(nd)):
    points.InsertNextPoint(nd[loopA,0], nd[loopA,1], nd[loopA,2])

  # insert the quad cells
  for loopA in range(len(el)):
    # create a new quad cell
    tri = vtk.vtkTriangle()
    tri.GetPointIds().SetId(0, el[loopA,0])
    tri.GetPointIds().SetId(1, el[loopA,1])
    tri.GetPointIds().SetId(2, el[loopA,2])
    cells.InsertNextCell(tri)

  # create the point data
  data = numpy_to_vtk(gID)
  data.SetNumberOfComponents(1)
  data.SetNumberOfTuples(len(gID))
  data.SetName('GlobalNodeID')

  # create the grid from the points and cells
  grid = vtk.vtkPolyData()
  grid.SetPoints(points)
  grid.SetPolys(cells)
  grid.GetPointData().SetScalars(data)

  return grid

def write_polydata(file_out,grid):

  # write the grid  
  writer = vtk.vtkXMLPolyDataWriter()
  writer.SetFileName(file_out)
  writer.SetInputData(grid)
  writer.Write()

def apply_velocity(vtk_face):

  num_nodes = vtk_face.GetNumberOfPoints()
  face_vels = np.zeros((num_nodes,3))
  face_vels[:,0] = 1.0

  data = numpy_to_vtk(face_vels)
  data.SetNumberOfComponents(3)
  data.SetNumberOfTuples(num_nodes)
  # Add velocity at time 0.0
  data.SetName('velocity_0.00')
  vtk_face.GetPointData().AddArray(data)
  # Add velocity at time 100.0
  data2 = vtk.vtkDoubleArray()
  data2.DeepCopy(data)
  data2.SetName('velocity_100.00')
  vtk_face.GetPointData().AddArray(data2)

  return vtk_face

# =========
# MAIN CODE
# =========
if __name__ == "__main__":

  # Set file name
  vtu_file  = 'cube.vtu'
  il_file   = 'top_surface.vtp'
  wall_file = 'walls.vtp'
  # output files
  wall_file_out = 'walls_bc.vtp'
  inlet_file_out = 'top_bc.vtp'

  # Read VTU
  print("Reading VTU model...")
  model_nd,model_el = read_nd_el(vtu_file)

  # Read INLET
  print("Reading inlet...")
  il_nd,il_el,il_gID = proc_gID(il_file,model_nd,model_el)

  # Read WALL
  print("Reading wall...")
  wall_nd,wall_el,wall_gID = proc_gID(wall_file,model_nd,model_el)

  # Write INLET
  print("Writing inlet...")
  il_faces = process_faces(il_nd,il_el,il_gID)
  il_faces = apply_velocity(il_faces)
  write_polydata(inlet_file_out,il_faces)

  # Write WALLS
  print("Writing wall...")
  wall_faces = process_faces(wall_nd,wall_el,wall_gID)
  write_polydata(wall_file_out,wall_faces)

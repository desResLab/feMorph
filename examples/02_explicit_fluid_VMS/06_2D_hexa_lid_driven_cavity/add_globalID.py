import vtk
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy
from vtk.util.numpy_support import numpy_to_vtk

# For HEXA8
num_vol = 8
num_sur = 4
# For TET4
# num_vol = 4
# num_sur = 3

def create_node_id(data):
  # Get the number of points in the dataset
  num_points = data.GetNumberOfPoints()

  # Create a new scalar array for the node numbers
  node_numbers = vtk.vtkIntArray()
  node_numbers.SetNumberOfComponents(1)
  node_numbers.SetName("GlobalNodeID")

  # Generate node numbers and assign them to the scalar array
  for i in range(num_points):
      # Global Node ID is equal to the node number PLUS 1!!!
      node_numbers.InsertNextValue(i+1)

  # Add the scalar array to the point data
  data.GetPointData().AddArray(node_numbers)

def find_tet(conn,el_at_nodes):
  res1 = list(set(el_at_nodes[conn[0]]) & set(el_at_nodes[conn[1]]))
  res2 = list(set(res1) & set(el_at_nodes[conn[2]]))
  if(len(res2) != 1):
    return -1
  else:
    return res2[0]

def create_element_id(data,vol_id):
  num_points = data.GetNumberOfPoints()
  
  # Get entity ID
  cell_entity_id = data.GetCellData().GetArray('CellEntityIds')

  # init empty list of lists
  el_at_nodes = [[] for x in range(num_points)]
  
  # Create list of tet4 connected to nodes
  num_cells = data.GetNumberOfCells()
  for i in range(num_cells):
    cell = data.GetCell(i)  
    # Check if the cell is a tet4
    if((cell.GetNumberOfPoints() == num_vol)and(cell_entity_id.GetValue(i) == vol_id)):
      for j in range(cell.GetNumberOfPoints()):
        curr_node = cell.GetPointId(j)
        el_at_nodes[curr_node].append(i)

  cell_tet = []
  # Create list of tet4 connected to nodes
  num_cells = data.GetNumberOfCells()
  for i in range(num_cells):
    cell = data.GetCell(i)
    if(cell.GetNumberOfPoints() == num_sur):
      # Get connectivities for current cell
      conn = []
      for j in range(cell.GetNumberOfPoints()):
        conn.append(cell.GetPointId(j))
      # Find tet element number       
      cell_tet.append(find_tet(conn,el_at_nodes))
    else:
      cell_tet.append(-1)

  # Return GlobalElementID for current triangle mesh
  return cell_tet

def get_volume_id(data):
  num_cells = data.GetNumberOfCells()
  cell_entity_id = data.GetCellData().GetArray('CellEntityIds')

  vol_id_list = []
  for i in range(num_cells):
    cell = data.GetCell(i)
    if(cell.GetNumberOfPoints() == num_vol):
      vol_id_list.append(cell_entity_id.GetValue(i))
  # Find unique ids
  return np.unique(vol_id_list)

def export_volume_mesh(data,vol_id):
  # Find the list of tet4 with the vol_id index
  num_cells = data.GetNumberOfCells()
  cell_entity_id = data.GetCellData().GetArray('CellEntityIds')

  el_list = []
  for i in range(num_cells):
    cell = data.GetCell(i)
    if(cell.GetNumberOfPoints() == num_vol)and(cell_entity_id.GetValue(i) == vol_id):
      el_list.append(i)      

  # Select it
  selection_node = vtk.vtkSelectionNode()
  selection_node.SetFieldType(vtk.vtkSelectionNode.CELL)
  selection_node.SetContentType(vtk.vtkSelectionNode.INDICES)
  selection_node.SetSelectionList(numpy_to_vtk(el_list))
  selection = vtk.vtkSelection()
  selection.AddNode(selection_node)

  # Create the extract selection filter
  extract_selection = vtk.vtkExtractSelection()
  extract_selection.SetInputData(0, data)
  extract_selection.SetInputData(1, selection)
  extract_selection.Update()

  # Save to VTU
  fn = 'volume_'+str(vol_id)+'.vtu'
  writer = vtk.vtkXMLUnstructuredGridWriter()
  writer.SetFileName(fn)
  writer.SetInputData(extract_selection.GetOutput())
  writer.Write()

  return extract_selection.GetOutput()

def export_by_cell_id(data,vol_id,el_id):
  # Get total number of cells
  num_cells = data.GetNumberOfCells()

  # Create CellGlobalID 
  cell_id = vtk.vtkIntArray()
  cell_id.SetNumberOfComponents(1)
  cell_id.SetName("GlobalElementID")
  # Generate node numbers and assign them to the scalar array
  for i in range(num_cells):
    # The Global Element ID is equal to the element number PLUS 1!!!
    cell_id.InsertNextValue(el_id[i]+1)
  # Add the scalar array to the point data
  data.GetCellData().AddArray(cell_id)

  # Find unique ids
  entity_id = numpy_to_vtk(data.GetCellData().GetArray('CellEntityIds'))
  all_id = np.unique(entity_id)
  for i in all_id:
    # Get all elements with this entity ID
    sel_list = np.where(np.logical_and((entity_id == i),(el_id != -1)))[0]

    if(len(sel_list) > 0):
      # Set up the selection
      selection_node = vtk.vtkSelectionNode()
      selection_node.SetFieldType(vtk.vtkSelectionNode.CELL)
      selection_node.SetContentType(vtk.vtkSelectionNode.INDICES)
      selection_node.SetSelectionList(numpy_to_vtk(sel_list))

      selection = vtk.vtkSelection()
      selection.AddNode(selection_node)

      # Create the extract selection filter
      extract_selection = vtk.vtkExtractSelection()
      extract_selection.SetInputData(0, data)
      extract_selection.SetInputData(1, selection)
      extract_selection.Update()

      # Extract surface
      surf = vtk.vtkGeometryFilter()
      surf.SetInputData(extract_selection.GetOutput());
      surf.Update()

      # Save to VTU
      fn = 'surf_vol_'+str(vol_id)+'_entityId_'+str(i)+'.vtp'
      writer = vtk.vtkXMLPolyDataWriter()
      writer.SetFileName(fn)
      writer.SetInputData(surf.GetOutput())
      writer.Write()

def extract_vol_and_boundaries(data,vol_id,el_id_curr_vol):
  num_cells = data.GetNumberOfCells()
  cell_entity_id = data.GetCellData().GetArray('CellEntityIds')

  el_list = []
  # Add cells for current volume
  count_face = 0
  for i in range(num_cells):
    cell = data.GetCell(i)
    is_ok_volume = (cell.GetNumberOfPoints() == num_vol)and(cell_entity_id.GetValue(i) == vol_id)
    is_ok_face = (cell.GetNumberOfPoints() == num_sur)and(el_id_curr_vol[i] > -1)
    if(is_ok_face):
      count_face += 1
    if(is_ok_volume or is_ok_face):
      el_list.append(i)    

  # Create selection
  selection_node = vtk.vtkSelectionNode()
  selection_node.SetFieldType(vtk.vtkSelectionNode.CELL)
  selection_node.SetContentType(vtk.vtkSelectionNode.INDICES)
  selection_node.SetSelectionList(numpy_to_vtk(el_list))

  selection = vtk.vtkSelection()
  selection.AddNode(selection_node)

  # Create the extract selection filter
  extract_selection = vtk.vtkExtractSelection()
  extract_selection.SetInputData(0, data)
  extract_selection.SetInputData(1, selection)
  extract_selection.Update()

  return extract_selection.GetOutput(),count_face

 
# ============
# MAIN ROUTINE
# ============
if __name__ == "__main__":
  # Path to the input file
  input_file = 'complete.vtk'

  # Create a reader for the VTK file
  reader = vtk.vtkGenericDataObjectReader()
  reader.SetFileName(input_file)
  reader.Update()
  # Get the output data from the reader
  data = reader.GetOutput()

  # Get number of ID for solids 
  tet_volume_ids = get_volume_id(data)

  for vol_id in tet_volume_ids:

    # Get list of elements connected to volume vol_id
    el_id_curr_vol = create_element_id(data,vol_id)

    # Extract Volume with attached elements
    vol_and_faces,tot_faces = extract_vol_and_boundaries(data,vol_id,el_id_curr_vol)

    # Create GlobalNodeID for all nodes
    create_node_id(vol_and_faces)

    # Export VTP by CellEntityIds
    only_vol_mesh = export_volume_mesh(vol_and_faces,vol_id)

    # Create GlobalElementID for all surface elements
    el_id = np.array(create_element_id(vol_and_faces,vol_id))

    # Shift non negative components
    el_id[el_id > -1] = el_id[el_id > -1]-tot_faces

    # Export VTP by CellEntityIds
    export_by_cell_id(vol_and_faces,vol_id,el_id)


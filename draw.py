#!/usr/bin/env python

import vtk

# Create several unstructured grids each containing a cell of a
# different type.
tetraPoints = vtk.vtkPoints()
tetraPoints.SetNumberOfPoints(4)
tetraPoints.InsertPoint(0, 0, 0, 0)
tetraPoints.InsertPoint(1, 1, 0, 0)
tetraPoints.InsertPoint(2, .5, 1, 0)
tetraPoints.InsertPoint(3, .5, .5, 1)
aTetra = vtk.vtkTetra()
aTetra.GetPointIds().SetId(0, 0)
aTetra.GetPointIds().SetId(1, 1)
aTetra.GetPointIds().SetId(2, 2)
aTetra.GetPointIds().SetId(3, 3)
aTetraGrid = vtk.vtkUnstructuredGrid()
aTetraGrid.Allocate(1, 1)
aTetraGrid.InsertNextCell(aTetra.GetCellType(), aTetra.GetPointIds())
aTetraGrid.SetPoints(tetraPoints)
aTetraMapper = vtk.vtkDataSetMapper()
aTetraMapper.SetInputData(aTetraGrid)
aTetraActor = vtk.vtkActor()
aTetraActor.SetMapper(aTetraMapper)
aTetraActor.AddPosition(4, 0, 0)
aTetraActor.GetProperty().SetDiffuseColor(0, 1, 0)

# Create the usual rendering stuff.
ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
renWin.SetSize(300, 150)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

ren.SetBackground(.1, .2, .4)

ren.AddActor(aTetraActor)

ren.ResetCamera()
ren.GetActiveCamera().Azimuth(30)
ren.GetActiveCamera().Elevation(20)
ren.GetActiveCamera().Dolly(2.8)
ren.ResetCameraClippingRange()

# Render the scene and start interaction.
iren.Initialize()
renWin.Render()
iren.Start()

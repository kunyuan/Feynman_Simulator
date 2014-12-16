# Visualizing data from a molecular dynamics simulation
# Data provided by Daniel Spangberg

from vtk import *
from ReadAtoms import *

# Read the data into a vtkPolyData using the functions in ReadPoints.py
data = vtk.vtkPolyData()
data.SetPoints(readPoints("Coordinates.txt"))
data.GetPointData().SetScalars(readScalars("Radii.txt")) 
data.SetLines(readConnections("Connections.txt"))

# Put spheres at each point in the dataset.
# The size and color of a sphere is determined by the
# scalar value (radius) of the corresponding point.

ball = vtk.vtkSphereSource()
ball.SetRadius(0.25)
ball.SetThetaResolution(8)
ball.SetPhiResolution(8)

ballGlyph = vtk.vtkGlyph3D()
ballGlyph.SetInput(data)
ballGlyph.SetSourceConnection(ball.GetOutputPort())
ballGlyph.SetScaleModeToScaleByScalar()
ballGlyph.SetColorModeToColorByScalar()
ballGlyph.SetScaleFactor(2.0)

colorTransferFunction = vtk.vtkColorTransferFunction()
colorTransferFunction.AddRGBPoint(0.37, 1.0, 1.0, 1.0)
colorTransferFunction.AddRGBPoint(0.68, 1.0, 1.0, 0.0)
colorTransferFunction.AddRGBPoint(0.73, 0.0, 0.0, 1.0)
colorTransferFunction.AddRGBPoint(0.74, 1.0, 0.0, 0.0)
colorTransferFunction.AddRGBPoint(0.77, 0.0, 0.0, 1.0)
colorTransferFunction.AddRGBPoint(2.0, 0.0, 1.0, 0.0)

ballMapper = vtkPolyDataMapper()
ballMapper.SetInputConnection(ballGlyph.GetOutputPort())
ballMapper.SetLookupTable(colorTransferFunction)

ballActor = vtkActor()
ballActor.SetMapper(ballMapper)

#Create tubes between connected molecules
tubeFilter = vtkTubeFilter()
tubeFilter.SetInput(data)
tubeFilter.SetRadius(0.15)
tubeFilter.SetNumberOfSides(7)

tubeMapper = vtkPolyDataMapper()
tubeMapper.SetInputConnection(tubeFilter.GetOutputPort())

# We want to be able to set the color ourselves!
tubeMapper.ScalarVisibilityOff() 

tubeActor = vtkActor()
tubeActor.SetMapper(tubeMapper)
tubeActor.GetProperty().SetColor(0.8,0.8,0.8)

tubeActor.GetProperty().SetSpecularColor(1, 1, 1)
tubeActor.GetProperty().SetSpecular(0.3)
tubeActor.GetProperty().SetSpecularPower(20)
tubeActor.GetProperty().SetAmbient(0.2)
tubeActor.GetProperty().SetDiffuse(0.8)

# A Bounding box
outlineData = vtkOutlineFilter()
outlineData.SetInputConnection(ballGlyph.GetOutputPort())

outlineMapper = vtkPolyDataMapper()
outlineMapper.SetInputConnection(outlineData.GetOutputPort())

outlineActor = vtkActor()
outlineActor.SetMapper(outlineMapper)
outlineActor.GetProperty().SetColor(0, 0, 0)

# Create the Renderer, Window and Interator
ren = vtkRenderer()
ren.AddActor(ballActor)
ren.AddActor(outlineActor)
ren.AddActor(tubeActor)
ren.SetBackground(0.4, 0.4, 0.4)

renWin = vtkRenderWindow()
renWin.AddRenderer(ren)
renWin.SetWindowName("Atoms")
renWin.SetSize(600, 600)

iren = vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)
iren.Initialize()
iren.Start()

#!/usr/bin/env python
# Visualizing data from a molecular dynamics simulation
# Data provided by Daniel Spangberg
from vtk import *
import argparse, os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
#add parentdir into PYTHONPATH, where IO module can be found
import IO

version=vtkVersion().GetVTKMajorVersion()

def Read(filename):
    Dict=IO.LoadDict(filename)
    data=Dict["Points"]
    points = vtkPoints()
    color = vtkFloatArray()
    captions = []
    SubLat=set()
    for vec, coord, sub in data:
        x, y = vec[0], vec[1]
        if len(vec)==3:
            z=vec[2]
        else:
            z=0
        points.InsertNextPoint(x, y, z)
        color.InsertNextValue(sub)
        captions.append([(x,y,z),str(tuple(coord))])
        SubLat.add(sub)

    Edges=vtkCellArray()
    lines=Dict["Lines"]
    for line, sub in lines:
        Edges.InsertNextCell(2)
        Edges.InsertCellPoint(line[0])
        Edges.InsertCellPoint(line[1])

    Interaction=vtkCellArray()
    lines=Dict["Interaction"]
    for vec,coord, sub in lines:
        Interaction.InsertNextCell(2)
        Interaction.InsertCellPoint(coord[0])
        Interaction.InsertCellPoint(coord[1])
    return (SubLat, points, captions, color, Edges, Interaction)

def Plot(InputFile, HasCaption):
    # Read the data into a vtkPolyData using the functions in ReadPoints.py
    SubLat, points, captions, color, Edges, Interaction=Read(InputFile)

    data = vtk.vtkPolyData()
    data.SetPoints(points)
    data.GetPointData().SetScalars(color) 
    data.SetLines(Edges)
    data2 = vtk.vtkPolyData()
    data2.SetPoints(points)
    data2.SetLines(Interaction)

    ball = vtk.vtkSphereSource()
    ball.SetRadius(0.04)
    ball.SetThetaResolution(8)
    ball.SetPhiResolution(8)

    ballGlyph = vtk.vtkGlyph3D()
    if version<=5:
        ballGlyph.SetInput(data)
    else:
        ballGlyph.SetInputData(data)
    ballGlyph.SetSourceConnection(ball.GetOutputPort())
    ballGlyph.SetScaleModeToDataScalingOff()
    ballGlyph.SetColorModeToColorByScalar()
    ballGlyph.SetScaleFactor(1.0)

    colorTransferFunction = vtk.vtkColorTransferFunction()
    colorTransferFunction.AddRGBPoint(0, 1.0, 0.0, 0.0)
    colorTransferFunction.AddRGBPoint(1, 0.0, 0.0, 1.0)
    colorTransferFunction.AddRGBPoint(2, 0.0, 1.0, 0.0)
    colorTransferFunction.AddRGBPoint(3, 1.0, 1.0, 0.0)
    colorTransferFunction.AddRGBPoint(4, 0.0, 1.0, 1.0)
    colorTransferFunction.AddRGBPoint(5, 1.0, 1.0, 0.0)

    ballMapper = vtkPolyDataMapper()
    ballMapper.SetInputConnection(ballGlyph.GetOutputPort())
    ballMapper.SetLookupTable(colorTransferFunction)
    ballActor = vtkActor()
    ballActor.SetMapper(ballMapper)

#lines within unit cell
    tubeFilter = vtkTubeFilter()
    if version<=5:
        tubeFilter.SetInput(data)
    else:
        tubeFilter.SetInputData(data)
    tubeFilter.SetRadius(0.01)
    tubeFilter.SetNumberOfSides(7)
    tubeMapper = vtkPolyDataMapper()
    tubeMapper.SetInputConnection(tubeFilter.GetOutputPort())
# We want to be able to set the color ourselves!
    tubeMapper.ScalarVisibilityOff() 
    tubeActor = vtkActor()
    tubeActor.SetMapper(tubeMapper)
    tubeActor.GetProperty().SetColor(0.9,0.9,0.9)
    tubeActor.GetProperty().SetOpacity(0.3)

#Interaction lines
    tube2 = vtkTubeFilter()
    if version<=5:
        tube2.SetInput(data2)
    else:
        tube2.SetInputData(data2)
    tube2.SetRadius(0.02)
    tube2.SetNumberOfSides(7)
    mapper = vtkPolyDataMapper()
    mapper.SetInputConnection(tube2.GetOutputPort())
    mapper.ScalarVisibilityOff() 
    InteractionActor = vtkActor()
    InteractionActor.SetMapper(mapper)
    InteractionActor.GetProperty().SetColor(0.95,0.95,0.95)

# A Bounding box
    outlineData = vtkOutlineFilter()
    outlineData.SetInputConnection(ballGlyph.GetOutputPort())
    outlineMapper = vtkPolyDataMapper()
    outlineMapper.SetInputConnection(outlineData.GetOutputPort())
    outlineActor = vtkActor()
    outlineActor.SetMapper(outlineMapper)
    outlineActor.GetProperty().SetColor(0.0, 0.0, 0.0)

    txtProp=vtkTextProperty()
    txtProp.BoldOff()
    txtProp.ItalicOff()
    txtProp.ShadowOff()
    txtActor=[]
    for vec, content in captions:
        #TODO: different captions now have their own actor, 
        #it would be better to merge them into one actor to have better efficiency
        txt=vtkCaptionActor2D()
        #txt.LeaderOff()
        txt.BorderOff()
        txt.SetAttachmentPoint(*vec)
        txt.SetCaptionTextProperty(txtProp)
        txt.SetCaption(content)
        txt.SetWidth(0.05)
        txt.SetHeight(0.04)
        txt.GetProperty().SetColor(0.0, 0.0, 0.0)
        txt.GetProperty().SetOpacity(0.9)
        txtActor.append(txt)

    legend=vtkLegendBoxActor()
    legend.SetNumberOfEntries(4)
    for e in SubLat:
        legend.SetEntry(e, ball.GetOutput(), str(e),
                colorTransferFunction.GetColor(e)) 
    legend.BorderOff()
    legend.SetWidth(0.1)
    legend.SetHeight(0.1)
    legend.SetDisplayPosition(10,5)
    #txtProp=legend.GetEntryTextProperty()
# Create the Renderer, Window and Interator
    ren = vtkRenderer()
    ren.AddActor(ballActor)
    ren.AddActor(outlineActor)
    ren.AddActor(tubeActor)
    if HasCaption:
        for e in txtActor:
            ren.AddActor(e)
    ren.AddActor(InteractionActor)
    ren.AddActor(legend)
    ren.SetBackground(0.4, 0.4, 0.4)
    return ren

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="use file path to find the input file")
    parser.add_argument("-nc", "--nocaption", action="store_true", help="turn of caption")
    args = parser.parse_args()
    InputFile=os.path.abspath(args.file)

    ren=Plot(InputFile, not args.nocaption)
    renWin = vtkRenderWindow()
    renWin.AddRenderer(ren)
    renWin.SetWindowName("Lattice")
    renWin.SetSize(1000, 800)
    iren = vtkRenderWindowInteractor()
    #iren.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
    iren.SetRenderWindow(renWin)

    iren.Initialize()
    iren.Start()

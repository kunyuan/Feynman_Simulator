# Reading Atom data

import string
import vtk

def readPoints(filename):
    points = vtk.vtkPoints()
    scalars = vtk.vtkFloatArray()
    captions = []
    file = open(filename)
    lines = file.read()
    data = eval(lines)
    for vec, coord, sub in data:
        x, y, z = float(vec[0]), float(vec[1]), float(vec[2])
        points.InsertNextPoint(x, y, z)
        scalars.InsertNextValue(sub)
        captions.append([(x,y,z),str(tuple(coord))])
    return (points, captions, scalars)

def readConnections(filename):
    connections=vtk.vtkCellArray()
    #file = open(filename)

    #line = file.readline()
    #while line:
        #data = string.split(line)
        #if data and data[0] != '#':
            #a, b = int(data[0]), int(data[1])
            #connections.InsertNextCell(2)
            #connections.InsertCellPoint(a)
            #connections.InsertCellPoint(b)
        #line = file.readline()
    connections.InsertNextCell(2)
    connections.InsertCellPoint(0)
    connections.InsertCellPoint(1)
    return connections

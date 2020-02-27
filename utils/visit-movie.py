#!/usr/bin/env python

# Make a movie from a QC-0 simulation

# See
# <https://visit-sphinx-github-user-manual.readthedocs.io/en/develop/cli_manual/index.html>
# for documentation on VisIt's CLI.

# module load visit/3.1.0
# visit -cli -nosplash -nowin -s visit-movie.py
# ffmpeg -framerate 10 -i qc0%04d.png qc0.mp4

OpenDatabase("/gpfs/eschnetter/simulations/qc0/output-0000/qc0/qc0.visit")

AddPlot("Pseudocolor", "ADMBASE__alp") 
p = PseudocolorAttributes()
p.colorTableName = "plasma"     
p.min, p.max = 0, 1
p.minFlag = p.maxFlag = True
SetPlotOptions(p)

AddOperator("Slice")
a = SliceAttributes()
a.originType = a.Point
a.normal, a.upAxis = (0,0,1), (0,1,0)  
SetOperatorOptions(a)

DrawPlots()
v = GetView2D()
v.windowCoords = (-10, 10, -10, 10)
SetView2D(v)

ann = AnnotationAttributes()
ann.databaseInfoFlag = False
ann.userInfoFlag = False
ann.timeInfoFlag = False
SetAnnotationAttributes(ann)

s = SaveWindowAttributes()
s.format = s.PNG
s.fileName = "qc0"
s.width, s.height = 1024, 768
s.screenCapture = 0
SetSaveWindowAttributes(s)

for state in range(0, TimeSliderGetNStates(), 10):
    SetTimeSliderState(state)
    SaveWindow()

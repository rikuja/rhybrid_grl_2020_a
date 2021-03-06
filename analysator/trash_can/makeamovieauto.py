#!/usr/bin/python -i

import sys
import visit as vis
import loadvisitsettings as visSettings
import subprocess
import numpy as np
from sets import Set
from vlsvfile import *
from useroptions import *


def make_movie_auto( variableName, boundaryBox, vlsvFileName, inputDirectory, inputFileName, outputDirectory, outputFileName, colorTableName="hot_desaturated", startFrame=-1, endFrame=-1, thresholdCoefficient=0.6 ):
   '''
   Function for making a movie
   
   Arguments:
   :param variableName                  Name of the variable
   :param boundaryBox                   Box for collecting min and max threshold (The movie will focus on that area)
   :param vlsvFileName                  Name of a vlsv file where the function collects the threshold for the boundary box
   :param inputDirectory                Path to input vlsv/silo files
   :param inputFileName                 Name of the file(s) so for example if the filenames are bulk.0000.silo, bulk.0001.silo, .. then inputFileName=\"bulk.*.silo\""
   :param outputDirectory               Path to output directory
   :param outputFileName                Name of the output file
   :param colorTableName="hot_desaturated"  Color table for the plots
   :param thresholdCoefficient          Sets the coefficient for a covariant collected from the values from boundary box. The lower this is, the more focused the movie will be on the boundary box area
   :param startFrame=-1                 Starting frame of the movie (-1 equals 0)
   :param endFrame=-1                   Starting frame of the movie (-1 equals last frame)
   '''
   if thresholdCoefficient < 0:
      print "thresholdCoefficient must be non-negative!"
      return
   # OPTIONS
   #################################################################
   # Input frame properties
   _startFrame = startFrame # Note: if _startFrame is set to -1 the start frame gets set to 0
   _endFrame = endFrame # Note: if _endFrame is set to -1 the _endFrame is automatically the number of frames in the database
   
   # Input variable
   _variableName = variableName

   # Input directory and file names
   #_outputDir = "/home/hannukse/MOVINGFRAME_MOVIES/AAJ_BZ_REMAKE/" # Set the output directory (Where .png s are saved)
   _outputDir = outputDirectory
   #_outputFileName = "BZ_FORESHOCK_2_" # The file names for the png files. These for ex. will be saved visit0000.png, visit0001.png, ..
   _outputFileName = outputFileName # The file names for the png files.
   #databaseName = "localhost:/home/hannukse/meteo/stornext/field/vlasiator/2D/AAJ/silo_files/bulk.*.silo database" # For navigating to the silo files
   databaseName = "localhost:" + inputDirectory + inputFileName + " database" # For navigating to the silo files
   # Note: a slice of the plot in z-axis is taken automatically
   #################################################################


   # Get the min and max values:
   # Get all cell ids within the boundary box:
   vlsvReader = VlsvReader(vlsvFileName)
   # Get global boundaries:
   # Get xmax, xmin and xcells_ini
   xmax = vlsvReader.read_parameter(name="xmax")
   xmin = vlsvReader.read_parameter(name="xmin")
   xcells = vlsvReader.read_parameter(name="xcells_ini")
   # Do the same for y
   ymax = vlsvReader.read_parameter(name="ymax")
   ymin = vlsvReader.read_parameter(name="ymin")
   ycells = vlsvReader.read_parameter(name="ycells_ini")
   # And for z
   zmax = vlsvReader.read_parameter(name="zmax")
   zmin = vlsvReader.read_parameter(name="zmin")
   zcells = vlsvReader.read_parameter(name="zcells_ini")
   #Calculate cell lengths:
   cell_lengths = np.array([(xmax - xmin)/(float)(xcells), (ymax - ymin)/(float)(ycells), (zmax - zmin)/(float)(zcells)])
   # Get cell indices:
   cell_indice_bounds = np.array([(int)(((float)(boundaryBox[0]) - xmin) / (float)(cell_lengths[0])),(int)(((float)(boundaryBox[1]) - xmin) / (float)(cell_lengths[0])), (int)(((float)(boundaryBox[2]) - ymin) / (float)(cell_lengths[1])), (int)(((float)(boundaryBox[3]) - ymin) / (float)(cell_lengths[1])), (int)(((float)(boundaryBox[4]) - zmin) / (float)(cell_lengths[2])), (int)(((float)(boundaryBox[5]) - zmin) / (float)(cell_lengths[2]))])
   # Get every cell id within the boundary box:
   cellids = []
   cell_indice = np.array([cell_indice_bounds[0], cell_indice_bounds[2], cell_indice_bounds[4]])
   while True:
      cellids.append(cell_indice[0] + cell_indice[1] * xcells + cell_indice[2] * xcells * ycells + 1)
      if cell_indice[0] < cell_indice_bounds[1]:
         cell_indice[0] = cell_indice[0] + 1
      elif cell_indice[1] < cell_indice_bounds[3]:
         cell_indice[1] = cell_indice[1] + 1
         cell_indice[0] = cell_indice_bounds[0]
      elif cell_indice[2] < cell_indice_bounds[5]:
         cell_indice[2] = cell_indice[2] + 1
         cell_indice[1] = cell_indice_bounds[1]
         cell_indice[0] = cell_indice_bounds[0]
      else:
         # Indice out of bounds -- got all cell ids
         break
   # Convert cell ids into set:
   cellids = Set(cellids)
   cellidlocations = []
   # Get all of the cell ids locations:
   allcellids = vlsvReader.read(name="SpatialGrid",tag="MESH")
   for i in xrange(len(allcellids)):
      if allcellids[i] in cellids:
         #This cell id is within the user-given boundary
         cellidlocations.append(allcellids[i])
   # Get all of the values:
   allvalues = vlsvReader.read_variables(name=_variableName)
   values = []
   # Get the values of the cell ids within the boundary
   for i in cellidlocations:
      values.append(allvalues[i])
   # We now have all the cell ids (and their locations in the arrays) from the area, set min and max thresholds:
   meanValue = np.mean(values)
   standardDeviationValue = np.std(values)
   maxValue = meanValue + (float)(thresholdCoefficient) * standardDeviationValue
   minValue = meanValue - (float)(thresholdCoefficient) * standardDeviationValue
   # Put threshold values:
   minVariableValue = minValue
   maxVariableValue = maxValue

   # LaunchNowin(vdir=visitBinDirectory)
   #dx = speedX * frameInSeconds # Note: This is in meters per frame!
   #dy = speedY * frameInSeconds # Note: This is in meters per frame!
   #LaunchNowin(vdir="/usr/local/visit/bin")
   #Set up window and annotations
   #vis.LaunchNowin(vdir="/usr/local/visit/bin")
   vis.OpenDatabase(databaseName, 0)

   #Load settings
   visSettings.load_visit_settings()

   vis.AddPlot("Pseudocolor", _variableName, 1, 1) #CONTINUE
   vis.SetActivePlots(0)
   vis.PseudocolorAtts = vis.PseudocolorAttributes()
   vis.PseudocolorAtts.legendFlag = 1
   vis.PseudocolorAtts.lightingFlag = 1
   vis.PseudocolorAtts.minFlag = 1
   vis.PseudocolorAtts.maxFlag = 1
   vis.PseudocolorAtts.centering = vis.PseudocolorAtts.Natural  # Natural, Nodal, Zonal
   vis.PseudocolorAtts.scaling = vis.PseudocolorAtts.Linear  # Linear, Log, Skew
   vis.PseudocolorAtts.limitsMode = vis.PseudocolorAtts.CurrentPlot  # OriginalData, CurrentPlot
   vis.PseudocolorAtts.min = minVariableValue
   vis.PseudocolorAtts.max = maxVariableValue
   vis.PseudocolorAtts.pointSize = 0.05
   vis.PseudocolorAtts.pointType = vis.PseudocolorAtts.Point  # Box, Axis, Icosahedron, Point, Sphere
   vis.PseudocolorAtts.skewFactor = 1
   vis.PseudocolorAtts.opacity = 1
   vis.PseudocolorAtts.colorTableName = colorTableName
   vis.PseudocolorAtts.invertColorTable = 0
   vis.PseudocolorAtts.smoothingLevel = 0
   vis.PseudocolorAtts.pointSizeVarEnabled = 0
   vis.PseudocolorAtts.pointSizeVar = "default"
   vis.PseudocolorAtts.pointSizePixels = 2
   vis.PseudocolorAtts.lineStyle = vis.PseudocolorAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
   vis.PseudocolorAtts.lineWidth = 0
   vis.PseudocolorAtts.opacityType = vis.PseudocolorAtts.Explicit  # Explicit, ColorTable
   vis.SetPlotOptions(vis.PseudocolorAtts)

   
   vis.SetActivePlots(0)
   vis.AddOperator("Slice", 1)
   vis.SetActivePlots(0)
   vis.SliceAtts = vis.SliceAttributes()
   vis.SliceAtts.originType = vis.SliceAtts.Intercept  # Point, Intercept, Percent, Zone, Node
   vis.SliceAtts.originPoint = (0, 0, 0)
   vis.SliceAtts.originIntercept = 0
   vis.SliceAtts.originPercent = 0
   vis.SliceAtts.originZone = 0
   vis.SliceAtts.originNode = 0
   vis.SliceAtts.normal = (0, 0, 1)
   vis.SliceAtts.axisType = vis.SliceAtts.ZAxis  # XAxis, YAxis, ZAxis, Arbitrary, ThetaPhi
   vis.SliceAtts.upAxis = (0, 1, 0)
   vis.SliceAtts.project2d = 1
   vis.SliceAtts.interactive = 1
   vis.SliceAtts.flip = 0
   vis.SliceAtts.originZoneDomain = 0
   vis.SliceAtts.originNodeDomain = 0
   vis.SliceAtts.meshName = "SpatialGrid"
   vis.SliceAtts.theta = 0
   vis.SliceAtts.phi = 90
   vis.SetOperatorOptions(vis.SliceAtts, 1)
   vis.DrawPlots()
   
   if _endFrame == -1:
      _endFrame = vis.TimeSliderGetNStates() - 1
   
   if _startFrame == -1:
      _startFrame = 0
   
   # Iterate through frames
   for i in xrange(_startFrame, _endFrame+1):
      vis.SetTimeSliderState(i)
      frame = i - _startFrame
      vis.SaveWindowAtts = vis.SaveWindowAttributes()
      vis.SaveWindowAtts.outputToCurrentDirectory = 0
      vis.SaveWindowAtts.outputDirectory = _outputDir
      vis.SaveWindowAtts.fileName = _outputFileName
      vis.SaveWindowAtts.family = 1
      vis.SaveWindowAtts.format = vis.SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
      vis.SaveWindowAtts.width = 3000
      vis.SaveWindowAtts.height = 300
      vis.SaveWindowAtts.screenCapture = 0
      vis.SaveWindowAtts.saveTiled = 0
      vis.SaveWindowAtts.quality = 100
      vis.SaveWindowAtts.progressive = 0
      vis.SaveWindowAtts.binary = 0
      vis.SaveWindowAtts.stereo = 0
      vis.SaveWindowAtts.compression = vis.SaveWindowAtts.PackBits  # None, PackBits, Jpeg, Deflate
      vis.SaveWindowAtts.forceMerge = 0
      vis.SaveWindowAtts.resConstraint = vis.SaveWindowAtts.ScreenProportions  # NoConstraint, EqualWidthHeight, ScreenProportions
      vis.SaveWindowAtts.advancedMultiWindowSave = 0
      vis.SetSaveWindowAttributes(vis.SaveWindowAtts)
      vis.SaveWindow()
   vis.DeleteActivePlots()
   vis.CloseDatabase(databaseName)
   # Make the movie:
   #subprocess.call("./moviecompilescript.sh " + _outputDir + " " + _outputFileName)
   pyVisitPath = "pyVisit/"
   #subprocess.call(pythonLibDirectoryPath + pyVisitPath + "moviecompilescript.sh")
   #subprocess.call(pythonLibDirectoryPath + pyVisitPath + "moviecompilescript.sh " + _outputDir + " " + _outputFileName)
   frameRate = "10"
   subprocess.call([pythonLibDirectoryPath + pyVisitPath + "moviecompilescript.sh", _outputDir, _outputFileName, frameRate])

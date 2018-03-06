#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

## read the .out file with crevasse positions and create a list containing everything
#b0 = list()
#b1 = list()
#s0 = list()
#s1 = list()
#xxx = list()
#zzz = list()
#fhand = open('/Users/yuema/Documents/free_slip_800-700/Remesh/crevs920.out','r')
#for i, line in enumerate(fhand):
#    if i == 0:
#        xxx = line.strip().split()
#        xxx = map(float,xxx)
#    if i == 1:
#        zzz = line.strip().split()
#        zzz = map(float,zzz)
#        for k in range(len(xxx)):
#            b0.append(xxx[k])
#            b0.append(zzz[k])
#            b0.append(0.0)
#    if i == 2:
#        xxx = line.strip().split()
#        xxx = map(float,xxx)
#    if i == 3:
#        zzz = line.strip().split()
#        zzz = map(float,zzz)
#        for k in range(len(xxx)):
#            b1.append(xxx[k])
#            b1.append(zzz[k])
#            b1.append(0.0)
#    if i == 4:
#        xxx = line.strip().split()
#        xxx = map(float,xxx)
#    if i == 5:
#        zzz = line.strip().split()
#        zzz = map(float,zzz)
#        for k in range(len(xxx)):
#            s0.append(xxx[k])
#            s0.append(zzz[k])
#            s0.append(0.0)
#    if i == 6:
#        xxx = line.strip().split()
#        xxx = map(float,xxx)
#    if i == 7:
#        zzz = line.strip().split()
#        zzz = map(float,zzz)
#        for k in range(len(xxx)):
#            s1.append(xxx[k])
#            s1.append(zzz[k])
#            s1.append(0.0)
#fhand.close()

# create a new 'PVD Reader'
tensile920pvd = PVDReader(FileName='/Users/yuema/Documents/free_slip_800-700/Remesh/tensile920.pvd')

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [913, 479]

# get color transfer function/color map for 'f246935'
f246935LUT = GetColorTransferFunction('f246935')
f246935LUT.RGBPoints = [-6170090.235706781, 0.231373, 0.298039, 0.752941, -2694094.588410979, 0.865003, 0.865003, 0.865003, 781901.0588848228, 0.705882, 0.0156863, 0.14902]
f246935LUT.ScalarRangeInitialized = 1.0

# show data in view
tensile920pvdDisplay = Show(tensile920pvd, renderView1)
# trace defaults for the display properties.
tensile920pvdDisplay.ColorArrayName = ['POINTS', 'f_246935']
tensile920pvdDisplay.LookupTable = f246935LUT
tensile920pvdDisplay.ScalarOpacityUnitDistance = 89.01075630298526

# reset view to fit data
renderView1.ResetCamera()

#changing interaction mode based on data extents
renderView1.CameraPosition = [2965.0560215897563, 335.3273779553915, 10000.0]
renderView1.CameraFocalPoint = [2965.0560215897563, 335.3273779553915, 0.0]

# show color bar/color legend
tensile920pvdDisplay.SetScalarBarVisibility(renderView1, True)

# get opacity transfer function/opacity map for 'f246935'
f246935PWF = GetOpacityTransferFunction('f246935')
f246935PWF.Points = [-6170090.235706781, 0.0, 0.5, 0.0, 781901.0588848228, 1.0, 0.5, 0.0]
f246935PWF.ScalarRangeInitialized = 1

# Properties modified on tensile920pvdDisplay
tensile920pvdDisplay.InterpolateScalarsBeforeMapping = 0

# Rescale transfer function
f246935LUT.RescaleTransferFunction(-6e6, 1e6)

# Rescale transfer function
f246935PWF.RescaleTransferFunction(-6e6, 1e6)

# Properties modified on f246935LUT
f246935LUT.RGBPoints = [-6e6, 0.23137254902, 0.298039215686, 0.752941176471, 0, 1.0, 1.0, 1.0, 1e6, 0.705882352941, 0.0156862745098, 0.149019607843]


# create a new 'Calculator'
calculator1 = Calculator(Input=tensile920pvd)
calculator1.Function = ''

# Properties modified on calculator4
calculator1.Function = 'f_246935+1020*9.8*(700-coordsY)*(coordsY<700)'

# get color transfer function/color map for 'Result'
resultLUT = GetColorTransferFunction('Result')
resultLUT.LockDataRange = 1
resultLUT.RGBPoints = [-1258258.6198942522, 0.231373, 0.298039, 0.752941, 498277.8320558979, 1.0, 1.0, 1.0, 2254814.284006048, 0.705882, 0.0156863, 0.14902]
resultLUT.ScalarRangeInitialized = 1.0

# show data in view
calculator1Display = Show(calculator1, renderView1)
# trace defaults for the display properties.
calculator1Display.ColorArrayName = ['POINTS', 'Result']
calculator1Display.LookupTable = resultLUT
calculator1Display.ScalarOpacityUnitDistance = 117.53109761528958

# hide data in view
Hide(tensile920pvd, renderView1)

# show color bar/color legend
calculator1Display.SetScalarBarVisibility(renderView1, True)

# get opacity transfer function/opacity map for 'Result'
resultPWF = GetOpacityTransferFunction('Result')
resultPWF.Points = [-1258258.6198942522, 0.0, 0.5, 0.0, 2254814.284006048, 1.0, 0.5, 0.0]
resultPWF.ScalarRangeInitialized = 1

# Rescale transfer function
resultLUT.RescaleTransferFunction(-500000.0, 500000.0)

# Rescale transfer function
resultPWF.RescaleTransferFunction(-500000.0, 500000.0)


# get color legend/bar for resultLUT in view renderView1
resultLUTColorBar = GetScalarBar(resultLUT, renderView1)
resultLUTColorBar.ComponentTitle = ''

# Properties modified on resultLUTColorBar
resultLUTColorBar.AutoOrient = 0
resultLUTColorBar.Orientation = 'Horizontal'
resultLUTColorBar.Title = ''
resultLUTColorBar.TitleFontSize = 12
resultLUTColorBar.LabelFontSize = 12
resultLUTColorBar.AutomaticLabelFormat = 0
resultLUTColorBar.DrawTickMarks = 0
resultLUTColorBar.DrawTickLabels = 0
resultLUTColorBar.AddRangeLabels = 0
resultLUTColorBar.AspectRatio = 16.0
resultLUTColorBar.Position = [0.27,0.25]

# Properties modified on resultLUT
resultLUT.Annotations = ['-5e5', '-0.5', '0', '0', '5e5', '0.5']
resultLUT.IndexedColors = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

# create a new 'Text'
text1 = Text()

# show data in view
text1Display = Show(text1, renderView1)

# Properties modified on text1Display
text1Display.Color = [0.0, 0.0, 0.0]

# Properties modified on text1
text1.Text = 'Largest Principal Stress (MPa)'

# Properties modified on text1Display
text1Display.WindowLocation = 'AnyLocation'
text1Display.Position = [0.27, 0.67]

# set active source
SetActiveSource(tensile920pvd)

# Properties modified on resultLUTColorBar
resultLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
resultLUTColorBar.LabelColor = [0.0, 0.0, 0.0]

# Properties modified on renderView1
renderView1.UseLight = 0

# Properties modified on renderView1
renderView1.LightSwitch = 1

# Properties modified on renderView1
renderView1.Background = [0.8, 0.8, 0.8]

# create a new 'Box'
box1 = Box()

# Properties modified on box1
box1.XLength = 5100.0
box1.YLength = 700.0
box1.ZLength = 0.0
box1.Center = [2550.0, 350.0, -1.0]

# show data in view
box1Display = Show(box1, renderView1)
# trace defaults for the display properties.
box1Display.ColorArrayName = [None, '']

# change solid color
box1Display.DiffuseColor = [0.0, 0.0, 0.4980392156862745]

## create a new 'Poly Line Source'
#polyLineSource1 = PolyLineSource()
#
## toggle 3D widget visibility (only when running from the GUI)
#Hide3DWidgets(proxy=polyLineSource1)
#
## Properties modified on polyLineSource1
#polyLineSource1.Points = b0
#
## show data in view
#polyLineSource1Display = Show(polyLineSource1, renderView1)
## trace defaults for the display properties.
#polyLineSource1Display.ColorArrayName = [None, '']
#
## change solid color
#polyLineSource1Display.DiffuseColor = [0, 0, 0]
#
## Properties modified on polyLineSource1Display
#polyLineSource1Display.LineWidth = 4.0
#
## create a new 'Poly Line Source'
#polyLineSource2 = PolyLineSource()
#
## toggle 3D widget visibility (only when running from the GUI)
#Hide3DWidgets(proxy=polyLineSource2)
#
## Properties modified on polyLineSource1
#polyLineSource2.Points = b1
#
## show data in view
#polyLineSource2Display = Show(polyLineSource2, renderView1)
## trace defaults for the display properties.
#polyLineSource2Display.ColorArrayName = [None, '']
#
## change solid color
#polyLineSource2Display.DiffuseColor = [0, 0, 0]
#
## Properties modified on polyLineSource1Display
#polyLineSource2Display.LineWidth = 4.0
#
## create a new 'Poly Line Source'
#polyLineSource3 = PolyLineSource()
#
## toggle 3D widget visibility (only when running from the GUI)
#Hide3DWidgets(proxy=polyLineSource3)
#
## Properties modified on polyLineSource1
#polyLineSource3.Points = s0
#
## show data in view
#polyLineSource3Display = Show(polyLineSource3, renderView1)
## trace defaults for the display properties.
#polyLineSource3Display.ColorArrayName = [None, '']
#
## change solid color
#polyLineSource3Display.DiffuseColor = [0, 0, 0]
#
## Properties modified on polyLineSource1Display
#polyLineSource3Display.LineWidth = 4.0
#
## create a new 'Poly Line Source'
#polyLineSource4 = PolyLineSource()
#
## toggle 3D widget visibility (only when running from the GUI)
#Hide3DWidgets(proxy=polyLineSource4)
#
## Properties modified on polyLineSource1
#polyLineSource4.Points = s1
#
## show data in view
#polyLineSource4Display = Show(polyLineSource4, renderView1)
## trace defaults for the display properties.
#polyLineSource4Display.ColorArrayName = [None, '']
#
## change solid color
#polyLineSource4Display.DiffuseColor = [0, 0, 0]
#
## Properties modified on polyLineSource1Display
#polyLineSource4Display.LineWidth = 4.0

# set active source
SetActiveSource(tensile920pvd)

# set active view
SetActiveView(renderView1)

#### saving camera placements for all active views
# current camera placement for renderView4
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [2446.076664329161, 394.4538140333408, 10000.0]
renderView1.CameraFocalPoint = [2446.076664329161, 394.4538140333408, 0.0]
renderView1.CameraParallelScale = 1684.3661561010308

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
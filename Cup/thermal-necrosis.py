
# Python script to calculate thermal necrosis of bone due to heat released  
# from exothermic reaction during polymerization of bone cement

# Thermal damage parameter is based on the info given by:
# Quarini GL, Learmonth ID, and S Gheduzzi. Numerical predictions of the thermal
# behaviour and resultant effects of grouting cements while setting prosthetic 
# components in bone. Proc IMechE Part H: J Engineering in Medicine, 2006, 220:625-634

# Author  : Michael Hogg
# Created : 4 January 2007

# ------------------------------------------------------------------------------------------

from abaqus import *
from abaqusConstants import *
import visualization
from odbAccess import *

odbfilename = 'C:\work\Cement polymerization\cup\cup_6-post.odb'

# Thermal damage parameters
A1=1.0e+98
A2=7.5e+04

# ------------------------------------------------------------------------------------------

# Open odb file
myOdb = openOdb(name=odbfilename, readOnly=FALSE)

# Get step
this_step=myOdb.steps['Step-1']

# Get frame data
results_frames=[]
time=[]
for i in this_step.frames:
  results_frames.append(i)
  time.append(i.frameValue)

# Get instance and bone element set
myInstance=myOdb.rootAssembly.instances['PART-1-1']
bone = myInstance.elementSets['PART-1-1_BONE']

# Get element labels of elements in elset bone
first_frame=this_step.frames[0].fieldOutputs['TEMP'].getSubset(region=bone).values
elm_labels=[]
for i in first_frame:
  if i.integrationPoint==1:
    elm_labels.append(i.elementLabel)

# Intitialise damage parameter to zero (at all points for all time steps)
damage=[]
length1=len(results_frames)
length2=len(first_frame)

for i in range(length1):
  damage.append([])
  for j in range(length2):
    damage[i].append([0.0])

# Calculate damage parameter by integrating damage function (itself a function of temperature) over time
# The trapezoidal rule is used for numerical integration
for i in range(1,length1):

  values1=this_step.frames[i-1].fieldOutputs['TEMP'].getSubset(region=bone).values
  values2=this_step.frames[i].fieldOutputs['TEMP'].getSubset(region=bone).values
  dt = time[i]-time[i-1]

  for j in range(length2):

    temperature1=values1[j].data
    temperature2=values2[j].data

    val1 = A1*exp(-A2/temperature1)
    val2 = A1*exp(-A2/temperature2)
    current_damage=0.5*(val1+val2)*dt

    damage[i][j][0]=damage[i-1][j][0]+current_damage

# ------------------------------------------------------------------------------------------

# Create fieldoutput of thermal damage parameter at each frame
for i in range(length1):
  thisframe=results_frames[i]
  try:
    myField = thisframe.FieldOutput(name='TD', description='Thermal damage', type=SCALAR)
  except:
    print 'Error creating field', i, 'This field Will not be created' 
  else:
    myField.addData(position=INTEGRATION_POINT, instance=myInstance, labels=elm_labels, data=damage[i])

# ------------------------------------------------------------------------------------------

# Save and close odb file
myOdb.save()
myOdb.close()

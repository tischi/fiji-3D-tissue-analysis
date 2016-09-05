from ij.io import OpenDialog
from ij.io import Opener
from ij.gui import GenericDialog
from ij.plugin import ZProjector, RGBStackMerge, SubstackMaker, Concatenator
from ij import IJ, ImagePlus, ImageStack, WindowManager
from ij.plugin import Duplicator
from ij.process import StackStatistics
from ij.plugin import ImageCalculator
from ij.measure import ResultsTable
from ij.plugin.frame import RoiManager
import os, os.path, re, sys
from jarray import array
from ij.process import ImageConverter
import math
from math import sqrt
from ij.macro import MacroRunner

from loci.plugins import BF
from loci.common import Region
from loci.plugins.in import ImporterOptions

from automic.table import TableModel			# this class stores the data for the table
from automic.table import ManualControlFrame 	#this class visualises TableModel via GUI
from java.io import File
from ah.utils import ROIManipulator


# import my analysis function collection
'''import os, sys, inspect
this_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
if this_folder not in sys.path:
  print this_folder
  sys.path.insert(0, this_folder)
import ct_analysis_functions as af
reload(af)'''

def close_all_image_windows():
  # forcefully closes all open images windows
  ids = WindowManager.getIDList();
  if (ids==None):
    return
  for i in ids:
     imp = WindowManager.getImage(i)
     if (imp!=None):
       win = imp.getWindow()
       if (win!=None):
         imp.changes = False # avoids the "save changes" dialog
         win.close()
         

def mean(x):
  mean = sum(x) / len(x)
  return mean
  
def sd(x):
  mean = sum(x) / len(x)
  differences = [xx - mean for xx in x]
  sq_differences = [xx**2 for xx in differences]
  sd = sqrt(sum(sq_differences)/len(x))
  return sd


def extractChannel(imp, nChannel, nFrame):
  """ Extract a stack for a specific color channel and time frame """
  stack = imp.getImageStack()
  ch = ImageStack(imp.width, imp.height)
  for i in range(1, imp.getNSlices() + 1):
    index = imp.getStackIndex(nChannel, i, nFrame)
    ch.addSlice(str(i), stack.getProcessor(index))
  return ImagePlus("Channel " + str(nChannel), ch)

def measureSumIntensity3D(imp):
  stats = StackStatistics(imp)
  return stats.mean * stats.pixelCount

def autoThreshold(imp, method):
  impout = imp.duplicate() 
  IJ.run(impout, "Auto Threshold", "method=" + method + " white stack use_stack_histogram");
  impout.setTitle("Auto Threshold")
  return impout

def threshold(_img, lower_threshold, upper_threshold):
  imp = Duplicator().run(_img)
  #imp.show(); time.sleep(0.2)
  #IJ.setThreshold(imp, mpar['lthr'], mpar['uthr'])
  IJ.setThreshold(imp, lower_threshold, upper_threshold)
  IJ.run(imp, "Convert to Mask", "stack")
  imp.setTitle("Threshold");
  #IJ.run(imp, "Divide...", "value=255 stack");
  #IJ.setMinAndMax(imp, 0, 1);
  return imp

def compute_overlap(tbModel, iDataSet, impA, iChannelA, iChannelB):
  imp_bw = ImageCalculator().run("AND create stack", impA[iChannelA-1], impA[iChannelB-1])
  overlap_AandB = measureSumIntensity3D(imp_bw)/255
  tbModel.setNumVal(overlap_AandB, iDataSet, "PBT_"+str(iChannelA)+"AND"+str(iChannelB))
  return tbModel

# Structure for batch analysis:
#  
# - main 
#  - parameters = get_analysis_parameters()
#  - folder = get_folder()
#  - table = init_results_table()
#  - get_data_info(folder, table) 
#  - batch_analyze(parameters, table) 
#    - for row in table 
#       - analyze(row, table, parameters)
#         - imp = load_imp(table, row)
#         - write results to table(row)
#         - write segmentation overlay images (use from table(row))


def analyze(iDataSet, tbModel, p, output_folder):

  #
  # LOAD FILES
  #

  filepath = tbModel.getFileAPth(iDataSet, "RAW", "IMG")
  filename = tbModel.getFileName(iDataSet, "RAW", "IMG") 
  print("Analyzing: "+filepath)
  IJ.run("Bio-Formats Importer", "open=["+filepath+"] color_mode=Default view=Hyperstack stack_order=XYCZT");
  imp = IJ.getImage()
  
  #
  # INIT
  #
  IJ.run("Options...", "iterations=1 count=1"); 

  #
  # SHOW DATA
  # 
  
  #imp.show()

  #
  # SCALING
  #
  
  #IJ.run(imp, "Scale...", "x="+str(p["scale"])+" y="+str(p["scale"])+" z=1.0 interpolation=Bilinear average process create"); 
  
  #
  # CONVERSION
  #
  
  #IJ.run(imp, "8-bit", "");
 
  #
  # CROPPING
  #
  
  #imp.setRoi(392,386,750,762);
  #IJ.run(imp, "Crop", "");

  
  #
  # BACKGROUND SUBTRACTION
  #
  
  # IJ.run(imp, "Subtract...", "value=32768 stack");

  #
  # REGION SEGMENTATION
  #
  # general issues:
  # - intensity decreases along z
  	
  impA = []
  sumIntensities = []
  sumIntensitiesBW = []

  nChannels = 3
  for iChannel in range(1,nChannels+1):
  	
    # gating
    lower_threshold_value = p["lower_th_ch"+str(iChannel)]
    upper_threshold_value = p["upper_th_ch"+str(iChannel)]
    print("gating channel "+str(iChannel)+" with values from "+str(lower_threshold_value)+" to "+str(upper_threshold_value))
    imp_c = extractChannel(imp, iChannel, 1)
    imp_bw = threshold(imp_c, lower_threshold_value, upper_threshold_value) 
    impA.append(imp_bw)  # store the binary mask (0, 255) for OR operations later

    # convert from 0,255 to 0,1 for the following calculations
    imp_01 = imp_bw.duplicate()
    IJ.run(imp_01, "Divide...", "value=255 stack"); 
    
    # set all pixels in the intensity image to zero that are not within the gate
    imp_c_gated = ImageCalculator().run("Multiply create stack", imp_c, imp_01)

    # now measure in imp_c_gated and imp_bw(0,1)
    tbModel.setNumVal(round(measureSumIntensity3D(imp_01),0), iDataSet, "PBT_"+str(iChannel))
    tbModel.setNumVal(round(measureSumIntensity3D(imp_c_gated),0), iDataSet, "SumIntensity_"+str(iChannel))
    tbModel.setNumVal(lower_threshold_value, iDataSet, "LOWER_TH_CH"+str(iChannel)) # also record the threshold used
    tbModel.setNumVal(upper_threshold_value, iDataSet, "UPPER_TH_CH"+str(iChannel)) # also record the threshold used


  # Compute overlaps of binary images
  tbModel = compute_overlap(tbModel, iDataSet, impA, 1, 2)
  tbModel = compute_overlap(tbModel, iDataSet, impA, 2, 3) 
  tbModel = compute_overlap(tbModel, iDataSet, impA, 1, 3)

  # Compute total volume, i.e. pixels that are above threshold in any of the channels
  imp_combined = ImageCalculator().run("OR create stack", impA[0], impA[1])
  imp_combined = ImageCalculator().run("OR create stack", imp_combined, impA[2])
  tbModel.setNumVal(round(measureSumIntensity3D(imp_combined)/255,0), iDataSet, "PBT_1OR2OR3")

  #impbw = ImageCalculator().run("OR create stack", impbw, impA[2])


#
# ANALYZE INPUT FILES
#
def determine_input_files(foldername, tbModel):

  print("Determine input files in:",foldername)
  PBTtern = re.compile('(.*).tif') 
  
  #PBTtern = re.compile('(.*)--beats.tif') 
   
  i = 0
  for root, directories, filenames in os.walk(foldername):
	for filename in filenames:
	   print("Checking:", filename)
	   if (filename == "Thumbs.db") or filename.startswith('.'):
	     continue
	   match = re.search(PBTtern, filename)
	   if (match == None) or (match.group(1) == None):
	     continue
	   tbModel.addRow()
	   tbModel.setFileAPth(foldername, filename, i, "RAW","IMG")
	   print("Accepted:", filename)
	   
	   i += 1
    
  return(tbModel)

#
# GET PARAMETERS
#
def get_parameters(p, keys, num_data_sets):
  gd = GenericDialog("Please enter parameters")

  gd.addMessage("Found "+str(num_data_sets)+" data sets")
  gd.addStringField("analyse", "all");

  gd.addMessage("Image analysis parameters:")
  gd.addMessage("Please note: the threshold values are inclusive!\nThus, to exlude pixels with value 255 the upper threshold needs to be 254")
  
  
  for k in keys:
    gd.addNumericField(k, p[k], 2);
  gd.showDialog()
  if gd.wasCanceled():
    return

  to_be_analyzed = gd.getNextString()
  for k in keys:
    p[k] = gd.getNextNumber()
    
  return to_be_analyzed, p

    
if __name__ == '__main__':

  print("")
  print("#")
  print("# Tissue analysis")
  print("#")
  print("")
  
  #
  # GET INPUT FOLDER
  #
  od = OpenDialog("Click on one of the image files in the folder to be analysed", None)
  input_folder = od.getDirectory()
  if input_folder is None:
    sys.exit("No folder selected!")
    

  #
  # MAKE OUTPUT FOLDER
  #
  output_folder = input_folder[:-1]+"--fiji"
  if not os.path.isdir(output_folder):
    os.mkdir(output_folder)

  #
  # DETERMINE INPUT FILES
  #
  tbModel = TableModel(input_folder)
  tbModel.addFileColumns('RAW','IMG')
  tbModel = determine_input_files(input_folder, tbModel)
  
  #
  # CHECK FIRST IMAGE
  #
  filepath = tbModel.getFileAPth(0, "RAW", "IMG")
  print("Analyzing: "+filepath)
  IJ.run("Bio-Formats Importer", "open=["+filepath+"] color_mode=Default view=Hyperstack stack_order=XYCZT");
  imp = IJ.getImage()
  bit_depth = imp.getBitDepth()
  imp.close()
  
  
  #
  # GET PARAMETERS
  #
  p = dict()
  keys = list() # necessary to have the correct order in the gui
  nChannels = 3
  for i in range(1,nChannels+1):
    p["lower_th_ch"+str(i)] = 30; keys.append("lower_th_ch"+str(i))
    if bit_depth==8:
      p["upper_th_ch"+str(i)] = 254; keys.append("upper_th_ch"+str(i))
    elif bit_depth==16:
      p["upper_th_ch"+str(i)] = 65534; keys.append("upper_th_ch"+str(i))
    else:
      IJ.log("ERROR: unknown bit depth!")
      sys.exit(0)
    
  to_be_analyzed, p = get_parameters(p, keys, tbModel.getRowCount())

  #
  # POPULATE AND SHOW INTERACTIVE TABLE
  #

  # thresholds
  for i in range(1, nChannels+1):
    tbModel.addValColumn("LOWER_TH_CH"+str(i), "NUM")
    tbModel.addValColumn("UPPER_TH_CH"+str(i), "NUM")
   
  # pixels above threshold (PBT)
  for i in range(1, nChannels+1):
    tbModel.addValColumn("PBT_"+str(i), "NUM")

  # total intensity in whole stack (no masking and no bg-subtraction)
  for i in range(1, nChannels+1):
    tbModel.addValColumn("SumIntensity_"+str(i), "NUM")

  # overlapping PBTs in different channels  
  tbModel.addValColumn("PBT_1AND3", "NUM")
  tbModel.addValColumn("PBT_1AND2", "NUM")
  tbModel.addValColumn("PBT_2AND3", "NUM")
   
  tbModel.addValColumn("PBT_1OR2OR3", "NUM")

  
  '''
  for i in range(p["n_rois"]):
    tbModel.addFileColumns('FFT_R'+str(i+1),'IMG')
  '''
  
  frame=ManualControlFrame(tbModel)
  frame.setVisible(True)
  
  #
  # ANALYZE
  #

  if not to_be_analyzed=="all":
    close_all_image_windows()
    analyze(int(to_be_analyzed)-1, tbModel, p, output_folder)
  else:
    for i in range(tbModel.getRowCount()):
      close_all_image_windows()
      analyze(i, tbModel, p, output_folder)

  close_all_image_windows()
from ij.io import OpenDialog, DirectoryChooser
from ij.gui import GenericDialog
from ij.io import Opener
from ij.plugin import ZProjector, RGBStackMerge
from ij import IJ, ImagePlus, ImageStack
from ij.plugin import Duplicator
from ij.process import StackStatistics
from ij.plugin import ImageCalculator
import os
from jarray import array

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

def threshold(_img, threshold):
  imp = Duplicator().run(_img)
  #imp.show(); time.sleep(0.2)
  #IJ.setThreshold(imp, mpar['lthr'], mpar['uthr'])
  IJ.setThreshold(imp, threshold, 1000000000)
  IJ.run(imp, "Convert to Mask", "stack")
  imp.setTitle("Threshold");
  #IJ.run(imp, "Divide...", "value=255 stack");
  #IJ.setMinAndMax(imp, 0, 1);
  return imp

def analyse(imp, root, filename, thresholds):
  imp.show()
  print ""
  print(imp.getTitle())
  print "**********************************"
  
  channels=(1,2,3)
  backgrounds=thresholds #(175,0,275)

  #IJ.run(imp, "Gaussian Blur...", "sigma=1 stack");
  
  impA = []
  sumIntensities = []
  sumIntensitiesBW = []
  
  for iChannel in channels:

  	# general issues:
  	# - intensity decreases along z
  	
    # measure total intensity
    # - bg subtraction?
    impc = extractChannel(imp, iChannel, 1)
    impbw = threshold(impc, thresholds[iChannel-1]) 
    
    impcBGcorr = Duplicator().run(impc)
    IJ.run(impcBGcorr, "Subtract...", "value="+str(backgrounds[iChannel-1])+" stack");
    sumIntensities.append(measureSumIntensity3D(impcBGcorr))

    # measure number of pixels containing a signification signal
    # - what is good here for thresholding and preprocessing?
    #impbw = autoThreshold(impc, "Default")
    
    #impbw.setTitle("bw"+str(iChannel))
    #impbw.show()
    impA.append(impbw)
    sumIntensitiesBW.append(measureSumIntensity3D(impbw)/255)
      
  # measure how many pixels are above threshold in either of the images (2015-09-21)
  impbw = ImageCalculator().run("OR create stack", impA[0], impA[1])
  impbw = ImageCalculator().run("OR create stack", impbw, impA[2])
  impA.append(impbw)
  sumIntensitiesBW.append(measureSumIntensity3D(impbw)/255)  # volume of pixels that are above threshold in either of the images
  sumIntensities.append(0.0)
 
  print "intensities,",str(sumIntensities)[1:-1]
  print "volumes,",str(sumIntensitiesBW)[1:-1]

  # normalise to the volume of the 2nd channel (which is the AF channel) (removed on 2015-09-21)
  #sumIntensitiesNorm = [x / sumIntensitiesBW[1] for x in sumIntensities]
  #sumIntensitiesBWNorm = [x / sumIntensitiesBW[1] for x in sumIntensitiesBW]

  # normalise to the volume of the 4th channel (which is the OR channel) (added on 2015-09-21)
  sumIntensitiesNorm = [x / sumIntensitiesBW[3] for x in sumIntensities]
  sumIntensitiesBWNorm = [x / sumIntensitiesBW[3] for x in sumIntensitiesBW]
  
  # print all the results without the brackets (that's why the 1:-1 and not :)
  print "intensities normalised,",str(sumIntensitiesNorm)[1:-1] 
  print "volumes normalised,",str(sumIntensitiesBWNorm)[1:-1] 
  print "thresholds,",str(thresholds)[1:-1]
  
  impMerge = RGBStackMerge.mergeChannels(array(impA, ImagePlus), False)
  impMerge.setDisplayMode(IJ.COLOR);
  impMerge.setC(1);
  IJ.run(impMerge, "Green", "");
  impMerge.setC(2);
  IJ.run(impMerge, "Grays", "");
  impMerge.setC(3);
  IJ.run(impMerge, "Blue", "");
  impMerge.setC(4);
  IJ.run(impMerge, "Red", "");
  impMerge.show()
  IJ.saveAs(impMerge, "Tiff", os.path.join(root, filename + "--segmented.tif"));
  #imp = IJ.getImage();
  #imp.setDisplayMode(IJ.COLOR);

# obtains the files and sends them to the analyse function
# todo: give the name of the analysis function as an argument!
def getFilesAndAnalyze(folder, thresholds):
  for root, directories, filenames in os.walk(folder):
    print ""
    print "directories:",directories
    for filename in filenames:
      if ( ("segmented.tif" == filename.split("--")[-1]) or (filename=="Thumbs.db") ):
        print "\nwill not analyze file:",filename
        continue
      else:
        #print root
        print "\nanalyzing file:",filename
        imp = Opener().openImage(root, filename)
        analyse(imp, root, filename, thresholds)

          
if __name__ == '__main__':
  srcDir = DirectoryChooser("Choose directory").getDirectory()
  #od = OpenDialog("Click on one of the image files in the folder to be analysed", None	)
  #filename = od.getFileName()
  #foldername = od.getDirectory()
  #foldername = "Z:/Tischi/projects/CellZome/2015-04-27--3D_Tissue/data"
  gd = GenericDialog("3D Tissue analysis");
  gd.addNumericField("Threshold Ch1: ", 30, 0);
  gd.addNumericField("Threshold Ch2: ", 30, 0);
  gd.addNumericField("Threshold Ch3: ", 30, 0);
  gd.showDialog();
  #if (gd.wasCanceled()):
  #  return
  thresholds = []
  thresholds.append(gd.getNextNumber())
  thresholds.append(gd.getNextNumber())
  thresholds.append(gd.getNextNumber())
  #foldername = "Z:/Tischi/projects/CellZome/2015-04-29--3D_Tissue/subset"
  
  if None != srcDir:
    print "\n\n\n\nfoldername = "+srcDir
    getFilesAndAnalyze(srcDir, thresholds)
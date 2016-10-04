# fiji-3D-tissue-analysis

## installation and running

- install Fiji: http://fiji.sc/#download
- dowload and unzip: https://github.com/tischi/fiji-3D-tissue-analysis/archive/master.zip
- copy __AutoMic_JavaTools-1.0.0-08032016.jar__ to your Fiji plugins folder 
- copy __3D-tissue-analysis.py__ into an arbiratry folder

Drag&drop __3D-tissue-analysis.py__ onto Fiji. The script editor window will open. Click [Run] at the bottom of this window.

## input data

A folder with tif files. Each tif file must be an ImageJ Composite image (also called Hyperstack in ImageJ). Currently, each file must contain exactly 3 channels, but can contain an arbitrary number of z-slices. If the folder contains multiple such files they will all be analyze automatically.

## input parameters

The user must specify an lower and an upper threshold intensity for each channel. The same threshold intensities will be used for all images in the folder.

## computation

...

## output

There will be an interactive table at the bottom of your screen with the analysis results. __PBT__ means "Pixels Between Thresholds" in the respective channel. __AND__ means that the pixels were selected in both channels. __OR__ means that the pixels were only selected in either of the channels. 

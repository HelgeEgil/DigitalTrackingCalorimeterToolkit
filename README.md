# Digital Tracking Calorimeter Toolkit
Toolkit for the Monte Carlo and data analysis for the Digital Tracking Calorimeter project.
For more details see the publication <i>Pettersen, H. E. S. et al., Proton tracking in a high-granularity Digital Tracking Calorimeter for proton CT purposes. Nuclear Instruments and Methods in Physics Research A 860 (C): 51 - 61 (2017)</i>
https://doi.org/10.1016/j.nima.2017.02.007 and the PhD thesis <i>A Digital Tracking Calorimeter for Proton Computed Tomography: University of Bergen 2018. http://bora.uib.no/handle/1956/17757</i>

## Overview
It is developed for the management, pre-processing, reconstruction, analysis and presentation of data from the beam test and from MC simulations.
It is of a modular and object-oriented design, such that it should be simple to extend
the software with the following purposes in mind: to include analysis data from multiple
sources such as the next-generation readout electronics and different MC simulation software;
to give the user a broad selection of different geometrical, physical and reconstruction
models; to include a more extensive proton CT simulation with complex phantoms and positional
trackers before and after the phantom or patient; and to facilitate further development
and usage of the software by making it available as a analysis library.

The software is written in C++
and ROOT, with some auxiliary tools written in Python and user scripts written in bash.
Several hands-on teaching workshops have been held at the University of Bergen in order to
demonstrate the usage of the software, and summaries of these are available at the project
documentation website https://wiki.uib.no/pct/. In total the framework with tools consists
of approximately 20 000 lines of code.

## Functionality:
* Generating GATE simulation files
* Perfoming GATE simulations
* Making range-energy tables, finding the straggling, etc.
* Tracking analysis: This can be done both simplified and full
* Simplified: No double-modelling of the pixel diffusion process (use MC provded energy loss), no track reconstruction (use eventID tag to connect tracks from same primary).
* The 3D reconstruction of phantoms using tracker planes has not yet been implemented
* Range estimation

The analysis toolchain has the following components:

![alt text](https://raw.githubusercontent.com/HelgeEgil/focal/master/projects/Readme/Analysis_chain.PNG "Logo Title Text 1")

The full tracking workflow is implemented in the function <code>DTCToolkit/HelperFunctions/getTracks.C::getTracks()</code>, and the tracking and range estimation workflow is found in <code>DTCToolkit/Analysis/Analysis.C::drawBraggPeakGraphFit()</code>.


## GATE simulations 
### Geometry scheme
The simplified simulation geometry for the future DTC simulations has been proposed as:

![alt text](https://raw.githubusercontent.com/HelgeEgil/focal/master/projects/Readme/Layer_schematics.PNG "Logo Title Text 1")

It is partly based on the ALPIDE design, and the FoCal design. The GATE geometry corresponding to this scheme is based on the following hierarchy:
```
World -> Scanner1 -> Layer -> Module + Absorber + Air gap
                              Module = Active sensor + Passive sensor + Glue + PCB + Glue
      -> Scanner2 -> [Layer] * Number Of Layers
```

The idea is that Scanner1 represents the first layer (where e.g. there is no absorber, only air), and that Scanner2 represents all the following (similar) layers which are repeated.

#### Generating the macro files
To generate the geometry files to run in Gate, a Python script is supplied.
It is located within the ''gate/python'' subfolder.
```
[gate/python] $ python gate/python/makeGeometryDTC.py
```

![alt text](https://raw.githubusercontent.com/HelgeEgil/focal/master/projects/Readme/GATE_geometry_builder.PNG "Logo Title Text 1")

Choose the wanted characteristics of the detector, and use ''write files'' in order to create the geometry file Module.mac, which is automatically included in Main.mac.
Note that the option "Use water degrader phantom" should be checked (as is the default behavior)!

### Creating the full simulations files for a range-energy look-up-table
In this step, 5000-10000 particles are usually sufficient in order to get accurate results.
To loop through different energy degrader thicknesses, run the script ''runDegraderFull.sh'':
```
[gate/python] $ ./runDegraderFull.sh <absorber thickness> <degraderthickness from> <degraderthickness stepsize> <degraderthickness to>
```
The brackets indicate the folder in the Github repository to run the code from. Please note that the program should not be executed using the <code>sh</code> command, as this refers do different shells in different Linux distribtions, and not all shells support the conditional bash expressions used in the script.


For example, with a 3 mm degrader, and simulating a 250 MeV beam passing through a phantom of 50, 55, 60, 65 and 70 mm water:
```
[gate/python] $ ./runDegraderFull.sh 3 50 5 70
```
Please note that there is a variable NCORES in this script, which ensures that NCORES versions of the Gate executable are run in parallel, and then waits for the last background process to complete before a new set of NCORES executables are run. So if you set NCORES=8, and run <code>sh runDegraderFull.sh 3 50 1 70</code>, first 50-57 will run in parallel, and when they're done, 58-65 will start, etc. The default value is NCORES=4.

### Creating the chip-readout simulations files for resolution calculation
In this step a higher number of particles is desired. I usually use 25000 since we need O(100) simulations. A sub 1-mm step size will really tell us if we manage to detect such small changes in a beam energy.

And loop through the different absorber thicknesses:
```
[gate/python] $ ./runDegrader.sh <absorber thickness> <degraderthickness from> <degraderthickness stepsize> <degraderthickness to>
```
The same parallel-in-sequential run mode has been configured here.

### Creating the basis for range-energy calculations
#### The range-energy look-up-table
Now we have ROOT output files from Gate, all degraded differently through a varying water phantom and therefore stopping at different places in the DTC.
We want to follow all the tracks to see where they end, and make a histogram over their stopping positions. This is of course performed from a looped script, but to give a small recipe:
* Retrieve the first interaction of the first particle. Note its event ID (history number) and edep (energy loss for that particular interaction)
* Repeat until the particle is outside the phantom. This can be found from the volume ID or the z position (the first interaction with {math|z>0}). Sum all the found edep values, and this is the energy loss inside the phantom. Now we have the "initial" energy of the proton before it hits the DTC
* Follow the particle, noting its z position. When the event ID changes, the next particle is followed, and save the last z position of where the proton stopped in a histogram
* Do a Gaussian fit of the histogram after all the particles have been followed. The mean value is the range of the beam with that particular "initial" energy. The spread is the range straggling. Note that the range straggling is more or less constant, but the contributions to the range straggling from the phantom and DTC, respectively, are varying linearly. 

This recipe has been implemented in <code>DTCToolkit/Scripts/findRange.C</code>. Test run the code on a few of the cases (smallest and biggest phantom size ++) to see that
* The correct start- and end points of the histogram looks sane. If not, this can be corrected for by looking how <code>xfrom</code> and <code>xto</code> is calculated and playing with the calculation.
* The mean value and straggling is calculated correctly
* The energy loss is calculated correctly
You can run <code>findRange.C</code> in root by compiling and giving it three arguments; Energy of the protons, absorber thickness, and the degrader thickness you wish to inspect. 

```
[DTCToolkit/Scripts] $ root 
ROOT [1] .L findRange.C+
// void findRange(Int_t energy, Int_t absorberThickness, Int_t degraderThickness)
ROOT [2] findRange f(250, 3, 50); f.Run();
```

The output should look like this: Correctly places Gaussian fits is a good sign.

![alt text](https://raw.githubusercontent.com/HelgeEgil/focal/master/projects/Readme/FindRanges.JPG "Logo Title Text 1")

If you're happy with this, then a new script will run <code>findRange.C</code> on all the different ROOT files generated earlier.
```
[DTCToolkit/Scripts] $ root 
ROOT [1] .L findManyRangesDegrader.C
// void findManyRanges(Int_t degraderFrom, Int_t degraderIncrement, Int_t degraderTo, Int_t absorberThicknessMm)
ROOT [2] findManyRanges(50, 5, 70, 3)
```

This is a serial process, so don't worry about your CPU.
The output is stored in <code>DTCToolkit/Output/findManyRangesDegrader.csv</code>.
It is a good idea to look through this file, to check that the values are not very jumpy (Gaussian fits gone wrong).

We need the initial energy and range in ascending order. The findManyRangesDegrader.csv files contains more rows such as initial energy straggling and range straggling for other calcualations. This is a bit tricky, but do (assuming a 3 mm absorber geometry):
```
[DTCToolkit] $ cat OutputFiles/findManyRangesDegrader.csv | awk '{print ($6 " " $3)}' | sort -n > Data/Ranges/3mm_Al.csv
```

NB: If there are many different absorber geometries in findManyRangesDegrader, either copy the interesting ones or use <code>| grep " X " |</code> to only keep X mm geometry

When this is performed, the range-energy table for that particular geometry has been created, and is ready to use in the analysis. Note that since the calculation is based on cubic spline interpolations, it cannot extrapolate -- so have a larger span in the full Monte Carlo simulation data than with the chip readout. For more information about that process, see this document: 

![alt text](https://raw.githubusercontent.com/HelgeEgil/focal/master/projects/Readme/HowCalculateResolution.PNG "Logo Title Text 1")

### Range straggling parameterization and <math>R_0 = \alpha E^p</math>
It is important to know the amount of range straggling in the detector, and the amount of energy straggling after the degrader. In addition, to calculate the parameters <math>\alpha, p</math> from the somewhat inaccurate Bragg-Kleeman equation <math>R_0 = \alpha E ^ p</math>, in order to correctly model the "depth-dose curve" <math>dE / dz = p^{-1} \alpha^{-1/p} (R_0 - z)^{1/p-1}</math>. This happens automatically in <code>DTCToolkit/GlobalConstants/MaterialConstants.C</code>, so there is no need to do this by hand. This step should be followed anyhow, since it is a check of the data produced in the last step: Outliers should be removed or "fixed", eg. by manually fitting the conflicting datapoints using <code>findRange.C</code>.

To find all this, run the script <code>DTCToolkit/Scripts/findAPAndStraggling.C</code>. This script will loop through all available data lines in the <code>DTCToolkit/OutputFiles/findManyRangesDegrader.csv</code> file that has the correct absorber thickness, so you need to clean the file first (or just delete it before running <code>findManyRangesDegrader.C</code>).

```
[DTCToolkit/Scripts] $ root
ROOT [0] .L findAPAndStraggling.C+
// void findAPAndStraggling(int absorberthickness)
ROOT [1] findAPAndStraggling(3)
```

The output from this function should be something like this:

![alt text](https://raw.githubusercontent.com/HelgeEgil/focal/master/projects/Readme/FindAPAndStraggling.JPG "Logo Title Text 1")

### Configuring the DTC Toolkit to run with correct geometry
The values from <code>findManyRanges.C</code> should already be in <code>DTCToolkit/Data/Ranges/3mm_Al.csv</code> (or the corresponding material / thickness).

Look in the file <code>DTCToolkit/GlobalConstants/Constants.h</code> and check that the correct absorber thickness values etc. are set:
```
...
39 Bool_t useDegrader = true;
...
52 const Float_t kAbsorberThickness = 3;
...
59 Int_t kEventsPerRun = 100000;
...
66 const Int_t kMaterial = kAluminum;
```

Since we don't use tracking but only MC truth in the optimization, the number kEventsPerRun (<math>n_p</math> in the NIMA article) should be higher than the number of primaries per energy.
If tracking is to be performed anyhow, turn on the line
```
const kDoTracking = true;
```

## Running the DTC Toolkit 

The following section will detail how to perform these separate steps. A quick review of the classes available:
* <code>Hit</code>: A (int x,int y,int layer, float edep) object from a pixel hit. edep information only from MC
* <code>Hits</code>: A <code>TClonesArray</code> collection of Hit objects
* <code>Cluster</code>: A (float x, float y, int layer, float clustersize) object from a cluster of <code>Hit</code>s The (x,y) position is the mean position of all involved hits.
* <code>Clusters</code>: A <code>TClonesArray</code> collection of <code>Cluster</code> objects.
* <code>Track</code>: A <code>TClonesArray</code> collection of <code>Cluster</code> objects... But only one per layer, and is connected through a physical proton track. Many helpful member functions to calculate track properties.
* <code>Tracks</code>: A <code>TClonesArray</code> collection of <code>Track</code> objects.
* <code>Layer</code>: The contents of a single detector layer. Is stored as a <code>TH2F</code> histogram, and has a <code>Layer::findHits</code> function to find hits, as well as the cluster diffusion model <code>Layer::diffuseLayer</code>. It is controlled from a <code>CalorimeterFrame</code> object.
* <code>CalorimeterFrame</code>: The collection of all <code>Layer</code>s in the detector.
* <code>DataInterface</code>: The class to talk to DTC data, either through semi-<code>Hit</code> objects as retrieved from Utrecht from the Groningen beam test, or from ROOT files as generated in Gate.

'''Important''': To load all the required files / your own code, include your C++ sources files in the <code>DTCToolkit/Load.C</code> file, after Analysis.C has loaded:
```
...
gROOT->LoadMacro("Analysis/Analysis.C+");
gROOT->LoadMacro("Analysis/YourFile.C+"); // Remember to add a + to compile your code
}
```

### Data readout: MC, MC + truth, experimental
In the class <code>DataInterface</code> there are several functions to read data in ROOT format.

```
int   getMCFrame(int runNumber, CalorimeterFrame *calorimeterFrameToFill, [..]) <- MC to 2D hit histograms
void  getMCClusters(int runNumber, Clusters *clustersToFill); <-- MC directly to clusters w/edep and eventID
void  getDataFrame(int runNumber, CalorimeterFrame *calorimeterFrameToFill, int energy); <- experimental data to 2D hit histograms
```

To e.g. obtain the experimental data, use
```
DataInterface *di = new DataInterface();
CalorimeterFrame *cf = new CalorimeterFrame();
  
for (int i=0; i<numberOfRuns; i++) { // One run is "readout + track reconstruction
   di->getDataFrame(i, cf, energy);
   // From here the object cf will contain one 2D hit histogram for each of the layers
   // The number of events to readout in one run: kEventsPerRun (in GlobalConstants/Constants.h)
}
```

Examples of the usage of these functions are located in <code>DTCToolkit/HelperFunctions/getTracks.C</code>.
Please note the phenomenological difference between experimental data and MC:
* Exp. data has some noise, represented as "hot" pixels and 1-pixel clusters
* Exp. data has diffused, spread-out, clusters from physics processes
* Monte Carlo data has no such noise, and proton hits are represented as 1-pixel clusters (with edep information)

### Pixel diffusion modelling (MC only)
To model the pixel diffusion process, i.e. the the diffusion of the electron-hole pair charges generated from the proton track towards nearby pixels, an empirical model has been implemented. It is described in the NIMA article [[http://dx.doi.org/10.1016/j.nima.2017.02.007]], and also in the source code in  <code>DTCToolkit/Classes/Layer/Layer.C::diffuseLayer</code>.

To perform this operation on a filled <code>CalorimeterFrame *cf</code>, use
```
TRandom3 *gRandom = new TRandom3(0); // use #import <TRandom3.h>
cf->diffuseFrame(gRandom);
```

#### Inverse pixel diffusion calculation (MC and exp. data)
This process has been inversed in a Python script, and performed with a large number of input cluster sizes. The result is a parameterization between the proton's energy loss in a layer, and the number of activated pixels:

![alt text](https://raw.githubusercontent.com/HelgeEgil/focal/master/projects/Readme/edep-and-cs.png "Logo Title Text 1")

The function <code>DTCToolkit/HelperFunctions/Tools.C::getEdepFromCS(n)</code> contains the parameterization:
```
Float_t getEdepFromCS(Int_t cs) {
   return -3.92 + 3.9 * cs - 0.0149 * pow(cs,2) + 0.00122 * pow(cs,3) - 1.4998e-5 * pow(cs,4);
}
```

### Cluster identification 
Cluster identification is the process to find all connected hits (activated pixels) from a single proton in a single layer. It can be done by several algorithms, simple looped neighboring, DBSCAN, ...
The process is such:
* All hits are found from the diffused 2D histograms and stored as <code>Hit</code> objects with <math>(x,y,layer)</math> in a TClonesArray list.
* This list is indexed by layer number (a new list with the index the first Hit in each layer) to optimize any search
* The cluster finding algorithm is applied. For every Hit, the Hit list is looped through to find any connected hits. The search is optimized by use of another index list on the vertical position of the Hits. All connected hits (vertical, horizontal and diagonal) are collected in a single Cluster object with <math>(x,y,layer,cluster size)</math>, where the cluster size is the number of its connected pixels.

This task is simply performed on a diffused <code>CalorimeterFrame *cf</code>:
```
Hits *hits = cf->findHits();
Clusters *clusters = hits->findClustersFromHits();
```

### Proton track reconstruction 
The process of track reconstruction is described fully in [[http://dx.doi.org/10.1016/j.nima.2017.02.007]].

From a collection of cluster objects, <code>Clusters * clusters</code>, use the following code to get a collection of the Track objects connecting them across the layers.
```
Tracks * tracks = clusters->findCalorimeterTracks();
```

Some optimization schemes can be applied to the tracks in order to increase their accuracy:
```
tracks->extrapolateToLayer0(); // If a track was found starting from the second layer, we want to know the extrapolated vector in the first layer
tracks->splitSharedClusters(); // If two tracks meet at the same position in a layer, and they share a single cluster, split the cluster into two and give each part to each of the tracks
tracks->removeTracksLeavingDetector(); // If a track exits laterally from the detector before coming to a stop, remove it
tracks->removeTracksEndingInBadChannnels(); // ONLY EXP DATA: Use a mask containing all the bad chips to see if a track ends in there. Remove it if it does.
```

### Putting it all together so far 
It is not easy to track a large number of proton histories simultaneously, so one may want to loop this analysis, appending the result (the tracks) to a larger Tracks list. This can be done with the code below:
```
DataInterface *di = new DataInterface();
CalorimeterFrame *cf = new CalorimeterFrame();
Tracks * allTracks = new Tracks();
 
for (int i=0; i<numberOfRuns; i++) { // One run is "readout + track reconstruction
   di->getDataFrame(i, cf, energy);
   TRandom3 *gRandom = new TRandom3(0); // use #import <TRandom3.h>
   cf->diffuseFrame(gRandom);
   Hits *hits = cf->findHits();
   Clusters *clusters = hits->findClustersFromHits();
   Tracks * tracks = clusters->findCalorimeterTracks();
   tracks->extrapolateToLayer0();
   tracks->splitSharedClusters();
   tracks->removeTracksLeavingDetector();
   tracks->removeTracksEndingInBadChannnels();
    
   for (int j=0; j<tracks->GetEntriesFast(); j++) {
      if (!tracks->At(j)) continue;
      allTracks->appendTrack(tracks->At(j));
   }

   delete tracks;
   delete hits;
   delete clusters;
}
```

### Individual tracks: Energy loss fitting 
To obtain the most likely residual range / stopping range from a Track object, use
```
track->doRangeFit();
float residualRange = track->getFitParameterRange();
```

What happens here is that a TGraph with the ranges and in-layer energy losses of all the Cluster objects is constructed. A differentiated Bragg Curve is fitted to this TGraph:

<math> f(z) = p^{-1} \alpha^{-1/p} (R_0 - z)^{1/p-1} </math>

With <math>p,\alpha</math> being the parameters found during the full-scoring MC simulations. The value <math>R_0</math>, or <code>track::getFitParameterRange</code> is stored.

![alt text](https://raw.githubusercontent.com/HelgeEgil/focal/master/projects/Readme/EnergyLossFit.PNG "Logo Title Text 1")

### (3D reconstruction / MLP estimation) 
When the volume reconstruction is implemented, it is to be put here:
* Calculate the residual range and incoming vectors of all protons
* Find the Most Likely Path (MLP) of each proton
* Divide the proton's average energy loss along the MLP
* Then, with a measure of a number of energy loss values in each voxel, perform some kind of average scheme to find the best value.

Instead, we now treat the complete detector as a single unit / voxel, and find the best SUM of all energy loss values (translated into range). The average scheme used in this case is described below, however this might be different than the best one for the above case.

# Residual range calculation 
To calculate the most likely residual range from a collection of individual residual ranges is not a simple task!
It depends on the average scheme, the distance between the layers, the range straggling etc. Different solutions have been attempted:
* In cases where the distance between the layers is large compared to the straggling, a histogram bin sum based on the depth of the first layer identified as containing a certain number of proton track endpoints is used. It is the method detailed in the NIMA article [[http://dx.doi.org/10.1016/j.nima.2017.02.007]], and it is implemented in <code>DTCToolkit/Analysis/Analysis.C::doNGaussianFit(*histogram, *means, *sigmas)</code>.
* In cases where the distance between the layers is small compared to the straggling, a single Gaussian function is fitted on top of all the proton track endpoints, and the histogram bin sum average value is calculated from minus 4 sigma to plus 4 sigma. This code is located in <code>DTCToolkit/Analysis/Analysis.C::doSimpleGaussianFit(*histogram, *means, *sigmas)</code>. This is the version used for the geometry optimization project.

With a histogram <code>hRanges</code> containing all the different proton track end points, use
```
float means[10] = {};
float sigmas[10] = {};
TF1 *gaussFit = doSimpleGaussianFit(hRanges, means, sigmas);
printf("The resulting range of the proton beam if %.2f +- %.2f mm.\n", means[9], sigmas[9]);
```

![alt text](https://raw.githubusercontent.com/HelgeEgil/focal/master/projects/Readme/ResidualRangeHistogram.JPG "Logo Title Text 1")

## Geometry optimization: How does the DTC Toolkit calculate resolution? 
The resolution in this case is defined as the width of the final range histogram for all protons.
The goal is to match the range straggling which manifests itself in the Gaussian distribution of the range of all protons in the DTC, from the full Monte Carlo simulations:

To characterize the resolution, a realistic analysis is performed. Instead of scoring the complete detector volume, including the massive energy absorbers, only the sensor chips placed at intervals (<math>\Delta z = 0.375\ \textrm{mm} + d_{\textrm{absorber}}</math>) are scored. Tracks are compiled by using the eventID tag from GATE, so that the track reconstruction efficiency is 100%. Each track is then put in a depth / edep graph, and a Bragg curve is fitted on the data:

The distribution of all fitted ranges (simple to calculate from fitted energy) should match the distribution above - with a perfect system. All degradations during analysis, sampling error, sparse sampling, mis-fitting etc. will ensure that the peak is broadened.

![alt text](https://raw.githubusercontent.com/HelgeEgil/focal/master/projects/Readme/HowCalculateResolution.PNG "Logo Title Text 1")

### Finding the resolution 
To find this resolution, or degradation in the straggling width, for a single energy, run the DTC toolkit analysis.
```
[DTCToolkit] $ root Load.C
// drawBraggPeakGraphFit(Int_t Runs, Int_t dataType = kMC, Bool_t recreate = 0, Float_t energy = 188, Float_t degraderThickness = 0)
ROOT [0] drawBraggPeakGraphFit(1, 0, 1, 250, 34)
```
This is a serial process, so don't worry about your CPU when analysing all ROOT files in one go.
With the result

![alt text](https://raw.githubusercontent.com/HelgeEgil/focal/master/projects/Readme/Distribution_after_analysis2.JPG "Logo Title Text 1")

The following parameters are then stored in <code>DTCToolkit/OutputFiles/results_makebraggpeakfit.csv</code>:



| Absorber thickness | Degrader thickness | Nominal WEPL range | Calculated WEPL range | Nominal WEPL straggling | Calculated WEPL straggling |
| ---------------    | -----------------  | ------------------ | --------------------- | --------------------    | -------------------------  | 
| 3 (mm)             | 34 (mm)            | 345 (mm WEPL)      | 345.382 (mm WEPL)     | 2.9 (mm WEPL)           | 6.78 (mm WEPL)             |


To perform the analysis on all different degrader thicknesses, use the script <code>DTCToolkit/makeFitResultPlotsDegrader.sh</code> (arguments: degrader from, degrader step and degrader to):
 ```
 [DTCToolkit] $ sh makeFitResultsPlotsDegrader.sh 1 1 380
 ```
This may take a few minutes...
When it's finished, it's important to look through the file results_makebraggpeakfit.csv to identify all problem energies, as this is a more complicated analysis than the range finder above.
If any is identified, run the drawBraggPeakGraphFit at that specific degrader thickness to see where the problems are.

### Displaying the results
If there are no problems, use the script <code>DTCToolkit/Scripts/makePlots.C</code> to plot the contents of the file <code>DTCToolkit/OutputFiles/results_makebraggpeakfit.csv</code>:
   [DTCToolkit/Scripts/optimization] $ root plotRangesAndStraggling.C
The output is a map of the accuracy of the range determination, and a comparison between the range resolution (#sigma of the range determination) and its lower limit, the range straggling.

![alt text](https://raw.githubusercontent.com/HelgeEgil/focal/master/projects/Readme/MakePlots_accuracy.JPG "Logo Title Text 1")

![alt text](https://raw.githubusercontent.com/HelgeEgil/focal/master/projects/Readme/MakePlots_resolution.JPG "Logo Title Text 1")

### A review of the different modules in the code 
To clone the project, run
```
git clone https://github.com/HelgeEgil/focal
```
in a new folder to contain the project. The folder structure will be
```
    DTCToolkit/                 <- the reconstruction and analysis code
    DTCToolkit/Analysis         <- User programs for running the code
    DTCToolkit/Classes          <- All the classes needed for the project
    DTCToolkit/Data             <- Data files: Range-energy look up tables, Monte Carlo code, LET data from experiments, the beam data from Groningen, ...
    DTCToolkit/GlobalConstants  <- Constants to adjust how the programs are run. Material parameters, geometry, ...
    DTCToolkit/HelperFunctions  <- Small programs to help running the code.
    DTCToolkit/OutputFiles      <- All output files (csv, jpg, ...) should be put here
    DTCToolkit/RootFiles        <- ROOT specific configuration files.
    DTCToolkit/Scripts          <- Independent scripts for helping the analysis. E.g. to create Range-energy look up tables from Monte Carlo data
    gate/                       <- All Gate-related files
    gate/python                 <- The DTC geometry builder
    projects/                   <- Other projects related to WP1
```

The best way to learn how to use the code is to look at the user programs, e.g. Analysis.C::DrawBraggPeakGraphFit which is the function used to create the Bragg Peak model fits and beam range estimation used in the 2017 NIMA article. From here it is possible to follow what the code does.
It is also a good idea to read through what the different classes are and how they interact:
* <code>Hit</code>: A (int x,int y,int layer, float edep) object from a pixel hit. edep information only from MC
* <code>Hits</code>: A <code>TClonesArray</code> collection of Hit objects
* <code>Cluster</code>: A (float x, float y, int layer, float clustersize) object from a cluster of <code>Hit</code>s The (x,y) position is the mean position of all involved hits.
* <code>Clusters</code>: A <code>TClonesArray</code> collection of <code>Cluster</code> objects.
* <code>Track</code>: A <code>TClonesArray</code> collection of <code>Cluster</code> objects... But only one per layer, and is connected through a physical proton track. Many helpful member functions to calculate track properties.
* <code>Tracks</code>: A <code>TClonesArray</code> collection of <code>Track</code> objects.
* <code>Layer</code>: The contents of a single detector layer. Is stored as a <code>TH2F</code> histogram, and has a <code>Layer::findHits</code> function to find hits, as well as the cluster diffusion model <code>Layer::diffuseLayer</code>. It is controlled from a <code>CalorimeterFrame</code> object.
* <code>CalorimeterFrame</code>: The collection of all <code>Layer</code>s in the detector.
* <code>DataInterface</code>: The class to talk to DTC data, either through semi-<code>Hit</code> objects as retrieved from Utrecht from the Groningen beam test, or from ROOT files as generated in Gate.

To run the code, do
```
[DTCToolkit] $ root Load.C
```
and ROOT will run the script <code>Load.C</code> which loads all code and starts the interpreter. From here it is possible to directly run scripts as defined in the <code>Analysis.C</code> file:
```
ROOT [1] drawBraggPeakGraphFit(...)
```

# MeltMigrator
MeltMigrator is a MATLAB®-based melt migration software developed to process three-dimensional mantle temperature and velocity data from user-supplied numerical models of mid-ocean ridges, calculate melt production and melt migration trajectories in the mantle, estimate melt flux along plate boundaries and predict crustal thickness distribution on the seafloor.

The software assumes the existance of a pre-constructed 3D numerical model of mantle flow and thermal structure. We mainly use finite element models developed with COMSOL Multiphysics, although the code could work with a variety of other models (by default, when not using COMSOL Multiphysics, the model result should be stored in a text file in columns with the sequence: x, y, z, vx, vy, vz, T).

The software first calculates the trajectory of melts according to a 3-step simplified melt extraction model:

  1) Melt is generated in the asthenosphere; it rises vertically due to its buoyancy and under the assumption that permeability is high enough for the background mantle flow ot be ignored.
  
  2) Melt accumulates in a decompaction channel lining a permeability barrier, typically at the base of the thermal lithosphere. Melt travels upward parallel to the barrier where it is slightly inclined.
  
  3) Melt reaches the surface when it encounters a Melt Extraction Zone, or MEZ, that is defined as a rectangular box around a plate boundary.

These steps were described in *Montési et al.* [2011] and *Gregg et al.* [2012]. This part of the software should work as long as the thermal structure is monotonic with respect to the vertical direction.
The end-product of this part of the softare is a prediciton of melt flux along the axis expressed as crustal thickness. Geochemical indicators are also computed.

Second, the software integrates the history of crustal accretion and predicts crustal thickness over the entire surface of the computational domain. This step is described by *Bai and Montési* [2015].

We have applied this code extensively to different mid-ocean ridge settings, and the same principles should be applicable to a broader range of tectonic enviroments. 

## Getting started:

To run the demo, launch MATLAB, change to MeltMigrator directory, and in command window, type `meltMain` and enter, the software will run and process the demo model automatically.

To use the software for their own models, users first need to convert their model results into MeltMigrator-compatible format (text file with x, y, z, vx, vy, vz, T stored in columns), modify the script **setParameters.m** to set key parameters based on their models, and then run `meltMain`. For debugging, we advise running meltMain initially cell by cell to examine the topology of the permeability barrier and ensure that the input parameters are appropriate for the calculation.

## User Workflow:

  1) Edit setParameters.m: Specifiy external model name, geometry, review choice of melting function and calibration.
  
  2) Run meltMain.m at least including the cell labelled lid sampling; The lid information is saved in matlab binary format.
  
    * Examine the graphical output of lidSample to understand the topology of the lid. 
 
    * Specify starting points for the saddle (variable SaddleInitialPoint in setParameters.m)
    
    * If desired, adjust the remaining parameters in setParameters.m
  
  3) Run meltMain.m to completion. Note that the lidSample.m should be skipped and the saved lid information be used instead.
  
Also verify from the graphical output that the tesselation of the lid does not include obvious artefacts such as overlapping tiles or gaps. If necessary, adjust the discretization parameters in setParameters.m to produce a better discretization.

## License and Citation:
MeltMigrator was developed by Laurent G.J. Montési, Mark D. Behn, Laura B. Hebert, and Hailong Bai.

MeltMigrator is released under the MIT license.

Please cite MeltMigrator in your publications if it helps your research:

> Bai, H., L. G. J. Montési, and M. D. Behn, 2017. MeltMigrator: a MATLAB-based Software for Modeling Three-dimensional Melt Migration and Crustal Thickness Variations at Mid-Ocean Ridges Following a Rules-Based Approach. Geochem. Geophys. Geosyst., in review.


## References:
> Montési, L. G. J, M. D. Behn, L. B. Hebert, J. Lin, and J. L. Barry, 2011. Controls on melt migration and extraction at the ultraslow Southwest Indian Ridge 10–16°E. J. Geophys. Res., 116, B10102. doi:10.1029/2011JB008259.
  
> Gregg, P. M., L. B. Hebert, L. G. J. Montési, and R. F. Katz, 2012. Geodynamic models of melt generation and extraction at mid-ocean ridges, Oceanography, 25, 8–88, doi:10.5670/oceanog.2012.05.
  
> Bai, H., and L. G. J. Montési, 2015. Slip-rate-dependent melt extraction at oceanic transform faults. Geochem. Geophys. Geosyst., 16, 401-419. doi:10.1002/2014GC005579.

# MeltMigrator
Matlab routines for solving melt migration out of the mantle and predicting crustal thickness.

The context for this work is the generation of partial melts in mid-ocean ridges, in geology/geophysics/geodynamics. 

The software assumes that the existance of a preexisting numerical model of mantle flow and associated structure. We use mainly Finite Element Models developed with COMSOL Multiphysics, although the code could work with a variety of other models (by default, when not using COMSOL, the model result should be stored in a text file in columns with the sequence: x, y, z, vx, vy, vz, T). We have applied this code extensively to different Mid-ocean ridge settings, although the same principles should be applicable to a broader range of tectonic enviroments. 

The software first calculates the trajectory of melts according to a 3-step simplified melt extraction model
  1) Melt is generated in the asthenosphere; it rises vertically due to its buoyancy and under the assumption that permeability is high enough for the background mantle flow ot be ignored.
  2) Melt accumulates in a decompaction channel lining a permeability barrier, typically at the base of the thermal lithosphere. Melt travels upward parallel to the barrier where it is slightly inclined. 
  3) Melt reaches the surface when it encounters a Melt Extraction Zone, or MEZ, that is defined as a rectangular box around a plate boundary
These steps were described in Montési et al., [2011] and Gregg et al. [2012]. This part of the software should work as long as the thermal structure is monotonic with respect to the vertical direction.
The end-product of this part of the softare is a prediciton of melt flux along the axis expressed as crustal thickness. Geochemical indicators are also computed.

Second, the software integrates the history of crustal accretion and predicts crustal thickness over the entire surface of the computational domain. This step is described by Bai and Montési [2015].

To run the demo, launch MATLAB, change to MeltMigrator directory, and in command window, type meltMain and enter, the software will run and process the demo model automatically.
To use the software for their own models, users first need to modify the script setParameters.m to set key parameters based on their models, and then run meltMain. For debugging, we advise running meltMain initially cell by cell to examine the topology of the permeability barrier and ensure that the input parameters are appropriate for the calculation.

A paper describing the code and its assumption is currently under review.

MeltMigrator was developed by Laurent G.J. Montési, Mark D. Behn, Laura B. Hebert, and Hailong Bai.

References: 
  Montési, L. G. J, M. D. Behn, L. B. Hebert, J. Lin, and J. L. Barry, 2011. Controls on melt migration and extraction at the ultraslow Southwest Indian Ridge 10–16°E. J. Geophys. Res., 116, B10102. doi:10.1029/2011JB008259.
  Gregg, P. M., L. B. Hebert, L. G. J. Montési, and R. F. Katz, 2012. Geodynamic models of melt generation and extraction at mid-ocean ridges, Oceanography, 25, 8–88, doi:10.5670/oceanog.2012.05.
  Bai, H., and L. G. J. Montési, 2015. Slip-rate-dependent melt extraction at oceanic transform faults. Geochem. Geophys. Geosyst., 16, 401-419. doi:10.1002/2014GC005579.




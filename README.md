Gridplates is an attempt to simulate the movements, sources, and sinks of sedimentary rocks through geologic time. Also the fluxes of carbon, and the climate state that would be constrained by that.  

The model is driven by a simulation within GPlates.  Any simulation should work but currently the only available simulation that has the required "dynamic polygon" representation of the plate outcrops is from Merdith.  Shapefiles were saved every million years, and converted to csv text files using pygplates and the python scripts provided in utils.  

Gridplates also uses the rotation file from GPlates, which specifies the rotations of everything on arbitrary time points.  Plates are rotated relative to other plates ultimately one reaches the non-rotating Earth.  This heirarcy of rotations is navigated by recursion.

The earth is represented as a 180 x 360 grid in latitude and longitude.  Each plate that is currently represented at the earth's surface is also allocated a global grid at the same resolution.  

There are many changes of plateID number in the reconstruction.  When this happens the information from the old ID has to be interpolated to the grid of the new ID.  The transitions were identified by eye and are read in from a file.  These transitions happen on a sub-time step of 1 myr.  

Each time step the plate rotation matrices are updated and plate information is interpolated to the world grid, in grid points where the GPlates plateID field says that a given plateID outcrops.  

Orogeny is controlled by giving the model global fields of 0, 1's to indicate the footprint of an uplift event, and settings in the code to control the time period of the event as well as scaling its magnitude.  The crust magically thickens and lifts up isostatically.  

Land sediment transport is diffusive with a transport coefficient set in params.jl . Parts of the domain are completely denuded of sediment in the solution of the time step.  These grid cells are found by iteration.  First they are presumed to be sediment-covered at the end of the step.  Then the step is taken.  If a cell's projected sediment thickness at the end of the step is less than zero, then it is excluded from the calculation in a subsequent iteration.  Adjacent cells then go negative, so multiple passes are required.

Ocean sediment deposition is governed by the constraint that we mustn't overfill the accomodation space (water depth).  Runoff from the continent is tabulated at coastal ocean grid cells, along with coastal CaCO3 depositon and whatever.  First the flux fields are smoothed based on a diffusion constant.  Then the coastal grid cells are queried for whether they overflow or not, and if they do, their mass is transferred inward, away from the coast.  This repeats until all the mass is accomodated.  

The file main.jl contains the master plan and major operations, including running a time series, creating plots and animations from the simulation, and creating an html directory with copies of all the files.  The model can be invoked at a Julia command line as 

include("main.jl")

and where you might want to go from there is basically in that file.  

Parameters tunable and otherwise are kept in params.jl .  Here you control what the simulation will be called, and how long it will span.  There is an option there for accelerated simulation of just the orography, without all the mud.  



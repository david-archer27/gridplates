Gridplates is an attempt to simulate the movements, sources, and sinks of sedimentary rocks through the Phanerozoic. Also the fluxes of carbon, and the climate state that would be constrained by them.  Latest results are presented at http://geosci.uchicago.edu/~archer/gridplates .  

The model is driven by a reconstruction within GPlates.  Any reconstruction should work but currently the only available simulation that has the required "dynamic polygon" representation of the plate outcrops is from Merdith.  Shapefiles were saved every million years, and converted to csv text files using pygplates and the python scripts provided in utils.  

Gridplates also uses the rotation file from GPlates, which specifies the rotations of everything on arbitrary time points.  Plates are rotated relative to other plates until ultimately one reaches the non-rotating Earth.  This heirarcy of rotations is navigated by recursion.

The earth is represented as a 180 x 360 grid in latitude and longitude.  Each plate at the earth's surface is also allocated a global grid at the same resolution.  

There are many renamings in the reconstrucrtion: changes of plateID number.  When this happens the information from the old ID has to be interpolated to the grid of the new ID.  The transitions were identified by eye and are read in from a file.  These transitions happen on a sub-time step of 1 myr.  

Each time step the plate rotation matrices are updated and plate information is interpolated to the world grid, in grid points where the GPlates plateID field says that a given plateID outcrops.  The model also reads a file containing the locations of continents in the GPlates reconstruction.  

Orogeny is controlled by giving the model global fields of 0, 1's to indicate the present-day footprint of an uplift event, and settings in the code to control the time period of the event as well as scaling its magnitude.  The crust magically thickens and lifts up isostatically at rates specified in create_orogenic_events(), in grid cells which rotate to the present-day locations of the orogenic footprint field.  

Land sediment transport is diffusive with a transport coefficient set in params.jl . Parts of the domain are completely denuded of sediment in the solution of the time step.  These grid cells are found by iteration.  First they are presumed to be sediment-covered at the end of the step.  Then the step is taken.  If a cell's projected sediment thickness at the end of the step is less than zero, then it is excluded from the calculation in a subsequent iteration.  Adjacent cells then go negative, so multiple passes are required.  

Ocean sediment deposition is governed by the constraint that we mustn't overfill the accomodation space (water depth).  Runoff from the continent is tabulated at coastal ocean grid cells, along with coastal CaCO3 depositon and whatever.  First the flux fields are smoothed into ocean grid points according to a diffusion constant.  Then the coastal grid cells are queried for whether they overflowed or not, and if they did, their mass is transferred inward, away from the coast.  This repeats until all the mass is accomodated. 

CaCO3 deposition in the ocean is based on three mechanisms.  Pelagic deposition depends on latitude, water depth, and an ocean CO3 parameter.  CaCO3 rolls off of shallow waters surrounding continents into coastal points at rates that depend on latitude and ocean CO3.  When sea level floods continental crust, CaCO3 deposits more-or-less up to sea level, although attenuated in high latitudes.  The global total burial rate of CaCO3 is specified, and the ocean CO3 parameter is iterated to find the value which gives the right flux.  

A "CaO" fraction of sediment (alongside unreactive "clays" and CaCO3) is intended to represent the CO2-consuming weathering flux as fresh clays weather.  Continental dissolution fluxes of both CaCO3 and CaO are based on a reconstructed latitudinal dependence of runoff, and observed (simplified) runoff / dissolution flux relationships.  Ultimately it should be possible to constrain an atmospheric CO2 parameter based on the constraint that the weathering of exposed bedrock plus clay "CaO" must provide enough Ca to convert a CO2 degassing flux to the atmoshere into CaCO3.   

The file main.jl contains the master plan and major operations, including running a time series, creating plots and animations from the simulation, and creating an html directory with copies of all the files.  The model can be invoked at a Julia command line as 

include("main.jl")

and where you might want to go from there is basically in that file.  Notable routines are run_timeseries(), animate_all(), create_timeseries_charts(), and create_output_html_directory( ), or just run_and_plot_simulation() and run_and_plot_orogeny_simulation().  

Parameters tunable and otherwise are kept in params.jl .  Here you control what the simulation will be called, and how long it will span.  The output files from the simulation are set to be placed alongside the "gridplates" git directory.  

The world state is saved every time step for diagnostics and graphics.  Plates files are saved periodically for use in restarting a simulation.  There has to be a separate file for each plate to avoid a 2 GB restriction in the BSON format, and you can't save the plates every step or they will eat your computer.  

The graphics are based on Makie and Plots, and can be used interactively at the command line as well as for making the saved output files. 



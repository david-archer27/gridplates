Gridplates is an attempt to simulate the movements, sources, and sinks of sedimentary rocks through the Phanerozoic. Also the fluxes of carbon, and the climate state that would be constrained by them.  Latest results are presented at http://geosci.uchicago.edu/~archer/gridplates .  

The model is driven by a reconstruction within GPlates.  Any reconstruction should work but currently the only available simulation that has the required "dynamic polygon" representation of the plate outcrops is from Merdith.  Shapefiles were saved every million years, and converted to csv text files using pygplates and the python scripts provided in utils.  

Gridplates also uses the rotation file from GPlates, which specifies the rotations of everything on arbitrary time points.  Plates are rotated relative to other plates until ultimately one reaches the non-rotating Earth.  This heirarcy of rotations is navigated by recursion.

The earth is represented as a 180 x 360 grid in latitude and longitude.  Each plate at the earth's surface is also allocated a global grid at the same resolution.  

There are many renamings in the reconstruction: changes of plateID number.  When this happens the information from the old ID has to be interpolated to the grid of the new ID.  The transitions were identified by eye and are read in from a file.  These transitions happen on a sub-time step of 1 myr.  There is an animation of all the locations of plate ID transplants in the 
gallery page.  

Each time step the plate rotation matrices are updated and plate information is interpolated to the world grid, in grid points where the GPlates plateID field says that a given plateID outcrops.  The model also reads a file containing the locations of continents in the GPlates reconstruction.  

Orogeny is currently controlled by comparison to the elevations from Scotese which have been rotated
by Chris into their present-day locations, from which they are rotated to match the Merdith 
reconstruction in gridplates.  Only elevations above 1500 meters are uplifted, which is done by 
thickening the crust.  

Inputs of large igneous provinces and ophiolites are similarly done in modern-day locations, twisted
back to their paleo locations at the appropriate times, which are the numbers in the input grid. A 
crust_composition variable has been added but not yet implemented for weathering.  Ice sheets are 
imposed based on CESM results using the Scotese continental configurations, rotated to present-day by
Chris, then rotated to Merdith positions at run time.  Ice is imposed where the temperature was less
than -10 degrees C.  

Land sediment transport is diffusive with a transport coefficient set in params.jl . Parts of the domain are completely denuded of sediment in the solution of the time step.  These grid cells are found by iteration.  First they are presumed to be sediment-covered at the end of the step.  Then the step is taken.  If a cell's projected sediment thickness at the end of the step is less than zero, then it is excluded from the calculation in a subsequent iteration.  Adjacent cells then go negative, so multiple passes are required.  Sediment that is denuded by transport in this way, or denuded by the presence of
ice sheets, is imposed as a boundary condition to the land sediment transport domain, by summing the 
initial sediment cover in a denuded zone, and distributing it evenly into grid cells surrounding the
denuded zone.  CaCO3 deposited on land surface points is considered fixed to bedrock, and is not subject to 
the downhill sediment transport, only dissolution.  

Ocean sediment deposition is governed by the constraint that we mustn't overfill the accomodation space (water depth).  Runoff from the continent is tabulated at coastal ocean grid cells, along with coastal CaCO3 depositon and whatever.  First the flux fields are smoothed into ocean grid points according to a diffusion constant.  Then the coastal grid cells are queried for whether they overflowed or not, and if they did, their mass is transferred inward, away from the coast.  This repeats until all the mass is accomodated. 

CaCO3 deposition in the ocean is based on three mechanisms.  Pelagic deposition depends on latitude, water depth, and an ocean CO3 parameter.  CaCO3 rolls off of shallow waters surrounding continents into coastal points at rates that depend on latitude and ocean CO3.  When sea level floods continental crust, CaCO3 deposits more-or-less up to sea level, although attenuated in high latitudes.  The global total burial rate of CaCO3 is specified, and the ocean CO3 parameter is iterated to find the value which gives the right flux.  There is currently some bug in the CaCO3 budget.  

A "CaO" fraction of sediment (alongside unreactive "clays" and CaCO3) is intended to represent the CO2-consuming weathering flux as fresh clays weather.  Continental dissolution fluxes of both CaCO3 and CaO are based on a rudimentary runoff field, and observed (simplified) runoff / dissolution flux relationships.  Ultimately it should be possible to constrain an atmospheric CO2 parameter based on the constraint that the weathering of exposed bedrock plus clay "CaO" must provide enough Ca to convert a CO2 degassing flux to the atmoshere into CaCO3. 

The runoff applied to the weathering parameterization is an eyeball fit to the results of Baum et al. 2022.  A comparison of the fit to the original is shown as an animation output.  The fraction denuded is compared with the estimation from Holland's book that igneous and metamorphic comprise 25% of the Earth's surface. The CaO weathering fluxes chart shows the relative fluxes from hard rock vs. sediments, compared with the results of BLAG table 1 that 2/3 of "CaO" or calcium silicate weathering is primary rather than from clays.  Maps of sediment distribution are compared with present-day as one of the charts also; currently the model is still somewhat low.  Finally, the elevations are compared with the Scotese reconstruction as an animation of maps, and a time series plot of the mean elevation in a chart. 

The file main.jl contains the master plan and major operations, including running a time series, creating plots and animations from the simulation, and creating an html directory with copies of all the files.  The model can be invoked at a Julia command line as 

include("main.jl")

and where you might want to go from there is basically in that file.  Notable routines are run_timeseries(), animate_all(), create_timeseries_charts(), and create_output_html_directory( ), or just run_and_plot_simulation() and run_and_plot_orogeny_simulation().  

Parameters tunable and otherwise are kept in params.jl .  Here you control what the simulation will be called, and how long it will span.  The output files from the simulation are set to be placed alongside the "gridplates" git directory as specified here.  Model diagnostic parameters are stored in two giant arrays, the contents of which are specified in this file.  One array is specific for fields that apply to the individual sediment fractions (clay, CaCO3, and CaO).  The other is for fields like runoff that do not have separate values for the different sediment fractions.  These arrays are saved in the world output files, and all the post-processing and plotting is done by reading the world files.  

The world state is saved every time step for diagnostics and graphics.  Plates files are saved periodically for use in restarting a simulation.  There has to be a separate file for each plate to avoid a 2 GB restriction in the BSON format, and you can't save the plates every step or they will eat your computer.  

The graphics are based on Makie and Plots, and can be used interactively at the command line as well as for making the saved output files. Makie is beautiful but sometimes has perspective problems or other inexplicible bugs.  On my apple M1 I have to use v0.4.7 ; more recent versions don't work.  Also Julia freezes at some point in making animations.  Hence two routines, animate_until_you_almost_puke(), which should be followed by restarting Julia, then doing animate_the_rest().

It is worth using multiple cores by setting that preference in VSC or invoking the julia executable using --threads=16 (or however many).  Because cores in CPUs can multi-thread, it is worthwhile to set
maybe twice as many threads as your CPU has cores.  








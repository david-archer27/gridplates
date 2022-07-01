# data structures
mutable struct rotation_struct
    id
    age
    latitudeaxis
    longitudeaxis
    angle
    parent
end
mutable struct plate_struct
    plateID         # scalar, from gplates rotation file
    crust_type       # 2d arrays
    crust_age
    crust_thickness
    crust_density
    tectonics
    sediment_thickness  # geologic record, dimensions of x,y
    sediment_surface_fractions # x,y,sed type 
    sediment_layer_thickness # x,y,time_bins, not used on land points
    sediment_layer_fractions  # x, y, sed type, time_bin
    rotationmatrix  # 3x3 heirarchy-resolved rotation matrix
    resolvetime
    parentstack     # list of parents at last resolve time
    lastparentstack
    firstappearance
    lastappearance
end
mutable struct orogenic_event_struct
    name
    footprint
    onset
    finish
    altitude
end
mutable struct world_struct
    # create and manipulate a global world for snapshot composite view of rotated plates
    age
    sealevel
    plateID         # pertaining to plate tectonics
    plateIDlist
    crust_type      # beginning of the 9 state variables
    crust_age
    crust_thickness # pertaining to geomorphology
    crust_density
    geomorphology    # for sediment transport, subaereal or submarine
    tectonics
    sediment_thickness  # total thickness
    sediment_surface_fractions # x,y,sed type 
    sediment_layer_thickness # x,y,time_bins, not used on land points
    sediment_layer_fractions  # x, y, sed type, time_bin
    elevation_offset
    surface_elevation
    freeboard  # end of 9 state variables
    diags
    frac_diags
end

# diagnostics variable system
world_diag_names = ["ocean_created_plate_area",
    "continent_2_ocean_plate_area",
    "continent_2_ocean_sediment_leak",
    "continent_created_plate_area",
    "ocean_2_continent_plate_area",
    "ocean_2_continent_sediment_leak",
    "ocean_subduct_plate_area",
    "ocean_subduct_age_plate_area",
    "ocean_subduct_sediment_volume",
    "continent_subduct_plate_area",
    "continent_subduct_sediment_volume",
    "ocean_2_continent_world_area",
    "continent_2_ocean_world_area",
    "continent_initialization_sediment_volume",
    "IDchange_plate_area",
    "continent_orogenic_uplift_rate",   # set in outer time loop, uplift +, meters / Myr
    "subduction_orogenic_uplift_rate",  
    "crust_erosion_rate",               # units m/Myr, calc on substep
    "crust_clay_source_rate",
    "aolean_clay_erosion_rate", 
    "aolean_clay_deposition_rate", 
    "land_orogenic_clay_flux", # in orogeny-neighboring land grid cells
    "land_CaCO3_dissolution_rate", # subaereal erosion
    "continental_CaCO3_deposition_rate", # when flooded
    "land_sediment_deposition_rate",
    "seafloor_sediment_deposition_rate",
    "global_sediment_deposition_rate",
    "world_grid_sediment_change",
    "coastal_orogenic_clay_flux", # in coastal ocean points, boundary fluxes
    "coastal_CaCO3_flux",
    "pelagic_CaCO3_deposition_rate",
    "seafloor_delta_CO3"]
world_frac_diag_names = [
    "ocean_subduct_sediment_fraction_volume",
    "continent_subduct_sediment_fraction_volume",
    "continent_initialization_sediment_fraction_volume",
    "continent_2_ocean_sediment_fraction_leak",
    "ocean_2_continent_sediment_fraction_leak",
    "land_sediment_fraction_deposition_rate",
    "land_trapped_sediment_rate",
    "coastal_sediment_fraction_runoff_flux",
    "ocean_sediment_fraction_influx",
    "seafloor_sediment_fraction_deposition_rate", 
    "seafloor_sediment_fraction_overflow",
    "global_sediment_fraction_deposition_rate",
    "global_sediment_fraction_deposition_ratio"]


# grid 
nx = 360; ny = 180
zoomplot_ix = 190; zoomplot_iy = 130; zoomplot_size = 20
xcoords = collect(-180. + 180/nx: 360/nx: 180. - 180/nx)
ycoords = collect(-90. + 180/nx: 360 / nx: 90. - 180/nx)
areabox = fill(0.,ny)
delta_x = fill(0.,ny); delta_y = 110.e3
for iy in 1:ny
    areabox[iy] = 110. * 110. * cos(ycoords[iy] / 180. * pi)  * 1e6 # m2
    delta_x[iy] = 110.e3 * cos(ycoords[iy] / 180. * pi) # m
end
shelf_depth_clastics = -100.; shelf_depth_CaCO3 = -10.

log_IO = 0 # dummy file handle for log output file

# world.crust_type[] values
ocean_crust = 1
continent_crust = 2 

# world.geomorphology
new_ocean_crust = 1
pelagic_seafloor = 2
coastal_depocenter = 3 
ocean_shelf = 4 # non-depositing coastal zone
filled_basin = 5
sedimented_land = 6 # diffusive zone
exposed_basement = 7 

# world.tectonics
new_ocean_crust = 1
new_continent_crust = 2
ocean_to_continent = 3 # set in update_world_continents_from_file()
continent_to_ocean = 4
subducting_ocean_crust = 5
subducting_continent_crust = 6

# plate.tectonics, maps of plate transformations
not_yet_defined = -1
not_at_surface = 0
#new_ocean_crust = 1
#new_continent_crust = 2
#ocean_to_continent = 3 # set in update_world_continents_from_file()
#continent_to_ocean = 4
#subducting_ocean_crust = 5
#subducting_continent_crust = 6
new_plateID = 7 # set in copy_plate_point_plate12coord!
deleted_plateID = 8

# sediment types
sediment_type_names = ["Clay","CaCO3"]; clay_sediment = 1; CaCO3_sediment = 2
n_sediment_types = length(sediment_type_names)
initial_sediment_type = clay_sediment 
initial_land_sediment_thickness = 1.; initial_ocean_sediment_thickness = 1.

# time 
earliesttime = 300.
time_step = 2. # Myr
sediment_record_resolution = 10.
sediment_time_bins = Float32[  ]
for atime in range( earliesttime, step = -sediment_record_resolution, stop = 0.)
    push!( sediment_time_bins, atime )
end
#Float32[635,541,250,65,0] # resolution of the sediment record on the plates
n_sediment_time_bins = length(sediment_time_bins)
geo_interval_time_bins = [ 540, 485, 444, 419,  # these are the beginnings
    359, 299, 251, 201, 145,
    66, 56, 34, 23, 5, 2, 0 ]
geo_interval_time_bins = convert(Array{Float32},geo_interval_time_bins) # keeps netcdf library from puking
geo_interval_names = [ "Cambrian", "Ordovician", "Silurian", "Devonian",
    "Carboniferous","Permian","Triassic","Jurassic","Cretaceous",
    "Paleocene","Eocene","Oligocene","Miocene","Pliocene","Pleistocene"]

code_base_directory = "/Users/archer/Synched/papers/spongeball/codebase"
animation_output_directory = "/Users/archer/Synched/papers/spongeball/animations/"
data_output_directory = "/Users/archer/Synched/papers/spongeball/outfiles"

# geophysics parameters
continent_crust_h0 = 16500.; continent_crust_h0_max = 35000.
ocean_crust_h0 = 4000.
rho_ocean_crust = 3.; rho_continent_crust = 2.5; rho_mantle = 3.5
rho_sediment = 2.; rho_seawater = 1.024
sediment_freeboard_expression = 1. - rho_sediment / rho_mantle
crust_freeboard_expression = 1. - rho_continent_crust / rho_mantle
mantle_T0 = 2000.
ocean_T0 = 0.

# geomorphology parameters
land_base_diffcoeff = 1.e5
seafloor_base_diffcoeff = 3.e4 # m2 / yr
crust_smooth_coeff = 1.e4
land_CaCO3_dissolving_altitude = 10. # meters
land_CaCO3_dissolution_rate = 0. # meters/Myr
# 1000 meters dissolves in 100 Myr
orogenic_base_uplift_rate = 7.5 # m / Myr
subduction_orogeny_parameter = 2.e-8 # small because mult by area
aolean_erosion_rate_constant = 1. / 2000. # Myr
orogenic_erosion_tau_apparent = 200. # Myr

# CaCO3 parameters
global_CO2_degassing_rate = 7.5E12 # mol / yr
global_terrestrial_CaCO3_dissolution_rate = 100.E15 / 100. / 1000.
#                                       g CaCO3/kyr mol/yr   mol/yr
CaCO3_deposition_lat_scale = 45.

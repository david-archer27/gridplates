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
    surface_type
    sediment_thickness  # geologic record, dimensions of x,y,time_bin
    sediment_fractions  # x, y, sed type, time_bin
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
    surface_type    # for sediment transport, subaereal or submarine
    sediment_thickness  # nx,ny,n_sediment_time_bins
    sediment_fractions  # nx,ny,n_sediment_types,n_sediment_time_bins
    elevation_offset
    surface_elevation
    freeboard  # end of 9 state variables
    diags
    frac_diags
end

# diagnostics variable system
world_diag_names = ["ocean_created_plate_area",
    "continent_2_ocean_plate_area",
    "continent_created_plate_area",
    "ocean_2_continent_plate_area",
    "ocean_subduct_plate_area",
    "ocean_subduct_age_plate_area",
    "ocean_subduct_sediment_plate_volume",
    "continent_subduct_plate_area",
    "ocean_2_continent_world_area",
    "continent_2_ocean_world_area",
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
    "sediment_deposition_rate",
    "total_sediment_thickness",   
    "coastal_orogenic_clay_flux", # in coastal ocean points, boundary fluxes
    "coastal_CaCO3_flux",
    "pelagic_CaCO3_deposition_rate",
    "seafloor_delta_CO3",
    "surface_type_changes"]
world_frac_diag_names = [
    "ocean_subduct_sediment_fraction_volume",
    "land_sediment_fraction_deposition_rate",
    "coastal_sediment_fraction_runoff_flux",
    "ocean_sediment_fraction_influx",
    "seafloor_sediment_fraction_deposition_rate", 
    "seafloor_sediment_fraction_overflow",]


# grid 
nx = 360; ny = 180
xcoords = collect(-180. + 180/nx: 360/nx: 180. - 180/nx)
ycoords = collect(-90. + 180/nx: 360 / nx: 90. - 180/nx)
areabox = fill(0.,ny)
delta_x = fill(0.,ny); delta_y = 110.e3
for iy in 1:ny
    areabox[iy] = 110. * 110. * cos(ycoords[iy] / 180. * pi)  * 1e6 # m2
    delta_x[iy] = 110.e3 * cos(ycoords[iy] / 180. * pi) # m
end
shelf_depth_clastics = -100.; shelf_depth_CaCO3 = -10.

# world.crust_type[] values
notyetdefined = -2; notatsurface = -1; ocean_crust = 1
continent_crust = 2 

# world.surface_type[] values
pelagic_seafloor = 1 
coastal_depocenter = 2 
ocean_shelf = 3 # non-depositing coastal zone
sedimented_land = 4 # diffusive zone
exposed_basement = 5 

# world.surface_type_change[] values
land_soil_denuded = 5 # highest value first sets color map
burying_land = 4
accomodation_filled = 3
incoming_turbidite = 2
basin_filled = 1

# sediment types
sediment_type_names = ["Clay","CaCO3"]; clay_sediment = 1; CaCO3_sediment = 2
n_sediment_types = length(sediment_type_names)
initial_sediment_type = clay_sediment; initial_sediment_thickness = 1.

# time 

earliesttime = 150.
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
land_CaCO3_dissolving_altitude = 10. # meters
land_CaCO3_dissolution_rate = 10. # meters/Myr
# 1000 meters dissolves in 100 Myr
aolean_erosion_rate_constant = 0. # 1. / 100. # Myr

# CaCO3 parameters
global_CO2_degassing_rate = 7.5E12 # mol / yr
global_terrestrial_CaCO3_dissolution_rate = 100.E15 / 100. / 1000.
#                                       g CaCO3/kyr mol/yr   mol/yr
CaCO3_deposition_lat_scale = 45.

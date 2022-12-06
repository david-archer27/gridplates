# data structures
mutable struct rotation_struct
    id
    age
    latitudeaxis
    longitudeaxis
    angle
    parent
end

mutable struct orogenic_event_struct
    name
    footprint
    onset
    finish
    uplift_scale
end
mutable struct world_struct
    # create and manipulate a global world for snapshot composite view of rotated plates
    age
    sealevel
    atmCO2
    plateID         # pertaining to plate tectonics
    plateIDlist
    crust_type      # beginning of the 10 state variables
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
    freeboard  
    subducted_land_sediment_volumes # [land,ocean], not maps
    subducted_ocean_sediment_volumes
    initial_ocean_sediment_inventories
    initial_land_sediment_inventories
    diags
    frac_diags
end
mutable struct plate_struct
    plateID         # scalar, from gplates rotation file
    crust_type       # beginning of the 10 state variables
    crust_age
    crust_thickness
    crust_density
    geomorphology
    tectonics
    #potential_uplift 
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
# diagnostics variable system
world_diag_names = ["ocean_created_plate_area",
    "continent_2_ocean_plate_area",
    "continent_2_ocean_sediment_leak",
    "continent_created_plate_area",
    "ocean_2_continent_plate_area",
    "ocean_subduct_plate_thickness",
    "ocean_subduct_age_plate_area",
    "ocean_subduct_age_thickness",
    #"ocean_subduct_sediment_thickness",
    "continent_subduct_plate_thickness",
    #"continent_subduct_sediment_thickness",
    "ocean_2_continent_world_area",
    "continent_2_ocean_world_area",
    "continent_initialization_sediment_volume",
    "IDchange_plate_area",
    "continent_orogenic_uplift_rate",   # set in outer time loop, uplift +, meters / Myr
    "subduction_orogenic_uplift_rate",  
    "crust_erosion_rate",               # units m/Myr, calc on substep
    "land_orogenic_Ca_source_rates",
    #"crust_clay_source_rate",
    #"aolean_clay_erosion_rate", 
    #"aolean_clay_deposition_rate", 
    #"land_orogenic_clay_flux", # in orogeny-neighboring land grid cells
    #"land_CaCO3_dissolution_rate", # subaereal erosion
    "continental_CaCO3_deposition_rate", # when flooded
    "land_sediment_deposition_rate",
    "land_Q_runoff_field",
    "land_sediment_weathering_index",
    "seafloor_sediment_deposition_rate",
    "global_sediment_deposition_rate",
    "world_grid_sediment_change",
    #"coastal_orogenic_clay_flux", # in coastal ocean points, boundary fluxes
    "coastal_CaCO3_flux",
    "pelagic_CaCO3_deposition_rate",
    "seafloor_delta_CO3"]
world_frac_diag_names = [
    "ocean_subduct_sediment_fraction_thickness",
    "continent_subduct_sediment_fraction_thickness",
    "continent_initialization_sediment_fraction_volume",
    "continent_2_ocean_sediment_fraction_leak",
    "ocean_2_continent_sediment_fraction_displaced",
    "continent_2_ocean_sediment_fraction_displaced",
    "denuded_sediment_source_fraction_flux",
    "denuded_land_boundary_fraction_flux",
    "denuded_coastal_boundary_fraction_flux",
    "aolean_erosion_fraction_flux",
    "aolean_deposition_fraction_flux",
    "crust_orogenic_fraction_flux",
    "land_orogenic_fraction_flux",
    "coastal_orogenic_fraction_flux",
    "land_sediment_fraction_dissolution_rate",
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
zoomplot_ix = 60; zoomplot_iy = 120; zoomplot_size = 20
#plot_resolution_x = 3000; plot_resolution_y = 1800
plot_resolution_x = 1600; plot_resolution_y = 1000
xcoords = collect(-180. + 180/nx: 360/nx: 180. - 180/nx)
ycoords = collect(-90. + 180/nx: 360 / nx: 90. - 180/nx)
areabox = fill(0.,ny)
delta_x = fill(0.,ny); delta_y = 110.e3
n_equiv_boxes = fill(1,ny)
for iy in 1:ny
    areabox[iy] = 110. * 110. * cos(ycoords[iy] / 180. * pi)  * 1e6 # m2
    delta_x[iy] = 110.e3 * cos(ycoords[iy] / 180. * pi) # m
    n_equiv_boxes[iy] = Int(floor(delta_y / delta_x[iy]))
end
shelf_depth_clastics = -100.; shelf_depth_CaCO3 = -0.1

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
burying_land = 7

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
sediment_type_names = ["Clay","CaCO3","CaO"]
clay_sediment = 1; CaCO3_sediment = 2; CaO_sediment = 3
n_sediment_types = length(sediment_type_names)
initial_sediment_fractions = [ 0.975,0.,0.025 ] # adjusted bc not pure CaO [ 0.85,0.,0.15 ] # present-day sed avg: Holland
orogenic_sediment_source_fractions = [ 0.95,0.,0.05 ] # a bit higher for fresh clay?
initial_land_sediment_thickness = 1.; initial_ocean_sediment_thickness = 1.

# Time 
earliesttime = 600.
main_time_step = 5. # Myr
sub_time_step = 1.
sediment_record_resolution = 20.
sediment_time_bins = Float32[  ]
for atime in range( earliesttime, step = -sediment_record_resolution, stop = 0.)
    push!( sediment_time_bins, atime )
end
#Float32[635,541,250,65,0] # resolution of the sediment record on the plates
n_sediment_time_bins = length(sediment_time_bins)
geo_interval_time_bins = Float64[ 540, 485, 444, 419,  # these are the beginnings
    359, 299, 251, 201, 145,
    66, 56, 34, 23, 5, 2, 0 ]
geo_interval_names = [ "Cambrian", "Ordovician", "Silurian", "Devonian",
    "Carboniferous","Permian","Triassic","Jurassic","Cretaceous",
    "Paleocene","Eocene","Oligocene","Miocene","Pliocene","Pleistocene"]
sealevel_timepoints = [550.,500.,450.,250., 75.,0.]
sealevel_values =     [  0.,  0.,200.,  0.,250.,0.]
sealevel_base = - 4714.7
atmCO2_timepoints = [550.,0.]
atmCO2_values = [400.,400.]
atmCO2_base = 400.

#sealevel_timepoints = [100.,50.,0.]
#sealevel_values = [0.,2.,0.]
#sealevel_timepoints = [100.,0.]
#sealevel_values = [0.,0.]

output_tag = "tuesday"

code_base_directory = pwd() # "gridplates"
plateID_input_directory = code_base_directory * "/plates"
continent_input_directory = code_base_directory * "/continents"
utils_directory = code_base_directory * "/utils"
data_directory = code_base_directory * "/data"

output_location = code_base_directory[1:findlast("/", code_base_directory )[1]-1]
output_directory = output_location * "/" * "outfiles." * output_tag
html_directory = output_location * "/html." * output_tag

animation_directory = output_directory * "/animations"
world_directory = output_directory * "/world"
plate_directory = output_directory * "/plates"
continents_directory = output_directory * "/continents"
charts_directory = output_directory * "/charts"

code_backup_directory = output_directory * "/code_bak"
scotese_data_directory = code_base_directory * "/data/scotese_elevation_files"
animation_n_step = 1
animation_initial_age = earliesttime; animation_final_age = 0
save_plate_transplants_image_number = 0

# configuration
enable_land_runoff = true
enable_aolean_transport = true
enable_cont_update_from_files = true
enable_crust_cont2ocn = true    
enable_crust_ocn2cont = true
enable_orogeny = true
enable_subduction_orogeny = false
eliminate_regrid_drift = true
enable_eyeball_changing_plateIDs = true
enable_sealevel_change = true
enable_geomorph = true
enable_atmCO2_change = false
enable_continental_CaCO3_deposition = true
enable_CaCO3_land_diss_ocean_repcp = true
#enable_watch_plate_transplant_misses = false
enable_watch_plate_transplants = false
enable_save_plate_transplant_images = false
enable_watch_denuding_landscape = false
enable_remask_plate_diagnostics = false
enable_changing_plateID_diagnostics = false
enable_step_tectonics_world_diagnostics = true
enable_step_tectonics_plate_diagnostics = false
enable_substep_tectonics_diagnostics = false
enable_watch_orogeny_substeps = false
enable_check_sediment_layer_inventory = false
enable_step_everything_diagnostics = true
enable_step_geomorph_diagnostics = true
enable_extraverbose_step_geomorph_diagnostics = true
enable_distribute_ocean_fluxes_diagnostics = true
enable_watch_ocean_offshore_transport = false

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
orogeny_smooth_coeff = 1000. 
crust_smooth_coeff = 1000.
orogenic_base_uplift_meters = 2.e4 # from 16 -> 35 km through the event
orogenic_uplift_parameter = 600. # m / Myr base
subduction_orogeny_parameter = 7500. 
subduction_orogeny_smooth_coeff = 1.e4
subduction_orogeny_hor_offset = 5 * delta_y
mountain_max_altitude_target = 10.e3 # m
max_uplift_rate_target = 600. # m / Myr
mean_elevation_land_target = 800. # m
orogenic_area_fraction_target = 0.1
orogenic_area_width = 1.e6 # m
orogenic_erosion_tau_apparent = mountain_max_altitude_target / 
    max_uplift_rate_target # Myr
land_base_diffcoeff = max_uplift_rate_target * 
    orogenic_area_width * orogenic_area_fraction_target * 
    orogenic_area_width / ( 2. * mean_elevation_land_target ) / 1.e6 # m/yr

land_base_diffcoeff *= 1.0
orogenic_erosion_tau_apparent *= 45.
orogenic_uplift_parameter *= 2.
#subduction_orogeny_smooth_coeff = 0.

runoff_tropical_max_rate = 1. # m / yr
runoff_tropical_max_width = 12.
runoff_east_coast_penetration_scale = 3.e6
runoff_subpolar_max_rate = 1.5
runoff_subpolar_max_lat = 45.
runoff_subpolar_max_width = 7.5
runoff_west_coast_penetration_scale = 3.e6
minimum_runoff = 0.1
runoff_smoothing = 100.


exposed_basement_CO2uptake_coeff = 0.25 # Suchet 2003 combined shield,basalt,acid volc
sediment_CaCO3_CO2uptake_coeff = 1.6 
sediment_CaO_CO2uptake_coeff = 0.6


function create_orogenies()
    orogenic_events = Dict()
    orogenic_events["Pan_African"] =
        create_orogenic_event("pan_african",650.,550.,0.5)  # africa, south america
    orogenic_events["Avalonian"] =
        create_orogenic_event("taconian",650.,500.,0.5)
    orogenic_events["East African"] =
        create_orogenic_event("e_african",550.,400.,0.25)
    orogenic_events["Taconian"] =
        create_orogenic_event("taconian",490.,440.,3.)
    orogenic_events["Calcedonian/Acadian"] =
        create_orogenic_event("calcedonian",460.,330.,0.5)
    orogenic_events["Borchgrevink"] =
        create_orogenic_event("borchgrevink",440.,400.,0.5)
    orogenic_events["Hercynian/Alleghenian/Uralian"] =
        create_orogenic_event("hercynian",300.,240.,0.25) # = variscan
        # changed end from 250 so something would be uplifting at 0 sl 250
    # Amurian (Japan) should be covered by subduction_orogeny
    orogenic_events["Indo Sinean"] =
        create_orogenic_event("indo_sinean",200.,180.,0.5)
    orogenic_events["Laramide"] =
        create_orogenic_event("laramide",400.,0.,0.2)
    orogenic_events["aolean_deposition_fraction_flux"] =
        create_orogenic_event("andean",400.,0.,0.5)
    orogenic_events["Cimmerian"] =
        create_orogenic_event("cimmerian",180.,150.,1.)
    orogenic_events["Mongol Okhotsk"] =
        create_orogenic_event("mongol_okhotsk",140.,130.,1.)
    orogenic_events["Wrangellian"] =
        create_orogenic_event("wrangellian",110.,80.,0.75)
    orogenic_events["Verkhoyansk"] =
        create_orogenic_event("verkhoyansk",100.,70.,0.75)
    orogenic_events["Alpine"] =
        create_orogenic_event("alpine",50.,0.,0.75)
    orogenic_events["Himalayan"] =
        create_orogenic_event("himalayan",40.,0.,0.75)
    return orogenic_events
end

if enable_aolean_transport
    aolean_erosion_rate_constant = 1.e-2 # Myr
else
    aolean_erosion_rate_constant = 0.
end

#land_CaCO3_dissolving_altitude = 10. # meters
land_sediment_dissolution_rate_constants = [0.0,0.01,0.01] # per Myr
land_sediment_dissolution_2xCO2 = [0.,0.,0.]
global_CaCO3_net_burial_flux = 1.e15 # today fluxes
seafloor_base_diffcoeff = 3.e5 # m2 / yr


# CaCO3 parameters
global_CO2_degassing_rate = 7.5E12 # mol / yr
global_terrestrial_CaCO3_dissolution_rate = 100.E15 / 100. / 1000.
#                                       g CaCO3/kyr mol/yr   mol/yr
CaCO3_deposition_lat_scale = 30.

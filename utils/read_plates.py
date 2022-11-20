import sys
sys.path.insert(1,'pygplates_rev18_python27_MacOS64/')
import pygplates
import pickle
#import matplotlib.pyplot as pyplot

#input_feature_collection = pygplates.FeatureCollection('Matthews_etal_GPC_2016_ContinentalPolygons.gpmlz')
#input_feature_collection = pygplates.FeatureCollection('Matthews_etal_GPC_2016_TopologyBuildingBlocks.gpmlz')
input_feature_collection = pygplates.FeatureCollection('ScoteseLab/myModernPlateBoundaries.gpml')

#print(input_feature_collection)

is_continent = []

for feature in input_feature_collection:
	    # Print the feature type (Coastline) and the name of the coastline.

    # Print the description of the coastline.
 #   print '  description: %s' % feature.get_description()

    # Print the plate ID of the coastline.
     # Print the length of the coastline geometry(s).
    # There could be more than one geometry per feature.
    for geometry in feature.get_geometries():
        try:
            polygon = pygplates.PolygonOnSphere(geometry, allow_one_or_two_points=False)
     #       print '  plate ID: %d' % feature.get_reconstruction_plate_id()
     #       print '  %f -> %f' % feature.get_valid_time()
     #       print '%s: %s' % (feature.get_feature_type().get_name(), feature.get_name())
            for i_grid in range(-180,180):
                is_continent_row = []
                for j_grid in range(-90,90):
                    if polygon.is_point_in_polygon((j_grid,i_grid)):
                        print i_grid, j_grid
        except:
            x=2

 
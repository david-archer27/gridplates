import sys
sys.path.insert(1,'pygplates_rev18_python27_MacOS64/')
import pygplates
import csv

def create_grid(collection):
    gridded = []
    for irow in range(90):
        gridded_row = []
        for jrow in range(180):
            gridded_row.append(-1)
        gridded.append(gridded_row)
    for feature in collection:
        for geometry in feature.get_geometries():
            id = feature.get_reconstruction_plate_id()
            feattype = feature.get_feature_type().get_name()
            featname = feature.get_name()
#            print '%d %s %s' % (id, feattype, featname)
#            if id == 304:
#                points = list(geometry)
#                print points
            polygon = pygplates.PolygonOnSphere(geometry, allow_one_or_two_points=False)
            for irow in range(90):
                lat = (irow*2. - 90.) + 1. 
                for icol in range(180):
                    longt = (icol*2 - 180.) + 1.
                    if polygon.is_point_in_polygon((lat,longt)):
                        gridded[irow][icol] = id
#                        if irow == 60:
#                        if icol == 90:
#                            print '%d %d %d ' % (lat, longt, id)
#                        print '%d ' % (id)
    return gridded

def list_plates(collection):
    for feature in collection:
        for geometry in feature.get_geometries():
            id = feature.get_reconstruction_plate_id()
            feattype = feature.get_feature_type().get_name()
            featname = feature.get_name()
           # print '%d %s %s' % (id, feattype, featname)

def write_grid(grid, filename):
    with open(filename,"wb") as csvfile:
        writer = csv.writer(csvfile)
        for row in grid:
            writer.writerow(row)

def filter_grid(grid, deleteIDs, goesintoID):
    for irow in range(90):
        for icol in range(180):
            for deleteID in deleteIDs:
                if grid[irow][icol] == deleteID:
                    grid[irow][icol] = goesintoID
            if grid[irow][icol] == -1:
                grid[irow][icol] = goesintoID
    return grid


timeSlices = range(110, 126, 1)
deleteIDs = []
goesintoID = 715

def process_slices():
  for iSlice in timeSlices:
#      print "time slice beginning", iSlice
      platesFile = 'Matthews_etal_GPC_2016_MesozoicCenozoic_PlateTopologies/topology_' + str(iSlice) + '.00Ma.shp'
      plates_collection = pygplates.FeatureCollection(platesFile)
      plateIDs = create_grid(plates_collection)
      plateIDs = filter_grid(plateIDs,deleteIDs,goesintoID)
      plateOutFile = "plateIDs." + str(iSlice) + ".csv"
      write_grid(plateIDs,plateOutFile)
      print "time slice ", iSlice

process_slices()







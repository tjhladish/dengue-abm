PROJ_HOME = /home/tjhladish/work/dengue/gis/extract_raster_data_by_region
#GDAL_LIB = $(PROJ_HOME)/gdal/lib
#GDAL_INC = $(PROJ_HOME)/gdal/include
GDAL_LIB = /usr/local/lib
GDAL_INC = /usr/local/include

default: clip_raster_by_region 
	export LD_LIBRARY_PATH=$(GDAL_LIB)

ogr_test: Makefile ogr_test.cpp
	g++ ogr_test.cpp -o ogr_test -I$(GDAL_INC) -L$(GDAL_LIB) -lgdal

gdal_test: Makefile gdal_test.cpp
	g++ gdal_test.cpp -o gdal_test -I$(GDAL_INC) -L$(GDAL_LIB) -lgdal

clip_raster_by_region: Makefile clip_raster_by_region.cpp
	g++ -O2 clip_raster_by_region.cpp -o extract_yucatan -I$(GDAL_INC) -L$(GDAL_LIB) -lgdal

pixel_failure_test: Makefile pixel_failure_test.cpp
	g++ pixel_failure_test.cpp -o pixel_failure_test -I$(GDAL_INC) -L$(GDAL_LIB) -lgdal

clean:
	rm -f test

# export LD_LIBRARY_PATH=$HOME/work/dengue/gis/eugene_work/gdal/lib


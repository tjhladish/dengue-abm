#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()
#include "math.h"
#include <cmath>
#include "ogrsf_frmts.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include "assert.h"
 
using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::stringstream;
using std::setprecision;
using std::vector;

struct Raster{
    GDALRasterBand* band;
    OGRGeometry* rasterBBox;
    double pixel_width;
    double pixel_height; 
    double raster_width; 
    double raster_height;
    double upper_left_x;
    double upper_left_y;
    double lower_right_x;
    double lower_right_y;
};

Raster* extract_raster_attributes( GDALDataset* poDataset ) {
    Raster* ra = new Raster;
    OGRGeometry* rasterBBox_tmp;
    double adfGeoTransform[6];
    if( poDataset->GetGeoTransform( adfGeoTransform ) == CE_None ) {
        ra->pixel_width  = fabs(adfGeoTransform[1]);
        ra->pixel_height = fabs(adfGeoTransform[5]);
        ra->raster_width  = poDataset->GetRasterXSize();
        ra->raster_height = poDataset->GetRasterYSize();
        ra->upper_left_x = adfGeoTransform[0];
        ra->upper_left_y = adfGeoTransform[3];
        ra->lower_right_x = adfGeoTransform[0] + ra->raster_width * adfGeoTransform[1];
        ra->lower_right_y = adfGeoTransform[3] + ra->raster_height * adfGeoTransform[5];
        fprintf(stderr, "Origin (upper left corner of upper left pixel) = (%.6f,%.6f)\n", ra->upper_left_x, ra->upper_left_y );
        fprintf(stderr, "Opposite corner (LR of LR) = (%.6f,%.6f)\n", ra->lower_right_x, ra->lower_right_y );
        fprintf(stderr, "Pixel Size = (%.6f,%.6f)\n", adfGeoTransform[1], adfGeoTransform[5] );
        
        stringstream ss;
        ss << "POLYGON((" << ra->upper_left_x  << " " << ra->upper_left_y  << ", "
                          << ra->lower_right_x << " " << ra->upper_left_y  << ", " 
                          << ra->lower_right_x << " " << ra->lower_right_y << ", " 
                          << ra->upper_left_x  << " " << ra->lower_right_y << ", "
                          << ra->upper_left_x  << " " << ra->upper_left_y  << "))";

        const string& wkt_s = ss.str();
        const char* wkt = wkt_s.c_str();
        // cast because OGR_G_CreateFromWkt will move the pointer 
        char* pszWkt = (char*) wkt;
        OGRSpatialReference ref = OGRSpatialReference(NULL);
        OGRErr err = OGRGeometryFactory::createFromWkt(&pszWkt, &ref, &rasterBBox_tmp);
    } else {
        cerr << "Could not extract geometry of raster data\n";
        exit(2);
    }

    ra->rasterBBox = (OGRPolygon*) rasterBBox_tmp;
    return ra;
}

Raster* import_raster(string raster_filename, int band_number) {
    GDALAllRegister();
    GDALDataset*  poDataset = (GDALDataset *) GDALOpen( raster_filename.c_str(), GA_ReadOnly );
    if( poDataset == NULL ) { 
        cerr << "Error: Could not open raster data file" << endl;
        exit(1); 
    }

    fprintf(stderr, "Driver: %s/%s\n", poDataset->GetDriver()->GetDescription(), poDataset->GetDriver()->GetMetadataItem( GDAL_DMD_LONGNAME ) );
    fprintf(stderr, "Size is %dx%dx%d\n", poDataset->GetRasterXSize(), poDataset->GetRasterYSize(), poDataset->GetRasterCount() );

    if( poDataset->GetProjectionRef()  != NULL ) cerr << "Projection is `" << poDataset->GetProjectionRef() << "'" << endl;;

    
    GDALRasterBand* poBand = poDataset->GetRasterBand( band_number );
    int nBlockXSize, nBlockYSize;
    poBand->GetBlockSize( &nBlockXSize, &nBlockYSize );
    fprintf(stderr, "Block=%dx%d Type=%s, ColorInterp=%s\n",
            nBlockXSize, nBlockYSize,
            GDALGetDataTypeName(poBand->GetRasterDataType()),
            GDALGetColorInterpretationName( poBand->GetColorInterpretation()) );

    Raster* raster = extract_raster_attributes( poDataset );
    raster->band = poBand;
    return raster;
}

void report_raster_data_within_polygon(Raster* raster, OGRPolygon* poPolygon, int feature_num, string feature_name) {
    if (raster->rasterBBox->Intersects( poPolygon ) ) {
        OGRPolygon *poIntersection = (OGRPolygon*) raster->rasterBBox->Intersection( poPolygon );
        OGREnvelope* envl = new OGREnvelope;
        poIntersection->getEnvelope(envl);
        //printf ("Bounding rect (minX, maxX; minY, maxY): %.3f, %.3f; %.3f, %.3f\n", envl->MinX, envl->MaxX, envl->MinY, envl->MaxY);
        const double pw = raster->pixel_width;
        const double ph = raster->pixel_height;
        int pixelOffset_X = round((envl->MinX - raster->upper_left_x) / pw);
        int pixelOffset_Y = round((raster->upper_left_y - envl->MaxY) / ph);
        int nCellsX = round((envl->MaxX - raster->upper_left_x) / pw) - pixelOffset_X;
        int nCellsY = round((raster->upper_left_y - envl->MinY) / ph) - pixelOffset_Y;
        OGRPoint* pixelPoint = new OGRPoint(0,0);

        float* pPixelValue = (float *) CPLMalloc(sizeof(float));

        cerr << "Feature size: " << nCellsX << "x" << nCellsY << endl;
        for ( int row = 0; row < nCellsY; ++row ) {
            for ( int col = 0; col < nCellsX; ++col ) {
                double xcenter = (col + pixelOffset_X + 0.5)*pw + raster->upper_left_x;
                double ycenter = raster->upper_left_y - (row + pixelOffset_Y + 0.5)*ph;

                pixelPoint->setX(xcenter);
                pixelPoint->setY(ycenter);

                if ( pixelPoint->Intersects(poPolygon)) {
                    raster->band->RasterIO( GF_Read, pixelOffset_X + col, pixelOffset_Y + row, 1,1,  pPixelValue, 1,1, GDT_Float32, 0, 0 );
                    cout << setprecision(10) << xcenter << "," << setprecision(10) <<  ycenter << "," <<  *pPixelValue << "," << feature_num << ","  << feature_name << endl;
                } 
            }
        }
        delete envl;
        CPLFree(pPixelValue);
    } else {
        cerr << "Feature does not intersect raster" << endl;
    }
    return;
}

OGRLayer* import_shapefile(string shapefile_name, int layer_number) {
    OGRRegisterAll();
    OGRDataSource *poDS;

    poDS = OGRSFDriverRegistrar::Open( shapefile_name.c_str(), FALSE );
    if( poDS == NULL ) { cerr << "Failed to open shapefile" << endl; exit( 1 ); }

    OGRLayer* poLayer = poDS->GetLayer(layer_number);
    cerr << "Shapefile layer name: " << poLayer->GetName() << endl;

//    char* proj[255];
//    poLayer->GetSpatialRef()->exportToWkt(proj);
//    cerr << "Shapefile projection: " << *proj << endl; 

    return poLayer;
}

int main() {
    // Read in raster data for night time lights
    int band_number = 1; // only one band, starts with one
    Raster* raster = import_raster("raster.tif", band_number);

    // Read in shapefile data containing municipality administrative regions
    int layer_number = 0; // only one layer, starts with zero
    OGRLayer* shapelayer = import_shapefile("MEX_adm2.shp", layer_number);

    shapelayer->SetAttributeFilter("ID_1 = 1834");       // Filter for Yucatan
    const int idx_of_number_field = 5;             // Column number of municipality number
    const int idx_of_name_field = 6;               // Column number of municipality name

    OGRFeature* poFeature;
    int feature_ctr = 0;
    while( (poFeature = shapelayer->GetNextFeature()) != NULL ) {
        cerr << "Feature: " << feature_ctr++ << "\t";
        int feature_num = poFeature->GetFieldAsInteger(idx_of_number_field);
        string feature_name = poFeature->GetFieldAsString(idx_of_name_field);

        OGRFeatureDefn *poFDefn = shapelayer->GetLayerDefn();
        for( int iField = 0; iField < poFDefn->GetFieldCount(); iField++ ) {
            OGRFieldDefn *poFieldDefn = poFDefn->GetFieldDefn( iField );
            if( poFieldDefn->GetType() == OFTString )  cerr << poFeature->GetFieldAsString(iField) << ",";
        }

        OGRGeometry* poGeometry = poFeature->GetGeometryRef();
        if( poGeometry != NULL) {
            // For contiguous regions
            if ( wkbFlatten(poGeometry->getGeometryType()) == wkbPolygon ) {
                cerr << " polygon" << endl;
                report_raster_data_within_polygon(raster, (OGRPolygon *) poGeometry, feature_num, feature_name);
            // For disjoint regions
            } else if ( wkbFlatten(poGeometry->getGeometryType()) == wkbMultiPolygon ) {
                cerr << " multipolygon" << endl;
                OGRMultiPolygon *multipolygon = (OGRMultiPolygon *) poGeometry;
                for (int i = 0; i<multipolygon->getNumGeometries(); i++) {
                    report_raster_data_within_polygon(raster, (OGRPolygon*) multipolygon->getGeometryRef(i), feature_num, feature_name);
                }
            // Is this really the right shapefile?
            } else {
                cerr << "No polygon or multipolygon geometry for this feature: " << poGeometry->getGeometryName() << endl;
            }
        } else {
            cerr << "No geometry for this feature" << endl;
        }
    }
    OGRFeature::DestroyFeature( poFeature );
}


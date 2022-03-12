// qh5_projlambert.cpp : This file contains the 'main' function. Program execution begins and ends there.
// 对极轨卫星数据，使用经纬度查找表投影变换到 兰伯特方位等积投影 缩写laea


#include <iostream>
#include "gdal_priv.h"
#include <vector>
#include "wGdalRaster.h"
#define WGDALRASTER_H_IMPLEMENTATION
#include "wGdalRaster.h"
#include <memory>
#include "ogr_spatialref.h"
#include <cassert>

#define RES_NEAR 0
#define RES_BIL  1

using namespace std;


struct DDPoint {
    inline DDPoint() :x(0), y(0) {}
    inline DDPoint(double x1, double y1) :x(x1), y(y1) {}
    double x, y;
};



/// <summary>
/// 四边形
/// </summary>
struct DDPointQuad {
    inline DDPointQuad() :qsize(0) {}
    int qsize;//total = qsize * qsize BOX
    vector<DDPoint> points;
    void appendPoints(vector<DDPoint>& pts);
};
void DDPointQuad::appendPoints(vector<DDPoint>& pts)
{
    for (int i = 0; i < pts.size(); ++i)points.push_back(pts[i]);
}


/// <summary>
/// qsize数量包含了起始点，中间点和结束点
/// </summary>
/// <param name="start"></param>
/// <param name="end"></param>
/// <param name="qsize"></param>
/// <returns></returns>
vector<DDPoint> buildEdge(DDPoint start, DDPoint end, int qsize)
{
    assert(qsize >= 2);
    vector<DDPoint> temp(qsize);
    temp[0] =start ;
    temp[qsize - 1] = end;
    double stepx = (end.x - start.x) / (qsize-1);
    double stepy = (end.y - start.y) / (qsize-1);
    for (int i = 1; i < qsize-1 ; ++i)
    {
        double midx = start.x + i * stepx;
        double midy = start.y + i * stepy;
        temp[i] = DDPoint(midx, midy);
    }
    return temp;
}

/// <summary>
/// p0-p3 have a clockwise or c-clockwise order.
/// </summary>
/// <param name="p0"></param>
/// <param name="p1"></param>
/// <param name="p2"></param>
/// <param name="p3"></param>
/// P0 ---- P1
/// |        |
/// |        |
/// P3 ---- P2
DDPointQuad buildQuad(DDPoint p0, DDPoint p1, DDPoint p2, DDPoint p3,int qsize)
{
    assert(qsize >= 2);
    DDPointQuad quad;
    quad.qsize = qsize;
    //四条边使用相同的等分，方便后续计算
    vector<DDPoint> edge03 = buildEdge(p0, p3, qsize);
    vector<DDPoint> edge12 = buildEdge(p1, p2, qsize);   
    for (int irow = 0; irow < qsize; ++irow)
    {
        vector<DDPoint> tempEdge = buildEdge(edge03[irow], edge12[irow], qsize);
        quad.appendPoints(tempEdge);
    }
    return quad;
}

void fillRasMultiBandData(wGdalRaster* outRas, int ix,int iy,vector<double> data)
{
    for (int ib = 0; ib < outRas->getNBand(); ++ib) {
        outRas->setValued(ix, iy, ib, data[ib]);
    }
}

vector<double> bilinearInterpolMultiBandData(
    int ixto0 , int iyto0, int qsize,
    vector<double>& data0,
    vector<double>& data1,
    vector<double>& data2,
    vector<double>& data3
    )
{
    vector<double> outdata(data0.size());
    double x0weight = 1.0 - (ixto0 + 1.0) / qsize;
    double y0weight = 1.0 - (iyto0 + 1.0) / qsize;

    for (int ib = 0; ib < data0.size(); ++ib)
    {
        double data01 = x0weight * data0[ib] + (1.0 - x0weight) * data1[ib];
        double data32 = x0weight * data3[ib] + (1.0 - x0weight) * data2[ib];
        double datafinal = data01 * y0weight + (1.0 - y0weight) * data32;
        outdata[ib] = datafinal;
    }
    return outdata;
}

/// <summary>
/// 填充到结果栅格数据
/// </summary>
/// <param name="quad"></param>
/// <param name="data0"></param>
/// <param name="data1"></param>
/// <param name="data2"></param>
/// <param name="data3"></param>
/// <param name="outRas"></param>
void fillOutputByQuad(DDPointQuad& quad,
    vector<double>& data0,
    vector<double>& data1,
    vector<double>& data2,
    vector<double>& data3,
    const int resample ,//RES_NEAR, RES_BIL
    const double outX0,const double outY1,
    double outReso,
    wGdalRaster* outRas
)
{
    for (int iy = 0; iy < quad.qsize; ++iy)
    {
        for (int ix = 0; ix < quad.qsize; ++ix)
        {
            DDPoint& pt = quad.points[iy * quad.qsize + ix];
            int ioutx = (pt.x - outX0) / outReso ;
            int iouty = (outY1 - pt.y) / outReso ;
            if (ioutx < 0 || ioutx >= outRas->getXSize() || iouty < 0 || iouty >= outRas->getYSize() ) {
                cout << "fillOutputByQuad Exception : " << ioutx << "," << iouty << endl;
                cout << "pt.x:" << pt.x << endl;
                cout << "pt.y:" << pt.y << endl;
                exit(20);
                return;
            }
            
            if( resample==RES_NEAR ){
                if (iy < quad.qsize / 2) {
                    if (ix < quad.qsize / 2) {
                        fillRasMultiBandData(outRas, ioutx, iouty, data0);
                    }
                    else {
                        fillRasMultiBandData(outRas, ioutx, iouty, data1);
                    }
                }
                else {
                    if (ix < quad.qsize / 2) {
                        fillRasMultiBandData(outRas, ioutx, iouty, data3);
                    }
                    else {
                        fillRasMultiBandData(outRas, ioutx, iouty, data2);
                    }
                }
            }
            else {
                //bilinear interpol
                vector<double> bildata = bilinearInterpolMultiBandData(ix, iy, quad.qsize, data0, data1, data2, data3);
                fillRasMultiBandData(outRas, ioutx, iouty, bildata);
            }
        }
    }

}


/// <summary>
/// dataRas必须是至少2x2的图像
/// </summary>
/// <param name="dataRas"></param>
/// <param name="dataScanRow0"></param>
/// <param name="dataScanRow1"></param>
/// <param name="lonRas"></param>
/// <param name="latRas"></param>
/// <param name="outX0"></param>
/// <param name="outY0"></param>
/// <param name="outReso"></param>
/// <param name="outRas"></param>
void processOneScan(wGdalRaster* dataRas, 
    int dataScanRow0, int dataScanRow1,
    wGdalRasterFloat& lonRas,
    wGdalRasterFloat& latRas,
    const double outX0,const double outY1, const double outReso,
    OGRCoordinateTransformation* pTrans,
    const int resample , //RES_NEAR, RES_BIL
    const double noDataValue ,//如果有无效值，那么该Quad自动使用最临近插值
    wGdalRaster* outRas)
{
    if (dataRas->getXSize() < 2 || dataRas->getYSize() < 2) {
        cout << "processOneScan exception: data size smaller 2x2" << endl;
        return;
    }
    const int nband = dataRas->getNBand();
    const int xsize = dataRas->getXSize();
    int per0 = -1;
    for (int irow = dataScanRow0; irow < dataScanRow1-1; ++irow)
    {
        for (int ix = 0; ix < xsize-1; ++ix) {
            double lonarr[4];
            double latarr[4];

            lonarr[0] = lonRas.getValuef(ix, irow, 0);
            lonarr[1] = lonRas.getValuef(ix + 1, irow, 0);
            lonarr[2] = lonRas.getValuef(ix + 1, irow + 1, 0);
            lonarr[3] = lonRas.getValuef(ix, irow + 1, 0);

            latarr[0] = latRas.getValuef(ix, irow, 0);
            latarr[1] = latRas.getValuef(ix + 1, irow, 0);
            latarr[2] = latRas.getValuef(ix + 1, irow + 1, 0);
            latarr[3] = latRas.getValuef(ix, irow + 1, 0);

            vector<double> dataPt0, dataPt1, dataPt2, dataPt3;
            bool hasNoData = false;
            int quadResample = resample;

            //获取Quad四个点的值，并检查是否有无效值，如果有的话，该quad使用最临近插值
            for (int ib = 0; ib < nband; ++ib) {
                dataPt0.push_back(dataRas->getValued(ix, irow, ib));
                dataPt1.push_back(dataRas->getValued(ix+1, irow, ib));
                dataPt2.push_back(dataRas->getValued(ix+1, irow+1, ib));
                dataPt3.push_back(dataRas->getValued(ix, irow+1, ib));
                if (dataPt0[ib] == noDataValue) hasNoData = true;
                if (dataPt1[ib] == noDataValue) hasNoData = true;
                if (dataPt2[ib] == noDataValue) hasNoData = true;
                if (dataPt3[ib] == noDataValue) hasNoData = true;
            }
            if (hasNoData == true) quadResample = RES_NEAR;

            //windows
            //pTrans->Transform(4, latarr, lonarr);//latarr-x[...] , lonarr-y[...]
            // double wsize0 = fabs(latarr[1] - latarr[0]);
            // double hsize0 = fabs(lonarr[1] - lonarr[0]);

            // double wsize1 = fabs(latarr[2] - latarr[1]);
            // double hsize1 = fabs(lonarr[2] - lonarr[1]);

            // double wsize2 = fabs(latarr[3] - latarr[2]);
            // double hsize2 = fabs(lonarr[3] - lonarr[2]);

            // double wsize3 = fabs(latarr[0] - latarr[3]);
            // double hsize3 = fabs(lonarr[0] - lonarr[3]);
            //windows
            // DDPointQuad quad = buildQuad(
            //     DDPoint(latarr[0], lonarr[0]),
            //     DDPoint(latarr[1], lonarr[1]),
            //     DDPoint(latarr[2], lonarr[2]),
            //     DDPoint(latarr[3], lonarr[3]),
            //     quadsize
            // );

            //linux
            pTrans->Transform(4, lonarr, latarr);
            double wsize0 = fabs(lonarr[1] - lonarr[0]);
            double hsize0 = fabs(latarr[1] - latarr[0]);

            double wsize1 = fabs(lonarr[2] - lonarr[1]);
            double hsize1 = fabs(latarr[2] - latarr[1]);

            double wsize2 = fabs(lonarr[3] - lonarr[2]);
            double hsize2 = fabs(latarr[3] - latarr[2]);

            double wsize3 = fabs(lonarr[0] - lonarr[3]);
            double hsize3 = fabs(latarr[0] - latarr[3]);

            double maxquadwid = max(wsize0, max(wsize1, max(wsize2, wsize3)));
            double maxquadhei = max(hsize0, max(hsize1, max(hsize2, hsize3)));

            int xcellnum = (maxquadwid / outReso) + 2;//这里+2强制缩小Quad的像素间隔小于输出分辨率一点点
            int ycellnum = (maxquadhei / outReso) + 2;

            int quadsize = max(2,max(xcellnum, ycellnum));
            //linux
            DDPointQuad quad = buildQuad(
                DDPoint(lonarr[0], latarr[0]),
                DDPoint(lonarr[1], latarr[1]),
                DDPoint(lonarr[2], latarr[2]),
                DDPoint(lonarr[3], latarr[3]),
                quadsize
            );

            
            fillOutputByQuad(quad, dataPt0, dataPt1, dataPt2, dataPt3, quadResample, outX0, outY1, outReso, outRas);
        }
        int per1 = (irow - dataScanRow0)*100.f / (dataScanRow1-dataScanRow0);
        if (per1 != per0) {
            per0 = per1;
            cout << per0 << " ";
        }
    }
    cout << " 100. " << endl;
}

void fillScanIntervalGap(wGdalRaster* ras, int ix, int iy, int ib, double nodatavalue)
{
    if (ix < 1 || iy < 1 || ix >= ras->getXSize() - 1 || iy >= ras->getYSize() - 1) {
        return;
    }
    //up down
    if( ras->getValued(ix,iy-1,ib) != nodatavalue && ras->getValued(ix,iy+1,ib) != nodatavalue)
    {
        ras->setValued(ix, iy, ib, ras->getValued(ix, iy - 1, ib)*0.5+ ras->getValued(ix, iy + 1, ib)*0.5 );
    }
    else if (ras->getValued(ix - 1, iy, ib) != nodatavalue && ras->getValued(ix + 1, iy, ib) != nodatavalue)
    {
        ras->setValued(ix, iy, ib, ras->getValued(ix-1, iy, ib) * 0.5 + ras->getValued(ix+1, iy, ib) * 0.5);
    }
}



int main(int argc ,char* argv[])
{
    std::cout << "A program to proj PolarScanData into lambert azi equal area projection. 2022-3-12" << endl;
    cout << "usage : qh5_projlambert lonfile latfile datafile nrowsPerScan outResoMeter inOutFilldata res outputfile" << endl;
    cout << "nrowsPerScan : number of rows per scan, if do not known use 0 for using full height as one scan." << endl;
    cout << "outResoMeter : output resolution in meter." << endl;
    cout << "inOutFilldata : in and out filldata." << endl;
    cout << "res : resample method near/bil ." << endl;
    cout << "lon/lat file : should be float32 data." << endl;
    cout << "v1.0.0 created 2022-3-12" << endl;
    cout << "v1.1.0 working version in linxu 2022-3-12." << endl;


    GDALAllRegister();

    if (argc != 9) {
        cout << "argc not 9." << endl;
        return 11;
    }

    string lonfile = argv[1];
    string latfile = argv[2];
    string datafile = argv[3];
    int numRowsPerScan = atof(argv[4]);
    double outReso = atof(argv[5]);
    double filldata = atof(argv[6]);
    string resample = argv[7];
    string outfile = argv[8];

    cout << "--------------------------------------- " << endl;
    cout << "lonfile " << lonfile << endl;
    cout << "latfile " << latfile << endl;
    cout << "datafile " << datafile << endl;
    cout << "numRowsPerScan " << numRowsPerScan << endl;
    cout << "outReso " << outReso << endl;
    cout << "filldata " << filldata << endl;
    cout << "resample " << resample << endl;
    cout << "outfile " << outfile << endl;
    cout << "--------------------------------------- " << endl;

    int ires = RES_NEAR;
    if (resample.compare("bil") == 0) ires = RES_BIL;

    if (outReso < 1) {
        cout << "bad outReso , too small for meter(1.0)." << endl;
        return 12;
    }

    wGdalRasterFloat lonRas, latRas;
    bool ok = lonRas.open(lonfile);
    if (ok == false) {
        cout << "bad lonfile" << endl;
        return 12;
    }

    ok = latRas.open(latfile);
    if (ok == false) {
        cout << "bad latfile" << endl;
        return 13;
    }

    shared_ptr<wGdalRaster> dataRasPtr(wGdalRasterFactory::OpenFile(datafile));
    if (dataRasPtr.get() == 0) {
        cout << "bad datafile" << endl; 
        return 14;
    }

    int xsize = dataRasPtr->getXSize();
    int ysize = dataRasPtr->getYSize();
    int nband = dataRasPtr->getNBand();

    if (xsize < 2 || ysize < 2) {
        cout << "bad data size " << xsize << "," << ysize << endl;
        return 14;
    }

    if (xsize != lonRas.getXSize() || xsize != latRas.getXSize() ||
        ysize != lonRas.getYSize() || ysize != lonRas.getYSize())
    {
        cout << "lon, lat, data raster sizes are not same." << endl;
        return 15;
    }


    //compute lon/lat center
    double centerlon = lonRas.getValuef(xsize / 2, ysize / 2, 0);
    double centerlat = latRas.getValuef(xsize / 2, ysize / 2, 0);

    double fourCornerLon[4];//00,01,10,11
    double fourCornerLat[4];
    {
        fourCornerLon[0] = lonRas.getValuef(0, 0, 0);
        fourCornerLon[1] = lonRas.getValuef(0, ysize-1, 0);
        fourCornerLon[2] = lonRas.getValuef(xsize-1, 0, 0);
        fourCornerLon[3] = lonRas.getValuef(xsize-1, ysize-1, 0);

        fourCornerLat[0] = latRas.getValuef(0, 0, 0);
        fourCornerLat[1] = latRas.getValuef(0, ysize-1, 0);
        fourCornerLat[2] = latRas.getValuef(xsize-1, 0, 0);
        fourCornerLat[3] = latRas.getValuef(xsize-1, ysize-1, 0);
    }

    OGRSpatialReference wgs84;
    wgs84.SetWellKnownGeogCS("WGS84");
    char laeaProj4Str[1024];
    sprintf(laeaProj4Str, "+proj=laea +lon_0=%.2f +lat_0=%.2f +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs", centerlon, centerlat);
    cout << "laeaProj4Str:" << laeaProj4Str << endl;
    OGRSpatialReference lambertRef ;
    lambertRef.importFromProj4(laeaProj4Str);
    lambertRef.Validate();
    OGRCoordinateTransformation* pConverter = OGRCreateCoordinateTransformation(&wgs84,
        &lambertRef);
    {
        //windows
        //double x = centerlat;//注意这里应该是反的，也就是说x对lat，y对lon
        //double y = centerlon;

        //linux 很奇怪，这里linux的gdal库和windows的Qgis GDAL库是反的
        double x = centerlon;
        double y = centerlat;
        pConverter->Transform(1, &x, &y);
        cout << "lon,lat:" << centerlon << "," << centerlat << endl;
        cout << "x,y:" << x << "," << y << endl;
    }

    double minLambertX(0), maxLambertX(0);
    double minLambertY(0), maxLambertY(0);
    double fourCornerLambertX[4];
    double fourCornerLambertY[4];
    for (int iii = 0; iii < 4; ++iii) {
        //win
        //double tx = fourCornerLat[iii];
        //double ty = fourCornerLon[iii];

        //linux
        double tx = fourCornerLon[iii];
        double ty = fourCornerLat[iii];

        pConverter->Transform(1, &tx, &ty);
        fourCornerLambertX[iii] = tx;
        fourCornerLambertY[iii] = ty;
        if (iii == 0) {
            minLambertX = maxLambertX = tx;
            minLambertY = maxLambertY = ty;
        }
        else {
            minLambertX = min(minLambertX, tx);
            maxLambertX = max(maxLambertX, tx);
            minLambertY = min(minLambertY, ty);
            maxLambertY = max(maxLambertY, ty);
        }
        cout << "lon,lat:" << fourCornerLon[iii] << "," << fourCornerLat[iii] << endl;
        cout << "lambert x,y:" << tx << "," << ty << endl;
    }
    cout << "lambert minx,maxx:" << minLambertX << "," << maxLambertX << endl;
    cout << "lambert miny,maxy:" << minLambertY << "," << maxLambertY << endl;

    minLambertX = minLambertX - outReso / 2; 
    maxLambertY = maxLambertY + outReso / 2;
    cout << "ajust lambert min x:" << minLambertX << endl;
    cout << "ajust lambert max y:" << maxLambertY << endl;

    int outxsize = (maxLambertX - minLambertX) / outReso+1.5;
    int outysize = (maxLambertY - minLambertY) / outReso+1.5;

    cout << "outsize:" << outxsize << "," << outysize << endl;
    if (outxsize > 30000 || outysize > 30000) {
        cout << "can not support too large outsize " << outxsize << "," << outysize << endl;
        return 15;
    }
    if (outxsize <= 0 || outysize <= 0) {
        cout << "bad outsize " << outxsize << "," << outysize << endl;
        return 16;
    }

    shared_ptr<wGdalRaster> outRas(wGdalRasterFactory::Create(outxsize, outysize, nband, dataRasPtr->getDataType()));
    if (outRas.get() == 0) {
        cout << "bad outRas create." << endl;
        return 17;
    }
    //fill no data value
    for (int ib = 0; ib < nband; ++ib) {
        outRas->fill(ib, filldata);
    }

    if (numRowsPerScan == 0 || numRowsPerScan> ysize) numRowsPerScan = ysize;
    if (numRowsPerScan == 1) numRowsPerScan = 2;
    cout << "adjusted numRowsPerScan:" << numRowsPerScan << endl;



    const int numScan = ceil(ysize*1.0 / numRowsPerScan) ;
    for (int iscan = 0; iscan < numScan; ++iscan)
    {
        int scanRow0 = iscan * numRowsPerScan;
        int scanRow1 = min(ysize,scanRow0 + numRowsPerScan) ;
        cout << "processing for iscan:" << iscan << " ; row range:" << scanRow0 << "-" << scanRow1 << endl;
        processOneScan(dataRasPtr.get(), scanRow0, scanRow1, lonRas, latRas, minLambertX, maxLambertY, outReso,
            pConverter, ires, filldata, outRas.get() );
        cout << "scan done." << endl;
    }
    cout << endl<<"all scan done." << endl;

    //填充由于四舍五入造成的Scan间空隙，这类空隙的特点是只有一个像素
    cout << "fill gap between scans..." << endl;
    for (int iy = 1; iy < outysize-1; ++iy) {
        for (int ix = 1; ix < outxsize-1; ++ix) {
            if (outRas->getValued(ix, iy, 0) == filldata) {
                //如果上下有值就上下取平均，否则左右取平均，再没有跳出
                for (int ib = 0; ib < nband; ++ib) {
                    fillScanIntervalGap(outRas.get(), ix, iy, ib, filldata);
                }
            }
        }
    }
    cout << "fill gap done." << endl;
    
    double outTrans[6];
    outTrans[0] = minLambertX ;
    outTrans[1] = outReso;
    outTrans[2] = 0.0;
    outTrans[3] = maxLambertY ;
    outTrans[4] = 0;
    outTrans[5] = -outReso;

    {
        char* tempWktBuff = 0;
        lambertRef.exportToWkt(&tempWktBuff);
        outRas->copyTrans(outTrans);
        outRas->copyProj(tempWktBuff);
        CPLFree(tempWktBuff);
    }
    
    OCTDestroyCoordinateTransformation(pConverter);
    
    cout << "saving..." << endl;
    bool saveok = outRas->save(outfile);
    if (saveok == false) {
        cout << "failed to save " << outfile << endl;
        return 20;
    }
    cout << "all_done" << endl;

    return 0;
}



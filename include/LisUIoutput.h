/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 1992, 2003, 2016, 2024  Victor Jetten
**  contact: v.g.jetten AD utwente DOT nl
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with this program. If not, see <http://www.gnu.org/licenses/>.
**
**  Authors: Victor Jetten, Bastian van de Bout, Meindert Commelin
**  Developed in: MingW/Qt/, GDAL, PCRaster
**  website, information and code: https://github.com/vjetten/openlisem
**
*************************************************************************/

/*!
  \file LisUIoutput.h
  \brief structure to pass variables form the model to the interface, visible by both
  */

#include <CsfMap.h>
#include <CsfRGBMap.h>
#include <QList>

#ifndef LISUIOUTPUT_H_
#define LISUIOUTPUT_H_

/// structure to pass variables form the model to the interface.
/// This tsructure is the link, visible by both

//typedef struct LDD_COOR {
//    int r;
//    int c;
//}  LDD_COOR;

struct output{
    int nrRunsDone; // nr runs without closing interface, needed to destroyd old data before start of a new run
    int runstep;
    int printstep;
    int maxstep;
    int cores;

    QList<int> OutletIndices;
    QList<int> OutletLocationX;
    QList<int> OutletLocationY;
    QList<QVector<double>*> OutletQ;
    QList<QVector<double>*> Wavein;
    QList<QVector<double>*> OutletQs;  //current kg/s
    QList<QVector<double>*> OutletC;   // avg concetration
    QList<QVector<double>*> OutletChannelWH;
    QVector<double> OutletQpeak;
    QVector<double> OutletQpeaktime;
    QVector<double> OutletQtot;
    QVector<double> OutletQstot;  // sum in kg
    QVector<double> Pmm;
    QVector<double> Time;
    QVector <double> Qtile;

 //   QVector <double> CulvertX;
 //   QVector <double> CulvertY;
    QVector <double> EndPointX;
    QVector <double> EndPointY;
    QVector <double> ObsPointX;
    QVector <double> ObsPointY;
    QVector <LDD_COORIN> lddch_;

    double timestep, CatchmentArea, t,time, maxtime, EndTime, BeginTime;
    double _llx, _lly, _dx;
    int _nrRows, _nrCols;

    double
    // water
    MB, Qtot,  Qtiletot, RunoffFraction, RainpeakTime, Rainpeak,
    Qtotmm,  IntercTotmm, IntercHouseTotmm, WaterVolTotmm,InfilTotmm,StormDrainTotmm, Qboundtotmm,
    RainTotmm, ETaTotmm, SurfStormm, InfilKWTotmm,  IntercLitterTotmm, //WaterVolTotchannelmm,
    floodBoundaryTot, floodBoundarySedTot, Theta1, Theta2, GWlevel, BaseFlowTotmm, PeakFlowTotmm,
    // channel
    ChannelVolTotmm, ChannelSedTot, ChannelDepTot, ChannelDetTot, ChannelWH,
    // flood
    FloodTotMax, FloodAreaMax, FloodArea, WHflood, Qflood, volFloodmm,
    FloodDetTot, FloodDepTot, FloodSedTot,
    // sediment
    MBs, DetTot, DetTotSplash, DetTotFlow, DepTot, SoilLossTot, SedTot, maxRainaxis;

    // map pointers for display
    cTMap *baseMap;
    cTMap *baseMapDEM;
    cTMap *channelMap;
    cTMap *outletMap;
    cTMap *roadMap;
    cTMap *houseMap;
    cTMap *hardsurfaceMap;
    cTRGBMap *Image;

    QList<double> graindiameters;

    QList<int> ComboLists;
    QList<cTMap *> ComboMaps;
    QList<cTMap *> ComboMapsSafe;
    QList<QList<double>> ComboColorMap;
    QList<QList<QString>> ComboColors;
    QList<bool> ComboLogaritmic;
    QList<bool> ComboSymColor;
    QStringList ComboMapNames;
    QStringList ComboUnits;
    QList<double> ComboScaling;
    QList<double> userMinV;
    QList<double> userMaxV;
    QList<double> comboStep;

    bool comboboxset;
    bool has_image;
    bool SwitchCorrectMB_WH;

    QString runfilename;
    QString userAppDir;
    QString format;
    QString timeStartRun;
    QString datestamp;

    bool doBatchmode;
  //  bool hasrunonce;
  //  int nrMapsCreated;
};



#endif /* LISUIOUTPUT_H_ */

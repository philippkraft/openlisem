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
  \file lisModel.cpp
  \brief Central model file with the main loop. From here all processes are called.

  functions: \n
  - void TWorld::DoModel() the main model function with the timeloop. It is a 'slot' linked to a signal.\n
  - void TWorld::run() Run is called from the interface to activate DoModel() \n
  - void TWorld::stop() Stops the loop on user request.\n

*/

#include <QtGui>
#include "lisemqt.h"
#include "model.h"
#include "global.h"


//---------------------------------------------------------------------------
TWorld::TWorld(QObject *parent) :
    QThread(parent)
{
   moveToThread(this);
}
//---------------------------------------------------------------------------
TWorld::~TWorld()
{
}
//---------------------------------------------------------------------------
void TWorld::run()
{
    QTimer::singleShot(0, this, SLOT(DoModel()));
    exec();
}
//---------------------------------------------------------------------------
void TWorld::stop()
{
    QMutexLocker locker(&mutex);
    stopRequested = true;
}
//---------------------------------------------------------------------------
void TWorld::saveMBerror2file( bool start) //bool doError,
{
    if (start) {
        //create error file
        QFile efout(resultDir+errorFileName);
        efout.open(QIODevice::WriteOnly | QIODevice::Text);
        QTextStream eout(&efout);
        eout << "#mass balance error (%)\n";

        efout.flush();
        efout.close();
    }


 //   if (doError) {
        QFile efout(resultDir+errorFileName);
        efout.open(QIODevice::Append | QIODevice::Text);
        QTextStream eout(&efout);
        eout << " " << runstep << "," << MB << "," << (SwitchErosion ? MBs : 0.0) << "\n";
        efout.flush();
        efout.close();
  //  }

}
//---------------------------------------------------------------------------
// the actual model with the main loop
void TWorld::DoModel()
{
    if (!op.doBatchmode)
        temprunname = QString(op.userAppDir+"openlisemtmp.run");
    else
        temprunname = op.runfilename;

    mapFormat = "PCRaster";

    errorFileName = QString(resultDir + "error-"+ op.timeStartRun +".csv");
    //errorSedFileName = QString(resultDir + "errorsed-"+ op.timeStartRun +".txt");
    time_ms.start();
    // get time to calc run length
    startTime=omp_get_wtime()/60.0;

    try
    {
        DestroyData();

        DEBUG("reading and initializing data");

        IntializeOptions(); // reset all options

        DEBUG("GetRunFile()");
        GetRunFile();
        DEBUG("Parse Runfile");
        ParseRunfileData();
        // get and parse runfile

        QString S = resultDir + QFileInfo(op.runfilename).fileName();
        QFile::copy(op.runfilename, S);


        //time vraiables in sec
        double btd = getvaluedouble("Begin time day");
        double btm = getvaluedouble("Begin time");
        double etd = getvaluedouble("End time day");
        double etm = getvaluedouble("End time");

        btd -= 1.0; // because day 1, minute 10 is in fact minute 10 in the first day
        etd -= 1.0;

        if (SwitchEventbased) {
            DEBUG("Day in start and end time is ignored.");
        }

        _dt = getvaluedouble("Timestep");

        if (SwitchEventbased) {
            BeginTime = (btm)*60; //for running in sec
            EndTime = (etm)*60;   //in sec
            op.BeginTime = BeginTime/60; // for graph drawing in min
            op.EndTime = EndTime/60;
        } else {
            BeginTime = (btd*1440+btm)*60; //for running in sec
            EndTime = (etd*1440+etm)*60;   //in sec
            op.BeginTime = BeginTime/60;// for graph drawing in min
            op.EndTime = EndTime/60;
        }

        //get all maps
        DEBUG("Get Input Maps");
        GetInputData();
        DEBUG("Intialize Database");
        IntializeData();

        // MC - no_ui probalbly this can be skipped for noInterface??
        setupDisplayMaps();
        // reset all display output maps for new job
        // must be done after Initialize Data because then we know how large the map is
        // clear() calls the destruction of all elements in the sturcture

        if (SwitchRainfall)
        {
            RainfallSeries.clear();
            RainfallSeriesMaps.clear();
            raintime.clear();

            DEBUG("Get Rainfall Data");
            if (SwitchRainfallSatellite) {
                GetSpatialMeteoData(rainSatFileName, 0);
            } else {
                GetRainfallStationData(rainFileName);
            }
        }

        if (SwitchIncludeET)
        {
            ETSeries.clear();
            ETSeriesMaps.clear();
            ETtime.clear();

            DEBUG("Get EvapoTranspiaration Data");
            if (SwitchETSatellite) {
                GetSpatialMeteoData(ETSatFileName, 1);
            } else {
                GetETStationData(ETFileName);
            }
        }

        // SwitchSnowmelt = false;
        // if (SwitchSnowmelt)
        // {
        //     SnowmeltSeries.clear();
        //     SnowmeltSeriesMaps.clear();
        //     snowmelttime.clear();
        //
        //     DEBUG("Get Snowmelt Data Information");
        //     if (SwitchSnowmeltSatellite) {
        //         GetSpatialMeteoData(snowmeltSatFileName, 2);
        //     } else {
        //         GetSnowmeltData(snowmeltFileName);
        //     }
        // }

        if (SwitchDischargeUser)
        {
            DischargeSeries.clear();
            dischargetime.clear();

            DEBUG("GetUserDischargeData()");
            GetUserDischargeData(dischargeinFileName);
        }

        if (SwitchWaveUser)
        {
            WHSeries.clear();
            WHtime.clear();

            DEBUG("GetWHboundaryData()");
            GetWHboundaryData(WaveinFileName);
        }

        // get all input data and create and initialize all maps and variables

        CountLandunits();
        //VJ 110110 for output totals per landunit

        runstep = 0; //  runstep is used to initialize graph!
        printstep = 1; // printstep determines report frequency in report()

      //  DEBUG("setupHydrographData()");
        setupHydrographData(); // reset hydrograph display

        //bool saveMBerror = true;
        //saveMBerror2file(true); //saveMBerror,

        SetFlowBarriers();     // update the presence of flow barriers, static for now, unless breakthrough
        GridCell();            // static for now

        _dt_user = _dt;

        DEBUG(" ");

        GetComboMaps(); // moved to outside timeloop!

        if (SwitchInfiltration && InfilMethod != INFIL_SWATRE )
            InfilEffectiveKsat();

        for (time = BeginTime; time < EndTime; time += _dt)
        {            
            // printstep determines report frequency in #define report(...)
            if (runstep > 0 && runstep % printinterval == 0)
                printstep++;

            runstep++;

            if(stopRequested) {
                mutex.lock();
                DEBUG("User interrupt... finishing time step");
                mutex.unlock();
            }

            if (waitRequested) {
                mutex.lock();
                DEBUG("User pause...");
                condition.wait(&mutex);
                mutex.unlock();
            }
            // check if user wants to quit or pause

            GetInputTimeseries(); // get rainfall, ET, snowmelt, discharge

            if (SwitchInfilCrust)
                InfilEffectiveKsat(); // if crusting recalc Ksateff and Poreff becuase of crusting effect

            HydrologyProcesses();  // hydrological processes in one loop, incl splash

            OverlandFlow(); // overland flow 1D (non threaded), 2Ddyn (threaded), if 2Ddyn then also SWOFsediment!

            // these are all non-threaded
            ChannelFlowandErosion();    // do ordered LDD solutions channel, tiles, drains, non threaded

            TileFlow();          // tile drain flow kin wave
                                 // storm drain flow kin wave
            //StormDrainFlow();

            TotalsHydro();       // calculate all totals and cumulative values
            TotalsFlow();
            TotalsSediment();

            MassBalance();       // check water and sed mass balance

            OutputUI();          // fill the "op" structure for screen output and calc some output maps

            reportAll();         // report maps and files to screen and disk

            emit show(noInterface); // send the 'op' structure with data to function worldShow in LisUIModel.cpp

            //saveMBerror2file(false); //saveMBerror

            if(stopRequested)
                time = EndTime;

            // show progress in console without GUI
            if (op.doBatchmode) {
                int x;
                x = std::round((op.t / op.maxtime) * 100) ;
                printf("\rprogress: %d %%                     ", x);
                // or use qDebug()
            }
             // MC - maybe not the most sophisticated solution but noInterface works again
        }

        if (SwitchEndRun)
            ReportMaps();

        emit done("finished");

        if (op.doBatchmode)
        {
            // delete all maps
            qDeleteAll(maplistCTMap.begin(),maplistCTMap.end());
            maplistCTMap.clear();

            //delete swatre 3D soil layer structure if exists
            if (initSwatreStructure)
                FreeSwatreInfo();

            qDebug() << "\nfinished after "<< op.maxtime << "minutes\n";
            if (noInterface)
                QCoreApplication::quit();
            else
                QApplication::quit();
            // close the world model
        }
    }
    catch(...)  // if an error occurred
    {
        emit done("ERROR STOP: "+ErrorString);
        if (op.doBatchmode) {qDebug() << "ERROR STOP "<< ErrorString;
            if (noInterface)
                QCoreApplication::quit();
            else
                QApplication::quit();
        }
    }
}
//---------------------------------------------------------------------------
void TWorld::GetInputTimeseries()
{
    // get meteo data
    if (SwitchRainfallSatellite)
        GetRainfallMapfromSat(time);         // get rainfall from maps
    else
        GetRainfallMapfromStations(time);  // get rainfall from stations

    if (SwitchIncludeET) {
        if (SwitchETSatellite)
            GetETSatMap(time); // get rainfall from maps
        else
            GetETMapfromStations(time);   // get rainfall from stations
    }

    if (SwitchDischargeUser) {
        GetDischargeMapfromStations(time);
    }

    if (SwitchWaveUser) {
        GetWHboundaryMap(time);
    }

//    if (SwitchSnowmelt) {
//        if (SwitchSnowmeltSatellite)
//            ; //TODO snowmelt satellite
//        else
//            GetSnowmeltMap(time);  // get snowmelt from stations
//    }

}
//---------------------------------------------------------------------------
// all hydrologuical processes in one big parallel loop for speed
void TWorld::HydrologyProcesses()
{
   // double soiltot1 = SoilWaterMass();

    if (SwitchIncludeET) {
        if (SwitchDailyET)
            ETafactor = getETaFactor();
        else
            ETafactor = 1.0;
    }

    // Do all hydrology in one big loop. Not sure if this is faster then a loop per process!
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        cell_Interception(r,c);
        // all interception on plants, houses, litter
        // result is rainnet (and leafdrip for erosion)

//        if (SwitchFloodInitial  && hmxInit->Drc > 0)
//            hmxInit->Drc += RainNet->Drc;

        if (FloodDomain->Drc > 0) {
            hmx->Drc += RainNet->Drc;// + Snowmeltc->Drc; // only used in kin wave pluf flood from channel, hmx is flood water
        } else {
            WH->Drc += RainNet->Drc;// + Snowmeltc->Drc;  // used in 2D flow and kin wave
        }

        if (SwitchWaveUser) {
            WHboundRain->Drc += RainNet->Drc;
            if (WHboundarea->Drc > 0) {
                // WHbound is the forced water level in area with value '1', ples cum rainfall
                WH->Drc = WHbound->Drc + WHboundRain->Drc;
            }
        }
        // if(std::isnan(Thetaeff->Drc))
        //     qDebug() << QString("A nan 1 %1 %2").arg(r).arg(c);

        if (SwitchInfiltration  && InfilMethod != INFIL_SWATRE) {
            switch (InfilMethod) {
                case INFIL_SOAP : cell_Soilwater(i_); break;
                case INFIL_GREENAMPT:
                case INFIL_SMITH:
                    // Green and Ampt + redistribution
                    cell_InfilMethods(r, c);

                    if (SwitchIncludeET)
                        cell_ETa(r,c);

                    if (SwitchTwoLayer) {
                        cell_Redistribution2(r, c);
                        //cell_Channelinfow2(r, c);
                    } else {
                        cell_Redistribution1(r, c);
                        //cell_Channelinfow1(r, c);
                    }

                    if (!SwitchImpermeable)
                        Perc->Drc = cell_Percolation(r, c, 1.0);

                    break;
            }
        }
    }}


    if (SwitchInfiltration && InfilMethod == INFIL_SWATRE) {
        // #pragma omp parallel for num_threads(userCores)
        // FOR_ROW_COL_MV_L {
        //     cell_InfilSwatre(i_, r,c);
        // }}
        InfilSwatre(); // does the same
    }


    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        // do not do this!
        //  cell_depositInfil(r,c);
        // deposit all sediment still in flow when infiltration causes WH to become minimum
        // gives huge MBs errors!

        cell_SurfaceStorage(r, c);
        //calc surf storage and total watervol and WHrunoff

        if (SwitchErosion)
            cell_SplashDetachment(r, c);

        if (SwitchSlopeStability)
            cell_SlopeStability(r, c);
    }}

    //MoistureContent();
    // double soiltot2 = SoilWaterMass();
    // if (InfilMethod != INFIL_SOAP)
    //     SoilMoistDiff = soiltot2 - soiltot1;

}
//---------------------------------------------------------------------------



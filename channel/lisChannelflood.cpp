/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 1992, 2003, 2016, 2024  Victor Jetten
**  contact:
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
 \file lisChannelflood.cpp
 \brief Channel flood using a various solutions of St Venant equations: \n
        and more stable 1st and 2nd order st Venant following the fullSWOF2D code (univ Orleans)\n
        called before ChannelFlow(), takes old channel overflow height and spreads it out, puts new channelWH \n
        back into channel before kin wave of channel is done in ChannelFlow()
        
functions: \n
- void TWorld::ChannelOverflow(void) Mixing of flood and overflow in channel cells, source of overflow
- void TWorld::ChannelFlood(void) Calculate channelflood height maps (hmx, QFlood, UVFlood) and FloodDomain
*/

#include <algorithm>
#include "lisemqt.h"
#include "model.h"
#include "operation.h"
#include "global.h"

#define HMIN3 0.001

void TWorld::ChannelOverflow(cTMap *_h, cTMap *V)
{
    if (!SwitchIncludeChannel)
         return;

     //#pragma omp parallel for num_threads(userCores)
     FOR_ROW_COL_MV_CHL {
         double maxQ = SwitchCulverts ? ChannelMaxQ->Drc  : -1;

         if (ChannelWidth->Drc > 0 && maxQ <= 0)
         {
             double chdepth = ChannelDepth->Drc;
             double dH = std::max(0.0, (ChannelWH->Drc-chdepth));

             if (dH <= HMIN3 && _h->Drc <= HMIN3)
                 continue;
             // no flow activity then continue

             if (fabs(dH - _h->Drc) < HMIN3)
                 continue;
             // no diff in water level, no flow, continue

             // VELOCITIES
             double VtoChan = V->Drc;
             double fracA = std::min(1.0, _dt*VtoChan/(0.5*ChannelAdj->Drc));
             // fraction from _h to channel based on average flood velocity
             double VfromChan = sqrt(2*9.804*dH); //Bernoulli
             //see https://www.engineeringtoolbox.com/velocity-head-d_916.html
             double fracC = std::min(1.0, _dt*VfromChan/(0.5*ChannelAdj->Drc));
             // fraction from channel to surrounding

             double cwa = ChannelWidth->Drc/ChannelAdj->Drc;

             bool dosimpel = false;

             if (dH > _h->Drc)   // flow from channel
             {
                 double dwh = fracC * (dH-_h->Drc);
                 // amount flowing from channel
                 if (_h->Drc + dwh*cwa > dH-dwh) {
                     // if flow causes situation to reverse (channel dips below _h)
                     dosimpel = true;
                 } else {

                     _h->Drc += dwh*cwa;
                     ChannelWH->Drc -= dwh;

                     if(SwitchErosion) {
                         double sed = ChannelSSConc->Drc * dwh*ChannelWidth->Drc*ChannelDX->Drc;
                         ChannelSSSed->Drc -= sed;
                         SSFlood->Drc += sed;
                     }
                 }
             }
             else   // flow to channel
             {
                 double dwh = fracA * (_h->Drc-dH);
                 // amount flowing to channel
                 if (dH + dwh/cwa > _h->Drc-dwh) {
                     // if too much flow
                     dosimpel = true;
                 } else {
                     _h->Drc -= dwh;
                     ChannelWH->Drc += (dwh/cwa);
                     if(SwitchErosion) {
                         double sed = fracA*SSFlood->Drc;
                         ChannelSSSed->Drc += sed;
                         SSFlood->Drc -= sed;
                     }
                 }
             }

             // instantaneous waterlevel exquilibrium acccross channel and adjacent
             if (dosimpel)
             {
                 double fc = ChannelWidth->Drc/_dx;
                 // fraction of the channel in the gridcell, 1-fc = (dx-chw)/dx = chanadj/dx
                 double whlevel = (ChannelWH->Drc-chdepth)*fc + _h->Drc*(1-fc);
                 // equilibrium water level = weighed values of channel surplus level + _h
                 // can be negative if channelwh is below channel depth and low _h level
                 if(whlevel > 0)
                 {
                     double sedch = 0;
                     double sed = 0;
                     if (SwitchErosion) {
                        sedch = ChannelSSSed->Drc;
                        sed = SSFlood->Drc;
                     }
                     double oldchwh = ChannelWH->Drc;
                     double oldwh = _h->Drc;
                     ChannelWH->Drc = whlevel + chdepth;
                     _h->Drc = whlevel;

                     // new equilibrium levels
                     if (SwitchErosion) {
                        // double sed_ = SSFlood->Drc + ChannelSSSed->Drc;
                         if (oldchwh > ChannelWH->Drc) {
                             double sed = (oldchwh-ChannelWH->Drc)*ChannelWidth->Drc*ChannelDX->Drc * ChannelSSConc->Drc;
                             ChannelSSSed->Drc -= sed;
                             SSFlood->Drc += sed;
                         } else {
                             double sed = (oldwh-_h->Drc)*CHAdjDX->Drc * SSCFlood->Drc;
                             SSFlood->Drc -=sed;
                             ChannelSSSed->Drc += sed;
                         }
                     }

                 }
                 else
                 {
                     //DO NOTHING
                     // this happens if there is very little flood water (< 5cm) and the channelWH is below the channeldepth
                     // we assume that there is no more flow towards the channel.
                 }
             }
         }
     }}
}

//---------------------------------------------------------------------------
//! Get flood level in channel from 1D kin wave channel
//! Instantaneous mixing of flood water and channel water in channel cells
//! note: ChannelDepth lets you also control which channels flood:
//! those that are 0 react as usual (infinite capacity)

void TWorld::ChannelOverflowIteration(cTMap *_h, cTMap *V)
{
    if (!SwitchIncludeChannel)
        return;

    double nrsteps = 0;
    bool go = false;
    Fill(*tma, 0);
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_CHL {
        double maxQ = SwitchCulverts ? ChannelMaxQ->Drc  : -1;
        if (ChannelWidth->Drc > 0 && maxQ <= 0)
        {
            double dH = std::max(0.0, (ChannelWH->Drc-ChannelDepth->Drc));

            if (dH <= HMIN3 && _h->Drc <= HMIN3)
                continue;
            // no flow activity then continue

            if (fabs(dH - _h->Drc) < HMIN3)
                continue;
            // no diff in water level, no flow, continue

            tma->Drc = 1;
            go = true;
            nrsteps += 1.0;
            double Vb = sqrt(2*GRAV*dH);

            double Vavg;
            //if (dH*Vb > _h->Drc*V->Drc)
            if (dH > V->Drc)
                Vavg = Vb; //Bernouilli
            else
                Vavg = V->Drc;
            // V from channel or reverse

            nrsteps = qMax(nrsteps,  qMax(1.0, (Vavg * _dt/(0.5*ChannelAdj->Drc))));
        }
    }}

    if (!go)
        return;

    // if every cell has its own step the result is an unstable hydrograph
    int step = (int)sqrt(nrsteps); // limit nr of steps, solution is fine anyway

    step = std::min(10,step);
    //qDebug() << step;

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_CHL {
        if (tma->Drc > 0)//ChannelWidth->Drc > 0 && ChannelMaxQ->Drc <= 0)
        {
            double chdepth = ChannelDepth->Drc;
            double dH = std::max(0.0, (ChannelWH->Drc-chdepth));
            double cwa = ChannelWidth->Drc/ChannelAdj->Drc;
            double fr = 1.0/(double)step * _dt/(0.5*ChannelAdj->Drc);

            for (int i = 0; i < step; i++) // do the flow twice as a kind of iteration
            {
                dH = std::max(0.0, (ChannelWH->Drc-chdepth));

                if (dH > _h->Drc) {
                    // flow from channel

                    double frac = std::min(1.0, fr * sqrt(2*GRAV*dH));  // start Vb here gives unstable hydrograph!
                    double dwh =  dH * frac;
                   // double vol = dwh * ChannelWidth->Drc * DX->Drc;
                   // double newh = _h->Drc + vol/CHAdjDX->Drc;

                    if (_h->Drc + dwh*cwa > dH-dwh) {
                        // if flow causes situation to reverse (channel dips below _h)
                        // find a smaller frac
                        while (_h->Drc + dwh*cwa > dH-dwh) {
                            frac = 0.9*(dH-_h->Drc)/((1-cwa)*dH);
                            dwh = dH*frac;
                        }
                    } else {
                        _h->Drc += dwh*cwa;
                        ChannelWH->Drc -= dwh;

                        if(SwitchErosion) {
                            double sed = frac * ChannelSSConc->Drc * dwh*ChannelWidth->Drc*DX->Drc;
                            ChannelSSSed->Drc -= sed;
                            SSFlood->Drc += sed;
                        }
                    }
                } else {
                    // flow to channel, dH can be 0 = channel wh below edge

                    double frac = std::min(1.0, fr * V->Drc);
                    double dwh = _h->Drc * frac;

                    if (dH + dwh/cwa > _h->Drc-dwh) {
                        while (dH + dwh/cwa > _h->Drc-dwh) {
                            frac = 0.9*(dH-_h->Drc)/_h->Drc * cwa/(1+cwa);
                            dwh = dH*frac;
                        }
                    } else {
                        _h->Drc -= dwh;
                        ChannelWH->Drc += (dwh/cwa);
                        if(SwitchErosion) {
                            double sed = frac * SSFlood->Drc;
                            ChannelSSSed->Drc += sed;
                            SSFlood->Drc -= sed;
                        }
                    }
                }

            } // step
        } // tma > 0
    }}
//qDebug() <<  MB << MBs;
}
//---------------------------------------------------------------------------
/**
 * @fn void TWorld::ToFlood(void)
 * @brief Calculates overland flow that flows into flooding water
 *
 * Calculates overland flow of water and sediment that flows into flooding water
 * based on the runoff partitioning factor. Depending on the parameter, water
 * is either transformed quickly or slowly. This imitates the effect that overland
 * flow would have on the velocity of the flood water.
 *
 * @return void
 * @see runoff_partitioning
 */
void TWorld::ToFlood()
{
    #pragma omp parallel for  num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if (hmx->Drc > HMIN && WHrunoff->Drc > HMIN) // && (WHrunoff->Drc > hmx->Drc))
        {
            double frac = 1.0;//1-exp(-2.0*hmx->Drc/(WHrunoff->Drc+HMIN));

            frac = std::max(std::min(frac, 1.0),0.0);
            double dwh = frac * WHrunoff->Drc;

            hmx->Drc += dwh;
            WH->Drc = WHstore->Drc;
            WHrunoff->Drc = 0;

            hmxWH->Drc = hmx->Drc + WH->Drc;
            WaterVolall->Drc = CHAdjDX->Drc*(WHrunoff->Drc + hmx->Drc) + MicroStoreVol->Drc;

            if(SwitchErosion)
            {
                double dsed = frac*Sed->Drc;
                SSFlood->Drc += dsed;
                Sed->Drc = 0;
                Conc->Drc = 0;

                SWOFSedimentLayerDepth(r,c,hmx->Drc, V->Drc);
                SWOFSedimentSetConcentration(r,c,hmx);
               // Conc->Drc = MaxConcentration(WaterVolall->Drc, Sed->Drc);
            }
        }
    }}
}
//---------------------------------------------------------------------------
// DO NOT MAKE PARALLEL
void TWorld::FloodMaxandTiming()
{
    // floodwater volume and max flood map
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if (hmxWH->Drc > minReportFloodHeight) {
            floodTime->Drc += _dt/60;
            floodHmxMax->Drc = std::max(floodHmxMax->Drc, hmxWH->Drc);
            // for output
        }

        floodVMax->Drc = std::max(floodVMax->Drc, V->Drc);
        floodVHMax->Drc = std::max(floodVHMax->Drc, V->Drc*hmxWH->Drc);
        // max velocity
        WHmax->Drc = std::max(WHmax->Drc, hmxWH->Drc);
    }}
    floodVolTotMax = 0;
    floodArea = 0;
    double area = _dx*_dx;

   // #pragma omp parallel for reduction(+:floodVolTotMax,floodArea) num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if (floodHmxMax->Drc > minReportFloodHeight) {
            floodVolTotMax += floodHmxMax->Drc*area;
        }
        if (hmxWH->Drc > minReportFloodHeight && floodTimeStart->Drc == 0)  {
            if (SwitchWaveUser)
                floodTimeStart->Drc = (time - BeginTime)/60.0;
            else
                floodTimeStart->Drc = (time - RainstartTime)/60.0;
            // time since first pixel received rainfall
        }
        if (hmxWH->Drc > minReportFloodHeight) {
            floodArea += area;
        }
    }}

    floodAreaMax = std::max(floodArea,floodAreaMax);
}
//---------------------------------------------------------------------------
// NOTE DEM has barriers included, done in shade map calculation !!!!
void TWorld::ChannelFlood(void)
{
    if (!SwitchIncludeChannel)
        return;

    ToFlood();
    // mix HWrunoff with hmx
    // if toflood before channeloverflow then MB error in sed

    if (SwitchChannel2DflowConnect)
        ChannelOverflowIteration(hmx, V);
    else
        ChannelOverflow(hmx, V);
    // determine overflow water => hmx      
    // hmx is flood water, WH is overlandflow, WHrunoff etc

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        ChannelWaterVol->Drc = ChannelWH->Drc * ChannelDX->Drc * ChannelWidth->Drc;
        WaterVolall->Drc = CHAdjDX->Drc*(WHrunoff->Drc + hmx->Drc) + MicroStoreVol->Drc;
        // do not recalc floodvol, MB errors

        // recalc channel water vol else big MB error
        if(SwitchErosion)
        {
            SWOFSedimentLayerDepth(r,c,hmx->Drc, V->Drc);
            SWOFSedimentSetConcentration(r,c, hmx);

            RiverSedimentLayerDepth(r,c);
            RiverSedimentMaxC(r, c);
            // all concentrations, possible ChannelDep when surplus
        }
    }}

    double dtflood = 0;

    startFlood = false;
    FOR_ROW_COL_MV {
        if (hmx->Drc > 0) {
            startFlood = true;
            break;
        }
    }

    dtflood = fullSWOF2open(hmx, Uflood, Vflood, DEM);
    // in kindyn hmx is the channel overflow/flood part of the surface water, the rest is kinwave WHrunoff

    //new flood domain
    nrFloodedCells = 0;
    FOR_ROW_COL_MV {
        if (hmx->Drc > 0) {
            FloodDomain->Drc = 1;
            nrFloodedCells += 1.0;
        }
        else
            FloodDomain->Drc = 0;
    }

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        Qflood->Drc = 0;
        if (FloodDomain->Drc > 0) {
            V->Drc = sqrt(Uflood->Drc*Uflood->Drc+Vflood->Drc*Vflood->Drc);
            Qflood->Drc = V->Drc * hmx->Drc * ChannelAdj->Drc;
            Qn->Drc = 0;//V->Drc * WHrunoff->Drc * ChannelAdj->Drc;
            // ??????????? why, wh
        }
    }}

    Boundary2Ddyn();
    // 2D boundary flow

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        WH->Drc = WHrunoff->Drc+ WHstore->Drc;
        // add new average waterlevel (A/dx) to stored water

        WaterVolall->Drc = CHAdjDX->Drc*(WHrunoff->Drc + hmx->Drc) + MicroStoreVol->Drc;

        hmxWH->Drc = WH->Drc + hmx->Drc;
        // all water on surface

        hmxflood->Drc = std::max(0.0, WHrunoff->Drc + hmx->Drc - minReportFloodHeight);

        FloodWaterVol->Drc = hmxflood->Drc * CHAdjDX->Drc;
        double WHrunoffOutput = std::min(WHrunoff->Drc + hmx->Drc, minReportFloodHeight);
        RunoffWaterVol->Drc = WHrunoffOutput * CHAdjDX->Drc;
        // these are only used for reporting totals on screen and in file

        if(SwitchErosion) {
            Conc->Drc = MaxConcentration(WaterVolall->Drc, Sed->Drc);
            if (FloodDomain->Drc  > 0) {
                double sed = SSFlood->Drc + BLFlood->Drc;
                Conc->Drc =  MaxConcentration(FloodWaterVol->Drc, sed);
                Qsn->Drc += Conc->Drc*Qflood->Drc;
            }
        }
     }}

    FloodMaxandTiming();

    double area = nrFloodedCells*_dx*_dx;
    if (area > 0)
        debug(QString("Flooding (dt %1 sec, n %2): area %3 m2, %4 cells").arg(dtflood,6,'f',3).arg(iter_n,4).arg(area,8,'f',1).arg(nrFloodedCells));//.arg(K2DQOutBoun));
    // some screen error reporting

}

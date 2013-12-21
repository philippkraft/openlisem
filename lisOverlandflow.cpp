
/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011  Victor Jetten
**  contact:
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with this program.  If not, see <http://www.gnu.org/licenses/>.
**
**  Author: Victor Jetten
**  Developed in: MingW/Qt/
**  website, information and code: http://lisem.sourceforge.net
**
*************************************************************************/

/*!
  \file lisOverlandflow.cpp
  \brief calculate fraction flowing in the channel, Q, V and call kin wave

functions: \n
- void TWorld::ToChannel(void) \n
- void TWorld::CalcVelDisch(void) \n
- void TWorld::OverlandFlow(void) \n
 */

#include "model.h"
#define tiny 1e-8
//---------------------------------------------------------------------------
//fraction of water and sediment flowing into the channel
void TWorld::ToChannel(void)
{
    if (SwitchIncludeChannel)
    {
        FOR_ROW_COL_MV_CH
        {
            double fractiontochannel = min(_dt*V->Drc/(0.5*max(0.01,ChannelAdj->Drc)), 1.0);
            double Volume = WHrunoff->Drc * FlowWidth->Drc * DX->Drc;

            // note ChanelAdj can be 0 then fraction is 1

            if (SwitchAllinChannel)
                if (Outlet->Drc == 1)
                    fractiontochannel = 1.0;
            // in catchment outlet cell, throw everything in channel

            if (SwitchBuffers)
                if (BufferID->Drc > 0)
                    fractiontochannel = 1.0;
            // where there is a buffer in the channel, all goes in the channel

            if (SwitchChannelFlood)
            {
                if (ChannelWH->Drc >= ChannelDepth->Drc+ChannelLevee->Drc)
                    fractiontochannel = 0;
                // no inflow when flooded
                if (ChannelMaxQ->Drc > 0)
                    fractiontochannel = 0;
                // no surface inflow when culverts and bridges
            }

            RunoffVolinToChannel->Drc = fractiontochannel*Volume;
            // water diverted to the channel
            WHrunoff->Drc *= (1-fractiontochannel);

            WH->Drc = WHrunoff->Drc + WHstore->Drc;
            //VJ 130425

            // adjust water height
            if (SwitchErosion)
            {
                SedToChannel->Drc = fractiontochannel*Sed->Drc;
                //sediment diverted to the channel
                Sed->Drc -= SedToChannel->Drc;
                // adjust sediment in suspension
            }
        }
        CalcVelDisch();
        // recalc velocity and discharge
    }
}
//---------------------------------------------------------------------------
void TWorld::CalcVelDisch(void)
{
  //  tm->fill(0);
 //   tma->fill(0);
    FOR_ROW_COL_MV
    {
        double Perim;
        const double beta = 0.6;
        const double _23 = 2.0/3.0;
        double beta1 = 1/beta;
        //double kinvisc = 1.1e-6; // 15 degrees celcius water
        double NN = N->Drc;

        if (SwitchChannelFlood)
            NN = N->Drc * qExp(mixing_coefficient*hmx->Drc);
        // slow down water in flood zone
        //NOTE CALCULATE REYNOLDS AND SEE IF TURBULENCE MAKES COEFFICIENT HIGHER
        //LAMINAR MEANS NON-MIXING SO NN = N
        //    tma->Drc = hmx->Drc * UVflood->Drc/1.1e-6;

        // avg WH from soil surface and roads, over width FlowWidth
        Perim = 2*WHrunoff->Drc+FlowWidth->Drc;

        if (Perim > 0)
            R->Drc = WHrunoff->Drc*FlowWidth->Drc/Perim;
        else
            R->Drc = 0;

        Alpha->Drc = pow(NN/sqrt(Grad->Drc) * pow(Perim, _23),beta);

        if (Alpha->Drc > 0)
            Q->Drc = pow((FlowWidth->Drc*WHrunoff->Drc)/Alpha->Drc, beta1);
        else
            Q->Drc = 0;

        V->Drc = pow(R->Drc, _23)*sqrt(Grad->Drc)/NN;
        //tm->Drc = V->Drc * R->Drc/1.1e-6;
    }
    // tm->report("reyn");
    //  tma->report("reynF");

}
//---------------------------------------------------------------------------
void TWorld::OverlandFlow(void)
{

    /*---- Water ----*/
    // recalculate water vars after subtractions in "to channel"
    FOR_ROW_COL_MV
    {
        WaterVolin->Drc = DX->Drc * (WHrunoff->Drc*FlowWidth->Drc + WHstore->Drc*SoilWidthDX->Drc);

        // WaterVolin total water volume in m3 before kin wave, WHrunoff may be adjusted in tochannel
        q->Drc = FSurplus->Drc*_dx/_dt;
        // infil flux in kin wave <= 0, in m2/s, use _dx bexcause in kiv wave DX is used
    }

    //NOTE if buffers: all water into channel

    /*---- Sediment ----*/
    if (SwitchErosion)
    {
        // calc seediment flux going in kin wave as Qs = Q*C
        FOR_ROW_COL_MV
        {
            Qs->Drc =  Q->Drc * Conc->Drc;
            // calc sed flux as water flux * conc m3/s * kg/m3 = kg/s
        }
    }

    Qn->setMV();
    Qsn->fill(0);
    // flag all new flux as missing value, needed in kin wave and replaced by new flux
    FOR_ROW_COL_MV
    {
        if (LDD->Drc == 5) // if outflow point, pit
        {
            /* TODO: WHEN MORE PITS QPEAK IS FIRST INSTEAD OF MAIN PIT? */
            Kinematic(r,c, LDD, Q, Qn, Qs, Qsn, q, Alpha, DX, WaterVolin, Sed, BufferVol, BufferSed);
            //VJ 110429 q contains additionally infiltrated water volume after kin wave in m3
        }
    }
//    KinWave(LDD, Q, Qn, q, Alpha, DX);
    /*
      routing of substances add here!
      do after kin wave so that the new flux Qn out of a cell is known
      you need to have the ingoing substance flux QS (mass/s)
      and it will give outgoing flux QSn (mass/s)
      and the current amount Subs (mass) in suspension+solution
    */

    if (SwitchPesticide)
    {
        // calc pesticide flux going in kin wave as Qp = Q*C
        FOR_ROW_COL_MV
        {
            Qp->Drc =  Q->Drc * C->Drc;
            // calc sed flux as water flux * conc m3/s * kg/m3 = kg/s
        }


        Qpn->fill(0);
        FOR_ROW_COL_MV
        {
            if (LDD->Drc == 5) // if outflow point, pit
            {
                routeSubstance(r,c, LDD, Q, Qn, Qp, Qpn, Alpha, DX, WaterVolin, Pest);
            }
        }
//        Qpn->report("qpn");
//        Qp->report("qp");
//        Pest->report("pesta");

    }
    // calculate resulting flux Qn back to water height on surface
    FOR_ROW_COL_MV
    {
        double WHoutavg = (Alpha->Drc*pow(Qn->Drc, 0.6))/ChannelAdj->Drc;
        // WH based on A/dx = alpha Q^beta / dx
        /* TODO _dx also needs to be corrected for wheeltracks and gullies */

        WHroad->Drc = WHoutavg;
        // set road to average outflowing wh, no surface storage.

        WH->Drc = WHoutavg + WHstore->Drc;
        // add new average waterlevel (A/dx) to stored water

        if (ChannelAdj->Drc > 0 && WHoutavg > 0)
            V->Drc = Qn->Drc/(WHoutavg*ChannelAdj->Drc);
        else
            V->Drc = 0;
        // recalc velocity for output to map ????

        WaterVolall->Drc = DX->Drc*(WH->Drc*SoilWidthDX->Drc + WHroad->Drc*RoadWidthDX->Drc);
        // new water volume after kin wave, all water incl depr storage
        //     WaterVolall->Drc = DX->Drc * (WHrunoff->Drc*FlowWidth->Drc + WHstore->Drc*SoilWidthDX->Drc);

        double diff = q->Drc*_dt + WaterVolin->Drc - WaterVolall->Drc - Qn->Drc*_dt;
        //diff volume is sum of incoming fluxes+volume before - outgoing flux - volume after

        double difff = diff;
        diff = min(diff, -FSurplus->Drc*ChannelAdj->Drc*_dx/_dt);
        if (FFull->Drc == 1)
            diff = 0;
        difkin->Drc = difff - diff;

        if (SwitchBuffers && BufferVol->Drc > 0)
        {
            //qDebug() << "slope" << BufferVol->Drc << q->Drc*_dt << WaterVolin->Drc << WaterVolall->Drc << Qn->Drc*_dt << diff;
            //NOTE: buffervolume is affected by sedimentation, this causes a water volume loss that is corrected in the
            // totals and mass balance functions
        }
        else
            InfilVolKinWave->Drc = diff;
        // correct infiltration in m3
        // TODO what if infiltration == none then correct watervolume out, but then Qn and watervolout do not match?

        if (SwitchErosion)
        {
            Conc->Drc = (Qn->Drc > 1e-6 ? Qs->Drc/Qn->Drc : 0);
            //MaxConcentration(WaterVolall->Drc, Sed->Drc);
            // CHANGED, MORE STABLE CONC 19/9/13
            // correct for very high concentrations, 850 after Govers et al
            // recalc sediment volume

            if (SwitchPesticide)
            {
                //C->Drc = ConcentrationP(WaterVolall->Drc, Pest->Drc);
                C->Drc = Qn->Drc > 1e-10 ? Qpn->Drc/Qn->Drc : 0;//ConcentrationP(WaterVolall->Drc, Pest->Drc);
                C_N->Drc = C->Drc;
                //qDebug()<< "ds overlandflow"<< C->Drc;
                //qDebug()<< "ds overlandflow"<< Pest->Drc;
            }

        }
    }
//    difkin->report("df");
}
//---------------------------------------------------------------------------

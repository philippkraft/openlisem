
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
  \file lisTotalsMB.cpp
  \brief calculate water and sediment flux totals and mass balance

functions: \n
- void TWorld::Totals() calculate the totals for all water and fluxes for mss balance and output\n
- void TWorld::MassBalance() water and sediment mass balance\n
 */


#include "model.h"


//---------------------------------------------------------------------------
// totals for screen and file output and mass balance
void TWorld::Totals(void)
{
    double rainfall, snowmelt;
    double oldrainpeak, oldsnowpeak;
    double catchmentAreaFlatMM = 1000.0/(_dx*_dx*nrCells);

    //    QFile fout("massbalancenew.txt");
    //    fout.open(QIODevice::WriteOnly | QIODevice::Text);
    //    fout.close();

    /***** WATER *****/

    if (SwitchRainfall)
    {
        RainAvgmm = Rain->mapAverage()*1000.0;
        RainTotmm += RainAvgmm;
        // spatial avg area rainfall in mm

//        tm->calcMapValue(Rain, (_dx*_dx), MUL); //in m3
//        rainfall = tm->mapTotal();
//        RainTot += rainfall; // in m3

        rainfall = 0;
        oldrainpeak = Rainpeak;
        FOR_ROW_COL_MV
        {
            rainfall += (Rain->Drc*_dx*_dx);
            Rainpeak = max(rainfall, Rainpeak);
        }

//        Rainpeak = max(Rainpeak, rainfall);
        if (oldrainpeak  < Rainpeak)
            RainpeakTime = time;
    }

    if (SwitchSnowmelt)
    {
        SnowAvgmm += Snowmelt->mapAverage()*1000;
        SnowTotmm += SnowAvgmm;

        tm->calcMapValue(Snowmelt, (_dx*_dx), MUL); //in m3
        snowmelt = tm->mapTotal();
        SnowTot += snowmelt; // in m3

        oldsnowpeak = Snowpeak;
        Snowpeak = max(Snowpeak, snowmelt);
        if (oldsnowpeak < Snowpeak)
            SnowpeakTime = time;
    }

    IntercTot = Interc->mapTotal();
    IntercTotmm = IntercTot*catchmentAreaFlatMM;
    // interception in mm and m3

    //houses
    IntercHouseTot = IntercHouse->mapTotal();
    IntercHouseTotmm = IntercHouseTot*catchmentAreaFlatMM;
    // interception in mm and m3

    InfilTot += InfilVol->mapTotal() + InfilVolKinWave->mapTotal() + InfilVolFlood->mapTotal(); //m3
    difkinTot =0;//+=  difkin->mapTotal();
    InfilKWTot += InfilVolKinWave->mapTotal(); // not really used, available for output when needed
    InfilTotmm = max(0,(InfilTot)*catchmentAreaFlatMM);
    // infiltration mm and m3

    tm->calcMapValue(WHstore, 1000, MUL); //mm
    SurfStoremm = tm->mapAverage();
    // surface storage CHECK THIS
    // does not go to MB, is already in tot water vol

    WaterVolTot = WaterVolall->mapTotal();//m3
    WaterVolTotmm = WaterVolTot*catchmentAreaFlatMM; //mm
    // water on the surface in runoff in m3 and mm
    //NOTE: surface storage is already in here so does not need to be accounted for in MB


    // runoff fracyion per cell calc as in-out/rainfall, indication of sinks and sources of runoff
    // exclude channel cells
    FOR_ROW_COL_MV
            if(ChannelWidthUpDX->Drc == 0)
    {
        runoffTotalCell->Drc += Qn->Drc * _dt;
    }
    upstream(LDD, runoffTotalCell, tm);

    FOR_ROW_COL_MV
    {
        runoffFractionCell->Drc = RainCumFlat->Drc > 0 ? (runoffTotalCell->Drc-tm->Drc)/(RainCumFlat->Drc*_dx*_dx) : 0;
    }

    // sum outflow m3 for all timesteps for the outlet
    FOR_ROW_COL_MV
    {
        if (LDD->Drc == 5)
            Qtot += Qn->Drc*_dt;
    }
    // sum outflow m3 for all timesteps for all outlets, in m3
    // needed for mass balance

    QtotOutlet += Qn->DrcOutlet*_dt;
    // for screen output, total main outlet in m3

    QtotPlot += Qn->DrcPlot * _dt;
    //QPlot = Qn->DrcPlot;
    //VJ 110701 for screen output, total in hydrograph point n in m3

    if (SwitchIncludeChannel)
    {
        WaterVolTot += ChannelWaterVol->mapTotal(); //m3
        // add channel vol to total
        WaterVolTotmm = WaterVolTot*catchmentAreaFlatMM; //mm
        // recalc in mm for screen output
        FOR_ROW_COL_MV_CH
        {
            if (LDDChannel->Drc == 5)
                Qtot += ChannelQn->Drc*_dt;
            ChannelQntot->Drc += ChannelQn->Drc*_dt;  //m3 spatial for output

        }
        // add channel outflow (in m3) to total for all pits
        //Qtotmm = Qtot*catchmentAreaFlatMM;
        // recalc in mm for screen output

        QtotOutlet += ChannelQn->DrcOutlet * _dt;
        // add channel outflow (in m3) to total for main outlet
        QtotPlot += ChannelQn->DrcPlot * _dt;
        // add channel outflow (in m3) to total for main outlet
        // TotalWatervol->calcMap(ChannelWaterVol,ADD);
        // add channel volume to total for sed conc calc

        if (SwitchChannelFlood)
        {
            floodVolTot = FloodWaterVol->mapTotal();
            floodTotmm = floodVolTot*catchmentAreaFlatMM;
        }
        if (runstep == 1)
            floodVolTotInit = floodVolTot;

    }

    if (SwitchIncludeTile)
    {
        WaterVolSoilTot = TileWaterVolSoil->mapTotal();
        // input for mass balance, is the water seeping from the soil, input
        // this is the water before the kin wave
        tm->calc2Maps(TileDrainSoil, TileWidth, MUL); //in m3
        tm->calcMap(TileDX, MUL); //in m3
        // tm->calcV(_dx, MUL); //in m3 ??? or DX?
        TileVolTot += tm->mapTotal(); // in m3

        // water after kin wave
        WaterVolTot += TileWaterVol->mapTotal(); //m3
        // add tile vol to total
        WaterVolTotmm = WaterVolTot*catchmentAreaFlatMM; //mm
        // recalc in mm for screen output

        //Qtot += TileQoutflow->DrcOutlet;
        FOR_ROW_COL_MV_TILE
                if (LDDTile->Drc == 5)
                Qtot += TileQn->Drc * _dt;
        // add tile outflow (in m3) to total for all pits
        //Qtotmm = Qtot*catchmentAreaFlatMM;
        // recalc in mm for screen output

        QtotOutlet += TileQn->DrcOutlet * _dt;
        // add channel outflow (in m3) to total for main outlet
        QtotPlot += TileQn->DrcPlot * _dt;
        // add channel outflow (in m3) to total for subcatch outlet

    }

    if (SwitchBuffers)
    {
        BufferVolTot = BufferVol->mapTotal(); // in m3
        if (SwitchIncludeChannel)
            BufferVolTot += ChannelBufferVol->mapTotal();
        //sum up all volume remaining in all buffers (so the non-water!)
        BufferVolin = BufferVolTotInit - BufferVolTot;
        //subtract this from the initial volume to get the total water inflow in the buffers
    }

    // output fluxes for reporting to file and screen in l/s!
    FOR_ROW_COL_MV
    {
        Qoutput->Drc = 1000*(Qn->Drc + ChannelQn->Drc + TileQn->Drc); // in l/s
        if (Qoutput->Drc < 0.0001)
            Qoutput->Drc = 0.0001;
        // added minimum here to avoid strange maps
    }
    QPlot = 1000*(Qn->DrcPlot + ChannelQn->DrcPlot + TileQn->DrcPlot);
    // plot point output in l/s

    Qtotmm = Qtot*catchmentAreaFlatMM;
    // recalc to mm for screen output

    oldrainpeak = Qpeak;
    Qpeak = max(Qpeak, Qoutput->DrcOutlet);
    if (oldrainpeak < Qpeak)
        QpeakTime = time;
    // peak flow and peak time calculation, based on sum channel and runoff

    QpeakPlot = max(QpeakPlot, Qoutput->DrcPlot);

    // do this last because of possible flood inf volume
    FOR_ROW_COL_MV
    {
        InfilVolCum->Drc += InfilVol->Drc + InfilVolKinWave->Drc + InfilVolFlood->Drc ;
        InfilmmCum->Drc = max(0, InfilVolCum->Drc*1000.0/CellArea->Drc);
    }


    /***** SEDIMENT *****/
    // note DETFLOW, DETSPLASH AND DEP ARE IN KG/CELL
    if (SwitchErosion)
    {
        DetSplashTot += DETSplash->mapTotal();
        DetFlowTot += DETFlow->mapTotal();
        DepTot += DEP->mapTotal();
        DetTot += DETSplash->mapTotal() + DETFlow->mapTotal();
        SedTot = Sed->mapTotal();
        // all in kg/cell

        //SoilLossTot += Qsoutflow->DrcOutlet;
        FOR_ROW_COL_MV
                if (LDD->Drc == 5)
                SoilLossTot += Qsn->Drc * _dt;
        // sum all sed in all pits (in kg), needed for mass balance

        SoilLossTotOutlet += Qsn->DrcOutlet * _dt;
        // for screen output, total main outlet sed loss in kg
        TotalSed->copy(Sed);
        // for sed conc

        if (SwitchIncludeChannel)
        {
            // units here in kg, conversion to ton in report functions
            ChannelDetTot += ChannelDetFlow->mapTotal();
            ChannelDepTot += ChannelDep->mapTotal();
            ChannelSedTot = ChannelSed->mapTotal();

            FOR_ROW_COL_MV_CH
                    if (LDDChannel->Drc == 5)
                    SoilLossTot += ChannelQsn->Drc * _dt;
            // add sed outflow for all pits to total soil loss

            SoilLossTotOutlet += ChannelQsn->DrcOutlet * _dt;
            // add channel outflow (in kg) to total for main outlet

            TotalSed->calcMap(ChannelSed, ADD);
            // needed for sed conc in file output
        }

        if (SwitchBuffers || SwitchSedtrap)
        {
            BufferSedTot = BufferSed->mapTotal();
            if (SwitchIncludeChannel)
                BufferSedTot += ChannelBufferSed->mapTotal();

        }
        /** TODO add gully, wheeltracks etc */

        // spatial totals for output all in kg/cell
        FOR_ROW_COL_MV
        {
            Qsoutput->Drc = Qsn->Drc + ChannelQsn->Drc;  // sum channel and OF sed output in kg/s

            TotalDetMap->Drc += DETSplash->Drc + DETFlow->Drc;
            TotalDepMap->Drc += DEP->Drc;
            if (SwitchIncludeChannel)
            {
                TotalDetMap->Drc += ChannelDetFlow->Drc;
                TotalDepMap->Drc += ChannelDep->Drc;
            }
            TotalSoillossMap->Drc = TotalDetMap->Drc + TotalDepMap->Drc;
        }

        SoilLossTotPlot += Qsoutput->DrcPlot * _dt;

        FOR_ROW_COL_MV
        {
            TotalConc->Drc = (Qoutput->Drc > 1e-6 ? Qsoutput->Drc/Qoutput->Drc : 0);
        }
    }

    if (SwitchPesticide)
    {
        FOR_ROW_COL_MV
        {
            // = WHoutavg->Drc*_dx*DX->Drc*C->Drc*1000*1000*1000; //µg
            PDisMixing->Drc = CM->Drc*epsil->Drc*poro->Drc*_dx*_dx*1000*1000*1000; //µg
            PSorMixing->Drc = CS->Drc*epsil->Drc*rhob->Drc*_dx*_dx*1000*1000*1000; //µg
            PInfilt->Drc = pestiinf->Drc*CM->Drc*_dx*_dx*_dt*1000*1000*1000; //µg
            PStorage->Drc= WHstore->Drc*_dx*_dx*C->Drc*1000*1000*1000; //µg
            PRunoffSpatial->Drc = Pest->Drc*1000*1000*1000; //µg

            //            PRunoffSpatialex->Drc= WHoutavg->Drc*_dx*DX->Drc*C_Kexplicit->Drc*1000*1000*1000; //µg
            //            PDisMixingex->Drc = CM_Kexplicit->Drc*epsil->Drc*poro->Drc*_dx*DX->Drc*1000*1000*1000; //µg
            //            PSorMixingex->Drc = CS_Kexplicit->Drc*epsil->Drc*rhob->Drc*_dx*DX->Drc*1000*1000*1000; //µg
            //            PInfiltex->Drc = pestiinf->Drc*CM_Kexplicit->Drc*_dx*DX->Drc*_dt*1000*1000*1000; //µg

        }

        Pestdetach += Pdetach ->mapTotal(); //KCM
        PestCinfilt += PCinfilt->mapTotal(); //fc
        PestCfilmexit += PCfilmexit->mapTotal(); //KC
        PestLossTotOutlet += Qn->DrcOutlet*C->DrcOutlet*_dt*1000*1000*1000; //µg
        PestRunoffSpatial = PRunoffSpatial->mapTotal();
        PestDisMixing = PDisMixing->mapTotal();
        PestSorMixing = PSorMixing->mapTotal();
        PestInfilt += PInfilt->mapTotal();
        PestStorage = PStorage->mapTotal();

        double MBtest=0.0;

        // if (PestLossTotOutlet > 1e-9)
        // MBtest = (Pestdetach-PestCinfilt-PestRunoffSpatial-PestLossTotOutlet)*100/Pestdetach;

        //if (Pestdetach > 1e-9)
        MBtest = Pestdetach-PestCinfilt-PestCfilmexit-PestLossTotOutlet;
        // qDebug()<< "pestdetach" << Pestdetach << "pestCinfilt"<< PestCinfilt << "pestCfilmexit"<< PestCfilmexit<< "pestlosstotoutlet"<<PestLossTotOutlet;
        // qDebug()<< "MBtest" << MBtest;
        double test=0.0;
        test += InfilVolKinWave->mapTotal();

        //        PestLossTotOutletex += Qn->DrcOutlet*C_Kexplicit->DrcOutlet*_dt*1000*1000*1000; //µg
        //        PestRunoffSpatialex = PRunoffSpatialex->mapTotal();
        //        PestDisMixingex = PDisMixingex->mapTotal();
        //        PestSorMixingex = PSorMixingex->mapTotal();
        //        PestInfiltex += PInfiltex->mapTotal();

        // flux en µg
        //        double flux1=epsil->DrcOutlet*rhob->DrcOutlet*kr->DrcOutlet*KD->DrcOutlet*CM->DrcOutlet*_dx*DX->DrcOutlet*_dt*1000*1000*1000;
        //        double flux2=kr->DrcOutlet*CS->DrcOutlet*rhob->DrcOutlet*epsil->DrcOutlet*_dx*DX->DrcOutlet*_dt*1000*1000*1000;
        //        double flux3=pestiinf->DrcOutlet*CM->DrcOutlet*_dx*DX->DrcOutlet*_dt*1000*1000*1000;
        //        double flux4=Kfilm->DrcOutlet*CM->DrcOutlet*_dx*DX->DrcOutlet*_dt*1000*1000*1000;
        //        double flux5=(Kfilm->DrcOutlet+pestiinf->DrcOutlet)*C->DrcOutlet*_dx*DX->DrcOutlet*_dt*1000*1000*1000;
        //        double flux6=(Kfilm->DrcOutlet+RainNet->DrcOutlet/_dt)*C->DrcOutlet*_dx*DX->DrcOutlet*_dt*1000*1000*1000;


        //            QFile fout("massbalancenew.txt");
        //            fout.open(QIODevice::Append | QIODevice::Text);
        //            QTextStream out(&fout);
        //            out.setRealNumberPrecision(3);
        //            out.setFieldWidth(0);
        //            out.setRealNumberNotation(QTextStream::FixedNotation);

        //            out << time/60 << " " << PestMassApplied << " " << PestDisMixing << " " << PestSorMixing << " " << PestLossTotOutlet << " " << PestRunoffSpatial
        //                 << " " << PestInfilt << " " << (PestMassApplied-PestLossTotOutlet-PestRunoffSpatial-PestDisMixing-PestSorMixing-PestInfilt-PestStorage)*100/PestMassApplied << " "
        //                 << RainTot << " " << WaterVolSoilTot << " " << IntercTot << " " << InfilTot << " " << Qtot*1000*1000 << " "
        //                 << MBtest << " " << test << " "<< flux3 << " "<< flux4 << " "<< flux5 << " "<< flux6 <<" "<< pestiinf->DrcOutlet*pow(10.0,9)<< " "<<CM->DrcOutlet*pow(10.0,6)<<" "
        //                 << CS->DrcOutlet*pow(10.0,6)<<" "<< fact->DrcOutlet*1000<< " "<< InfilVol->DrcOutlet*1000*1000<<" "<<Qn->DrcOutlet*pow(10.0,6) << " "<< PDisMixing->DrcOutlet << " "<< poro->DrcOutlet
        //                 << " "<< epsil->DrcOutlet<< " "<< DX->DrcOutlet << " " << switchrunoff << " "<< K1->DrcOutlet << " "<< Q->DrcOutlet*pow(10.0,6)<< " "<< C->DrcOutlet*pow(10.0,10)
        //                 << " "<< WHoutavg->DrcOutlet << " "<< WHoutavgold->DrcOutlet<<" "<< (PestMassApplied-PestLossTotOutletex-PestRunoffSpatialex-PestDisMixingex-PestSorMixingex-PestInfiltex)*100/PestMassApplied
        //                 << " " << InfilVol->DrcOutlet*1000*1000 << " " << InfilVolold->DrcOutlet*1000*1000<< " " << Vup->DrcOutlet << " " << Vup_old->DrcOutlet << " "<< Cold->DrcOutlet*pow(10.0,10);
        //            out << "\n";
        //            out << MBp << "\n";
        //            fout.close();

    }

}
//---------------------------------------------------------------------------
void TWorld::MassBalance()
{
    // Mass Balance water, all in m3
    // VJ 110420 added tile volume here, this is the input volume coming from the soil after swatre
    if (RainTot + SnowTot > 0)
        MB = (RainTot + SnowTot + WaterVolSoilTot + floodVolTotInit
              - IntercTot - IntercHouseTot - InfilTot - WaterVolTot - floodVolTot - Qtot - BufferVolin - difkinTot)/
                (RainTot + SnowTot + WaterVolSoilTot + floodVolTotInit)*100;
    //watervoltot includes channel and tile

    //qDebug() << MB << RainTot << IntercTot << IntercHouseTot << InfilTot << WaterVolTot << floodVolTot << BufferVolin << Qtot<< InfilKWTot;

    // Mass Balance sediment, all in kg
    //   if (SwitchErosion && (DetTot + ChannelDetTot) > 0)
    //      MBs = (DetTot + ChannelDetTot - SoilLossTot - SedTot - ChannelSedTot +
    //             DepTot + ChannelDepTot - BufferSedTot)/(DetTot + ChannelDetTot)*100;
    //VJ 110825 forgot to include channeldettot in denominator in MBs!
    if (SwitchErosion && SoilLossTot > 1e-9)
        MBs = (1-(DetTot + ChannelDetTot - SedTot - ChannelSedTot +
                  DepTot + ChannelDepTot - BufferSedTot)/(SoilLossTot))*100;
    //VJ 121212 changed to mass balance relative to soil loss

    if (SwitchPesticide)
    {
        MBp = (PestMassApplied-PestLossTotOutlet-PestRunoffSpatial-PestDisMixing-PestSorMixing-PestInfilt-PestStorage)*100/PestMassApplied;
        //MBpex = (PestMassApplied-PestLossTotOutletex-PestRunoffSpatialex-PestDisMixingex-PestSorMixingex-PestInfiltex)*100/PestMassApplied;
        //(PestMassApplied-PestLossTotOutlet-PestRunoffSpatial-PestDisMixing-PestSorMixing-PestInfilt-PestStorage)*100/PestMassApplied
        debug(QString("mbp: %1").arg(MBp));
    }
}
//---------------------------------------------------------------------------

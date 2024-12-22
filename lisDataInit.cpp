﻿/*************************************************************************
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


#include <algorithm>
#include <qstring.h>
#include "io.h"
#include "lisemqt.h"
#include "global.h"

#include "model.h"
#include "operation.h"
#include "CsfRGBMap.h"
//---------------------------------------------------------------------------
void TWorld::GetInputData(void)
{
    InitParameters();

    InitStandardInput();
    //## Basic data start of map list etc.

    InitMeteoInput();

    InitLULCInput();
    //## surface related variables

    InitSoilInput();
    //## soil/infiltration data

    InitNewSoilProfile();
    // fin element soil init

    InitErosion();
    //extended sediment stuff

    InitChannel();
    //## read and initialize all channel maps and variables

    InitFlood();
    // vars for dyn wave

    InitBoundary();
    //find domain boundaries

    InitShade();
    InitImages();

    //## read and initialize all tile drain system maps and variables
    InitTiledrains();

    //## get flow barriers;
    InitFlowBarriers();

    InitScreenChanNetwork();

}
//---------------------------------------------------------------------------
void TWorld::InitParameters(void)
{       
    PBiasCorrection = getvaluedouble("Rainfall Bias Correction");
    ETBiasCorrection = getvaluedouble("ET Bias Correction");
    rainfallETa_threshold = getvaluedouble("Rainfall ET threshold");
    rainIDIfactor = getvaluedouble("IDI factor");

    HinitValue = getvaluedouble("Initial matrix potential");

    GW_recharge = getvaluedouble("GW recharge factor");
    GW_flow = getvaluedouble("GW flow factor");
    //GW_inflow = getvaluedouble("GW river inflow factor");
    GW_slope = getvaluedouble("GW slope factor");
    GW_deep = getvaluedouble("GW deep percolation"); // in mm/day
    GW_deep *= 0.001/3600*_dt; //mm/h to m/s
    GW_threshold = getvaluedouble("GW threshold factor");

    // get calibration parameters
    gsizeCalibrationD50 = getvaluedouble("Grain Size calibration D50");
    gsizeCalibrationD90 = getvaluedouble("Grain Size calibration D90");

    ksatCalibration = getvaluedouble("Ksat calibration");
    ksat2Calibration = getvaluedouble("Ksat2 calibration");
    ksat3Calibration = getvaluedouble("Ksat3 calibration");

    SmaxCalibration = getvaluedouble("Smax calibration");
    RRCalibration = getvaluedouble("RR calibration");

    nCalibration = getvaluedouble("N calibration");
    if (nCalibration == 0)
    {
        ErrorString = QString("Calibration: the calibration factor for Mannings n for slopes cannot be zero.");
        throw 1;
    }

    thetaCalibration = getvaluedouble("Theta calibration");
    psiCalibration = getvaluedouble("Psi calibration");
    SD1Calibration = getvaluedouble("SoilDepth1 calibration");
    SD2Calibration = getvaluedouble("SoilDepth2 calibration");

    ChnCalibration = getvaluedouble("Channel N calibration");

    WaveCalibration = getvaluedouble("Boundary water level calibration");


    //ChnTortuosity = 1.0;
    ChnTortuosity = getvaluedouble("Channel tortuosity");
    if (ChnCalibration == 0)
    {
        ErrorString = QString("Calibration: the calibration factor for Mannings n for channels cannot be zero.");
        throw 1;
    }

    ChKsatCalibration = getvaluedouble("Channel Ksat calibration");
    SplashDelivery =getvaluedouble("Splash Delivery Ratio");
    DepositedCohesion = 1.0;//getvaluedouble("Particle Cohesion of Deposited Layer");
    BulkDens =1350;//getvaluedouble("Sediment bulk density");
    //StemflowFraction = getvaluedouble("Stemflow fraction");
    CanopyOpeness = 0.45;//getvaluedouble("Canopy Openess");

    // VJ 170923 moved all 2D switches here
    minReportFloodHeight = getvaluedouble("Minimum reported flood height");
    courant_factor = getvaluedouble("Flooding courant factor");
    courant_factorSed = std::min(0.2,courant_factor);
    // courant_factor_sed = getvaluedouble("Flooding courant factor diffusive");
    TimestepfloodMin = getvaluedouble("Timestep flood");
    SwitchMUSCL = getvalueint("Flood Solution") == 1;
    SwitchHeun = getvalueint("Flood Heun 2nd order") == 1;
    F_pitValue = getvaluedouble("Pit Value");

    SwitchCorrectMB_WH = getvalueint("Correct MB with WH") == 1;
    op.SwitchCorrectMB_WH = SwitchCorrectMB_WH;

    if (SwitchAdvancedOptions) {
        F_MaxIter = getvalueint("Flood max iterations");
        F_fluxLimiter = getvalueint("Flooding SWOF flux limiter"); //minmax, vanleer, albeda
        F_scheme = getvalueint("Flooding SWOF Reconstruction");   //HLL HLL2 Rusanov
        F_minWH = getvaluedouble("Min WH flow");   //HLL HLL2 Rusanov
        //SwitchErosionInsideLoop = getvalueint("Calculate erosion inside 2D loop") == 1;
        SwitchLinkedList = getvalueint("Use linked List") == 1;
        SwitchPerimeterKW = getvalueint("Use Perimeter KW") == 1;
        _dtCHkin = getvaluedouble("Channel Kinwave dt");
        SwitchChannel2DflowConnect = getvalueint("Channel 2D flow connect") == 1;
        SwitchChannelWFinflow = false;//getvalueint("Channel WF inflow") == 1;
    } else {
        F_MaxIter = 200;
        F_minWH = 0.00001;
        F_fluxLimiter = 1; //minmod, vanleer, albeda
        F_scheme = 4;   //Rusanov HLL HLL2 HLL2c
        F_pitValue = _dx/100;
        SwitchLinkedList = false;
        SwitchPerimeterKW = false;
        _dtCHkin = 60.0;
        SwitchChannel2DflowConnect = false;
        SwitchChannelWFinflow = false;

        nN1_ = 3;
        nN2_ = 3;
        nN3_ = 6;
        SoilWBdtfactor = 10;


        //SwitchGWChangeSD = true;
    }
    _CHMaxV = 20.0;
    if (SwitchChannelMaxV)
       _CHMaxV =  getvaluedouble("Channel Max V");

    int wave = getvalueint("Routing Kin Wave 2D");
    if (wave == 0) SwitchKinematic2D = K2D_METHOD_KIN;
    if (wave == 1) SwitchKinematic2D = K2D_METHOD_KINDYN;
    if (wave == 2) SwitchKinematic2D = K2D_METHOD_DYN;
    if (wave < 2) SwitchWaveUser = false; // waveuser is an incoming wave at the boundary (tsunami type)

    userCores = getvalueint("Nr user Cores");
    int cores = omp_get_max_threads();
    if (userCores == 0 || userCores > cores)
        userCores = cores;
    op.cores = userCores;

}
//---------------------------------------------------------------------------
void TWorld::InitStandardInput(void)
{   
    //## catchment data
    LDD = InitMask(getvaluename("ldd"));
    // THIS SHOULD BE THE FIRST MAP
    // LDD is also mask and reference file, everthing has to fit LDD
    // channels use channel LDD as mask

    FOR_ROW_COL_MV {
        if (LDD->Drc == 0)
            SET_MV_REAL8(&LDD->Drc);
    }


    tm = NewMap(0); // temp map for aux calculations
    tma = NewMap(0); // temp map for aux calculations
    tmb = NewMap(0); // temp map for aux calculations
    tmc = NewMap(0); // temp map for aux calculations
    tmd = NewMap(0); // temp map for aux calculations

    nrValidCells = 0;
    FOR_ROW_COL_MV {
        nrValidCells++;
    }

    FOR_ROW_COL_MV {
        // LDD_COOR *newcr = new LDD_COOR;
        // newcr->r = r;
        // newcr->c = c;
        LDD_COOR newcr;
        newcr.r = r;
        newcr.c = c;
        cr_ << newcr;
    }

    FOR_ROW_COL_MV {
        if (LDD->Drc == 5) {
        LDD_COOR newcr;
        newcr.r = r;
        newcr.c = c;
        crldd5_ << newcr;
        }
    }
    nrValidCellsLDD5 = crldd5_.size();

    crlinkedldd_ = MakeLinkedList(LDD);

    DEM = ReadMap(LDD, getvaluename("dem"));
    MBm = NewMap(0);

    Grad = ReadMap(LDD, getvaluename("grad"));  // must be SINE of the slope angle !!!
    checkMap(*Grad, LARGER, 1.0, "Gradient cannot be larger than 1: must be SINE of slope angle (not TANGENT)");
    sqrtGrad = NewMap(0);
    FOR_ROW_COL_MV {
        sqrtGrad->Drc = sqrt(Grad->Drc);
    }

    SwitchSlopeStability = false;
    if (SwitchSlopeStability) {
        tanGrad = NewMap(0);
        cosGrad = NewMap(0);
        BulkDensity = NewMap(1400);
        AngleFriction = NewMap(1.0);  // tan phi = 45 degrees

        FSlope = NewMap(0);

        FOR_ROW_COL_MV {
            // grad = grad = sin(atan(slope(DEMm)))
            double at = asin(Grad->Drc);
            tanGrad->Drc = tan(at);
            cosGrad->Drc = cos(at);
        }
        report(*cosGrad,"cosgrad.map");
        report(*tanGrad,"tangrad.map");
    }


    if (SwitchCorrectDEM)
        CorrectDEM(DEM, Grad);

    if (SwitchBuffers) {
        Buffers = ReadMap(LDD, getvaluename("buffers"));
        calcMap(*DEM, *Buffers, ADD);
    } 

    // int cnt = 0;
    // Outlet = NewMap(0);
    // FOR_ROW_COL_MV {
    //     if(LDD->Drc == 5) {
    //         cnt++;
    //         //qDebug() << "ldd" << r << c << cnt;
    //         Outlet->Drc = cnt;
    //     }
    // }
    bool found = false;
    Outlet = ReadMap(LDD,getvaluename("outlet"));
    FOR_ROW_COL_MV {
        if(Outlet->Drc > 0) {
            found = true;
        }
    }
    if (!found) {
        ErrorString="outlet.map does not have any outlets";
        throw 1;
    }

    // points are user observation points. they should include outlet points
    PointMap = ReadMap(LDD,getvaluename("outpoint"));
    //map with points for output data
    // VJ 110630 show hydrograph for selected output point
    found = false;
    FOR_ROW_COL_MV {
        if(PointMap->Drc > 0) {
            found = true;
        }
    }
    if (!found) {
        copy(*PointMap, *Outlet);
        found = true;
    }


    if (found) {
        crout_.clear();
        FOR_ROW_COL_MV {
            if(PointMap->Drc > 0) {
                LDD_COORout newcr;
                newcr.r = r;
                newcr.c = c;
                newcr.nr = (int)PointMap->Drc ;
                crout_ << newcr;
            }
        }
    } else {
        ErrorString = QString("Outpoint.map has no values above 0. Copy at least outlet(s).");
        throw 1;
    }

    ChannelAdj = NewMap(_dx);
    CHAdjDX = NewMap(0);

}
//---------------------------------------------------------------------------
void TWorld::InitMeteoInput(void)
{
    SwitchUseIDmap = true;

    if (!SwitchRainfallSatellite)
    {
        RainZone = ReadMap(LDD,getvaluename("ID"));
        if (SwitchIDinterpolation) {
            if (SwitchUseIDmap){
                IDRainPoints = ReadMap(LDD,getvaluename("IDGauges"));
            } else {
                IDRainPoints = NewMap(0);
            }
        }
    }

    //### rainfall and interception maps
    RainTot = 0;
    RainTotmm = 0;
    Rainpeak = 0;
    RainpeakTime = 0;
    RainstartTime = -1;
    rainStarted = false;
    ETStarted = false;
    RainAvgmm = 0;
    SnowAvgmm = 0;
    SnowTot = 0;
    SnowTotmm = 0;
    Snowpeak = 0;
    SnowpeakTime = 0;

    Rain = NewMap(0);
    Rainc = NewMap(0);
    RainCum = NewMap(0);
    RainCumFlat = NewMap(0);
    RainNet = NewMap(0);

    if (SwitchIncludeET) {
        ETa = NewMap(0);
        ETaCum = NewMap(0);
        ETp = NewMap(0);
        ETpCum = NewMap(0);

        if (!SwitchETSatellite){
            ETZone = ReadMap(LDD,getvaluename("ETID"));
        }
    }

    if (SwitchSnowmelt) {
        Snowcover = NewMap(0);
        Snowmelt = NewMap(0);
        Snowmeltc = NewMap(0);
        SnowmeltCum = NewMap(0);
        if (!SwitchSnowmeltSatellite)
        {
            SnowmeltZone = ReadMap(LDD,getvaluename("SnowID"));
            FOR_ROW_COL_MV {
                Snowcover->Drc = (SnowmeltZone->Drc == 0 ? 0 : 1.0);
            }
        }
    }

}
//---------------------------------------------------------------------------
//## landuse and surface data
void TWorld::InitLULCInput(void)
{
    //===== surface =====
    N = ReadMap(LDD,getvaluename("manning"));
    checkMap(*N, SMALLER, 1e-6, "Manning's N must be > 0.000001");
    calcValue(*N, nCalibration, MUL);

    Norg = NewMap(0);
    copy(*Norg, *N); //ed in sed trap... if trap is full go back to original N

    RR = ReadMap(LDD,getvaluename("RR"));
    checkMap(*RR, SMALLER, 0.0, "Random roughness RR must be >= 0");
    calcValue(*RR, RRCalibration, MUL);

    //===== interception =====
    LAI = ReadMap(LDD,getvaluename("lai"));
    checkMap(*LAI, SMALLER, 0.0, "LAI must be >= 0");
    Cover = ReadMap(LDD,getvaluename("cover"));
    checkMap(*Cover, SMALLER, 0.0, "Cover fraction must be >= 0");
    checkMap(*Cover, LARGER, 1.0, "Cover fraction must be <= 1.0");

    LeafDrain = NewMap(0);
    CStor = NewMap(0);
    Interc = NewMap(0);
    InterceptionmmCum = NewMap(0);

    if (SwitchIncludeET)
        IntercETa = NewMap(0);

    InterceptionLAIType = getvalueint("Canopy storage equation");

    if (InterceptionLAIType < 8) {
        CanopyStorage = NewMap(0); //in m !!!
        FOR_ROW_COL_MV_L {
            switch (InterceptionLAIType)
            {
                case 0: CanopyStorage->Drc = 0.4376 * LAI->Drc + 1.0356;break; // gives identical results
                    //0.935+0.498*LAI->Drc-0.00575*(LAI->Drc * LAI->Drc);break;
                case 1: CanopyStorage->Drc = 0.2331 * LAI->Drc; break;
                case 2: CanopyStorage->Drc = 0.3165 * LAI->Drc; break;
                case 3: CanopyStorage->Drc = 1.46 * pow(LAI->Drc,0.56); break;
                case 4: CanopyStorage->Drc = 0.0918 * pow(LAI->Drc,1.04); break;
                case 5: CanopyStorage->Drc = 0.2856 * LAI->Drc; break;
                case 6: CanopyStorage->Drc = 0.1713 * LAI->Drc; break;
                case 7: CanopyStorage->Drc = 0.59 * pow(LAI->Drc,0.88); break;

            }
        }}
    } else {
        CanopyStorage = ReadMap(LDD,getvaluename("smax"));
    }
    calcValue(*CanopyStorage, SmaxCalibration, MUL);
    calcValue(*CanopyStorage, 0.001, MUL); // from mm to m
    //NOTE: LAI is still needed for canopy openness


    if (SwitchRoadsystem)
    {
        RoadWidthDX  = ReadMap(LDD,getvaluename("road"));
        checkMap(*RoadWidthDX, LARGER, _dx, "road width cannot be larger than gridcell size");       
    }
    else
        RoadWidthDX = NewMap(0);

    if (SwitchHardsurface)
    {
        HardSurface = ReadMap(LDD,getvaluename("hardsurf"));
        calcValue(*HardSurface, 1.0, MIN);
        calcValue(*HardSurface, 0.0, MAX);
    }
    else
        HardSurface = NewMap(0);

    RoadWidthHSDX = NewMap(0);
    if (SwitchRoadsystem || SwitchHardsurface)
        FOR_ROW_COL_MV_L {
            //double frac = std::min(1.0,(HardSurface->Drc*_dx + RoadWidthDX->Drc)/_dx);
            RoadWidthHSDX->Drc = std::min(_dx, RoadWidthDX->Drc + HardSurface->Drc*_dx);
        }}

    if (SwitchHouses)
    {
        HStor = NewMap(0);
        IntercHouse = NewMap(0);
        DStor = NewMap(0);

        HouseCover = ReadMap(LDD,getvaluename("housecover"));
        RoofStore = ReadMap(LDD,getvaluename("roofstore"));
        calcValue(*RoofStore, 0.001, MUL); // from mm to m
        DrumStore = ReadMap(LDD,getvaluename("drumstore"));
    }
    else
        HouseCover = NewMap(0);

    if (SwitchLitter)
    {
        LCStor = NewMap(0);
        LInterc = NewMap(0);
        Litter = ReadMap(LDD,getvaluename("litter"));
        checkMap(*Litter, SMALLER, 0.0, "Litter cover fraction must be >= 0");
        checkMap(*Litter, LARGER, 1.0, "Litter cover fraction must be <= 1.0");
        LitterSmax = getvaluedouble("Litter interception storage");
    }

    fractionImperm = NewMap(0);
    FOR_ROW_COL_MV_L {
        double frac = 0;
        if (SwitchHouses && !SwitchRoadsystem)
            frac = HouseCover->Drc;
        if (SwitchRoadsystem && !SwitchHouses )
            frac = RoadWidthHSDX->Drc/_dx;
        if (SwitchRoadsystem && SwitchHouses )
            frac = RoadWidthHSDX->Drc/_dx + HouseCover->Drc;
        fractionImperm->Drc = std::min(std::max(0.0, frac), 1.0);
        // 0 is fully permeable, 1 = impermeable
    }}

    report(*fractionImperm,"fracimpermeable.map");

    GrassFraction = NewMap(0);
    if (SwitchGrassStrip)
    {
        KsatGrass = ReadMap(LDD,getvaluename("ksatgras"));
        PoreGrass = ReadMap(LDD,getvaluename("poregras"));
        CohGrass = ReadMap(LDD,getvaluename("cohgras"));
        GrassWidthDX = ReadMap(LDD,getvaluename("grasswidth"));
        copy(*GrassFraction, *GrassWidthDX);
        calcValue(*GrassFraction, _dx, DIV);
        StripN = getvaluedouble("Grassstrip Mannings n");
        FOR_ROW_COL_MV_L {
            if (GrassWidthDX->Drc != 0)
            {
                N->Drc = N->Drc*(1-GrassFraction->Drc)+StripN*GrassFraction->Drc;
                Cover->Drc = Cover->Drc*(1-GrassFraction->Drc) + 0.95*GrassFraction->Drc;
                LAI->Drc = LAI->Drc*(1-GrassFraction->Drc) + 5.0*GrassFraction->Drc;
            }
            //adjust mann N Cover and height
        }}
    }

    //## make shaded relief map for display.
    if (SwitchHouses && SwitchAddBuildingsDEM) {
        double AddBuildingFraction = getvaluedouble("Add Building fraction");
        double AddBuildingHeight = getvaluedouble("Add Building height");
        FOR_ROW_COL_MV_L {
            double dem = DEM->Drc;
            dem += HouseCover->Drc > AddBuildingFraction  ? AddBuildingHeight: 0.0;
            dem = RoadWidthDX->Drc > 0.1 ? DEM->Drc : dem;
            DEM->Drc = dem;
        }}
    }
}
//---------------------------------------------------------------------------
void TWorld::InitSoilInput(void)
{
    // safeguard for deleting, set to null pointer
    SwatreSoilModel = nullptr;
    SwatreSoilModelCrust = nullptr;
    SwatreSoilModelCompact = nullptr;
    SwatreSoilModelGrass = nullptr;

    LandUnit = ReadMap(LDD,getvaluename("landunit"));  //VJ 110107 added
    ThetaI1a = NewMap(0); // used for screen output
    ThetaI2a = NewMap(0); // for output, average soil layer 2

    if(SwitchOMCorrection)
        OMcorr = ReadMap(LDD,getvaluename("OMmap"));
    else
        OMcorr = NewMap(0);

    if(SwitchDensCorrection)
        DensFact = ReadMap(LDD,getvaluename("Densmap"));
    else
        DensFact  = NewMap(0);

    if (SwitchInfilCrust) {
        CrustFraction = ReadMap(LDD,getvaluename("crustfrc"));
        checkMap(*CrustFraction, LARGER, 1.0, "crust fraction cannot be more than 1");
    } else {
        CrustFraction = NewMap(0);
    }

    if (SwitchInfilCompact) {
        CompactFraction = ReadMap(LDD,getvaluename("compfrc"));
        checkMap(*CompactFraction, LARGER, 1.0, "compacted area fraction cannot be more than 1");
    } else {
        CompactFraction = NewMap(0);
    }
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if (CrustFraction->Drc + CompactFraction->Drc > 1.0) {
            CrustFraction->Drc = 1.0-CompactFraction->Drc;
        }
    }}

    //## infiltration data
    if(InfilMethod != INFIL_SWATRE)
    {
        nrSoilLayers = getvalueint("Nr input layers");
        //SwitchPsiUser = false; // MC - moved to the section with defaults

        SoilDepth1 = ReadMap(LDD,getvaluename("soildep1"));
        calcValue(*SoilDepth1, 1000, DIV);
        calcValue(*SoilDepth1, SD1Calibration, MUL);
        SoilDepth1init = NewMap(0);
        copy(*SoilDepth1init, *SoilDepth1);

        ThetaS1 = ReadMap(LDD,getvaluename("thetas1"));
        ThetaI1 = ReadMap(LDD,getvaluename("thetai1"));
        calcValue(*ThetaI1, thetaCalibration, MUL);
        calcMap(*ThetaI1, *ThetaS1, MIN);
        copy(*ThetaI1a, *ThetaI1);

        Ksat1 = ReadMap(LDD,getvaluename("ksat1"));

        ThetaR1 = NewMap(0);
        lambda1 = NewMap(0);
        psi1ae = NewMap(0);
        ThetaFC1 = NewMap(0);
        vgalpha1 = NewMap(0);
        vgn1 = NewMap(0);

        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            //bca1->Drc = 5.55*qPow(Ksat1->Drc,-0.114);  // old and untracable! and wrong
            // comes form CHARIM somehow
            //Saxton and Rawls 2006
            //  lambda1->Drc = 0.0384*log(Ksat1->Drc)+0.0626;
            //rawls et al., 1982

            double ks = log(std::min(1000.0,std::max(0.5,Ksat1->Drc))); //NOTE ln = log, log = log10
            // see the excel file with the regression equations in auxfiles
            // the regression fit has cm as output unit.
            lambda1->Drc = 0.0849*ks+0.159;
            lambda1->Drc = std::min(std::max(0.1,lambda1->Drc),0.7);

            psi1ae->Drc = exp( -0.3012*ks + 3.5164) * 0.01; // 0.01 to convert to m

            vgalpha1->Drc = 100*0.0119*exp(0.4657*ks);//(0.02*ks + 0.0095); // in m-1
            vgn1->Drc = 0.2656*ks + 1.1042;
            ThetaR1->Drc = 0.0673*exp(-0.238*log(ks));
            ThetaFC1->Drc = -0.0519*log(ks) + 0.3714;
            // NOTE alpha must have the reverse units of H. If H is in m, alpha must be in 1/m
        }}
        //report(*vgalpha1,"valpha.map");
        //report(*vgn1,"vn.map");

        if (SwitchPsiUser) {
            Psi1 = ReadMap(LDD,getvaluename("psi1"));
            calcValue(*Psi1, psiCalibration, MUL); //calibration of user input psi
            calcValue(*Psi1, 0.01, MUL); //cm to m
        } else {
            Psi1 = NewMap(0);
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                Psi1->Drc = exp(-0.3382*log(std::max(0.5,Ksat1->Drc)) + 3.3425)*0.01;
               // double ks = std::max(0.5,std::min(1000.0,log(Ksat1->Drc)));
                //double psiae = exp( -0.3012*ks + 3.5164) * 0.01; // 0.01 to convert to m
                Psi1->Drc = std::max(Psi1->Drc, psi1ae->Drc);
            }}
        // MC - note that this regression can cause very high psi values. In Rawls et al 1983 large ranges of possible values are reported!
        }
        calcValue(*Psi1, psiCalibration, MUL);
        calcValue(*Ksat1, ksatCalibration, MUL);
            // apply calibration after all empirical relations


        if (nrSoilLayers == 2) {
            SwitchTwoLayer = true;
            SwitchThreeLayer = false;
        }
        if (nrSoilLayers == 3) {
            SwitchTwoLayer = true;
            SwitchThreeLayer = true;
        }

        if (SwitchTwoLayer) {

            SoilDepth2 = ReadMap(LDD,getvaluename("soilDep2"));
            calcValue(*SoilDepth2, 1000, DIV);
            calcValue(*SoilDepth2, SD2Calibration, MUL);

            SoilDepth2init = NewMap(0);
            copy(*SoilDepth2init, *SoilDepth2);

            ThetaS2 = ReadMap(LDD,getvaluename("thetaS2"));
            ThetaI2 = ReadMap(LDD,getvaluename("thetaI2"));
            calcValue(*ThetaI2, thetaCalibration, MUL); //VJ 110712 calibration of theta
            calcMap(*ThetaI2, *ThetaS2, MIN); //VJ 110712 cannot be more than porosity
            copy(*ThetaI2a, *ThetaI2);

            Ksat2 = ReadMap(LDD,getvaluename("ksat2"));

            ThetaR2 = NewMap(0);
            lambda2 = NewMap(0);             // lambda brooks corey
            psi2ae = NewMap(0);

            ThetaFC2 = NewMap(0);
            vgalpha2 = NewMap(0);
            vgn2 = NewMap(0);
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                // regression eq from data from Saxton and rawls 2006, excel file
                double ks = log(std::min(1000.0,std::max(0.5,Ksat2->Drc)));
                vgalpha2->Drc = 100*0.0119*exp(0.4657*ks);//(0.02*ks + 0.0095); // in m-1
                vgn2->Drc = 0.2656*ks + 1.1042;

                lambda2->Drc = 0.0849*ks+0.159;
                lambda2->Drc = std::min(std::max(0.1,lambda2->Drc),0.7);
                psi2ae->Drc = exp( -0.3012*ks + 3.5164) * 0.01; // 0.01 to convert to m
                ThetaR2->Drc = 0.0673*exp(-0.238*log(ks));
                ThetaFC2->Drc = -0.0519*log(ks) + 0.3714;
            }}

            // wetting front psi
            if (SwitchPsiUser) {
                Psi2 = ReadMap(LDD,getvaluename("psi2"));
                //   calcValue(*Psi2, psiCalibration, MUL); //VJ 110712 calibration of psi
                calcValue(*Psi2, 0.01, MUL);
            } else {
                Psi2 = NewMap(0);
                #pragma omp parallel for num_threads(userCores)
                FOR_ROW_COL_MV_L {
                    Psi2->Drc = exp(-0.3382*log(std::max(0.5,Ksat2->Drc)) + 3.3425)*0.01;
                    //double ks = std::max(0.5,std::min(1000.0,log(Ksat2->Drc)));
                    //double psi2ae = exp( -0.3012*ks + 3.5164) * 0.01; // 0.01 to convert to m
                    Psi2->Drc = std::max(Psi2->Drc,psi2ae->Drc);
                }}
            }
            calcValue(*Psi2, psiCalibration, MUL);
            calcValue(*Ksat2, ksat2Calibration, MUL);
        }

        if (SwitchThreeLayer)
        {
            SoilDepth3 = ReadMap(LDD,getvaluename("soilDep3"));
            calcValue(*SoilDepth3, 1000, DIV);
            //calcValue(*SoilDepth2, SD2Calibration, MUL);

            SoilDepth3init = NewMap(0);
            copy(*SoilDepth3init, *SoilDepth3);

            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                if (SoilDepth3->Drc < 0)
                {
                    ErrorString = QString("SoilDepth3 values < 0 at row %1, col %2").arg(r).arg(c);
                    throw 1;
                }
            }}

            ThetaS3 = ReadMap(LDD,getvaluename("thetaS3"));
            ThetaI3 = ReadMap(LDD,getvaluename("thetaI3"));
            ThetaI3a = NewMap(0); // for output, average soil layer 2
            calcValue(*ThetaI3, thetaCalibration, MUL);
            calcMap(*ThetaI3, *ThetaS3, MIN);
            copy(*ThetaI3a, *ThetaI3);

            Ksat3 = ReadMap(LDD,getvaluename("ksat3"));

            ThetaR3 = NewMap(0);
            lambda3 = NewMap(0);             // lambda brooks corey
            psi3ae = NewMap(0);
            ThetaFC3 = NewMap(0);
            vgalpha3 = NewMap(0);
            vgn3 = NewMap(0);
            FOR_ROW_COL_MV_L {
                // regression eq from data from Saxton and rawls 2006, excel file
                double ks = log(std::min(1000.0,std::max(0.5,Ksat3->Drc)));
                vgalpha3->Drc = 100*0.015*exp(0.3816*ks);//(0.02*ks + 0.0095); // in m-1
                vgn3->Drc = 0.2656*ks + 1.1042;

                lambda3->Drc = 0.0849*ks+0.159;
                lambda3->Drc = std::min(std::max(0.1,lambda3->Drc),0.7);

                psi3ae->Drc = exp( -0.3012*ks + 3.5164) * 0.01; // 0.01 to convert to m
                ThetaR3->Drc = 0.0673*exp(-0.238*log(ks));
                ThetaFC3->Drc = -0.0519*log(ks) + 0.3714;
            }}

            // wetting front psi
            if (SwitchPsiUser) {
                Psi3 = ReadMap(LDD,getvaluename("psi3"));
                //   calcValue(*Psi3, psiCalibration, MUL); //VJ 110712 calibration of psi
                calcValue(*Psi3, 0.01, MUL);
            } else {
                Psi3 = NewMap(0);
                #pragma omp parallel for num_threads(userCores)
                FOR_ROW_COL_MV_L {
                    Psi3->Drc = exp(-0.3382*log(Ksat2->Drc) + 3.3425)*0.01;
                    //double ks = std::max(0.5,std::min(1000.0,log(Ksat3->Drc)));
                    //double psi3ae = exp( -0.3012*ks + 3.5164) * 0.01; // 0.01 to convert to m
                    Psi3->Drc = std::max(Psi2->Drc,psi3ae->Drc);
                }}
            }
            calcValue(*Psi3, psiCalibration, MUL);
            calcValue(*Ksat3, ksat3Calibration, MUL);
        }

        if (SwitchInfilCrust) {
            KsatCrust = ReadMap(LDD,getvaluename("ksatcrst"));
            calcValue(*KsatCrust, ksatCalibration, MUL);
			//DO THIS, else inconsistency, and Ksat can be smaller than ksatcrust

            PoreCrust = ReadMap(LDD,getvaluename("porecrst"));
        } else {
            KsatCrust = NewMap(0);
            PoreCrust = NewMap(0);
        }

        if (SwitchInfilCompact)
        {
            KsatCompact = ReadMap(LDD,getvaluename("ksatcomp"));
            calcValue(*KsatCompact, ksatCalibration, MUL);
            //DO THIS, else inconsistency, Ksat can be smaller than ksatcomp
            PoreCompact = ReadMap(LDD,getvaluename("porecomp"));                       
        } else {
            KsatCompact = NewMap(0);
            PoreCompact = NewMap(0);
        }

        // check if compacted porosity is smaller than initial theta1
        if (SwitchInfilCompact) {
            double cnt = 0;
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                if(PoreCompact->Drc*CompactFraction->Drc+(1-CompactFraction->Drc)*ThetaS1->Drc < ThetaI1->Drc)
                    cnt+=1.0;
            }}
            if (cnt > 0) {
                ErrorString = QString("WARNING: Compacted porosity is smaller than initial moisture content in %1% of the cells, these cells will be seen as impermeable.").arg(cnt/nrCells*100);
                DEBUG(ErrorString);
             }
        }

        // check if crusted porosity is smaller than initial theta1
        if (SwitchInfilCrust) {
            double cnt = 0;
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                if(PoreCrust->Drc*CrustFraction->Drc+(1-CrustFraction->Drc)*ThetaS1->Drc < ThetaI1->Drc)
                    cnt+=1.0;
            }}
            if (cnt > 0) {
                ErrorString = QString("WARNING: Porosity in crusted cells is smaller than initial moisture content in %1% of the cells, these cells will be seen as impermeable.").arg(cnt/nrCells*100);
                DEBUG(ErrorString);
            }
        }

        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            Ksat1->Drc *= _dt/3600000.0; // mm/h to m
            if (SwitchTwoLayer)
                Ksat2->Drc *= _dt/3600000.0;
            if (SwitchThreeLayer)
                Ksat3->Drc *= _dt/3600000.0;
            if (SwitchInfilCrust)
                KsatCrust->Drc *= _dt/3600000.0;
            if (SwitchInfilCompact)
                KsatCompact->Drc *= _dt/3600000.0;
        }}


    } // not swatre

    // SWATRE infiltration read maps and structures
    if (InfilMethod == INFIL_SWATRE) {

        inith = new QVector<cTMap*>();

        Hswatre = NewMap(0);
        // read all Swatre profile maps
        ProfileID = ReadMap(LDD,getvaluename("profmap"));

        // if (SwitchDumphead)
        //     SwatreOutput = ReadMap(LDD,getvaluename("swatreout"));

        if (SwitchGrassStrip)
            ProfileIDGrass = ReadMap(LDD,getvaluename("profgrass"));

        if (SwitchInfilCrust)
            ProfileIDCrust = ReadMap(LDD,getvaluename("profcrst"));

        if (SwitchInfilCompact)
            ProfileIDCompact = ReadMap(LDD,getvaluename("profcomp"));

        // read the swatre tables and make the information structure ZONE etc
        // this does not make the profile information
        ReadSwatreInputNew();
    }
}
//---------------------------------------------------------------------------
void TWorld::InitBoundary(void)
{
    BoundaryQ = 0;
    BoundaryQs = 0;

    // make a 1 cell edge around the domain, used to determine flood at the edge
    DomainEdge = NewMap(0);
    for (int r = 1; r < _nrRows-1; r++)
        for (int c = 1; c < _nrCols-1; c++)
            if(!pcr::isMV(LDD->data[r][c]))
            {
                if (DomainEdge->Drc == 0 && pcr::isMV(LDD->data[r-1][c  ])) DomainEdge->Drc = 8; // use ldd logic for clarity
                if (DomainEdge->Drc == 0 && pcr::isMV(LDD->data[r+1][c  ])) DomainEdge->Drc = 2;
                if (DomainEdge->Drc == 0 && pcr::isMV(LDD->data[r  ][c-1])) DomainEdge->Drc = 4;
                if (DomainEdge->Drc == 0 && pcr::isMV(LDD->data[r  ][c+1])) DomainEdge->Drc = 6;
            }
    // if ldd touches the edge
    FOR_ROW_COL_MV {
        if(r == 0)          DomainEdge->Drc = 8;
        if(r == _nrRows-1)  DomainEdge->Drc = 2;
        if(c == 0)          DomainEdge->Drc = 4;
        if(c == _nrCols-1)  DomainEdge->Drc = 6;
    }

    FlowBoundary = NewMap(0);

    if(FlowBoundaryType == 1) // potential outflow everywhere
    {
        // determine dynamically in function K2DDEMA
        // for flood DomainEdge is used
        copy( *FlowBoundary, *DomainEdge);
    }
    if (FlowBoundaryType == 2 ) // user defined outflow (0 close, >0 outflow)
    {
        FlowBoundary = ReadMap(LDD,getvaluename("flowboundary"));
        // use flowboundary for domainedge
        FOR_ROW_COL_MV {
            if (FlowBoundary->Drc > 0)
                FlowBoundary->Drc = DomainEdge->Drc;
            // this sets the cells to 1 row/col instead of a user sloppy digitizing
        }
    }
    // always set outlet to 1
    FOR_ROW_COL_MV {
        if(LDD->Drc == 5)
            FlowBoundary->Drc = 1;
    }
    if (SwitchIncludeChannel) {
        FOR_ROW_COL_MV_CH {
            if(LDDChannel->Drc == 5)
                FlowBoundary->Drc = 1;
        }
    }
 //  report(*FlowBoundary, "bound.map");
 //  report(*DomainEdge, "edge.map");

}
//---------------------------------------------------------------------------
// read and Intiialize all channel variables and maps
void TWorld::InitChannel(void)
{
    // channel vars and maps that must be there even if channel is switched off
    ChannelVolTotmm = 0;
    ChannelSedTot = 0;
    ChannelDepTot = 0;
    ChannelDetTot = 0;
    BaseFlowTotmm = 0;
    PeakFlowTotmm = 0;
    QuserInTot = 0;

    if(!SwitchIncludeChannel)
        return;

    //SedToChannel = NewMap(0);
    //ChannelFlowWidth = NewMap(0);
    ChannelWidthMax = NewMap(0);
    ChannelWaterVol = NewMap(0);
    ChannelQ = NewMap(0);
    ChannelQb = NewMap(0); //baseflow
    ChannelQn = NewMap(0);
    ChannelQntot = NewMap(0);
    ChannelSed = NewMap(0);
    ChannelQs = NewMap(0);
    ChannelQsn = NewMap(0);
    ChannelQsr = NewMap(0);
    ChannelV = NewMap(0);//
    ChannelU = NewMap(0);//
    ChannelWH = NewMap(0);
    //Channelq = NewMap(0);//
    ChannelAlpha = NewMap(0);//
    ChannelDX = NewMap(0);
    ChannelInfilVol = NewMap(0);

    maxChannelflow = NewMap(0);//
    maxChannelWH = NewMap(0);//

    //## channel maps
    LDDChannel = InitMaskChannel(getvaluename("lddchan"));
    // LDDChannel is the mask for channels

    FOR_ROW_COL_MV_CH {
        if (LDDChannel->Drc == 0)
            SET_MV_REAL8(&LDDChannel->Drc);
    }

    nrValidCellsCH = 0;
    FOR_ROW_COL_MV_CH {
        nrValidCellsCH++;
    }
    //crch_ = (LDD_COOR*) malloc(sizeof(LDD_COOR)*nrValidCellsCH);
    //long i = 0;
    FOR_ROW_COL_MV_CH {
        LDD_COOR newcr;
        newcr.r = r;
        newcr.c = c;
        // LDD_COOR *newcr = new LDD_COOR;
        // newcr->r = r;
        // newcr->c = c;
        crch_ << newcr;
    }
    crlinkedlddch_= MakeLinkedList(LDDChannel);

    crlddch5_.clear();
    FOR_ROW_COL_MV_CH {
        if (LDDChannel->Drc == 5) {
            LDD_COOR newcr;
            newcr.r = r;
            newcr.c = c;
            crlddch5_ << newcr;
        }
    }
    nrValidCellsLDDCH5 = crlddch5_.size();


    // for 1D or 2D overland flow: channel outlet points are checked, leading
    // FOR_ROW_COL_MV_CH {
    //     if(Outlet->Drc > 0 && LDDChannel->Drc != 5)
    //     {
    //         //qDebug() << r << c << LDDChannel->Drc << Outlet->Drc;
    //         ErrorString = QString("Outlet points (outlet.map) do not coincide with Channel LDD endpoints: %1 %2.").arg(Outlet->Drc).arg(LDDChannel->Drc);
    //         throw 1;
    //     }
    // }

    ChannelWidth = ReadMap(LDDChannel, getvaluename("chanwidth")); // bottom width in m

    //     ChannelWidth->checkMap(LARGER, _dx, "Channel width must be smaller than cell size");
    //ChannelWidth->checkMap(SMALLEREQUAL, 0, "Channel width must be larger than 0 in channel cells");
    ChannelDepth = ReadMap(LDDChannel, getvaluename("chandepth"));
    cover(*ChannelWidth, *LDD,0);
    cover(*ChannelDepth, *LDD,0);

    ChannelWidthO = NewMap(0);
 //   ChannelDepthO = NewMap(0);

    FOR_ROW_COL_MV_CH
    {
        ChannelWidthO->Drc = ChannelWidth->Drc;
      //  ChannelDepthO->Drc = ChannelDepth->Drc;

        SwitchChannelAdjustCHW = true;
        if (SwitchChannelAdjustCHW && ChannelWidth->Drc  > 0.95* _dx) {
            ChannelWidth->Drc = 0.95*_dx;
            ChannelDepth->Drc *= ChannelWidthO->Drc/ChannelWidth->Drc; //(0.95*_dx);
        }

        if (ChannelWidth->Drc <= 0)
        {
            ErrorString = QString("Map %1 contains channel cells with width = 0").arg(getvaluename("chanwidth"));
            throw 1;
        }
    }

    ChannelSide = ReadMap(LDDChannel, getvaluename("chanside"));
    ChannelGrad = ReadMap(LDDChannel, getvaluename("changrad"));
    checkMap(*ChannelGrad, LARGER, 1.0, "Channel Gradient must be SINE of slope angle (not tangent)");
    //calcValue(*ChannelGrad, 0.001, MAX);
    //VJ 171002 better to check and set Q to 0 in the code
    ChannelN = ReadMap(LDDChannel, getvaluename("chanman"));

    cover(*ChannelGrad, *LDD, 0);
    cover(*ChannelSide, *LDD, 0);
    cover(*ChannelN, *LDD, 0);

    ChannelNcul = NewMap(0);
    ChannelQSide = NewMap(0);

    calcValue(*ChannelN, ChnCalibration, MUL);
    copy(*ChannelNcul, *ChannelN);

    if (SwitchChannelInfil)
    {
        ChannelKsat = ReadMap(LDDChannel, getvaluename("chanksat"));
        cover(*ChannelKsat, *LDD, 0);
        calcValue(*ChannelKsat, ChKsatCalibration, MUL);
        // ChannelStore = NewMap(0.050); // 10 cm deep * 0.5 porosity
        // store not used?
    }

    if (SwitchCulverts) {

        ChannelMaxQ = ReadMap(LDDChannel, getvaluename("chanmaxq"));
        cover(*ChannelMaxQ, *LDD,0);
        ChannelMaxAlpha = NewMap(0);

        FOR_ROW_COL_MV_CHL {
            if (ChannelMaxQ->Drc > 0) {
                // ChannelWidth->Drc = 0.3;
                // ChannelWidthO->Drc = 0.3;
                // ChannelDepth->Drc = 0.3;
                ChannelN->Drc = 0.015;
                ChannelGrad->Drc = 0.002;

                ChannelMaxAlpha->Drc = (ChannelWidth->Drc*ChannelDepth->Drc)/std::pow(ChannelMaxQ->Drc, 0.6);
            }
        }}

        // there can be no side inflow in a culvert (which is actually not true!)
        for (int i = 0; i < crlinkedlddch_.size(); i++) {
            int c = crlinkedlddch_.at(i).c;
            int r = crlinkedlddch_.at(i).r;
            if (ChannelMaxQ->Drc > 0) {
                LDD_COORIN hoi = crlinkedlddch_.at(i);
                hoi.ldd *= -1;
                crlinkedlddch_.replace(i, hoi) ;
            }
        }

    } else {
        ChannelMaxQ = NewMap(-1);
        ChannelMaxAlpha = NewMap(-1);
    }


    FOR_ROW_COL_MV_CH
    {
        ChannelWidthMax->Drc = ChannelWidth->Drc; // not used!
        // make always a rectangular channel
        ChannelDX->Drc = _dx/cos(asin(Grad->Drc)); // same as DX else mass balance problems
    }

    if (SwitchGWflow) {

        LDDbaseflow = ReadMap(LDD, getvaluename("lddbase"));
        crlinkedlddbase_= MakeLinkedList(LDDbaseflow);

        BaseflowL = ReadMap(LDDChannel, getvaluename("basereach")); // bottom width in m
        FOR_ROW_COL_MV_L {
            BaseflowL->Drc = pow(_dx/BaseflowL->Drc,GW_slope*2);
        }}

        GWVol = NewMap(0); //ReadMap(LDD, getvaluename("gwlevel")); // bottom width in m
        Qbase = NewMap(0);
        GWWH = NewMap(0);
        GWU = NewMap(0);
        GWV = NewMap(0);
        GWN = NewMap(0);
        GWWHmax = NewMap(0);

        GWdeep = NewMap(0);
        GWrecharge = NewMap(0);
        GWout = NewMap(0);
        GWz = NewMap(0);
        GWgrad = NewMap(0);

        FOR_ROW_COL_MV_L {
            //GWz->Drc = DEM->Drc - SoilDepth1->Drc - (SwitchTwoLayer ? SoilDepth2->Drc : 0.0);
            if (SwitchTwoLayer)
                GWz->Drc = DEM->Drc - SoilDepth2->Drc;
            else
                GWz->Drc = DEM->Drc - SoilDepth1->Drc;
            tm->Drc = SoilDepth2->Drc;
        }}
        Average3x3(*GWz, *LDD, false);

        Average3x3(*tm, *LDD, false);
        FOR_ROW_COL_MV_L {
            GWN->Drc = 0.1+pow(tm->Drc,2.0/3.0)*qSqrt(0.1)/(Ksat2->Drc/3600000/_dt);
        }}

        Average3x3(*GWN, *LDD, false);
        report(*GWN,"gwn.map");

    }

    if(SwitchErosion) {
        TotalChanDetMap = NewMap(0);
        TotalChanDepMap = NewMap(0);
        ChannelDetFlow = NewMap(0);
        ChannelDep = NewMap(0);
        ChannelSSSed = NewMap(0);
        ChannelSSConc = NewMap(0);
        ChannelSSTC = NewMap(0);
        ChannelSSDepth = NewMap(0);
        ChannelQSSs = NewMap(0);
        ChannelQSSsn = NewMap(0);
        if (SwitchUse2Phase) {
            ChannelBLSed = NewMap(0);
            ChannelBLConc = NewMap(0);
            ChannelBLTC = NewMap(0);
            ChannelBLDepth = NewMap(0);
            ChannelQBLs = NewMap(0);
            ChannelQBLsn = NewMap(0);
        }

        ChannelConc = NewMap(0);
        ChannelTC = NewMap(0);
        ChannelY = NewMap(0);

        D50CH = NewMap(0);
        SwitchD50CHavg = false;
        if(SwitchD50CHavg) {
            double D50ch = mapAverage(*D50);
            FOR_ROW_COL_MV_CHL {
                D50CH->Drc = D50ch;
            }}
        } else {
            FOR_ROW_COL_MV_CHL {
                D50CH->Drc = D50->Drc;
            }}
        }

        if (SwitchUse2Phase) {
            D90CH = NewMap(0);
            if(SwitchD50CHavg) {
                double D90ch = mapAverage(*D90);
                FOR_ROW_COL_MV_CHL {
                    D90CH->Drc = D90ch;
                }}
            } else {
                FOR_ROW_COL_MV_CHL {
                    D90CH->Drc = D90->Drc;
                }}
            }
        }

        COHCHCalibration = 1.0;
        UcrCHCalibration = 1.0;

        ChannelCohesion = ReadMap(LDDChannel, getvaluename("chancoh"));
        COHCHCalibration = getvaluedouble("Cohesion Channel calibration");
        UcrCHCalibration = getvaluedouble("Ucr Channel calibration");
        DirectEfficiency = getvaluedouble("Direct efficiency channel");

        FOR_ROW_COL_MV_CHL {
            if (ChannelCohesion->Drc > 0)
                ChannelCohesion->Drc *= COHCHCalibration;

            if (ChannelCohesion->Drc == 0) {
                ChannelY->Drc = 1.0;
            } else {
                if (SwitchEfficiencyDETCH == 1)
                    ChannelY->Drc = std::min(1.0, 1.0/(0.89+0.56*fabs(ChannelCohesion->Drc)));
                else
                    if (SwitchEfficiencyDETCH == 2)
                        ChannelY->Drc = std::min(1.0, 0.79*exp(-0.85*fabs(ChannelCohesion->Drc)));
                    else
                        if (SwitchEfficiencyDETCH == 3)
                            ChannelY->Drc = std::min(1.0, 1.0/(2.0*fabs(ChannelCohesion->Drc)));
                        else
                            if (SwitchEfficiencyDETCH == 4)
                                ChannelY->Drc = DirectEfficiency;
               }
            if (ChannelCohesion->Drc < 0)
                ChannelY->Drc = 0;
        }}
    }
}
//---------------------------------------------------------------------------
void TWorld::InitFlood(void)
{
    FloodSedTot = 0;
    FloodDepTot = 0;
    FloodDetTot = 0;

    prepareFlood = true;
    Qflood = NewMap(0);
    hmxWH = NewMap(0);
    FloodWaterVol = NewMap(0);
    RunoffWaterVol = NewMap(0);
    floodTimeStart = NewMap(0);

    hs = NewMap(0);
    Uflood = NewMap(0);
    Vflood = NewMap(0);
    hmx = NewMap(0);
    hmxflood = NewMap(0);
    FloodDomain = NewMap(0);

    floodHmxMax = NewMap(0);//
    floodVMax = NewMap(0);//
    floodVHMax = NewMap(0);//
    floodTime = NewMap(0);//
    FloodDT = NewMap(0);

    iter_n = 0;

    dcr_.clear();
    if (Switch2DDiagonalFlow)
        DiagonalFlowDEM();

    // OBSOLETE
    // if (!SwitchSWOFopen) {
    //FloodT = NewMap(0);
    //     vs = NewMap(0);
    //     us = NewMap(0);
    //     z1r = NewMap(0);
    //     z1l = NewMap(0);
    //     z2r = NewMap(0);
    //     z2l = NewMap(0);
    //     h1r = NewMap(0);
    //     h1l = NewMap(0);
    //     h2r = NewMap(0);
    //     h2l = NewMap(0);
    //     v1r = NewMap(0);
    //     v1l = NewMap(0);
    //     v2r = NewMap(0);
    //     v2l = NewMap(0);
    //     u1r = NewMap(0);
    //     u1l = NewMap(0);
    //     u2r = NewMap(0);
    //     u2l = NewMap(0);

    //     delzc1 = NewMap(0);
    //     delzc2 = NewMap(0);

    //     f1 = NewMap(0);
    //     f2 = NewMap(0);
    //     f3 = NewMap(0);
    //     cflx = NewMap(0);
    //     cfly = NewMap(0);
    //     g1 = NewMap(0);
    //     g2 = NewMap(0);
    //     g3 = NewMap(0);
    //     f1o = NewMap(0);
    //     f2o = NewMap(0);
    //     f3o = NewMap(0);
    //     g1o = NewMap(0);
    //     g2o = NewMap(0);
    //     g3o = NewMap(0);
    //     h1d = NewMap(0);
    //     h1g = NewMap(0);
    //     h2d = NewMap(0);
    //     h2g = NewMap(0);
    //     delz1 = NewMap(0);
    //     delz2 = NewMap(0);
    //     prepareFloodZ(DEM);
    // }

    if (SwitchErosion) {
        BLDepthFlood = NewMap(0);
        SSDepthFlood = NewMap(0);
        BLFlood = NewMap(0);
        BLCFlood = NewMap(0);
        BLTCFlood = NewMap(0);
        BLDetFlood = NewMap(0);

        SSFlood = NewMap(0);
        SSCFlood = NewMap(0);
        SSTCFlood = NewMap(0);
        SSDetFlood = NewMap(0);
        DepFlood = NewMap(0);
    }
}
//---------------------------------------------------------------------------
void TWorld::DiagonalFlowDEM()
{
    Fill(*tma,0);
    Fill(*tmb,0);
    FOR_ROW_COL_MV_L {
        double Z = DEM->Drc;
        double z_x1 =  c > 0 && !MV(r,c-1)         ? DEM->data[r][c-1] : Z;
        double z_x2 =  c < _nrCols-1 && !MV(r,c+1) ? DEM->data[r][c+1] : Z;
        double z_y1 =  r > 0 && !MV(r-1,c)         ? DEM->data[r-1][c] : Z;
        double z_y2 =  r < _nrRows-1 && !MV(r+1,c) ? DEM->data[r+1][c] : Z;


        double z_x11 =  c > 0 && r > 0 && !MV(r-1,c-1)         ? DEM->data[r-1][c-1] : Z;
        double z_x21 =  c > 0 && r < _nrRows-1 && !MV(r+1,c-1) ? DEM->data[r+1][c-1] : Z;
        double z_y11 =  r > 0 && c < _nrCols-1 && !MV(r-1,c+1)         ? DEM->data[r-1][c+1] : Z;
        double z_y21 =  r < _nrRows-1 && c < _nrCols-1 && !MV(r+1,c+1) ? DEM->data[r+1][c+1] : Z;

/*
        //note: true blockage if the diagonal cells are higher than the centre cell will not be flagged
        // left blockage
        if (z_x1 > Z+F_pitValue && z_y1 > Z+F_pitValue && z_y2 > Z+F_pitValue) {
            bool z1 = z_x11 < Z+F_pitValue;
            bool z2 = z_y11 < Z+F_pitValue;
            if (z1 && z2) {
                if (z_x11 < z_y11)
                    z2 = false;
            }

            if(z1) tma->Drc = 7;
            if(z2) tma->Drc = 1;
        }
        // right blockage
        if (z_x2 > Z+F_pitValue && z_y1 > Z+F_pitValue && z_y2 > Z+F_pitValue) {
            bool z1 = z_y11 < Z+F_pitValue;
            bool z2 = z_y21 < Z+F_pitValue;
            if (z1 && z2) {
                if (z_y11 < z_y21)
                    z2 = false;
            }
            if(z1) tma->Drc = 9;
            if(z2) tma->Drc = 3;
        }
        // upper blockage
        if (z_y1 > Z+F_pitValue && z_x1 > Z+F_pitValue && z_x2 > Z+F_pitValue) {
            bool z1 = z_x11 < Z+F_pitValue;
            bool z2 = z_x21 < Z+F_pitValue;
            if (z1 && z2) {
                if (z_x11 < z_x21)
                    z2 = false;
            }
            if(z1) tma->Drc = 7;
            if(z2) tma->Drc = 9;
        }
        //lower blockage
        if (z_y2 > Z+F_pitValue && z_x1 > Z+F_pitValue && z_x2 > Z+F_pitValue) {
            bool z1 = z_x21 < Z+F_pitValue;
            bool z2 = z_y21 < Z+F_pitValue;
            if (z1 && z2) {
                if (z_x21 < z_y21)
                    z2 = false;
            }
            if(z1) tma->Drc = 1;
            if(z2) tma->Drc = 3;
        }

// ldd map based:
        int ldd = (int) LDD->Drc;
        if (z_x1 > Z+F_pitValue && z_y1 > Z+F_pitValue && z_y2 > Z+F_pitValue) {
            if (ldd == 1 || ldd == 7)
                tmb->Drc = ldd;
        }
        if (z_x2 > Z+F_pitValue && z_y1 > Z+F_pitValue && z_y2 > Z+F_pitValue) {
            if (ldd == 3 || ldd == 9)
                tmb->Drc = ldd;
        }
        if (z_y1 > Z+F_pitValue && z_x1 > Z+F_pitValue && z_x2 > Z+F_pitValue) {
            if (ldd == 7 || ldd == 9)
                tmb->Drc = ldd;
        }
        if (z_y2 > Z+F_pitValue && z_x1 > Z+F_pitValue && z_x2 > Z+F_pitValue) {
            if (ldd == 1 || ldd == 3)
                tmb->Drc = ldd;
        }
*/
        int ldd = (int) LDD->Drc;
        if (z_y1 > Z+F_pitValue && z_y2 > Z+F_pitValue && z_x1 > Z+F_pitValue && z_x2 > Z+F_pitValue) {
            if (ldd == 1 || ldd == 3 || ldd == 7 || ldd == 9)
                tma->Drc = ldd;
        }

//          DEMdz->Drc = tma->Drc;
        // do not include channels, channels will do the outflow
        if(SwitchIncludeChannel && ChannelWidth->Drc > 0) {
            tma->Drc = 0;
            tmb->Drc = 0;
        }

        // make a list of pits
        if (tma->Drc > 0) {
            LDD_COORi dclrc;
            dclrc.r = r;
            dclrc.c = c;
            dclrc.ldd = (int) tma->Drc;
            dcr_ << dclrc;
        }
    }}
    //report(*tma,"diagflow.map");
}
//---------------------------------------------------------------------------
void TWorld::CorrectDEM(cTMap *h, cTMap * g)
{
    QList <double> zmin;
    Fill(*tma,-9999);
    Fill(*tmb,0);
    FOR_ROW_COL_MV_L {
        double Z = h->Drc;
        double z_x1 =  c > 0 && !MV(r,c-1)         ? h->data[r][c-1] : Z;
        double z_x2 =  c < _nrCols-1 && !MV(r,c+1) ? h->data[r][c+1] : Z;
        double z_y1 =  r > 0 && !MV(r-1,c)         ? h->data[r-1][c] : Z;
        double z_y2 =  r < _nrRows-1 && !MV(r+1,c) ? h->data[r+1][c] : Z;
        double z_x11 =  c > 0 && r > 0 && !MV(r-1,c-1)         ? h->data[r-1][c-1] : Z;
        double z_x21 =  c > 0 && r < _nrRows-1 && !MV(r+1,c-1) ? h->data[r+1][c-1] : Z;
        double z_y11 =  r > 0 && c < _nrCols-1 && !MV(r-1,c+1)         ? h->data[r-1][c+1] : Z;
        double z_y21 =  r < _nrRows-1 && c < _nrCols-1 && !MV(r+1,c+1) ? h->data[r+1][c+1] : Z;

        zmin.clear();
        zmin << z_x1 << z_x2 << z_y1 << z_y2 << z_x11 << z_y11 << z_x21 << z_y21;
        std::sort(zmin.begin(), zmin.end());
        if (Z < zmin.at(0)) {
           tma->Drc = zmin.at(0);
        }
    }}
    FOR_ROW_COL_MV_L {
        if (tma->Drc > -9999) {
            tmb->Drc = tma->Drc - h->Drc + 0.001*_dx;
            h->Drc = tma->Drc-0.001*_dx;
            g->Drc = 0.001;
        }
    }}
    //report(*tmb, "dempits.map");
}
//---------------------------------------------------------------------------
void TWorld::InitErosion(void)
{
//qDebug() << "hoi"; //SwitchSlopeStability ||
//if (SwitchErosion) {
//        COHCalibration = getvaluedouble("Cohesion calibration");
//        Cohesion = ReadMap(LDD,getvaluename("coh"));
//        RootCohesion = ReadMap(LDD,getvaluename("cohadd"));
//        FOR_ROW_COL_MV_L {
//            if (RootCohesion->Drc < 0) // root cohesion can be used to avoid surface erosion base don land use
//                CohesionSoil->Drc = -1;
//            else
//                CohesionSoil->Drc = COHCalibration*(Cohesion->Drc + Cover->Drc*RootCohesion->Drc);
//            if (SwitchGrassStrip)
//                CohesionSoil->Drc = CohesionSoil->Drc  *(1-GrassFraction->Drc) + GrassFraction->Drc * CohGrass->Drc;
//        }}
//    }


    if (!SwitchErosion)
        return;

    PlantHeight = ReadMap(LDD,getvaluename("CH"));
    checkMap(*PlantHeight, SMALLER, 0.0, "Cover fraction must be >= 0");

    StoneFraction  = ReadMap(LDD,getvaluename("stonefrc"));

  //  LandUnit = ReadMap(LDD,getvaluename("landunit"));  //VJ 110107 added

    COHCalibration = getvaluedouble("Cohesion calibration");
    Cohesion = ReadMap(LDD,getvaluename("coh"));

    RootCohesion = ReadMap(LDD,getvaluename("cohadd"));

    ASCalibration = getvaluedouble("Aggregate stability calibration");
    AggrStab = ReadMap(LDD,getvaluename("AggrStab"));

    D50 = ReadMap(LDD,getvaluename("D50"));
    //SwitchNeedD90 = SwitchErosion && (SwitchChannelFlood || (SwitchUse2Phase && !R_BL_Method == RGOVERS) || (SwitchEstimateGrainSizeDistribution && SwitchUseGrainSizeDistribution);
    if(SwitchUse2Phase)// && !SwitchUseGrainSizeDistribution)
    {
        D90 = ReadMap(LDD,getvaluename("D90"));
    }

    FOR_ROW_COL_MV {
        D50->Drc = D50->Drc *gsizeCalibrationD50;
        if(SwitchUse2Phase)// && !SwitchUseGrainSizeDistribution)
        {
            D90->Drc = D90->Drc * gsizeCalibrationD90;
        }
    }

    cgovers = NewMap(0);
    dgovers = NewMap(0);
    FOR_ROW_COL_MV {
        cgovers->Drc = pow((D50->Drc+5)/0.32, -0.6);
        dgovers->Drc = pow((D50->Drc+5)/300, 0.25);
    }

    SedimentFilter = NewMap(0);
    if (SwitchSedtrap)
    {
        SedMaxVolume = ReadMap(LDD,getvaluename("sedretmax"));
        SedTrapN = getvaluedouble("Sediment Trap Mannings n");
        FOR_ROW_COL_MV {
            if (SedMaxVolume->Drc > 0)
                N->Drc = SedTrapN;
        }
    }
    else {
        SedTrapN = 0;
        SedMaxVolume = NewMap(0);
    }

    //default
    R_BL_Method = FSRIJN;
    R_SS_Method = FSRIJN;
    FS_BL_Method = FSRIJN;
    FS_SS_Method = FSGOVERS;

    FS_SS_Method = getvalueint("Flooding SS method")-1;
    qDebug() << FS_SS_Method ;
    FS_BL_Method = getvalueint("Flooding BL method")-1;
    R_SS_Method  = getvalueint("River SS method")-1;
    R_BL_Method  = getvalueint("River BL method")-1;

    FS_SigmaDiffusion = getvaluedouble("Sigma diffusion");
    R_SigmaDiffusion = getvaluedouble("Sigma diffusion"); // same diffusion for river and OF

    SVCHCalibration = 1.0;
    SVCHCalibration = getvaluedouble("SV calibration");


//    if (SwitchUse2Phase && SwitchUseGrainSizeDistribution) {
//        R_BL_Method = FSWUWANGJIA;
//        R_SS_Method = FSWUWANGJIA;  // ignore because it has to be 3 when 2 layer and graisizedist
//        FS_BL_Method = FSWUWANGJIA;
//        FS_SS_Method = FSWUWANGJIA;
//    }
//    else
//        if(!SwitchUse2Phase && !SwitchUseGrainSizeDistribution) {
//            R_BL_Method = FSRIJN;     // if single layer and no grainsize = simple erosion, then govers
//            R_SS_Method = FSGOVERS;
//            FS_BL_Method = FSRIJN;
//            FS_SS_Method = FSGOVERS;
//        }

    Qs = NewMap(0);
    Qsn = NewMap(0);

    DetSplashTot = 0;
    DetFlowTot = 0;
    DepTot = 0;
    DetTot = 0;
    SoilLossTot = 0;
    SoilLossTot_dt = 0;
    SedTot = 0;

    TotalSoillossMap = NewMap(0);
    TotalSed = NewMap(0);
    TotalConc = NewMap(0);

    DETFlow = NewMap(0);
    DETSplash = NewMap(0);
    DETSplashCum = NewMap(0);
    DETFlowCum = NewMap(0);
    DEP = NewMap(0);
    DEPCum = NewMap(0);
    //DEPBLCum = NewMap(0);
    Sed = NewMap(0);
    TC = NewMap(0);
    Conc = NewMap(0);

    SettlingVelocitySS = NewMap(0);
    SettlingVelocityBL = NewMap(0);
    CohesionSoil = NewMap(0);
    Y = NewMap(0);
    //splashb = NewMap(0);

    FOR_ROW_COL_MV {
        SettlingVelocitySS->Drc = GetSV(D50->Drc/gsizeCalibrationD50);
        if (SwitchUse2Phase)
            SettlingVelocityBL->Drc = GetSV(D90->Drc/gsizeCalibrationD90);
    }

    SplashStrength = NewMap(0);

   // qDebug() << "SwitchEfficiencyDET" <<SwitchEfficiencyDET;

    FOR_ROW_COL_MV {
        if (RootCohesion->Drc < 0) // root cohesion can be used to avoid surface erosion base don land use
            CohesionSoil->Drc = -1;

        if (CohesionSoil->Drc >= 0)
            CohesionSoil->Drc = COHCalibration*(Cohesion->Drc + Cover->Drc*RootCohesion->Drc);

        // soil cohesion everywhere, plantcohesion only where plants
        if (SwitchGrassStrip)
            CohesionSoil->Drc = CohesionSoil->Drc  *(1-GrassFraction->Drc) + GrassFraction->Drc * CohGrass->Drc;

        if (CohesionSoil->Drc == 0)
            Y->Drc = 1.0;
        else {
            if (SwitchEfficiencyDET == 1)
                Y->Drc = std::min(1.0, 1.0/(0.89+0.56*fabs(CohesionSoil->Drc)));
            else
                if (SwitchEfficiencyDET == 2)
                    Y->Drc = std::min(1.0, 0.79*exp(-0.85*fabs(CohesionSoil->Drc)));
                else
                    if (SwitchEfficiencyDET == 3)
                        Y->Drc = std::min(1.0, 1.0/(2.0*fabs(CohesionSoil->Drc)));
        }
        if (CohesionSoil->Drc < 0)
            Y->Drc = 0; // to force max strength

        // empirical analysis based on Limburg data, dating 1989
        // aggr stab is Lowe test median drops to halve an aggregate
        if (SwitchSplashEQ == 1) {
            if (AggrStab->Drc > 0)
                SplashStrength->Drc = 5.331*pow(std::max(ASCalibration*AggrStab->Drc, 1.0),-0.238);
                //SplashStrength->Drc = 2.82/std::max(ASCalibration*AggrStab->Drc, 1.0);
                //splashb = 2.96;
                //redone as y = 5.3361x^-0.238  excell
        }

        // Eurosem method, aggr stab is not strength but sed delivery so the opposite
        if (SwitchSplashEQ == 2) {                      
           SplashStrength->Drc = (1/ASCalibration)*AggrStab->Drc;
        }
        if (AggrStab->Drc < 0 || RootCohesion->Drc < 0)
                SplashStrength->Drc = -1;
        // negative values give no splash
    }
}

//---------------------------------------------------------------------------
/// called after get input data, initializes non-input maps and variables
void TWorld::IntializeData(void)
{
    //TO DO add units and descriptions --> TMmapVariables.h

    //totals for mass balance
    MB = 0;
    MBs = 0;
    nrCells = 0;
    FOR_ROW_COL_MV {
        nrCells+=1;
    }

    DX = NewMap(0);
    CellArea = NewMap(0);
    FOR_ROW_COL_MV
    {
        DX->Drc = _dx/cos(asin(Grad->Drc));
        CellArea->Drc = DX->Drc * _dx;
    }
    CatchmentArea = mapTotal(*CellArea);

    SoilWidthDX = NewMap(0);

    // surface storage
    MDS = NewMap(0);
    FOR_ROW_COL_MV {
        double RRmm = 10 * RR->Drc;
        MDS->Drc = std::max(0.0, 0.243*RRmm + 0.010*RRmm*RRmm - 0.012*RRmm*tan(asin(Grad->Drc))*100);
        MDS->Drc /= 1000; // convert to m
    }

    //combination display
    COMBO_SS = NewMap(0);
    COMBO_BL = NewMap(0);
    COMBO_TC = NewMap(0);
    COMBO_V = NewMap(0);

    SoilETMBcorrection = 0;
    //### infiltration maps
    InfilTot = 0;
    InfilTotmm = 0;
    InfilKWTot = 0;
    IntercTot = 0;
    IntercETaTot = 0;
    IntercTotmm = 0;
    IntercETaTotmm = 0;
    ETaTot = 0;
    ETaTotmm = 0;
    ETaTotVol = 0;
    GWlevel = 0;
    theta1tot = 0;
    theta2tot = 0;
    thetai1tot = 0;
    thetai2tot = 0;
    thetai1cur = 0;
    thetai2cur = 0;

    BaseFlowTot = 0;
    BaseFlowInit = 0;
    SoilMoistTot = 0;
    SoilMoistDiff = 0;

    //houses
    IntercHouseTot = 0;
    IntercHouseTotmm = 0;
    IntercLitterTot = 0;
    IntercLitterTotmm = 0;
    WaterVolTot = 0;
    WaterVolSoilTileTot = 0;
    WaterVolTotmm = 0;
    WaterVolRunoffmm = 0;
    StormDrainTotmm = 0;
    ChannelVolTot = 0;
    QSideVolTot = 0;
    StormDrainVolTot = 0;
    floodVolTotmm= 0;
    floodVolTot = 0;
    //floodVolTotInit = 0;
    floodVolTotMax = 0;
    floodAreaMax = 0;
    floodBoundaryTot = 0;
    floodBoundarySedTot = 0;

    // infiltration
    InfilVolFlood = NewMap(0);
    //InfilVolKinWave = NewMap(0);
    InfilVol = NewMap(0);
    InfilmmCum = NewMap(0);
    InfilVolCum = NewMap(0);
    fact = NewMap(0); // used in SWATRE
    Ksateff = NewMap(0);
    Poreeff = NewMap(0);
    Thetaeff = NewMap(0);
    FFull = NewMap(0);
    Perc = NewMap(0);
    PercmmCum = NewMap(0);
    Fcum = NewMap(0);
    Lw = NewMap(0);
    Lwmm = NewMap(0);

    //### runoff maps
    Qtot = 0;
    Qtot_dt = 0;
    QTiletot = 0;
    QfloodoutTot = 0;
    Qfloodout = 0;
    Qtotmm = 0;
    Qboundtotmm = 0;
    GWdeeptot = 0;
    Qpeak = 0;
    QpeakTime = 0;

    WH = NewMap(0);
    WHbef = NewMap(0);
    WHrunoff = NewMap(0);
    WHmax = NewMap(0);
    WHstore = NewMap(0);
    MicroStoreVol = NewMap(0);
    FlowWidth = NewMap(0);
    V = NewMap(0);
    VH = NewMap(0);
    Alpha = NewMap(0);
    Q = NewMap(0);
    Qn = NewMap(0);

    if (SwitchDischargeUser) {
        DischargeUserPoints = ReadMap(LDD,getvaluename("qinpoints"));
        QuserIn = NewMap(0);

        FOR_ROW_COL_MV_L {
            if (DischargeUserPoints->Drc > 0 && ChannelWidth->Drc == 0) {
                //message
                int p = (int) DischargeUserPoints->Drc ;
                ErrorString = QString("Discharge input point %1 is not in a channel!").arg(p);
                DEBUG(ErrorString);
                throw 1;
            }
        }}

    }

    if (SwitchWaveUser) {
        WHboundarea = ReadMap(LDD,getvaluename("whbound"));
        FOR_ROW_COL_MV_L {
            if (WHboundarea->Drc != 0.0)
                WHboundarea->Drc = 1.0;
        }}
        WHbound = NewMap(0);
        WHboundRain = NewMap(0);
        report(*WHboundarea,"b.map"); // test delete later?
    }

    //flowmask = NewMap(0);
    K2DOutlets = NewMap(0);
    //K2DQ = NewMap(0);

    if(SwitchPesticide)
    {
        K2DQP = NewMap(0);
        K2DQPX = NewMap(0);
        K2DQPY = NewMap(0);
        K2DP = NewMap(0);
        K2DPC = NewMap(0);
        K2DPCN = NewMap(0);
    }

    QinKW = NewMap(0);
    Qoutput = NewMap(0);
    Qm3total = NewMap(0);
    Qm3max = NewMap(0);
    FHI = NewMap(0);
    Qsoutput = NewMap(0);

    WaterVolin = NewMap(0);
    WaterVolall = NewMap(0);

    WHinitVolTot = 0;

    if (SwitchKinematic2D != K2D_METHOD_DYN)
        SwitchFloodInitial = false;

    if (SwitchFloodInitial) {
        hmxInit = ReadMap(LDD, getvaluename("whinit"));
        report(*hmxInit,"wh_init.map");
    } else {
        hmxInit = NewMap(0);
    }

    // needs to be done here because profile uses data like impermable fration, tiledrain etc
    if (InfilMethod == INFIL_SWATRE) {
        thetaTop = NewMap(0);

        // VJ 110420 added tiledrain depth for all profiles, is all used in infiltration
        SwatreSoilModel = InitSwatre(ProfileID);
        if (SwatreSoilModel == nullptr)
            throw 3;

        if (SwitchInfilCrust)
        {
            SwatreSoilModelCrust = InitSwatre(ProfileIDCrust);
            if (SwatreSoilModelCrust == nullptr)
                throw 3;
        }
        if (SwitchInfilCompact)
        {
            SwatreSoilModelCompact = InitSwatre(ProfileIDCompact);
            if (SwatreSoilModelCompact == nullptr)
                throw 3;
        }
        if (SwitchGrassStrip)
        {
            SwatreSoilModelGrass = InitSwatre(ProfileIDGrass);
            if (SwatreSoilModelGrass == nullptr)
                throw 3;
        }
        initSwatreStructure = true;
        // flag: structure is created and can be destroyed in function destroydata
    }

    // SwitchUseMaterialDepth not active!
    SwitchUseMaterialDepth = false;
    // if(SwitchErosion && SwitchUseMaterialDepth)
    // {
    //     Storage = ReadMap(LDD, getvaluename("detmat"));
    //     StorageDep = NewMap(0.0);
    //     SedimentMixingDepth = ReadMap(LDD, getvaluename("sedmixdepth"));
    //     FOR_ROW_COL_MV
    //     {
    //         if(Storage->Drc != -1)
    //         {
    //             Storage->Drc = Storage->Drc * ChannelAdj->Drc * DX->Drc;
    //         }else
    //         {
    //             Storage->Drc = -999999;
    //         }
    //         SedimentMixingDepth->Drc  = std::max(0.01, SedimentMixingDepth->Drc);
    //     }
    // }

    // if(SwitchIncludeChannel) {
    //     if(SwitchErosion && SwitchUseMaterialDepth)
    //     {
    //         RStorageDep = NewMap(0.0);
    //         RSedimentMixingDepth = ReadMap(LDD, getvaluename("chansedmixdepth"));
    //         RStorage = ReadMap(LDD, getvaluename("chandetmat"));
    //         FOR_ROW_COL_MV
    //         {
    //             if(RStorage->Drc != -1)
    //             {
    //                 RStorage->Drc = RStorage->Drc * ChannelWidth->Drc * DX->Drc;
    //             }else
    //             {
    //                 RStorage->Drc = -999999;
    //             }
    //             RSedimentMixingDepth->Drc = std::max(RSedimentMixingDepth->Drc, 0.01);
    //         }
    //     }
    // }

    if (SwitchChannelBaseflowStationary)
        FindStationaryBaseFlow();

}
//---------------------------------------------------------------------------
void TWorld::InitPesticide(void)
{

    if (SwitchPesticide)
    {
        //### pesticides maps
        PestMassApplied = 0.0;
        PestLossTotOutlet = 0.0;
        PestFluxTotOutlet = 0.0;
        PestRunoffSpatial = 0.0;
        PestDisMixing = 0.0;
        PestSorMixing = 0.0;
        PestInfilt = 0.0;
        PestStorage = 0.0;
        MaxVup = 0.0;
        PestRunoffSpatialex = 0.0;
        PestDisMixingex = 0.0;
        PestSorMixingex = 0.0;
        PestInfiltex = 0.0;
        PestLossTotOutletex = 0.0;
        Maxsolubility=530e-3; // max solubility kg/m3 metolachlor
        Pestdetach = 0.0;
        PestCinfilt=0.0;
        PestCfilmexit=0.0;

        KD=NewMap(0);
        kr=NewMap(0);
        rhob=NewMap(0);
        pestiinf=NewMap(0);
        pestiinfold=NewMap(0);
        poro=NewMap(0);
        PCA=NewMap(0);
        epsil=NewMap(0);
        Kfilm=NewMap(0);
        K1=NewMap(0);
        AX=NewMap(0);

        C=NewMap(0);
        C_Kn=NewMap(0);
        CS=NewMap(0);
        CM=NewMap(0);
        Qp=NewMap(0);
        Qpn=NewMap(0);
        Pest=NewMap(0);
        PCinfilt=NewMap(0);
        PCfilmexit=NewMap(0);

        C_N=NewMap(0);
        CM_N=NewMap(0);
        CS_N=NewMap(0);

        C_K=NewMap(0);
        C_Kold=NewMap(0);
        CM_K=NewMap(0);
        CS_K=NewMap(0);

        Fkold=NewMap(0);
        Fk=NewMap(0);
        Fmk=NewMap(0);
        flagpest=NewMap(0);

        PMassApplied=NewMap(0);
        PRunoffSpatial=NewMap(0);
        PDisMixing=NewMap(0);
        PSorMixing=NewMap(0);
        PInfilt=NewMap(0);
        PStorage=NewMap(0);

        PRunoffSpatialex=NewMap(0);
        PDisMixingex=NewMap(0);
        PSorMixingex=NewMap(0);
        PInfiltex=NewMap(0);

        Pdetach=NewMap(0);
    }
    if (SwitchPesticide)
    {
        N_SPK=1;

        //test Joyce papier
        //PCA=NewMap(0.000180); //kg/m²
        //epsil=NewMap(0.25E-2); //m
        //KD=NewMap(0.00941);//m3/kg
        //kr=NewMap(0.000833333); // /s
        //poro=NewMap(0.47);
        //Kfilm=NewMap(1.16667E-5); // m/s

        // test 5-22
        PCA=NewMap(0.0000174); //kg/m²
        epsil=NewMap(0.001); //m
        KD=NewMap(0.00617);//m3/kg
        //KD=NewMap(0.0);//m3/kg
        kr=NewMap(0.0012); // /s
        poro=NewMap(0.37);
        Kfilm=NewMap(1.16667E-5); // m/s

        // qDebug()<< "initial " ;

        FOR_ROW_COL_MV
        {
            PMassApplied->Drc = PCA->Drc*_dx*_dx*1000*1000*1000; //*SnowmeltZone->Drc; //µg for partial appli //DX
            rhob->Drc=2.65E3*(1.0-poro->Drc);// soil bulk density g/m3 rhob=NewMap(1404.5); // kg/m3
            C_N->Drc= 0.0; // initialisation for t=0 kg/m3
            // partial application
            //            CM_N->Drc= (PCA->Drc*SnowmeltZone->Drc)/(epsil->Drc*poro->Drc + rhob->Drc*epsil->Drc*KD->Drc); // initialisation for t=0 kg/kg

            //VJ             CM_N->Drc= PCA->Drc*poro->Drc/epsil->Drc + (1-poro->Drc)*PCA->Drc/epsil->Drc*KD->Drc*rhob->Drc; // initialisation for t=0 kg/kg
            CM_N->Drc= (PCA->Drc)/(epsil->Drc*poro->Drc + rhob->Drc*epsil->Drc*KD->Drc); // initialisation for t=0 kg/kg
            CS_N->Drc = CM_N->Drc*KD->Drc; // ! initialisation for t=0 kg/m3
            //     qDebug()<< "initial C:"<< C->Drc << "cm"<< CM->Drc << "CS"<< CS->Drc;

            // no sorption
            // CS_N->Drc=0.0;
            // CM_N->Drc=(PCA->Drc)/(epsil->Drc*poro->Drc);

            PDisMixing->Drc = CM_N->Drc*epsil->Drc*poro->Drc*_dx*_dx*1000*1000*1000; //µg
            PSorMixing->Drc = CS_N->Drc*epsil->Drc*rhob->Drc*_dx*_dx*1000*1000*1000; //µg
        }

        PestMassApplied = mapTotal(*PMassApplied);
        PestDisMixing = mapTotal(*PDisMixing);
        PestSorMixing = mapTotal(*PSorMixing);

        if(Switchheaderpest)
        {
            Switchheaderpest=false;
            QFile fout("massbalancenew.txt");
            fout.open(QIODevice::WriteOnly | QIODevice::Text);
            QTextStream out(&fout);
            out.setRealNumberPrecision(3);
            out.setFieldWidth(0);
            out.setRealNumberNotation(QTextStream::FixedNotation);
            out << "time" << " " << "PestMassApplied" << " " << "PestDisMixing" << " " << "PestSorMixing" << " " << "PestLossTotOutlet" << " " << "PestRunoffSpatial"
                << " " << "PestInfilt" << " " << "MBp" << " "
                << "RainTot" << " " << "WaterVolSoilTot" << " " << "IntercTot" << " " << "InfilTot" << " " << "Qtot*1000*1000" << " "
                << "flux1" << " " << "flux2" << " "<< "flux3" << " "<< "flux4" << " "<< "flux5" << " "<< "flux6" <<" "<< "pestiinf*pow(10.0,9)"<<" "<<"CM*pow(10.0,6)"<<" "
                << "CS*pow(10.0,6"<<" "<< "fact*1000"<< " "<< "InfilVol*1000*1000"<<" "<<"Qn*pow(10.0,6)" << " "<< "PDisMixing" << " "<< "poro"
                << " "<< "epsil"<< " "<< "DX" << " "<< "switchrunoff" << " "<< "K1"<< " "<< "Q*pow(10.0,6)"<< " "<< "C*pow(10.0,10)"<< " "<< "iterconv"
                << " "<< "WHoutavg" << " "<< "WHoutavgold"<< " " << "MBpex" << " " << "InfilVol"<< " " << "InfilVolold";
            out << "\n";

            out << "EI" << " " << PestMassApplied << " " << PestDisMixing << " " << PestSorMixing << " " << "PestLossTotOutlet" << " " << "PestRunoffSpatial"
                << " " << "PestInfilt" << " " << PestMassApplied-PestDisMixing-PestSorMixing << " "
                << "RainTot" << " " << "WaterVolSoilTot" << " " << "IntercTot" << " " << "InfilTot" << " " << "Qtot*1000*1000" << " "
                << "flux1" << " " << "flux2" << " "<< "flux3" << " "<< "flux4" << " "<< "flux5" << " "<< "flux6" <<" "<< "pestiinf"<< " "<<"CM"<<" "
                << "CS"<<" "<< "fact"<< " "<< "InfilVol"<<" "<<"Qn" << " "<< "PDisMixing" << " "<< "poro"
                << " "<< "epsil"<< " "<< "DX" << " "<< "switchrunoff" << " "<< "K1"<< " "<< "Q*pow(10.0,6)"<< " "<< "C*pow(10.0,10)" << " "<< "iterconv"
                << " "<< "WHoutavg" << " "<< "WHoutavgold" << " " << "MBpex"<< " " << "InfilVol"<< " " << "InfilVolold"<< " " << "Vup" << " " << "Vup_old" << " "<< "Cold";
            out << "\n";
        }
    }
}
//---------------------------------------------------------------------------
//TODO: are all switches and options initialised here?
//TODO: add calibration factors here to set to 1.0
void TWorld::IntializeOptions(void)
{
    nrRainfallseries = 0;
    nrSnowmeltseries = 0;

    //dirs and names
    resultDir.clear();
    inputDir.clear();
    resultFileName = QString("totals.csv");
    outflowFileName = QString("outlets.csv");
    totalSeriesFileName = QString("totalseries.csv");
    totalLandunitFileName = QString("totlandunit.csv");
    floodStatsFileName = QString("floodstats.csv");

    totalErosionFileName = QString("erosion.map");
    totalDepositionFileName = QString("deposition.map");
    totalChanErosionFileName = QString("chandet.map");
    totalChanDepositionFileName = QString("chandep.map");
    totalSoillossFileName = QString("soilloss.map");

    rainfallMapFileName = QString("rainfall.map");
    interceptionMapFileName = QString("interception.map");
    infiltrationMapFileName = QString("infiltration.map");
    runoffMapFileName = QString("Flowcumm3.map");
    channelDischargeMapFileName = QString("chandism3.map");
    floodMaxQFileName = QString("chanmaxq.map");
    floodMaxChanWHFileName = QString("chanmaxwh.map");

    floodTimeFileName = QString("floodtime.map");
    floodFEWFileName = QString("floodstart.map");
    floodMaxVFileName = QString("Vmax.map");
    floodMaxVHFileName = QString("VHmax.map");
    floodWHmaxFileName= QString("WHmax.map");
    tileWaterVolfilename= QString("drainvol.map");
    //tileQmaxfilename= QString("drainqmax.map");

    rainFileName.clear();
    rainFileDir.clear();
    rainSatFileName.clear();
    rainSatFileDir.clear();
    ETFileName.clear();
    ETFileDir.clear();
    ETSatFileName.clear();
    ETSatFileDir.clear();
    dischargeinFileDir.clear();
    dischargeinFileName.clear();
//    snowmeltFileName.clear();
//    snowmeltFileDir.clear();
//    snowmeltSatFileName.clear();
//    snowmeltSatFileDir.clear();
    SwatreTableDir.clear();
    SwatreTableName.clear();
    resultFileName.clear();
    outflowFileName.clear();
    totalSeriesFileName.clear();
    dischargeinFileName.clear();
    dischargeinFileDir.clear();
    WaveinFileName.clear();
    WaveinFileDir.clear();

    SwitchUserCores = false;

    SwitchImage = false;
    SwitchResultDatetime = false;
    SwitchOutputTimestamp = false;
    SwitchVariableTimestep = false;
    SwitchWriteCommaDelimited = true;
    SwitchWritePCRtimeplot = false;
    //SwitchOutputTimeUser = false;
    SwitchSeparateOutput = false;
    SwitchWriteHeaders = true; // write headers in output files in first timestep
    SwitchEndRun = false;

    SwitchAdvancedOptions = false;
    SwitchPsiUser = false;
    SwitchRainfall = true;
    SwitchSnowmelt = false;
    SwitchRoadsystem = false;
    SwitchHouses = false;
    SwitchRaindrum = false;

    SwitchLitter = false;

    SwitchLinkedList = false;
    SwitchChannelKinWave = true;
    SwitchTimeavgV = true;
    SwitchChannelKinwaveDt = false;
    SwitchChannelMaxV = true;
    Switch2DDiagonalFlow = true;
    SwitchSWOFopen = true;
    SwitchMUSCL = false;
    SwitchFloodInitial = false;
    SwitchFlowBarriers = false;
    SwitchBuffers = false;
    SwitchHeun = false;
    SwitchErosion = false;
    SwitchUse2Phase = false;
    SwitchUseGrainSizeDistribution = false;
    SwitchReadGrainSizeDistribution = false;
    SwitchSedtrap = false;
    SwitchEfficiencyDET = 1;
    SwitchEfficiencyDETCH = 2;
    SwitchKETimebased = false;
    SwitchIncludeDiffusion = false;
    SwitchIncludeRiverDiffusion = false;
    SwitchUseMaterialDepth = false;

    SwitchIncludeChannel = false;
    //SwitchChannelBaseflow = false;
    SwitchGWflow = false;
    SwitchGW2Dflow =  false;
    SwitchGWSWOFflow =  false;
    SwitchLDDGWflow = false;
    SwitchSWATGWflow = false;
    SwitchChannelBaseflowStationary = false;
    SwitchChannelInfil = false;
    SwitchCulverts = false;
    SwitchDischargeUser = false;
    SwitchWaveUser = false;
    SwitchIncludeTile = false;
    SwitchIncludeStormDrains = false;

    SwitchHardsurface = false;
    SwitchInfilCompact = false;
    SwitchInfilCrust = false;
    SwitchGrassStrip = false;
    SwitchDumphead = false;
    initSwatreStructure = false;  // check to flag when swatre 3D structure is created, needed to clean up data
    SwitchGeometric = true;
    SwitchImpermeable = false;
    SwitchTwoLayer = false;
    SwitchThreeLayer = false;

    SwitchPesticide = false;
    Switchheaderpest = true;

    addedbaseflow = false;
}
//---------------------------------------------------------------------------
void TWorld::FindStationaryBaseFlow()
{
    int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
    int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};

    BaseFlowDischarges = ReadMap(LDD, getvaluename("baseflow"));
    BaseFlowInflow = NewMap(0.0);
    BaseFlowInitialVolume = NewMap(0.0);

    FOR_ROW_COL_MV_CH
    {
        pcr::setMV(tma->Drc);
        pcr::setMV(tmb->Drc);
        tmc->Drc = 0;
        tmd->Drc = 0;
    }

    for (int  ro = 0; ro < _nrRows; ro++){
        for (int  co = 0; co < _nrCols; co++){
            if(!pcr::isMV(LDDChannel->data[ro][co]))
            {
                if(LDDChannel->data[ro][co] == 5)
                {

                    int ncells = 0;
                    double inflow = 0;
                    double baseflow = BaseFlowDischarges->data[ro][co];
                    // in m3/s
                    if (BaseFlowDischarges->data[ro][co] == 0)
                        break;

                    LDD_LINKEDLIST *list = nullptr, *temp = nullptr;
                    list = (LDD_LINKEDLIST *)malloc(sizeof(LDD_LINKEDLIST));

                    list->prev = nullptr;
                    list->rowNr = ro;
                    list->colNr = co;

                    while (list != nullptr)
                    {
                        int i = 0;
                        bool  subCachDone = true;
                        int rowNr = list->rowNr;
                        int colNr = list->colNr;

                        for (i=1; i<=9; i++)
                        {
                            int r, c;
                            int ldd = 0;

                            // this is the current cell
                            if (i==5)
                                continue;

                            r = rowNr+dy[i];
                            c = colNr+dx[i];

                            if (INSIDE(r, c) && !pcr::isMV(LDDChannel->Drc))
                                ldd = (int) LDDChannel->Drc;
                            else
                                continue;

                            // check if there are more cells upstream, if not subCatchDone remains true
                            if (pcr::isMV(tma->Drc) &&
                                    FLOWS_TO(ldd, r, c, rowNr, colNr) &&
                                    INSIDE(r, c))
                            {
                                temp = (LDD_LINKEDLIST *)malloc(sizeof(LDD_LINKEDLIST));
                                temp->prev = list;
                                list = temp;
                                list->rowNr = r;
                                list->colNr = c;
                                subCachDone = false;
                            }
                        }

                        // all cells above a cell are linked in a "sub-catchment or branch
                        // continue with water and sed calculations
                        // rowNr and colNr are the last upstream cell linked
                        if (subCachDone)
                        {
                            int r = rowNr;
                            int c = colNr;
                            tma->Drc = 0;
                            ncells ++;

                            temp=list;
                            list=list->prev;
                            free(temp);
                            // go to the previous cell in the list

                        }/* eof subcatchment done */
                    } /* eowhile list != nullptr */


                    inflow = baseflow/ ncells;

                    list = nullptr;
                    temp = nullptr;
                    list = (LDD_LINKEDLIST *)malloc(sizeof(LDD_LINKEDLIST));

                    list->prev = nullptr;
                    // start gridcell: outflow point of area
                    list->rowNr = ro;
                    list->colNr = co;

                    while (list != nullptr)
                    {
                        int i = 0;
                        bool  subCachDone = true; // are sub-catchment cells done ?
                        int rowNr = list->rowNr;
                        int colNr = list->colNr;

                        for (i=1; i<=9; i++)
                        {
                            int r, c;
                            int ldd = 0;

                            // this is the current cell
                            if (i==5)
                                continue;

                            r = rowNr+dy[i];
                            c = colNr+dx[i];

                            if (INSIDE(r, c) && !pcr::isMV(LDDChannel->Drc))
                                ldd = (int) LDDChannel->Drc;
                            else
                                continue;

                            // check if there are more cells upstream, if not subCatchDone remains true
                            if (pcr::isMV(tmb->Drc) &&
                                    FLOWS_TO(ldd, r, c, rowNr, colNr) &&
                                    INSIDE(r, c))
                            {
                                temp = (LDD_LINKEDLIST *)malloc(sizeof(LDD_LINKEDLIST));
                                temp->prev = list;
                                list = temp;
                                list->rowNr = r;
                                list->colNr = c;
                                subCachDone = false;
                            }
                        }

                        // all cells above a cell are linked in a "sub-catchment or branch
                        // continue with water and sed calculations
                        // rowNr and colNr are the last upstreM cell linked
                        if (subCachDone)
                        {
                            int r = list->rowNr;
                            int c = list->colNr;
                            tmb->Drc = 0;
                            BaseFlowInflow->Drc = inflow; // will be added in baseflow every timestep is in m3/s!

                            tmc->Drc += 1;

                            for (i=1;i<=9;i++)
                            {

                                int r, c, ldd = 0;

                                if (i==5)  // Skip current cell
                                    continue;

                                r = rowNr+dy[i];
                                c = colNr+dx[i];

                                if (INSIDE(r, c) && !pcr::isMV(LDDChannel->Drc))
                                    ldd = (int) LDDChannel->Drc;
                                else
                                    continue;

                                if (INSIDE(r, c) &&
                                        FLOWS_TO(ldd, r,c,rowNr, colNr) &&
                                        !pcr::isMV(LDDChannel->Drc) )
                                {
                                    tmc->data[list->rowNr][list->colNr] += tmc->Drc;
                                    tmd->data[list->rowNr][list->colNr] += tmd->Drc;
                                }
                            }

                            r = list->rowNr;
                            c = list->colNr;

                            double q = (tmc->Drc * inflow - tmd->Drc);

                            double h, h1;
                            h = 1;
                            // first guess new h with old alpha
                            h1 = h;
                            double A = 0;

                            // newton raphson iteration
                            if (q > 0)
                            {
                                double F, dF;
                                int count = 0;

                                do{
                                    h = h1;
                                    if (h < 1e-10)
                                        break;

                                    double P,R;
                                    double FW = ChannelWidth->Drc;
                                    P = FW + 2.0*h;
                                    A = FW*h;
                                    F = std::max(0.0, 1.0 - q/(sqrt(ChannelGrad->Drc)/ChannelN->Drc*A*pow(A/P,2.0/3.0)));
                                    dF = (5.0*FW+6.0*h)/(3.0*h*P);
                                    h1 = h - F/dF;
                                    // function divided by derivative
                                    count++;
                                }while(fabs(h1-h) > 1e-10 && count < 20);
                            }

                            if (h > ChannelDepth->data[list->rowNr][list->colNr]) {
                                h = ChannelDepth->data[list->rowNr][list->colNr];
                                A = ChannelWidth->Drc*h;
                            }
                            BaseFlowInitialVolume->data[list->rowNr][list->colNr] = A*DX->Drc;

                            temp=list;
                            list=list->prev;
                            free(temp);
                            // go to the previous cell in the list

                        }/* eof subcatchment done */
                    } /* eowhile list != nullptr */
                }
            }
        }
    }

    FOR_ROW_COL_MV_CH
    {
        tmc->Drc = 0;
        tmd->Drc = 0;
    }
    report(*BaseFlowInitialVolume,"baseflowinitm3s.map");
    report(*BaseFlowInflow,"baseinflow.map");


    BaseFlowInit = MapTotal(*BaseFlowInitialVolume);

}
//---------------------------------------------------------------------------
void TWorld::InitImages()
{
    if(SwitchImage && QFileInfo(satImageFileName).exists())

    {
        DEBUG("Reading background satellite inage");
        cTRGBMap *image = readRasterImage(satImageFileName);
        //        qDebug() << "sat image" <<  image->cellSize()  << image->nrCols() << image->nrRows();
        this->RGB_Image = image;
    }
}
//---------------------------------------------------------------------------
// read and Intiialize all Tile drain variables and maps
void TWorld::InitTiledrains(void)
{
    if (SwitchIncludeTile || SwitchIncludeStormDrains) {

        //switch tile is soil draimn and switch drain is urban drains

        // channel vars and maps that must be there even if channel is switched off
        TileVolTot = 0;
        TileWaterVol = NewMap(0);
        RunoffVolinToTile = NewMap(0);
        TileQ = NewMap(0);
        TileMaxQ = NewMap(0);
        TileQn = NewMap(0);
        Tileq = NewMap(0);
        TileAlpha = NewMap(0);
        TileMaxAlpha = NewMap(0);

        //##### Tile maps #####

        LDDTile = InitMaskTiledrain(getvaluename("lddtile"));
        // must be first LDDTile is the mask for tile drains
        FOR_ROW_COL_MV_TILE {
            if (LDDTile->Drc == 0)
                SET_MV_REAL8(&LDDTile->Drc);
        }


        nrValidCellsTile = 0;
        FOR_ROW_COL_MV_TILE {
            nrValidCellsTile++;
        }
        FOR_ROW_COL_MV_TILE {
            LDD_COOR newcr;
            newcr.r = r;
            newcr.c = c;
            crtile_ << newcr;
        }
        crlinkedlddtile_= MakeLinkedList(LDDTile);


        TileArea = NewMap(0);
        TileDiameter = NewMap(0);
        TileInlet = ReadMap(LDDTile, getvaluename("tilesink"));
        TileGrad = ReadMap(LDDTile, getvaluename("tilegrad"));
        checkMap(*TileGrad, LARGER, 1.0, "Tile drain gradient must be SINE of slope angle (not tangent)");
        calcValue(*TileGrad, 0.001, MAX);
        TileN = ReadMap(LDDTile, getvaluename("tileman"));
        cover(*TileGrad, *LDD, 0);
        cover(*TileN, *LDD, 0);
        cover(*TileInlet, *LDD, 0);
        TileWaterVolSoil = NewMap(0);
        TileWidth = ReadMap(LDDTile, getvaluename("tilewidth"));
        TileHeight = ReadMap(LDDTile, getvaluename("tileheight"));
        TileDiameter = ReadMap(LDDTile, getvaluename("tilediameter"));

        cover(*TileN, *LDD, 0);
        cover(*TileGrad, *LDD, 0);
        cover(*TileDiameter, *LDD, 0);
        cover(*TileWidth, *LDD, 0);
        cover(*TileHeight, *LDD, 0);

        FOR_ROW_COL_MV_TILE {
            TileInlet->Drc = std::min(TileInlet->Drc, 0.9*_dx*_dx);
        }

        if (SwitchIncludeTile) {
            TileDepth = ReadMap(LDDTile, getvaluename("tiledepth"));
            cover(*TileDepth, *LDD, -1); //VJ non tile cells flagged by -1 value, needed in swatre init
            TileDrainSoil = NewMap(0);
        }

        // dimensions rectangular or circular
        if (SwitchIncludeStormDrains && SwitchStormDrainCircular) {
            TileDiameter = ReadMap(LDDTile, getvaluename("tilediameter"));
            FOR_ROW_COL_MV_TILE {
                double area = (TileDiameter->Drc*0.5)*(TileDiameter->Drc*0.5)*PI;// PI r^2
                TileArea->Drc = area * 2; // two sides of the street
            }
            CalcMAXDischCircular();
        }
        if (SwitchIncludeStormDrains && !SwitchStormDrainCircular) {
            //rectangular drainage
            FOR_ROW_COL_MV_TILE {
                TileArea->Drc = 2 * TileWidth->Drc*TileHeight->Drc;
                // two sides of the street
            }
            CalcMAXDischRectangular();
        }
    }
}
//---------------------------------------------------------------------------
// Make a shaded relief map from the DEM for map display
//shade=cos(I)sin(S)cos(A-D)+sin(I)cos(S)
//barriers should be added to the DEM already

void TWorld::InitShade(void)
{
    ShadeBW = NewMap(0);
    Fill(*tma,0);
    double maxDem = -1e9;
    double minDem = 1e9;

    FOR_ROW_COL_MV
    {
        //        double Incl = 15.0/180.0*PI;
        //        double Decl = 300/180.0*PI;
        double mat[9];
        double dx, dy;//, aspect;
        double factor = 1.0;

        minDem = std::min(DEM->Drc, minDem);
        maxDem = std::max(DEM->Drc, maxDem);

        for (int i = 0; i < 9; i++)
            mat[i] = DEM->Drc;
        if (r > 0 && r < _nrRows-1 && c > 0 && c < _nrCols-1)
        {
            if(!pcr::isMV(LDD->data[r-1][c-1]))
                mat[0] = DEM->data[r-1][c-1];
            if(!pcr::isMV(LDD->data[r-1][c  ]))
                mat[1] = DEM->data[r-1][c  ];
            if(!pcr::isMV(LDD->data[r-1][c+1]))
                mat[2] = DEM->data[r-1][c+1];
            if(!pcr::isMV(LDD->data[r  ][c-1]))
                mat[3] = DEM->data[r  ][c-1];

            if(!pcr::isMV(LDD->data[r  ][c+1]))
                mat[5] = DEM->data[r  ][c+1];
            if(!pcr::isMV(LDD->data[r+1][c-1]))
                mat[6] = DEM->data[r+1][c-1];
            if(!pcr::isMV(LDD->data[r+1][c  ]))
                mat[7] = DEM->data[r+1][c  ];
            if(!pcr::isMV(LDD->data[r+1][c+1]))
                mat[8] = DEM->data[r+1][c+1];
        }
        for (int i = 0; i < 9; i++)
            mat[i] *= factor;

        dx = (mat[2] + 2*mat[5] + mat[8] - mat[0] -2*mat[3] - mat[6])/(8*_dx);
        dy = (mat[0] + 2*mat[1] + mat[2] - mat[6] -2*mat[7] - mat[8])/(8*_dx);

        //        if (dy < 0)
        //            aspect = atan(dx/dy)+2*PI;
        //        else
        //            if (dy > 0)
        //                aspect = atan(dx/dy)+PI;
        //            else
        //                aspect = 0;
        //Shade->Drc = cos(Incl)*Grad->Drc*cos(aspect-Decl) + sin(Incl)*cos(asin(Grad->Drc));


        //http://edndoc.esri.com/arcobjects/9.2/net/shared/geoprocessing/spatial_analyst_tools/how_hillshade_works.htm
        //Burrough, P. A. and McDonell, R.A., 1998. Principles of Geographical Information Systems (Oxford University Press, New York), p. 190.
        double z_factor = 2.0;
        double Slope_rad = atan( z_factor * sqrt ( dx*dx+dy*dy) );
        double Aspect_rad = 0;
        if( dx != 0)
        {
            Aspect_rad = atan2(dy, -dx);
            if (Aspect_rad < 0)
                Aspect_rad = 2*PI + Aspect_rad;
        }
        else
        {
            if(dy > 0)
                Aspect_rad = PI/2.0;
            else
                Aspect_rad = 2*PI - PI/2.0;
        }
        double Zenith_rad = 70.0 * PI / 180.0;
        double Azimuth_rad = 240 * PI / 180.0;
        tma->Drc = 255.0 * ( ( cos(Zenith_rad) * cos(Slope_rad) ) + ( sin(Zenith_rad) * sin(Slope_rad) * cos(Azimuth_rad - Aspect_rad) ) );
    }
    double MaxV = mapMaximum(*tma);
    double MinV = mapMinimum(*tma);

    FOR_ROW_COL_MV
    {
        tma->Drc = (tma->Drc-MinV)/(MaxV-MinV);
        // VJ add a bit of elevation for enhanced effect
        tma->Drc = 0.8*tma->Drc+0.2*(DEM->Drc - minDem)/(maxDem-minDem);
        //ShadeBW->Drc = Shade->Drc;
    }
    MaxV = mapMaximum(*tma);
    MinV = mapMinimum(*tma);
    FOR_ROW_COL_MV
    {
        ShadeBW->Drc = (tma->Drc-MinV)/(MaxV-MinV);
        // VJ add a bit of elevation for enhanced effect
    }

}
//---------------------------------------------------------------------------
// for drawing onscreen
void TWorld::InitScreenChanNetwork()
{
    if(!SwitchIncludeChannel)
        return;

    op.lddch_.clear();
    op.lddch_.append(crlinkedlddch_);

 //   op.CulvertX.clear();
 //   op.CulvertY.clear();
    op.EndPointX.clear();
    op.EndPointY.clear();
    op.ObsPointX.clear();
    op.ObsPointY.clear();

    FOR_ROW_COL_MV_CH {
        if (LDDChannel->Drc == 5){
            op.EndPointX << _llx + c*_dx + 0.5*_dx;
            op.EndPointY << _lly + (_nrRows-r-1)*_dx + 0.5*_dx;
        }
    }
    FOR_ROW_COL_MV_CH {
        if (PointMap->Drc > 0){
            op.ObsPointX << _llx + c*_dx + 0.5*_dx;
            op.ObsPointY << _lly + (_nrRows-r-1)*_dx + 0.5*_dx;

        }
    }
}
//---------------------------------------------------------------------------
void TWorld::InitNewSoilProfile()
{
    if(InfilMethod != INFIL_SOAP)
        return;

    nN1_ = getvalueint("SoilWB nodes 1");
    nN2_ = 0;
    if (SwitchTwoLayer)
        nN2_ = getvalueint("SoilWB nodes 2");
    nN3_ = 0;
    if (SwitchThreeLayer)
        nN3_ = getvalueint("SoilWB nodes 3");
    SoilWBdtfactor = getvaluedouble("SoilWB dt factor");
    swatreDT = std::min(SoilWBdtfactor, _dt);
    KavgType = getvalueint("Infil Kavg");
    int vg = getvalueint("Van Genuchten");
    SwitchBrooksCorey = bool(vg == 1);
    SwitchVanGenuchten = !SwitchBrooksCorey;

    nNodes = nN1_ + nN2_ + nN3_ + 1;
    qDebug() << SwitchThreeLayer << nN3_ << nNodes;

    FOR_ROW_COL_MV {
        SOIL_LIST sr;
        sr.r = r;
        sr.c = c;
        sr.ponded = false;
        sr.dts = _dt*SoilWBdtfactor;
        sr.dtsum = 0;
        sr.drain = 0;
        sr.Infact = 0;
        sr.InfPot = 0;
        sr.SD = 0;

        sr.h.clear();
        sr.hb.clear();
        sr.Ks.clear();
        sr.pore.clear();
        sr.theta.clear();
        sr.thetar.clear();
        sr.lambda.clear();
        sr.vg_n.clear();
        sr.vg_alpha.clear();
        sr.dz.clear();
        sr.z.clear();
        sr.rootz.clear();

        sr.pore.resize(nNodes);
        sr.Ks.resize(nNodes);
        sr.h.resize(nNodes);
        sr.hb.resize(nNodes);
        sr.theta.resize(nNodes);
        sr.thetar.resize(nNodes);
        sr.lambda.resize(nNodes);
        sr.vg_n.resize(nNodes);
        sr.vg_alpha.resize(nNodes);
        sr.dz.resize(nNodes);
        sr.z.resize(nNodes);
        sr.rootz.resize(nNodes);

        crSoil << sr;
    }

    FOR_ROW_COL_MV_L {
        // use replace first time, else array doesn't initialise???
        crSoil[i_].SD = SoilDepth1->Drc;
        if (SwitchTwoLayer)
            crSoil[i_].SD = SoilDepth2->Drc;
        if (SwitchThreeLayer)
            crSoil[i_].SD = SoilDepth3->Drc;

        double dz, dz2, dz3;
        dz = SoilDepth1->Drc / nN1_;
        double facta = 1;//e100;
        if (SwitchTwoLayer)
            dz2 = (SoilDepth2->Drc - SoilDepth1->Drc) / nN2_;
        if (SwitchThreeLayer)
            dz3 = (SoilDepth3->Drc - SoilDepth2->Drc) / nN3_;

        double z = 0;
        for (int j = 0; j < nN1_+1; j++) {
            if (j == 0)
                crSoil[i_].z.replace(j, 0);
            if (j >= 1)
                crSoil[i_].z.replace(j, 0.5*dz + (j-1)*dz);
            z = 0.5*dz + (j-1)*dz;
            crSoil[i_].dz.replace(j, dz);
            crSoil[i_].theta.replace(j, ThetaI1->Drc);
            crSoil[i_].pore.replace(j, ThetaS1->Drc);
            crSoil[i_].Ks.replace(j, Ksat1->Drc/3600000); // calibrated Ksat ! so do not use for lambda etc
            crSoil[i_].thetar.replace(j, ThetaR1->Drc);
            crSoil[i_].lambda.replace(j, lambda1->Drc);
            crSoil[i_].vg_alpha.replace(j, vgalpha1->Drc*facta);
            crSoil[i_].vg_n.replace(j, vgn1->Drc);
            crSoil[i_].hb.replace(j, -psi1ae->Drc);
        }
        crSoil[i_].dz[0] = dz/2;

        if (SwitchTwoLayer) {
            for (int j = nN1_+1; j < nN1_+nN2_+1; j++) {
                if (j == nN1_+1)
                    z += 0.5*dz + 0.5*dz2;
                else
                    z += dz2;

                crSoil[i_].z.replace(j, z);
                crSoil[i_].dz.replace(j, dz2);
                crSoil[i_].theta.replace(j, ThetaI2->Drc);
                crSoil[i_].pore.replace(j, ThetaS2->Drc);
                crSoil[i_].Ks.replace(j, Ksat2->Drc/3600000); // calibrated Ksat ! so do not use for lambda etc

                crSoil[i_].thetar.replace(j,  ThetaR2->Drc);
                crSoil[i_].lambda.replace(j,  lambda2->Drc);
                crSoil[i_].vg_alpha.replace(j, vgalpha2->Drc*facta);
                crSoil[i_].vg_n.replace(j, vgn2->Drc);

                crSoil[i_].hb.replace(j,  -psi2ae->Drc);
            }
        }
        if (SwitchThreeLayer) {
            for (int j = nN1_+nN2_+1; j < nN1_+nN2_+nN3_+1; j++) {
                if (j == nN2_+1)
                    z += 0.5*dz2 + 0.5*dz3;
                else
                    z += dz3;

                crSoil[i_].z.replace(j, z);

                crSoil[i_].dz.replace(j, dz3);


                crSoil[i_].theta.replace(j, ThetaI3->Drc);
                crSoil[i_].pore.replace(j, ThetaS3->Drc);
                crSoil[i_].Ks.replace(j, Ksat3->Drc/3600000); // calibrated Ksat ! so do not use for lambda etc

                crSoil[i_].thetar.replace(j,  ThetaR3->Drc);
                crSoil[i_].lambda.replace(j,  lambda3->Drc);
                crSoil[i_].vg_alpha.replace(j, vgalpha3->Drc*facta);
                crSoil[i_].vg_n.replace(j, vgn3->Drc);
                crSoil[i_].hb.replace(j,  -psi3ae->Drc);
            }
        }

        // calc h
        for (int j = 0; j < nNodes; j++) {
            double se = (crSoil[i_].theta[j] - crSoil[i_].thetar[j])/(crSoil[i_].pore[j]-crSoil[i_].thetar[j]);
            if (SwitchBrooksCorey) {
                double hh = std::pow(se, (1.0/crSoil[i_].lambda[j]));
                crSoil[i_].h.replace(j,crSoil[i_].hb[j]/hh);
            } else {
                double n = crSoil[i_].vg_n[j];
                double m = 1-1/n;
                crSoil[i_].h.replace(j, -std::pow((std::pow(1/se,1/m)-1),1/n)/crSoil[i_].vg_alpha[j]);
                // kPa to m water
            }

          //  qDebug() << j << crSoil[i_].h[j];
        }
        // for (int j = 0; j < nNodes; j++) {
        //     getHfromTheta(j,crSoil[i_]);
        // }



        double sum = 0;
        double rootmax = 0.8;
        // linear root distribution following dz, sum = 1
        for (int j = 0; j < nNodes; j++) {
            crSoil[i_].rootz[j] = (j+1)*crSoil[i_].dz[j];
            if (crSoil[i_].rootz[j] > rootmax)
                crSoil[i_].rootz[j] = 0;
            else
                crSoil[i_].rootz[j] = (rootmax - crSoil[i_].rootz[j])/rootmax;

            sum = sum + crSoil[i_].rootz[j];
        }
        for (int j = 0; j < nNodes; j++)
            crSoil[i_].rootz[j] = crSoil[i_].rootz[j]/sum;

    }}

}

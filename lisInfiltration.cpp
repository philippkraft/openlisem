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
  \file lisInfiltration.cpp
  \brief Simplified infiltraton processes: Green and Ampt, Smith and Parlanage, both 1 and 2 layer. SWATRE has separate files.

functions: \n
- void TWorld::InfilEffectiveKsat(void)
- void TWorld::InfilSwatre(cTMap *_WH)
- void TWorld::InfilMethods(cTMap * _Ksateff, cTMap *_WH, cTMap *_fpot, cTMap *_fact, cTMap *_L1, cTMap *_L2, cTMap *_FFull)
- double TWorld::IncreaseInfiltrationDepth(int r, int c, double fact, double *L1p, double *L2p, double *FFullp)
- void TWorld::Infiltration(void)
 */

#include <algorithm>
#include "lisemqt.h"
#include "global.h"
#include "model.h"
#include "operation.h"

//---------------------------------------------------------------------------
// Calculate effective Ksat based on surface structure, impermeable etc.
void TWorld::InfilEffectiveKsat(bool first)
{
    // moved to datainit!
    // if (first) {
    //     #pragma omp parallel for num_threads(userCores)
    //     FOR_ROW_COL_MV_L {
    //         Ksat1->Drc *= _dt/3600000.0; // mm/h to m
    //         if (SwitchTwoLayer)
    //             Ksat2->Drc *= _dt/3600000.0;
    //         if (SwitchThreeLayer)
    //             Ksat3->Drc *= _dt/3600000.0;
    //         if (SwitchInfilCrust)
    //             KsatCrust->Drc *= _dt/3600000.0;
    //         if (SwitchInfilCompact)
    //             KsatCompact->Drc *= _dt/3600000.0;
    //     }}
    // }

    if (SwitchInfiltration && InfilMethod != INFIL_SWATRE)
    {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            Ksateff->Drc = Ksat1->Drc;
            Poreeff->Drc = ThetaS1->Drc;

            // exponential crusting proces with cumulative rainfall
            if (SwitchInfilCrust) {
                double factor = 1.0-exp(-0.2*std::max(0.0, RainCum->Drc/1000-5.0));  //
                // exponential decline until from no crusting to full crusting at ~ 30 mm,
                //old research Jean Boiffin, multiple rainfall in a growing season, progressive crusting

                double ksatdiff = std::max(0.0,Ksat1->Drc - KsatCrust->Drc);
                Ksateff->Drc = KsatCrust->Drc + ksatdiff * factor;

                double porediff = std::max(0.0,ThetaS1->Drc - PoreCrust->Drc);
                Poreeff->Drc = PoreCrust->Drc + porediff * factor;
            }
            Thetaeff->Drc = std::max(ThetaR1->Drc,ThetaI1->Drc);

            // compacted surfaces
            if (SwitchInfilCompact) {
                Ksateff->Drc = Ksateff->Drc*(1-CompactFraction->Drc) + KsatCompact->Drc*CompactFraction->Drc;
                Poreeff->Drc = Poreeff->Drc*(1-CompactFraction->Drc) + PoreCompact->Drc*CompactFraction->Drc;
            }

            // grass strips? old concept?
            if (SwitchGrassStrip) {
                Ksateff->Drc = Ksateff->Drc*(1-GrassFraction->Drc) + KsatGrass->Drc*GrassFraction->Drc;
                Poreeff->Drc = ThetaS1->Drc*(1-GrassFraction->Drc) + PoreGrass->Drc*GrassFraction->Drc;
            }

            // density factor and OM corrections directly in LISEM (instead of dbase creator)
            // because SWATRE also needs this
            // these correction come from calculations based on Saxton and Rawls
            if (SwitchOMCorrection) {
                double OM2 = OMcorr->Drc*OMcorr->Drc;
                double corrKsOA = -0.0065*OM2 - 0.0415*OMcorr->Drc + 1.0001;
                double corrKsOB = -2.6319*OMcorr->Drc + 0.0197;
                double corrPOA =  -0.1065*OM2 - 0.0519*OMcorr->Drc + 0.9932;
                double corrPOB = 0.0532*OM2 + 0.008*OMcorr->Drc + 0.0037;
                Ksateff->Drc = corrKsOA*Ksateff->Drc + corrKsOB;
                Poreeff->Drc = corrPOA*Poreeff->Drc + corrPOB;
            }
            if (SwitchDensCorrection) {
                double corrKsDA =  -3.28*DensFact->Drc + 4.3;
                double corrKsDB = -96.65*DensFact->Drc + 98.28;
                double corrPDA = DensFact->Drc;
                double corrPDB = -1.0 * DensFact->Drc + 1.0;
                Ksateff->Drc = corrKsDA*Ksateff->Drc + corrKsDB;
                Poreeff->Drc = corrPDA*Poreeff->Drc + corrPDB;
            }

            Ksateff->Drc *= 1.0-fractionImperm->Drc;

            // to avoid pore is less than thetaR else nan in redistribution
            if (Poreeff->Drc < ThetaR1->Drc)
                ThetaR1->Drc = 0.5*Poreeff->Drc;
            Ksateff->Drc = std::max(0.0, Ksateff->Drc);

            // may be a problem in for instance redistribution            
            if (SwitchWaveUser) {
                // when incoming wave, no infil in that area
                if (WHboundarea->Drc > 0) {
                    Ksateff->Drc = 0;
                    Poreeff->Drc = 0;
                    Ksat1->Drc = 0;
                    Ksat2->Drc = 0;
                    ThetaS1->Drc = 0;
                    ThetaS2->Drc = 0;
                    ThetaI1->Drc = 0;
                    ThetaI2->Drc = 0;
                }
            }

        }}
    }

}
//---------------------------------------------------------------------------
// this function is not used!
void TWorld::Infiltration()
{
    if (!SwitchInfiltration)
        return;

    //NOTE fact and fpot have a unit of m (not m/s)
    if (InfilMethod == INFIL_SWATRE) {
        //InfilSwatre();
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            cell_InfilSwatre(i_, r, c);
        }}
    }
    else
    {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            cell_InfilMethods(r, c);
        }}
    }
}

//---------------------------------------------------------------------------
// Infiltration by Green and Ampt,Smith and Parlange
// All the same except for calculation of the potential infiltration fpot
// 1 layer and 2 layers
/*!
\brief function to calculate potential and actula infiltration rate according to
Green and Ampt, or Smith and Parlange.

This function calculates the potential infiltration according to G&A or S&P \n
then calls IncreaseInfiltrationDepth to increase the wetting front.
*/
void TWorld::cell_InfilMethods(int r, int c)
{
    // default vars are first layer vars
    double Ks = Ksateff->Drc;  //in m
    double Psi = Psi1->Drc; // in m
    double fwh = 0;
    double fpot_ = 0;
    double fact_ = 0;
    double SoilDep1 = SoilDepth1->Drc;
    double SoilDep2 = 0;

    if (Ksateff->Drc == 0)
        return;

    if (FloodDomain->Drc == 0) {
        fwh = WH->Drc; //runoff in kinwave or dyn wave
    } else {
        fwh = hmx->Drc; // flood in kin wave
    }
    // select the appropriate domain water height for overpressure

    fwh += MBm->Drc; // mass balance correction
    fwh = std::max(0.0,fwh);

    // only do infiltration on permeable soils, is now incorporated in ksateff
    //if (SoilWidthDX->Drc > 0 && fwh > 0) {
    if (fwh > 0) {
        //calculate potential infiltration rate fpot
        if (SwitchTwoLayer || SwitchThreeLayer) {
            SoilDep2 = SoilDepth2->Drc;
            // if wetting front in second layer set those vars
            if (Lw->Drc > SoilDep1 && Lw->Drc < SoilDep2) {
                //weighed harmonic mean:
                //https://corporatefinanceinstitute.com/resources/data-science/harmonic-mean/
                // sum (weights) / sum (weight/variable)
 //               Ks = Havg(Ksateff->Drc,Ksat2->Drc,SoilDep1,Lw->Drc-SoilDep1);
                Ks = Lw->Drc/(SoilDep1/Ksateff->Drc+(Lw->Drc-SoilDep1)/Ksat2->Drc);
                // if wetting front > layer 1 than ksat is determined weighted average (harmonic mean)
                Psi = Psi2->Drc; //in m
            }
        }

        if (InfilMethod == INFIL_GREENAMPT)
            fpot_ = Ks*(1.0+(Psi+fwh)/std::max(1e-3, Lw->Drc));
        else {
            // smith parlange, not really tested
            double space = Poreeff->Drc-Thetaeff->Drc;
            if (Lw->Drc > SoilDepth1->Drc)
                space = ThetaS2->Drc-ThetaI2->Drc;
            double B = (fwh + Psi)*space;
            if (B > 0.01) {
                fpot_ = Ks*exp(Fcum->Drc/B)/(exp(Fcum->Drc/B)-1);
            } else
                fpot_ = Ks;
        }

        fact_ = std::min(fpot_, fwh);
        if (fact_ < 1e-10)
            fact_ = 0;
        // actual infil in m, cannot have more infil than water on the surface, includes rainfall

        if (fact_ > 0) {
//            if (SwitchThreeLayer)
//                fact_ = IncreaseInfiltrationDepthNew3(fact_, r, c);
//            else
                if (SwitchTwoLayer)
                    fact_ = IncreaseInfiltrationDepthNew2(fact_, r, c);
                else
                    fact_ = IncreaseInfiltrationDepthNew1(fact_, r, c);
        }
        // adjust fact and increase Lw, for twolayer, impermeable etc



        if (fwh < fact_)
        {
            fact_ = fwh;
            fwh = 0;
        }
        else
            fwh -= fact_;

        if(FloodDomain->Drc == 0)
            WH->Drc = fwh;
        else
            hmx->Drc = fwh;
        // adjust the WH in the correct domain with new fact

        Fcum->Drc += fact_; // for Smith and Parlange
        // increase cumulative infil in m
       // fact->Drc = fact_;
        //InfilVol->Drc = fact_* SoilWidthDX->Drc * DX->Drc;
        InfilVol->Drc = fact_* FlowWidth->Drc * DX->Drc;
        // calc infiltrated volume for mass balance
        // use flowwidth because Ksateff included impermeable surfaces anyway
    } else {
       // fact->Drc = 0;
        InfilVol->Drc = 0;
    }
}

//---------------------------------------------------------------------------
/*!
\brief function to increase wetting front and deal with 2nd layer and impermeable subsoil
 returns actual infiltration rate.

 this function is called form all infiltration functions except Swatre:\n
 - one layer or two layers
  - returns depth of the wetting front (Lw)\n
 - returns actual infiltration in mm, NOT rate in mm/h
*/
double TWorld::IncreaseInfiltrationDepthNew1(double fact_in, int r, int c)
{
    double dtheta1 = std::max(0.0,Poreeff->Drc-Thetaeff->Drc); // space in the top layer
    double L = Lw->Drc;
    double SoilDep1 = SoilDepth1->Drc;
    double fact_out = 0;
    double space = 0;
    double Lnew = 0;

    // impermeable and L reached SD1, no more infil
    if (SwitchImpermeable && L > SoilDep1 - 0.001) {
        Lw->Drc = SoilDep1;
        return 0;
    }

    if (SwitchGWflow) {
        if (GWWH->Drc >= SoilDepth1init->Drc-HMIN) {
            return 0;
        }
    }

    Lnew = L + fact_in/std::max(dtheta1,0.01);
    // increase wetting front
    space = (SoilDep1 - L)*dtheta1;
    if(Lnew > SoilDep1 || space < fact_in) {
        if (SwitchImpermeable)
            // if impermeable remaining space is infiltration
            fact_out = space;
        else
            fact_out = Perc->Drc;
        Lnew = SoilDep1;
    } else {
        fact_out = fact_in;
    }

    Lw->Drc = std::min(SoilDep1,std::max(0.0, Lnew));
    return std::max(0.0, fact_out);
}
//---------------------------------------------------------------------------
double TWorld::IncreaseInfiltrationDepthNew2(double fact_in, int r, int c)
{
    double dtheta1 = std::max(0.0,Poreeff->Drc-Thetaeff->Drc); // space in the top layer
    double dtheta2 = std::max(0.0,ThetaS2->Drc-ThetaI2->Drc);
    double SoilDep1 = SoilDepth1->Drc;
    double SoilDep2 = SoilDepth2->Drc;
    double fact_out = 0;
    double space = 0;
    double Lnew = 0;
    double L = Lw->Drc;
    double dfact2 = 0;
    bool passing = false;
    double space2 = 0;

    // profile is full
    if (SwitchImpermeable && L > SoilDep2 - 0.001) {
        Lw->Drc = SoilDep2;
        return 0;
    }

    if (SwitchGWflow) {
       if (L >= SoilDep1 && GWWH->Drc >= SoilDepth2init->Drc-HMIN) {
           Lw->Drc = SoilDep1;
           return 0;
       }
       // when GWWH fills osildep2 than soildep2 is 0 anyway
    }


    // L is in layer 1
    if (L <= SoilDep1) {
        Lnew = L + fact_in/std::max(0.01,dtheta1);
        space = (SoilDep1-L)*dtheta1;

        if(fact_in > space || Lnew > SoilDep1) {
            // water is moving into layer 2
            passing = true;
            dfact2 = fact_in - space;
            // remaining water for layer 2
        } else {
            // all remains SD1
            fact_out = fact_in;
        }
    }

    // L is in layer 2
    if (L > SoilDep1) {
        //L already in layer 2

        Lnew = L + fact_in/std::max(0.01,dtheta2);
        space2 = (SoilDep2-L)*dtheta2;

        if (Lnew > SoilDep2 || fact_in > space2) {
            if (SwitchImpermeable)
                fact_out = space2;
            else
                fact_out = Perc->Drc;

            Lnew = SoilDep2;
            // L at bottom
        } else {
            fact_out = fact_in;
            // everything fitted
        }
    }
    // Lnew is now soildep2 or the actual depth

    // L is moving from layer 1 into 2 in this timestep
    if (passing) {
        // second layer still at initial
        space2 = (SoilDep2-SoilDep1)*dtheta2;
        Lnew = SoilDep1 + dfact2/std::max(0.01,dtheta2);
        dfact2 = std::min(dfact2, space2);

        if (dtheta2 < 0.01 || Lnew > SoilDep2) {
            if (SwitchImpermeable)
                fact_out = space+space2;
            else
                fact_out = Perc->Drc;
            Lnew = SoilDep2;
        } else
            fact_out = fact_in; // everything fitted
    }

    Lw->Drc = std::min(SoilDep2,std::max(0.0, Lnew));
    return std::max(0.0,fact_out);
}
//---------------------------------------------------------------------------
// 3 layer infiltration! not used yet
double TWorld::IncreaseInfiltrationDepthNew3(double fact_in, int r, int c)
{
    double dtheta1 = std::max(0.0,Poreeff->Drc-Thetaeff->Drc); // space in the top layer
    double dtheta2 = std::max(0.0,ThetaS2->Drc-ThetaI2->Drc);
    double dtheta3 = std::max(0.0,ThetaS3->Drc-ThetaI3->Drc);
    double SoilDep1 = SoilDepth1->Drc;
    double SoilDep2 = SoilDepth2->Drc;
    double SoilDep3 = SoilDepth3->Drc;
    double fact_out = 0;
    double space = 0;
    double Lnew = 0;
    double L = Lw->Drc;
    double dfact12 = 0;
    double dfact23 = 0;
    bool passing12 = false;
    bool passing23 = false;
    double space2 = 0;
    double space3 = 0;

    // profile is full
    if (SwitchImpermeable && L > SoilDep2 - 0.001) {
        Lw->Drc = SoilDep2;
        return 0;
    }

    // L is in layer 1
    if (L <= SoilDep1) {
        Lnew = L + fact_in/std::max(0.01,dtheta1);
        space = (SoilDep1-L)*dtheta1;

        if(fact_in > space || Lnew > SoilDep1) {
            // water is moving into layer 2
            passing12 = true;
            dfact12 = fact_in - space;
            // remaining water for layer 2
        } else {
            // all remains SD1
            fact_out = fact_in;
        }
    }

    // L is in layer 2
    if (L > SoilDep1 && L <= SoilDep2) {
        //L already in layer 2 but not in 3
        Lnew = L + fact_in/std::max(0.01,dtheta2);
        space2 = (SoilDep2-L)*dtheta2;

        if (fact_in > space2 || Lnew > SoilDep2) {
            passing23 = true;
            dfact23 = fact_in - space2;
        } else {
            // all remains SD2
            fact_out = fact_in;
        }
    }

    // L is in layer 3
    if (L > SoilDep2 && L <= SoilDep3) {
        //L already in layer 2 but not in 3
        Lnew = L + fact_in/std::max(0.01,dtheta3);
        space3 = (SoilDep3-L)*dtheta3;

        if (fact_in > space3 || Lnew > SoilDep3) {
            if (SwitchImpermeable)
                fact_out = space3;
            else
                fact_out = Perc->Drc;

            Lnew = SoilDep3;
            // L at bottom
        } else {
            // all remains SD3
            fact_out = fact_in;
        }
    }
    // Lnew is now soildep3 or the actual depth

    // L is moving from layer 1 into 2 in this timestep
    if (passing12) {
        // second layer still at initial
        space2 = (SoilDep2-SoilDep1)*dtheta2;
        Lnew = SoilDep1 + dfact12/std::max(0.01,dtheta2);
        dfact12 = std::min(dfact12, space2);

        if (dtheta2 < 0.01 || Lnew > SoilDep2) {
            passing23 = true;
            dfact23 = fact_in - space2;
            // also does not fit in SD2, passing to SD3
        } else {
            fact_out = fact_in;
            // everything fitted in SD2
        }
    }

    // L is moving from layer 2 into 3 in this timestep
    if (passing23) {
        // second layer still at initial
        space3 = (SoilDep3-SoilDep2)*dtheta3;
        Lnew = SoilDep2 + dfact23/std::max(0.01,dtheta3);
        dfact23 = std::min(dfact23, space3);

        if (dtheta3 < 0.01 || Lnew > SoilDep2) {
            if (SwitchImpermeable)
                fact_out = space+space3;
            else
                fact_out = Perc->Drc;
            Lnew = SoilDep3;
        } else
            fact_out = fact_in; // everything fitted
    }

    Lw->Drc = std::min(SoilDep3,std::max(0.0, Lnew));
    return std::max(0.0,fact_out);

}

//---------------------------------------------------------------------------
// NOT USED
void TWorld::cell_InfilSwatre(long i_, int r, int c)
{
    //profile number 0 is impeermeable so no need to do anything
    double frac = std::min(1.0, RoadWidthHSDX->Drc/_dx + HouseCover->Drc);

    if (ProfileID->Drc <= 0 || frac == 1.0) {
        fact->Drc = 0;
        InfilVol->Drc = 0;
        return;
    }

    if (FloodDomain->Drc == 0)
        tm->Drc = WH->Drc;
    else
        tm->Drc = hmx->Drc;

    WHbef->Drc = tm->Drc;
    tm->Drc *= frac;

    SwatreStep(i_, r, c, SwatreSoilModel, tm, TileDrainSoil, thetaTop);

    // actual infil is dif between WH before and after
    fact->Drc = std::max(0.0, WHbef->Drc - tm->Drc);

    if (FloodDomain->Drc == 0)
        WH->Drc = tm->Drc;
    else
        hmx->Drc = tm->Drc;

    if (CrustFraction->Drc > 0) {
        tm->Drc = WHbef->Drc;
        tma->Drc = 0;
        tmb->Drc = 0;

        SwatreStep(i_, r, c, SwatreSoilModelCrust, tm, tma, tmb);
        // calculate crust SWATRE and get the soil moisture of the top node

        double _wh;
        if (FloodDomain->Drc == 0)
            _wh = WH->Drc;
        else
            _wh = hmx->Drc;
        // water level on crusted areas

        double wha = tm->Drc*CrustFraction->Drc + _wh*(1-CrustFraction->Drc);
        // weighted average
        if (FloodDomain->Drc == 0)
            WH->Drc = wha;
        else
            hmx->Drc = wha;

        fact->Drc = (WHbef->Drc - wha);
        thetaTop->Drc = tmb->Drc*CrustFraction->Drc + thetaTop->Drc*(1-CrustFraction->Drc);
    }

    if (SwitchInfilCompact)
    {
        tm->Drc = WHbef->Drc;
        tma->Drc = 0;
        tmb->Drc = 0;

        SwatreStep(i_, r, c, SwatreSoilModelCompact, tm, tma, tmb);
        // calculate crust SWATRE and get the soil moisture of the top node

        double _wh;
        if (FloodDomain->Drc == 0)
            _wh = WH->Drc;
        else
            _wh = hmx->Drc;
        // water level on crusted areas

        double wha = tm->Drc*CompactFraction->Drc + _wh*(1-CompactFraction->Drc);
        // weighted average
        if (FloodDomain->Drc == 0)
            WH->Drc = wha;
        else
            hmx->Drc = wha;

        fact->Drc = (WHbef->Drc - wha);
        thetaTop->Drc = tmb->Drc*CompactFraction->Drc + thetaTop->Drc*(1-CompactFraction->Drc);
    }

    if (SwitchGrassStrip)
    {
        tm->Drc = WHbef->Drc;
        tma->Drc = 0;
        tmb->Drc = 0;

        SwatreStep(i_, r, c, SwatreSoilModelGrass, tm, tma, tmb);
        // calculate crust SWATRE and get the soil moisture of the top node

        double _wh;
        if (FloodDomain->Drc == 0)
            _wh = WH->Drc;
        else
            _wh = hmx->Drc;
        // water level on crusted areas

        double wha = tm->Drc*GrassFraction->Drc + _wh*(1-GrassFraction->Drc);
        // weighted average
        if (FloodDomain->Drc == 0)
            WH->Drc = wha;
        else
            hmx->Drc = wha;

        fact->Drc = (WHbef->Drc - wha);
        thetaTop->Drc = tmb->Drc*GrassFraction->Drc + thetaTop->Drc*(1-GrassFraction->Drc);
    }

    InfilVol->Drc = fact->Drc * FlowWidth->Drc * DX->Drc;

}

//---------------------------------------------------------------------------
/// SWATRE infiltration, takes WH and calculateds new WH and infiltration surplus for kin wave
/// OBSOLETE we use cell_infilswatre
void TWorld::InfilSwatre()
{
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {

        double frac = std::min(1.0, RoadWidthHSDX->Drc/_dx + HouseCover->Drc);

        if (ProfileID->Drc <= 0 || frac == 1.0) {
            fact->Drc = 0;
            InfilVol->Drc = 0;
            continue;
        }

        if (FloodDomain->Drc == 0)
            tm->Drc = WH->Drc;
        else
            tm->Drc = hmx->Drc;

        WHbef->Drc = tm->Drc;

        SwatreStep(i_, r, c, SwatreSoilModel, tm, TileDrainSoil, thetaTop);
        // tm = new water level after infiltration
        // thetatop is not used, meant for pesticides

        fact->Drc = std::max(0.0, WHbef->Drc - tm->Drc);
        // actual infil is dif between WH before and after

        if (FloodDomain->Drc == 0)
            WH->Drc = tm->Drc;
        else
            hmx->Drc = tm->Drc;

        //TODO test infil swatre for crusts and compaction
        if (CrustFraction->Drc > 0) {
            tm->Drc = WHbef->Drc;
            tma->Drc = 0;
            tmb->Drc = 0;

            SwatreStep(i_, r, c, SwatreSoilModelCrust, tm, tma, tmb);
            // calculate crust SWATRE and get the soil moisture of the top node

            double _wh;
            if (FloodDomain->Drc == 0)
                _wh = WH->Drc;
            else
                _wh = hmx->Drc;
            // new water level from regular swatre

            double wha = tm->Drc*CrustFraction->Drc + _wh*(1-CrustFraction->Drc);
            // weighed average
            if (FloodDomain->Drc == 0)
                WH->Drc = wha;
            else
                hmx->Drc = wha;

            fact->Drc = (WHbef->Drc - wha);
            thetaTop->Drc = tmb->Drc*CrustFraction->Drc + thetaTop->Drc*(1-CrustFraction->Drc);
        }

        if (SwitchInfilCompact) {
            tm->Drc = WHbef->Drc;
            tma->Drc = 0;
            tmb->Drc = 0;

            SwatreStep(i_, r, c, SwatreSoilModelCompact, tm, tma, tmb);
            // calculate crust SWATRE and get the soil moisture of the top node

            double _wh;
            if (FloodDomain->Drc == 0)
                _wh = WH->Drc;
            else
                _wh = hmx->Drc;
            // water level on crusted areas

            double wha = tm->Drc*CompactFraction->Drc + _wh*(1-CompactFraction->Drc);
            // weighted average
            if (FloodDomain->Drc == 0)
                WH->Drc = wha;
            else
                hmx->Drc = wha;

            fact->Drc = (WHbef->Drc - wha);
            thetaTop->Drc = tmb->Drc*CompactFraction->Drc + thetaTop->Drc*(1-CompactFraction->Drc);
        }

        if (SwitchGrassStrip) {
            tm->Drc = WHbef->Drc;
            tma->Drc = 0;
            tmb->Drc = 0;

            SwatreStep(i_, r, c, SwatreSoilModelGrass, tm, tma, tmb);
            // calculate crust SWATRE and get the soil moisture of the top node

            double _wh;
            if (FloodDomain->Drc == 0)
                _wh = WH->Drc;
            else
                _wh = hmx->Drc;
            // water level on crusted areas

            double wha = tm->Drc*GrassFraction->Drc + _wh*(1-GrassFraction->Drc);
            // weighted average
            if (FloodDomain->Drc == 0)
                WH->Drc = wha;
            else
                hmx->Drc = wha;

            fact->Drc = (WHbef->Drc - wha);
            thetaTop->Drc = tmb->Drc*GrassFraction->Drc + thetaTop->Drc*(1-GrassFraction->Drc);
        }

        InfilVol->Drc = fact->Drc * FlowWidth->Drc * DX->Drc;
        // calc infilvolume from fact

    }}

}

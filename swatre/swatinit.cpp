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
\file swatinit.cpp
\brief SWATRE: initialize soil profile with inithead maps data and clean up after run

functions:
- SOIL_MODEL * TWorld::InitSwatre(cTMap *profileMap); \n
- void TWorld::CloseSwatre(SOIL_MODEL *s); \n
*/

#include "lerror.h"
#include "model.h"

//--------------------------------------------------------------------------------
// make the 3D structure PIXEL_INFO, based on profile numbers in map
// needs zone info which needs to be done before in readswatreinput
// read optional Hinit maps
SOIL_MODEL *TWorld::InitSwatre(cTMap *profileMap)
{
    SOIL_MODEL *s = (SOIL_MODEL *)malloc(sizeof(SOIL_MODEL));

    // TODO check if this needs freeing when error;
    // why not new?

    s->minDt = swatreDT;
    s->pixel = new PIXEL_INFO[(long)nrCells];

    // set initial values
    for (long i = 0; i < (long)nrCells; i++) {
        s->pixel[i].profile = nullptr;
        s->pixel[i].dumpHid = 0;  //set to 1 for output of a pixel
        s->pixel[i].tiledrain = 0;
        s->pixel[i].wh = 0;
        s->pixel[i].percolation = 0;
        s->pixel[i].tilenode = -1;      // set tiledrain to 0, and tiledepth to -1 (above surface)
        // first node ksat adujtedd with this for partial impermeability
        s->pixel[i].impfrac = 0;
        s->pixel[i].corrKsOA = 1.0;
        s->pixel[i].corrKsOB = 0.0;
        s->pixel[i].corrKsDA = 1.0;
        s->pixel[i].corrKsDB = 0.0;
        s->pixel[i].corrPOA = 1.0;
        s->pixel[i].corrPOB = 0.0;
        s->pixel[i].corrPDA = 1.0;
        s->pixel[i].corrPDB = 0.0;
    }

    // give each pixel a profile
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        int profnr = swatreProfileNr.indexOf((int)profileMap->Drc);

        if (profnr > 0)
            s->pixel[i_].profile = profileList[profnr];  // pointer to profile
        // profile = <= 0 now set to impermeable
        s->pixel[i_].impfrac = fractionImperm->Drc;

        if (SwitchOMCorrection) {
            // these correction come from calculations based on Saxton and rawls
            s->pixel[i_].corrKsOA = -0.0065*(OMcorr->Drc*OMcorr->Drc) - 0.0415*OMcorr->Drc + 1.0001;
            s->pixel[i_].corrKsOB = -2.6319*OMcorr->Drc + 0.0197;
            s->pixel[i_].corrPOA =  -0.1065*(OMcorr->Drc*OMcorr->Drc) - 0.0519*OMcorr->Drc + 0.9932;
            s->pixel[i_].corrPOB = 0.0532*(OMcorr->Drc*OMcorr->Drc) + 0.008*OMcorr->Drc + 0.0037;
        }
        if (SwitchDensCorrection) {
            //A	-3.28	4.30
            //B	-96.65	98.28
            s->pixel[i_].corrKsDA =  -3.28*DensFact->Drc + 4.3;
            s->pixel[i_].corrKsDB = -96.65*DensFact->Drc + 98.28;
            s->pixel[i_].corrPDA =  DensFact->Drc;
            s->pixel[i_].corrPDB = -1.0 * DensFact->Drc + 1.0;
        }

        // TODO: does not work yet!
        if(SwitchDumphead) {
            s->pixel[i_].dumpHid = SwatreOutput->Drc;
        } else
            s->pixel[i_].dumpHid = 0;
    }}


    // fill the inithead structure of each pixel and set tiledrain depth if any
    double hi = HinitValue*psiCalibration;

    for (int k = 0; k < zone->nrNodes; k++) {

        if (!SwitchHinit4all) {
            QString fname = QString("%1.%2").arg(initheadName).arg(k+1, 3, 10, QLatin1Char('0'));
            // make inithead.001 to .00n name
            inith = ReadMap(LDD,fname);
        }

        // get inithead information
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            if (!SwitchHinit4all)
                hi = inith->Drc*psiCalibration;
            s->pixel[i_].h.append(hi);

            // find depth of tilenode
            if (SwitchIncludeTile) {
                if (!pcr::isMV(TileDepth->Drc) && TileDepth->Drc > 0) {
                    // NOTE depth is in m while node info is in cm, so *100
                    // endComp is the depth at the bottom of the compartment, so the tile is <= endcomp
                    if (s->pixel[i_].profile->zone->endComp[k] > TileDepth->Drc*100)
                        s->pixel[i_].tilenode = k-1;
                }
            }
        }}
    }

   //qDebug() << "DONE InitSwatre";
    return(s);
}
//--------------------------------------------------------------------------------
/// soil model instance to be freed
void TWorld::CloseSwatre(SOIL_MODEL *s)
{
    if (s == nullptr)
        return;

    swatreProfileDef.clear();
    swatreProfileNr.clear();

    delete[] s->pixel;

    free(s);
    s = nullptr;
    //qDebug() << "closed swatre";
}
//--------------------------------------------------------------------------------

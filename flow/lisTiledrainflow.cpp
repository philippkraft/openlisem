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
  \file lisTiledrainflow.cpp
  \brief calculate tile drain system flow as a kinematic wave, no sediment functions

functions: \n
- void TWorld::ToTiledrain(void) \n
- void TWorld::CalcVelDischTile() \n
- void TWorld::TileFlow(void)\n
 */

#include <algorithm>
#include "model.h"
#include "operation.h"

//TODO convert flow to linked list

//---------------------------------------------------------------------------
// flow in all road cells to tiledrain
//fraction of water and sediment flowing from the surface to the tiledrain system
void TWorld::ToTiledrainAll()
{

    if (SwitchIncludeStormDrains)  //SwitchIncludeTile ||
    {
        Fill(*RunoffVolinToTile,0);

        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_TILEL {
            //double fractiontotile = std::max(1.0, 2*0.09/CHAdjDX->Drc) * RoadWidthDX->Drc/_dx;
            //or:
            double fractiontotile  = std::min(1.0, _dt*V->Drc/(0.09/CHAdjDX->Drc)) * RoadWidthDX->Drc/_dx;;
            // 30x30cm = 0.09 m2, 2 is both sides of the street
            //Street inlet is assumed to be a hole in the street
            // fraction based on surface, simpel! or as velocity

            double MaxVol;
            if (SwitchStormDrainCircular)
                MaxVol = DX->Drc*PI*TileDiameter->Drc*TileDiameter->Drc*0.25; //(pi r^2)
            else
                MaxVol = DX->Drc*TileWidth->Drc*TileHeight->Drc;

            if (TileWaterVol->Drc >= MaxVol)
                fractiontotile = 0;
            else {
                double vol = fractiontotile*WaterVolall->Drc;
                double dh = vol/CHAdjDX->Drc;
                RunoffVolinToTile->Drc = vol;
                // adjust water height
                WHrunoff->Drc -= dh;
                WH->Drc -= dh;
                WaterVolall->Drc -= vol;
            }
        }}
    }
}
//---------------------------------------------------------------------------
//fraction of water and sediment flowing from the surface to the tiledrain system
void TWorld::ToTiledrain()
{
    if (SwitchIncludeStormDrains)  //SwitchIncludeTile ||
    {
        Fill(*RunoffVolinToTile,0);

        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_TILEL {
            if(TileInlet->Drc > 0) {// && WHrunoff->Drc > 1e-6) {
                double fractiontotile = 2 * std::max(1.0, std::min(0.0,TileInlet->Drc/CHAdjDX->Drc)) * RoadWidthDX->Drc/_dx;
                // 2 is both sides of the street
                // fraction based on surface, simpel!
                //Street inlet is assumed to be a hole in the street

                double MaxVol;
                if (SwitchStormDrainCircular)
                    MaxVol = DX->Drc*PI*TileDiameter->Drc*TileDiameter->Drc*0.25; //(pi r^2)
                else
                    MaxVol = DX->Drc*TileWidth->Drc*TileHeight->Drc;

                if (TileWaterVol->Drc >= MaxVol)
                    fractiontotile = 0;
                else {
                    double vol = fractiontotile*WaterVolall->Drc;//std::max(0.0,(WaterVolall->Drc-MicroStoreVol->Drc));
                    //  vol = std::min(vol, MaxVol - TileWaterVol->Drc);
                    double dh = vol/CHAdjDX->Drc;
                    RunoffVolinToTile->Drc = vol;
                    // adjust water height
                    WHrunoff->Drc -= dh;
                    //WHroad->Drc -= dh;
                    WH->Drc -= dh;
                    WaterVolall->Drc -= vol;
                }
            }
        }}
    }
}
//---------------------------------------------------------------------------
// V, alpha and Q in the Tile
void TWorld::CalcVelDischRectangular()
{
    double Perim, Area, TileV_;
    const double _23 = 2.0/3.0;

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_TILEL {
        double gradN = sqrt(TileGrad->Drc)/TileN->Drc;
        Area = TileWaterVol->Drc/DX->Drc;
        Perim = TileWidth->Drc + Area/TileWidth->Drc; //(=w+2*h)

        TileV_ = powl(Area/Perim,_23) * gradN;
        TileQ->Drc = Area*TileV_;
        TileAlpha->Drc  = Area/std::pow(TileQ->Drc, 0.6);
    }}
}
//---------------------------------------------------------------------------
// called from dataini, needed in kin wave
void TWorld::CalcMAXDischRectangular()
{    
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_TILEL {
        double Area = TileArea->Drc;
        double width = Area/TileHeight->Drc;
        double Perim = width+ 2*TileHeight->Drc; // factor 2 width for two drains in a strteet

        double TileV_ = powl(Area/Perim,2.0/3.0) * sqrt(TileGrad->Drc)/TileN->Drc;
        TileMaxQ->Drc = Area*TileV_;
        TileMaxAlpha->Drc  = Area/std::pow(TileMaxQ->Drc, 0.6);
    }}
}
//---------------------------------------------------------------------------
// V, alpha and Q in the Tile
// https://www.engineersedge.com/fluid_flow/partially_full_pipe_flow_calculation/partiallyfullpipeflow_calculation.htm
// https://www.ajdesigner.com/phphydraulicradius/hydraulic_radius_equation_pipe.php
// Neweton iteration to derive drain water height
void TWorld::CalcVelDischCircular()
{
   #pragma omp parallel for num_threads(userCores)
   FOR_ROW_COL_MV_TILEL {

      double gradN = sqrt(TileGrad->Drc)/TileN->Drc;
      double rr = TileDiameter->Drc/2;
      double Perim, K, theta;

      double Area = TileWaterVol->Drc / DX->Drc;

      theta = PI;
      double fx, Fx;
      double Ar = 2*Area/(rr*rr);
      for (int k = 0 ; k < 20; k++) {
            fx = 1-cos(theta);
            Fx = -sin(theta) + theta - Ar;
            theta = fx > 0 ? theta - Fx/fx : 0.0;
            if( Fx < 1e-6)
                break;
      }
      // newton rapson iteration to get the water height in a circular pipe for the wet perimeter

      K = rr*rr*(theta-sin(theta))*0.5;
      double TWH = rr - cos(theta/2.0)*rr;
      Perim = Area < 0.5*rr*rr*PI ? rr*theta : PI*TileDiameter->Drc-rr*theta;

      double minV = std::pow(TWH,2.0/3.0) * gradN;
      double TileV_ = Perim > 1e-10 ? std::min(minV, std::pow(Area/Perim,2.0/3.0) * gradN) : 0.0;
      //TileV_ = std::min(TileV_, 2.0);
      //limit velocity to 2 m/s?
      TileQ->Drc = Area*TileV_;
      TileAlpha->Drc  = Area/std::pow(TileQ->Drc, 0.6);
      //TileAlpha->Drc = std::pow(std::pow(Perim, 2.0/3.0)/gradN , 0.6);
   }}
}//---------------------------------------------------------------------------
// called from dataini, needed in kin wave
void TWorld::CalcMAXDischCircular()
{
   #pragma omp parallel for num_threads(userCores)
   FOR_ROW_COL_MV_TILEL {

      double Area = TileArea->Drc;
      double Perim = PI*TileDiameter->Drc;
      double TileV_ = std::pow(0.95*Area/Perim,2.0/3.0) * sqrt(TileGrad->Drc)/TileN->Drc;

      TileMaxQ->Drc = Area*TileV_;
      TileMaxAlpha->Drc  = Area/std::pow(TileMaxQ->Drc, 0.6);

   }}
}
//---------------------------------------------------------------------------
//- calc Tileflow, Tileheight, kin wave
void TWorld::TileFlow(void)
{
   if (!SwitchIncludeTile && !SwitchIncludeStormDrains)
      return;

   // get water from surface
   if (SwitchIncludeStormDrains) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_TILEL {
            TileWaterVol->Drc += RunoffVolinToTile->Drc;
                // add water from the surface
        }}
   }

    // get water from soil
   if (SwitchIncludeTile) {
       #pragma omp parallel for num_threads(userCores)
       FOR_ROW_COL_MV_TILEL {
          TileWaterVol->Drc += TileDrainSoil->Drc * TileDiameter->Drc * DX->Drc;
          // asume water can come from all sides!
          // add inflow to Tile in m3, tiledrainsoil is in m per timestep

          TileWaterVolSoil->Drc += TileDrainSoil->Drc * TileDiameter->Drc  * DX->Drc;
          // soil only used for MB correction

       }}
   }

   if (SwitchStormDrainCircular)
      CalcVelDischCircular();
   else
      CalcVelDischRectangular();

   TileQn->setAllMV();

   Fill(*QinKW, 0.0);
   // flag all new flux as missing value, needed in kin wave and replaced by new flux
   FOR_ROW_COL_MV_TILE {
      if (LDDTile->Drc == 5)
            Kinematic(r,c, LDDTile, TileQ, TileQn, TileAlpha, DX, TileMaxQ, TileMaxAlpha);
   }
//   KinematicExplicit(crlinkedlddtile_, TileQ, TileQn, TileAlpha, DX, TileMaxQ, TileMaxAlpha);

   cover(*TileQn, *LDD, 0); // avoid missing values around Tile for adding to Qn for output

   #pragma omp parallel for num_threads(userCores)
   FOR_ROW_COL_MV_TILEL {
        TileQn->Drc = std::min(QinKW->Drc+TileWaterVol->Drc/_dt, TileQn->Drc);
        TileWaterVol->Drc = TileWaterVol->Drc + _dt*(QinKW->Drc - TileQn->Drc);
        TileWaterVol->Drc = std::max(0.0, TileWaterVol->Drc);

// prob gives mass balance problems
        TileWaterVol->Drc = std::min(TileWaterVol->Drc, TileArea->Drc * DX->Drc);
        TileQ->Drc = TileQn->Drc;
   }}

}

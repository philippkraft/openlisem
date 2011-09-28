
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
 \file lisDatainit.cpp
 \brief All data handling functions, read, parse, initialize, make and delete the map structures etc.

 functions: \n
 - TMMap TWorld::NewMap(double value) Make a map filled with a value and LDD as mask \n
 - TMMap TWorld::ReadMap(cTMap *Mask, QString name) Make a map and fill it with data from a file\n
 - void TWorld::DestroyData(void) Destroy the list with all maps and pointers, no need to do this one by one.\n
 - TMMap TWorld::InitMask(QString name) Make and read the reference map which is the LDD \n
 - TMMap TWorld::InitMaskChannel(QString name) Make and read the channel map (channel LDD)\n
 - TMMap TWorld::InitMaskTiledrain(QString name) Make and read the tiell drain map (tile LDD)\n
 - void TWorld::InitTiledrains(void) read and intiialize all Tile drain variables and maps \n
 - void TWorld::InitBuffers(void)  read and intiialize all buffer and sediment fence variables and maps \n
 - void TWorld::InitChannel(void)  read and intiialize all channel variables and maps \n
 - void TWorld::GetInputData(void) Read all input maps \n
 - void TWorld::IntializeData(void) Initialize all auxillary maps and variables\n
 - void TWorld::IntializeOptions(void) Initialize all options (Switch...) before the run file.\n
*/

#include <qstring.h>

#include "model.h"


//---------------------------------------------------------------------------
/** \n void TWorld::InitMapList(void)
*  blabla
*/
void TWorld::InitMapList(void)
{

   maplistnr = 0;
   for (int i = 0; i < NUMNAMES; i++)
   {
      maplistTMMap[i].m = NULL;
   }
}
//---------------------------------------------------------------------------
TMMap *TWorld::NewMap(double value)
{
   TMMap *_M = new TMMap();

   _M->MakeMap(LDD, value);
   // changed to LDD instead of Mask

   if (_M)
   {
      maplistTMMap[maplistnr].m = _M;
      maplistnr++;
   }

   return(_M);
}
//---------------------------------------------------------------------------
TMMap *TWorld::ReadMap(cTMap *Mask, QString name)
{

   TMMap *_M = new TMMap();

   _M->PathName = /*inputdir + */name;

   bool res = _M->LoadFromFile();
   if (!res)
   {
      ErrorString = "Cannot find map " +_M->PathName;
      throw 1;
   }

   for (int r = 0; r < _nrRows; r++)
      for (int c = 0; c < _nrCols; c++)
         if (!IS_MV_REAL8(&Mask->Drc) && IS_MV_REAL8(&_M->Drc))
         {
            QString sr, sc;
            sr.setNum(r); sc.setNum(c);
            ErrorString = "Missing value at row="+sr+" and col="+sc+" in map: "+name;
            throw 1;
         }

   if (_M)
   {
      maplistTMMap[maplistnr].m = _M;
      maplistnr++;
   }

   return(_M);

}
//---------------------------------------------------------------------------
void TWorld::DestroyData(void)
{
   for (int i = 0; i < maplistnr; i++)
   {
      if (maplistTMMap[i].m != NULL)
      {
         maplistTMMap[i].m->KillMap();
         maplistTMMap[i].m = NULL;
      }
   }
   if (nrrainfallseries > 1)
   {
      for (int r=0; r < nrrainfallseries; r++)
         delete[] RainfallSeries[r];
      delete[] RainfallSeries;
   }

   if (InfilMethod == INFIL_SWATRE && initSwatreStructure)
   {
      FreeSwatreInfo();
      if (SwatreSoilModel)
         CloseSwatre(SwatreSoilModel);
      if (SwatreSoilModelCrust)
         CloseSwatre(SwatreSoilModelCrust);
      if (SwatreSoilModelCompact)
         CloseSwatre(SwatreSoilModelCompact);
      if (SwatreSoilModelGrass)
         CloseSwatre(SwatreSoilModelGrass);
   }
}
//---------------------------------------------------------------------------
TMMap *TWorld::InitMask(QString name)
{
   // read map and make a mask map

   TMMap *_M = new TMMap();

   _M->PathName = /*inputdir + */name;
   bool res = _M->LoadFromFile();
   if (!res)
   {
      ErrorString = "Cannot find map " +_M->PathName;
      throw 1;
   }

   if (_M)
   {
      maplistTMMap[maplistnr].m = _M;
      maplistnr++;
   }

   _dx = _M->MH.cellSizeX*1.0000000;
   _nrRows = _M->nrRows;
   _nrCols = _M->nrCols;

   return(_M);

}
//---------------------------------------------------------------------------
TMMap *TWorld::InitMaskChannel(QString name)
{

   TMMap *_M = new TMMap();

   _M->PathName = /*inputdir + */name;
   bool res = _M->LoadFromFile();
   if (!res)
   {
      ErrorString = "Cannot find map " +_M->PathName;
      throw 1;
   }

   if (_M)
   {
      maplistTMMap[maplistnr].m = _M;
      maplistnr++;
   }

   return(_M);

}
//---------------------------------------------------------------------------
TMMap *TWorld::InitMaskTiledrain(QString name)
{

   TMMap *_M = new TMMap();

   _M->PathName = /*inputdir + */name;
   bool res = _M->LoadFromFile();
   if (!res)
   {
      ErrorString = "Cannot find map " +_M->PathName;
      throw 1;
   }

   if (_M)
   {
      maplistTMMap[maplistnr].m = _M;
      maplistnr++;
   }

   return(_M);

}
//---------------------------------------------------------------------------
// read and Intiialize all Tile drain variables and maps
void TWorld::InitTiledrains(void)
{
   // channel vars and maps that must be there even if channel is switched off
   TileVolTot = 0;
   TileWaterVol = NewMap(0);
   TileWaterVolSoil = NewMap(0);
   //TileQoutflow = NewMap(0);
   RunoffVolinToTile = NewMap(0);
   TileQ = NewMap(0);
   TileQn = NewMap(0);
   TileQs = NewMap(0);
   TileQsn = NewMap(0);
   TileWH = NewMap(0);
   Tileq = NewMap(0);
   TileAlpha = NewMap(0);
   TileDrainSoil = NewMap(0);
   TileV = NewMap(0);
   TileDX = NewMap(0);

   if (!SwitchIncludeTile)
      TileDepth = NewMap(-1);
   // must exist for swatre

   // maybe needed later for erosion in tiledrain
   //TileSedTot = 0;
   //TileDepTot = 0;
   //TileDetTot = 0;
   //TileQsoutflow = NewMap(0);
   //TileDetFlow = NewMap(0);
   //TileDep = NewMap(0);
   //TileSed = NewMap(0);
   //TileConc = NewMap(0);
   //TileTC = NewMap(0);
   //TileY = NewMap(0);
   //SedToTile = NewMap(0);

   if (SwitchIncludeTile)
   {
      //## Tile maps
      LDDTile = InitMaskTiledrain(getvaluename("lddtile"));
      // must be first" LDDTile is the mask for tile drains


      TileSinkhole = ReadMap(LDDTile, getvaluename("tilesink"));
      TileWidth = ReadMap(LDDTile, getvaluename("tilewidth"));
      TileHeight = ReadMap(LDDTile, getvaluename("tileheight"));
      TileDepth = ReadMap(LDDTile, getvaluename("tiledepth"));
      TileGrad = ReadMap(LDDTile, getvaluename("tilegrad"));
      TileGrad->calcV(0.001, MAX);
      TileN = ReadMap(LDDTile, getvaluename("tileman"));
      //TileCohesion = ReadMap(LDDTile, getvaluename("chancoh"));

      TileGrad->cover(LDD, 0);
      TileWidth->cover(LDD, 0);
      TileHeight->cover(LDD, 0);
      TileDepth->cover(LDD, -1); //VJ non tile cells flaaged by -1 value, needed in swatre init
      TileN->cover(LDD, 0);
      TileSinkhole->cover(LDD, 0);

      /* TODO ? */

      FOR_ROW_COL_MV_TILE
      {
         TileDX->Drc = _dx/cos(asin(TileGrad->Drc));
         TileSinkhole->Drc = min(TileSinkhole->Drc, 0.9*_dx*_dx); //? why?
         //TileY->Drc = min(1.0, 1.0/(0.89+0.56*TileCohesion->Drc));
      }

      if (useSorted)
         lddlisttile = makeSortedNetwork(LDDTile, &lddlisttilenr);
      //VJ 110123 sorted tiledrain network
   }

}
//---------------------------------------------------------------------------
void TWorld::InitBuffers(void)
{

   BufferVolin = 0;
   BufferVolTot = 0;
   BufferVolTotInit = 0;
   BufferSedTot = 0;

   //## buffers read maps
   if (SwitchBuffers || SwitchSedtrap)
   {
      BufferID = ReadMap(LDD,getvaluename("bufferID"));
      BufferVol = ReadMap(LDD,getvaluename("bufferVolume"));
      BulkDens = getvaluedouble("Sediment bulk density");
      // also sed trap use bufffervol to calculate the max sed store

      FOR_ROW_COL_MV
      {
         if (SwitchBuffers && BufferID->Drc > 0)
         {
            Grad->Drc = 0.001;
            RR->Drc = 0.01;
            N->Drc = 0.25;
            // note ksateff in filtration is also set to 0

            // very arbitrary!!!
            /* TODO link tis to interface */
            Cover->Drc = 0;
            if (SwitchIncludeChannel && ChannelGrad->Drc > 0)
            {
               ChannelGrad->Drc = 0.001;
               ChannelN->Drc = 0.25;
               if (SwitchChannelInfil)
                  ChannelKsat->Drc = 0;
               //no infil in buffer
            }
         }
         // adjust soil and surface properties for buffercells, not sed traps
      }
   }
   else
   {
      BufferID = NewMap(0);
      BufferVol = NewMap(0);
   }

   //VJ 100514 buffer and sedtrap maps
   /* TODO how calculate max sed store with only sed traps? */
   // use slope of cell:        | /
   //                           |/
   // then max store is _dx/cos = DX*height fence * bulk dens?
   if (SwitchBuffers)
   {
      BufferVolInit = NewMap(0);
      ChannelBufferVolInit = NewMap(0);

      if (SwitchIncludeChannel)
      {
         ChannelBufferVol = NewMap(0);
         FOR_ROW_COL_MV_CH
               if (BufferID->Drc > 0)
         {
            ChannelBufferVol->Drc = BufferVol->Drc;
            BufferVol->Drc = 0;
         }
         //split buffers in channel buffers and slope buffers
         // in "ToCHannel" all flow in a buffer is dumped in the channel
         ChannelBufferVolInit->copy(ChannelBufferVol);
         // copy initial max volume of buffers in channels
      }
      BufferVolInit->copy(BufferVol);
      // copy initial max volume of remaining buffers on slopes
      BufferVolTotInit = BufferVolInit->MapTotal() + ChannelBufferVolInit->MapTotal();
      // sum up total initial volume available in buffers
      //BufferVol->fill(0);
      //ChannelBufferVol->fill(0);
      // rset to zero to fill up

   }

   BufferSed = NewMap(0);
   ChannelBufferSed = NewMap(0);

// no initial sediment assumed, volume must reflect empty status

//   if (SwitchBuffers || SwitchSedtrap)
//   {
//      BufferSedInit = NewMap(0);
//      ChannelBufferSedInit = NewMap(0);

//      BufferSed->calc2V(BufferVol, BulkDens, MUL);
//      //NOTE: buffer sed vol is maximum store in kg and will decrease while it
//      // fills up. It is assumed that the sedimented part contains a pore volume
//      // that can contain water, like  loose soil. Thsi is determined by the bulkdensity
//      if (SwitchIncludeChannel)
//      {
//         ChannelBufferSed = NewMap(0);
//         FOR_ROW_COL_MV_CH
//               if (BufferID->Drc > 0)
//         {
//            ChannelBufferSed->Drc = BufferSed->Drc;
//            BufferSed->Drc = 0;
//         }
//         //split buffers in channel buffers and slope buffers
//         // in "ToCHannel" all flow in a buffer is dumped in the channel
//         ChannelBufferSedInit->copy(ChannelBufferSed);
//         // copy initial max volume of buffers in channels
//      }
//      BufferSedInit->copy(BufferSed);
//      // copy initial max volume of remaining buffers on slopes
//      BufferSedTotInit = BufferSedInit->MapTotal() + ChannelBufferSedInit->MapTotal();
//      // sum up total initial volume available in buffers
      //BufferSed->fill(0);
      //ChannelBufferSed->fill(0);
      // rset to zero to fill up
//   }

}
//---------------------------------------------------------------------------
// read and Intiialize all channel variables and maps
void TWorld::InitChannel(void)
{
   // channel vars and maps that must be there even if channel is switched off
   ChannelVolTot = 0;
   ChannelSedTot = 0;
   ChannelDepTot = 0;
   ChannelDetTot = 0;
   SedToChannel = NewMap(0);
   ChannelWidthUpDX = NewMap(0);
   ChannelWaterVol = NewMap(0);
   //ChannelQoutflow = NewMap(0);
   RunoffVolinToChannel = NewMap(0);
   //ChannelQsoutflow = NewMap(0);
   ChannelQ = NewMap(0);
   ChannelQn = NewMap(0);
   ChannelQs = NewMap(0);
   ChannelQsn = NewMap(0);
   ChannelV = NewMap(0);
   ChannelWH = NewMap(0);
   Channelq = NewMap(0);
   ChannelAlpha = NewMap(0);
   ChannelPerimeter = NewMap(0); //VJ 110109 added for channel infil
   ChannelDX = NewMap(0);
   ChannelDetFlow = NewMap(0);
   ChannelDep = NewMap(0);
   ChannelSed = NewMap(0);
   ChannelConc = NewMap(0);
   ChannelTC = NewMap(0);
   ChannelY = NewMap(0);

   if (SwitchIncludeChannel)
   {
      //## channel maps
      LDDChannel = InitMaskChannel(getvaluename("lddchan"));
      // must be first" LDDChannel is the mask for channels

      ChannelWidth = ReadMap(LDDChannel, getvaluename("chanwidth"));
      ChannelSide = ReadMap(LDDChannel, getvaluename("chanside"));
      ChannelGrad = ReadMap(LDDChannel, getvaluename("changrad"));
      ChannelGrad->calcV(0.001, MAX);
      ChannelN = ReadMap(LDDChannel, getvaluename("chanman"));
      ChannelCohesion = ReadMap(LDDChannel, getvaluename("chancoh"));
      ChannelGrad->cover(LDD, 0);
      ChannelSide->cover(LDD, 0);
      ChannelWidth->cover(LDD, 0);
      ChannelN->cover(LDD, 0);

      ChannelN->calcV(ChnCalibration, MUL);
      if (SwitchChannelInfil)
      {
         ChannelKsat = ReadMap(LDDChannel, getvaluename("chanksat"));
         ChannelKsat->cover(LDD, 0);
         ChannelKsat->calcV(ChKsatCalibration, MUL);
         ChannelStore = NewMap(0.050); // 10 cm deep * 0.5 porosity
      }
      ChannelWidthUpDX->copy(ChannelWidth);
      ChannelWidthUpDX->cover(LDD, 0);
      FOR_ROW_COL_MV_CH
      {
         ChannelDX->Drc = _dx/cos(asin(ChannelGrad->Drc)); //VJ BUG fix, wrong calc here, div by gradient instead of cos
         ChannelY->Drc = min(1.0, 1.0/(0.89+0.56*ChannelCohesion->Drc));
      }

      if (useSorted)
         lddlistch = makeSortedNetwork(LDDChannel, &lddlistchnr);
      //VJ 110123 sorted channel network
   }
}
//---------------------------------------------------------------------------
void TWorld::InitMulticlass(void)
{
   if (!SwitchMulticlass)
      return;

   SubsMaps[0].m = ReadMap(LDD, getvaluename("fractionmu0"));
   SubsMaps[1].m = ReadMap(LDD, getvaluename("fractionmu1"));
   SubsMaps[2].m = ReadMap(LDD, getvaluename("fractionmu2"));
   SubsMaps[3].m = ReadMap(LDD, getvaluename("fractionmu3"));
   SubsMaps[4].m = ReadMap(LDD, getvaluename("fractionmu4"));
   SubsMaps[5].m = ReadMap(LDD, getvaluename("fractionmu5"));
   nrSubsMaps = 6;

}
//---------------------------------------------------------------------------
void TWorld::GetInputData(void)
{
   //## calibration factors
   ksatCalibration = getvaluedouble("Ksat calibration");
   nCalibration = getvaluedouble("N calibration");
   thetaCalibration = getvaluedouble("Theta calibration");
   psiCalibration = getvaluedouble("Psi calibration");
   ChnCalibration = getvaluedouble("Channel N calibration");
   ChKsatCalibration = getvaluedouble("Channel Ksat calibration");
   SplashDelivery = getvaluedouble("Splash Delivery Ratio");
   StemflowFraction = getvaluedouble("Stemflow fraction");
   CanopyOpeness = getvaluedouble("Canopy Openess");
   //VJ 110829 water repellency
   waterRep_a = getvaluedouble("Water Repellency A");
   waterRep_b = getvaluedouble("Water Repellency B");

   //## catchment data
   LDD = InitMask(getvaluename("ldd"));
   // THIS SHOULD BE THE FIRST MAP
   // LDD is also mask and reference file, everthing has to fit LDD
   // channels use channel LDD as mask

   tm = NewMap(0); // temp map for aux calculations
   tma = NewMap(0); // temp map for aux calculations
   tmb = NewMap(0); // temp map for aux calculations
   tmc = NewMap(0); // temp map for aux calculations
   for (int i = 0; i < 32; i++)
      SubsMaps[i].m = NULL;  // initialize substance structures

   Grad = ReadMap(LDD,getvaluename("grad"));  // must be SINE of the slope angle !!!
   Outlet = ReadMap(LDD,getvaluename("outlet"));
   Outlet->cover(LDD, 0);
   // fill outlet with zero, some users have MV where no outlet
   FOR_ROW_COL_MV
   {
      if (Outlet->Drc == 1)
      {
         if (LDD->Drc != 5)
         {
            ErrorString = "Main outlet gridcell does not coincide with pit in LDD";
            throw 1;
         }
         else
         {
            c_outlet = c;
            r_outlet = r;
            r_plot = r_outlet;
            c_plot = c_outlet;
         }
      }
   }

   PointMap = ReadMap(LDD,getvaluename("outpoint"));
   //map with points for output data


   if (SwitchRainfall)
   {
      RainZone = ReadMap(LDD,getvaluename("id"));
   }
   Snowcover = NewMap(0);
   if (SwitchSnowmelt)
   {
      SnowmeltZone = ReadMap(LDD,getvaluename("SnowID"));
      FOR_ROW_COL_MV
      {
         Snowcover->Drc = (SnowmeltZone->Drc == 0 ? 0 : 1.0);
      }
   }

   //## landuse and surface data
   N = ReadMap(LDD,getvaluename("manning"));
   N->calcV(nCalibration, MUL); //VJ 110112 moved

   RR = ReadMap(LDD,getvaluename("RR"));
   PlantHeight = ReadMap(LDD,getvaluename("CH"));
   LAI = ReadMap(LDD,getvaluename("lai"));
   Cover = ReadMap(LDD,getvaluename("cover"));
   LandUnit = ReadMap(LDD,getvaluename("landunit"));  //VJ 110107 added

   if (SwitchGrassStrip)
   {
      GrassWidthDX = ReadMap(LDD,getvaluename("grasswidth"));
      GrassFraction->copy(GrassWidthDX);
      GrassFraction->calcV(_dx, DIV);
      StripN = getvaluedouble("Grassstrip Mannings n");
      FOR_ROW_COL_MV
      {
         if (GrassWidthDX > 0)
         {
            N->Drc = N->Drc*(1-GrassFraction->Drc)+StripN*GrassFraction->Drc;
         }
      }
   }
   else
      GrassFraction = NewMap(0);
   StoneFraction  = ReadMap(LDD,getvaluename("stonefrc"));
   // WheelWidth  = ReadMap(LDD,getvaluename("wheelwidth"));
   RoadWidthDX  = ReadMap(LDD,getvaluename("road"));
   HardSurface = ReadMap(LDD,getvaluename("hardsurf"));

   //## infiltration data
   if(InfilMethod != INFIL_NONE && InfilMethod != INFIL_SWATRE)
   {
      Ksat1 = ReadMap(LDD,getvaluename("ksat1"));
      SoilDepth1 = ReadMap(LDD,getvaluename("soildep1"));
      SoilDepth1->calcV(1000, DIV);
      //VJ 101213 fixed bug: convert from mm to m
      // can be zero for outcrops
      FOR_ROW_COL_MV
      {
         if (SoilDepth1->Drc < 0)
         {
            ErrorString = QString("SoilDepth1 values < 0 at row %1, col %2").arg(r).arg(c);
            throw 1;
         }
      }

      ThetaS1 = ReadMap(LDD,getvaluename("thetas1"));
      ThetaI1 = ReadMap(LDD,getvaluename("thetai1"));

      ThetaI1->calcV(thetaCalibration, MUL); //VJ 110712 calibration of theta
      ThetaI1->calc(ThetaS1, MIN); //VJ 110712 cannot be more than porosity

      //VJ 101221 all infil maps are needed except psi
      if(InfilMethod != INFIL_KSAT)
      {
         Psi1 = ReadMap(LDD,getvaluename("psi1"));
         Psi1->calcV(psiCalibration, MUL); //VJ 110712 calibration of psi
         Psi1->calcV(0.01, MUL);
      }

      if (SwitchTwoLayer)
      {
         ThetaS2 = ReadMap(LDD,getvaluename("thetaS2"));
         ThetaI2 = ReadMap(LDD,getvaluename("thetaI2"));

         ThetaI2->calcV(thetaCalibration, MUL); //VJ 110712 calibration of theta
         ThetaI2->calc(ThetaS2, MIN); //VJ 110712 cannot be more than porosity

         //VJ 101221 all infil maps are needed except psi
         if(InfilMethod != INFIL_KSAT)
         {
            Psi2 = ReadMap(LDD,getvaluename("psi2"));
            Psi2->calcV(psiCalibration, MUL); //VJ 110712 calibration of psi
            Psi2->calcV(0.01, MUL);
         }

         Ksat2 = ReadMap(LDD,getvaluename("ksat2"));
         SoilDepth2 = ReadMap(LDD,getvaluename("soilDep2"));
         SoilDepth2->calcV(1000, DIV);
         //VJ 101213 fixed bug: convert from mm to m

         FOR_ROW_COL_MV
         {
            if (SoilDepth2->Drc < 0)
            {
               ErrorString = QString("SoilDepth2 values < 0 at row %1, col %2").arg(r).arg(c);
               throw 1;
            }
         }
      }
      if (SwitchInfilCrust)
      {
         CrustFraction = ReadMap(LDD,getvaluename("crustfrc"));
         KsatCrust = ReadMap(LDD,getvaluename("ksatcrst"));
      }
      else
         CrustFraction = NewMap(0);
      if (SwitchInfilCompact)
      {
         CompactFraction = ReadMap(LDD,getvaluename("compfrc"));
         KsatCompact = ReadMap(LDD,getvaluename("ksatcomp"));
      }
      else
         CompactFraction = NewMap(0);
   }

   if (InfilMethod == INFIL_SWATRE)
   {
      // read all Swatre profiles
      ProfileID = ReadMap(LDD,getvaluename("profmap"));

      if (SwitchGrassStrip)
         ProfileIDGrass = ReadMap(LDD,getvaluename("profgrass"));

      if (SwitchInfilCrust || SwitchWaterRepellency)
      {
         CrustFraction = ReadMap(LDD,getvaluename("crustfrc"));
         ProfileIDCrust = ReadMap(LDD,getvaluename("profcrst"));
      }
      else
         CrustFraction = NewMap(0);

      if (SwitchInfilCompact)
      {
         CompactFraction = ReadMap(LDD,getvaluename("compfrc"));
         ProfileIDCompact = ReadMap(LDD,getvaluename("profcomp"));
      }
      else
         CompactFraction = NewMap(0);

      // read the swatre tables and make the information structure ZONE etc
      int res = ReadSwatreInput(SwatreTableName, SwatreTableDir);
      if (res)
         throw res;


   }

   //## erosion maps
   if (SwitchErosion)
   {
      Cohesion = ReadMap(LDD,getvaluename("coh"));
      RootCohesion = ReadMap(LDD,getvaluename("cohadd"));
      AggrStab = ReadMap(LDD,getvaluename("AggrStab"));
      D50 = ReadMap(LDD,getvaluename("D50"));
   }

   //## read and initialize all channel maps and variables
   InitChannel();

   //## read and initialize all buffer maps and variables
   InitBuffers();

   //## read and initialize all tile drain system maps and variables
   InitTiledrains();

   InitMulticlass();

   // not more than 1.0
//   CrustFraction->calcV(1.0, MAX);
//   CompactFraction->calcV(1.0, MAX);

}
//---------------------------------------------------------------------------
void TWorld::IntializeData(void)
{

   //TO DO add units and descriptions --> TMmapVariables.h

   //totals for mass balance
   MB = 0;
   MBs = 0;
   nrCells = 0;
   FOR_ROW_COL_MV
   {
      nrCells+=1;
   }

   //### terrain maps
   DX = NewMap(0);
   CellArea = NewMap(0);
   FOR_ROW_COL_MV
   {
      Grad->Drc = max(0.0001, Grad->Drc);
      DX->Drc = _dx/cos(asin(Grad->Drc));
      CellArea->Drc = DX->Drc * _dx;
   }
   CatchmentArea = CellArea->MapTotal();

   WheelWidth = NewMap(0);
   WheelWidthDX = NewMap(0);
   SoilWidthDX = NewMap(0);
   GullyWidthDX = NewMap(0);
   MDS = NewMap(0);
   FOR_ROW_COL_MV
   {
      double RRmm = 10* RR->Drc;
      MDS->Drc = max(0, 0.243*RRmm + 0.010*RRmm*RRmm - 0.012*RRmm*tan(asin(Grad->Drc))*100);
      MDS->Drc /= 1000; // convert to m
   }

   //### rainfall and interception maps
   RainTot = 0;
   RainTotmm = 0;
   Rainpeak = 0;
   RainpeakTime = 0;
   RainAvgmm = 0;
   SnowAvgmm = 0;
   SnowTot = 0;
   SnowTotmm = 0;
   Snowpeak = 0;
   SnowpeakTime = 0;
   Rain = NewMap(0);
   Rainc = NewMap(0);
   RainCum = NewMap(0);
   RainNet = NewMap(0);
   LeafDrain = NewMap(0);
   //not used RainIntensity = NewMap(0);
   //not used RainM3 = NewMap(0);
   CStor = NewMap(0);
   Interc = NewMap(0);

   Snowmelt = NewMap(0);
   Snowmeltc = NewMap(0);
   SnowmeltCum = NewMap(0);

   InterceptionLAIType = getvalueint("Canopy storage equation");
   if (InterceptionLAIType == 8)
   {
      SwitchInterceptionLAI = false;
      CanopyStorage = ReadMap(LDD,getvaluename("smax"));
   }
   if (SwitchInterceptionLAI)
   {
      CanopyStorage = NewMap(0); //m
      FOR_ROW_COL_MV
      {
         switch (InterceptionLAIType)
         {
         case 0: CanopyStorage->Drc = 0.935+0.498*LAI->Drc-0.00575*(LAI->Drc * LAI->Drc);break;
         case 1: CanopyStorage->Drc = 0.2331 * LAI->Drc; break;
         case 2: CanopyStorage->Drc = 0.3165 * LAI->Drc; break;
         case 3: CanopyStorage->Drc = 1.46 * pow(LAI->Drc,0.56); break;
         case 4: CanopyStorage->Drc = 0.0918 * pow(LAI->Drc,1.04); break;
         case 5: CanopyStorage->Drc = 0.2856 * LAI->Drc; break;
         case 6: CanopyStorage->Drc = 0.1713 * LAI->Drc; break;
         case 7: CanopyStorage->Drc = 0.59 * pow(LAI->Drc,0.88); break;

         }
      }
   }
   CanopyStorage->calcV(0.001, MUL); // to m
   //NOTE: LAI is still needed for canopy openness, can be circumvented with cover


   //### infiltration maps
   InfilTot = 0;
   InfilTotmm = 0;
   InfilKWTot = 0;
   IntercTot = 0;
   IntercTotmm = 0;
   WaterVolTot = 0;
   WaterVolSoilTot = 0;
   WaterVolTotmm = 0;
   InfilVolKinWave = NewMap(0);
   InfilVol = NewMap(0);
   InfilVolCum = NewMap(0);
   fact = NewMap(0);
   fpot = NewMap(0);
   factgr = NewMap(0);
   fpotgr = NewMap(0);
   Ksateff = NewMap(0);
   FSurplus = NewMap(0);
   FFull = NewMap(0);

   if (InfilMethod != INFIL_SWATRE && InfilMethod != INFIL_NONE)
   {
      Fcum = NewMap(1e-10);
      L1 = NewMap(1e-10);
      L2 = NewMap(1e-10);
      Fcumgr = NewMap(1e-10);
      L1gr = NewMap(1e-10);
      L2gr = NewMap(1e-10);
      if (InfilMethod != INFIL_KSAT)
      {
         Soilwater = NewMap(0);
         Soilwater2 = NewMap(0);
         Soilwater->calc2(ThetaI1, SoilDepth1, MUL);
         if (SwitchTwoLayer)
         {
            Soilwater2->calc2(ThetaI2, SoilDepth2, MUL);
         }
      }
   }

   //### runoff maps
   Qtot = 0;
   QtotOutlet = 0;
   QtotPlot = 0;
   Qtotmm = 0;
   Qpeak = 0;
   QpeakTime = 0;
   WH = NewMap(0);
   WHbef = NewMap(0);
   WHrunoff = NewMap(0);
   WHrunoffCum = NewMap(0);
   WHstore = NewMap(0);
   WHroad = NewMap(0);
   WHGrass = NewMap(0);
   FlowWidth = NewMap(0);
   fpa = NewMap(0);
   V = NewMap(0);
   Alpha = NewMap(0);
   Q = NewMap(0);
   Qn = NewMap(0);
   Qoutput = NewMap(0);
   Qsoutput = NewMap(0);
   //Qoutflow = NewMap(0); // value of Qn*dt in pits only
   q = NewMap(0);
   R = NewMap(0);
   Perim = NewMap(0);
   WaterVolin = NewMap(0);
   WaterVolall = NewMap(0);

   SwatreSoilModel = NULL;
   SwatreSoilModelCrust = NULL;
   SwatreSoilModelCompact = NULL;
   SwatreSoilModelGrass = NULL;
   // swatre get input data is called before, ReadSwatreInput
   if (InfilMethod == INFIL_SWATRE)
   {
      thetaTop = NewMap(0);

      precision = 5.0;
      // note "5" is a precision factor dewtermining next timestep, set to 5 in old lisem

      // VJ 110420 added tiledrain depth for all profiles, is all used in infiltration
      SwatreSoilModel = InitSwatre(ProfileID, initheadName, TileDepth, swatreDT);
      if (SwatreSoilModel == NULL)
         throw 3;

      if (SwitchInfilCrust || SwitchWaterRepellency)
      {
         SwatreSoilModelCrust = InitSwatre(ProfileIDCrust, initheadName, TileDepth, swatreDT);
         if (SwatreSoilModelCrust == NULL)
            throw 3;
      }
      if (SwitchInfilCompact)
      {
         SwatreSoilModelCompact = InitSwatre(ProfileIDCompact, initheadName, TileDepth, swatreDT);
         if (SwatreSoilModelCompact == NULL)
            throw 3;
      }
      if (SwitchGrassStrip)
      {
         SwatreSoilModelGrass = InitSwatre(ProfileIDGrass, initheadName, TileDepth, swatreDT);
         if (SwatreSoilModelGrass == NULL)
            throw 3;
      }
      initSwatreStructure = true;
      // flag: structure is created and can be destroyed in function destroydata
   }

   //### erosion maps
   DetSplashTot = 0;
   DetFlowTot = 0;
   DepTot = 0;
   DetTot = 0;
   DepTot = 0;
   SoilLossTot = 0;
   SoilLossTotOutlet = 0;
   SedTot = 0;
   Qs = NewMap(0);
   Qsn = NewMap(0);
   //Qsoutflow = NewMap(0);
   DETSplash = NewMap(0);
   DETFlow = NewMap(0);
   DEP = NewMap(0);
   Sed = NewMap(0);
   TC = NewMap(0);
   Conc = NewMap(0);
   CG = NewMap(0);
   DG = NewMap(0);
   SettlingVelocity = NewMap(0);
   CohesionSoil = NewMap(0);
   Y = NewMap(0);

   TotalDetMap = NewMap(0);
   TotalDepMap = NewMap(0);
   TotalSoillossMap = NewMap(0);
   TotalSed = NewMap(0);
   TotalWatervol = NewMap(0);
   TotalConc = NewMap(0);


   if (SwitchErosion)
   {
      FOR_ROW_COL_MV
      {
         CG->Drc = pow((D50->Drc+5)/0.32, -0.6);
         DG->Drc = pow((D50->Drc+5)/300, 0.25);
         SettlingVelocity->Drc = 2*(2650-1000)*9.80*pow(D50->Drc/2000000, 2)/(9*0.001);
         CohesionSoil->Drc = Cohesion->Drc + Cover->Drc*RootCohesion->Drc;
         // soil cohesion everywhere, plantcohesion only where plants
         Y->Drc = min(1.0, 1.0/(0.89+0.56*CohesionSoil->Drc));
      }
   }

   if (useSorted)
      lddlist = makeSortedNetwork(LDD, &lddlistnr);

   //VJ 110113 all channel and buffer initialization moved to separate functions

}
//---------------------------------------------------------------------------
void TWorld::IntializeOptions(void)
{
   nrrainfallseries = 0;
   useSorted = false; // do not use alternative kin wave for now!

   //dirs and names
   resultDir.clear();
   inputDir.clear();
   outflowFileName = QString("totals.txt");//.clear();
   totalErosionFileName = QString("erososion.map");//.clear();
   totalDepositionFileName = QString("deposition.map");//.clear();
   totalSoillossFileName = QString("soilloss.map");//.clear();
   totalLandunitFileName = QString("totlandunit.txt");//.clear();
   outflowFileName = QString("hydrohgraph.csv");//.clear();
   rainFileName.clear();
   rainFileDir.clear();
   snowmeltFileName.clear();
   snowmeltFileDir.clear();
   SwatreTableDir.clear();
   SwatreTableName = QString("profile.inp");//.clear();
   resultFileName.clear();

   SwitchHardsurface = false;
   SwitchLimitTC = false;
   SwatreInitialized = false;
   SwitchInfilGA2 = false;
   SwitchWheelPresent = false;
   SwitchCompactPresent = false;
   SwitchIncludeChannel = false;
   SwitchChannelBaseflow = false;
   startbaseflowincrease = false;
   SwitchChannelInfil = false;
   SwitchAllinChannel = false;
   SwitchErosion = false;
   SwitchAltErosion = false;
   SwitchSimpleDepression = false;
   SwitchBuffers = false;
   SwitchSedtrap = false;
   SwitchRainfall = true; //VL 110103 add rainfall default true
   SwitchSnowmelt = false;
   SwitchRunoffPerM = false;
   SwitchInfilCompact = false;
   SwitchInfilCrust = false;
   SwitchGrassStrip = false;
   SwitchImpermeable = false;
   SwitchDumphead = false;
   SwitchWheelAsChannel = false;
   SwitchMulticlass = false;
   SwitchNutrients = false;
   SwitchGullies = false;
   SwitchGullyEqualWD = false;
   SwitchGullyInfil = false;
   SwitchGullyInit = false;
   SwitchOutputTimeStep = false;
   SwitchOutputTimeUser = false;
   SwitchMapoutRunoff = false;
   SwitchMapoutConc = false;
   SwitchMapoutWH = false;
   SwitchMapoutWHC = false;
   SwitchMapoutTC = false;
   SwitchMapoutEros = false;
   SwitchMapoutDepo = false;
   SwitchMapoutV = false;
   SwitchMapoutInf = false;
   SwitchMapoutSs = false;
   SwitchMapoutChvol = false;
   SwitchWritePCRnames = false;
   SwitchWriteCommaDelimited = true;
   SwitchWritePCRtimeplot = false;
   SwitchNoErosionOutlet = false;
   SwitchDrainage = false;
   SwitchPestout = false;
   SwitchSeparateOutput = false;
   SwitchSOBEKOutput = false;
   SwitchInterceptionLAI = false;
   SwitchTwoLayer = false;
   SwitchSimpleSedKinWave = false;
   SwitchSOBEKoutput = false;
   SwitchPCRoutput = false;
   SwitchSoilwater = false;
   SwitchGeometric = true;
   SwitchBacksubstitution = true;

   SwitchWriteHeaders = true; // write headers in output files in first timestep

   initSwatreStructure = false;
   // check to flag when swatre 3D structure is created, needed to clean up data
}
//---------------------------------------------------------------------------

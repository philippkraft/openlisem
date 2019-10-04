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
\file TMmapVariables.h
\brief List of maps with descriptions and units. Linked directly in the model class.
*/

// separte here for easier Doxygen comments

cTMap

//*_MASK,
*DEM,                        //!< DEM [m]
*Shade,                      //!< Shaded relief for display [0-1]
*ShadeBW,                      //!< Shaded relief for display [0-1]
*DX,                         //!< cell length divided by cosine slope (so corrected for terrain gradient) [m]
*CellArea,                   //!< cell area = DX * _dx [m^2]
*Grad,                       //!< sine of the DEM gradient [-]
*LDD,                        //!< local drain direction map [-]
*Outlet,                     //!< main outlet of the catchment, value 5 in LDD map [-]
*PointMap,                   //!< map with output points, values > 0 [-]
*FlowBoundary,               //!< map with open boundary fior diffusive runoff (1) or closed boundary (0)

*RainZone,                   //!< rainfall zone map (clasified map, numers corrspond to raingaug number in rainfall file) [-]
*Rain,                       //!< map with rain from tis time intervall [m]
*Rainc,                      //!< map with rain from tis time intervall, spread over the surface (corrected or slope) [m]
*RainCum,                    //!< cumulative rainfall, as spreadoutover slope [m]
*RainCumFlat,                //!< cumulative rainfall [m]
*RainNet,                    //!< net rainfall after interception [m]
*LeafDrain,                  //!< drainge from canopy, storage overflow [m]
*CStor,                      //!< actual canopy storage [m]
*Interc,                     //!< actual canopy storage volume, corrected for surfaces without vegetation (like roads) [m^3]
*LCStor,                     //!< actual Litter storage [m]
*LInterc,                    //!< actual Litter storage volume, corrected for surfaces without vegetation (like roads) [m^3]
*DStor,                      //!< actual drum storage of rainwater [m^3]
*HStor,                      //!< actual roof storage of rainwater [m]
*IntercHouse,                //!< actual roof storage volume [m^3]
*HouseCover,                 //!< fraction cover of house in pixel [-]
*HouseWidthDX,
*RoofStore,                  //!< Max storage of roof in [mm]
*DrumStore,                  //!< Max storage of rainwter drums [m^3]
*InterceptionmmCum,

*SnowmeltZone,               //!< snowmelt zone map, values corrspond to snowmelt gauge numbers [-]
*Snowcover,                  //!< snowmelt cover map, value 1.0 if there is snowcover, 0 without [-]
*Snowmelt,                   //!< snowmelt depth in water equivalent [m]
*Snowmeltc,                  //!< snowmelt depth in water equivalent, corrected for DEM gradient [m]
*SnowmeltCum,                //!< cumulative showmelt depth [m]

*WH,                         //!< water height on the surface [m]
*WHbef,                      //!< water height on the surface before infiltration [m]
*WHroad,                     //!< water height on the roads [m]
*WHhard,                     //!< water height on the roads [m]
*WHrunoff,                   //!< water height available for runoff [m]
*WHmax,                      //!< max runoff wh in m for reporting
*WHstore,                    //!< water heigth stored in micro depressions [m]
*WaterVolall,                //!< water volume total (incl surface storage) [m^3]
*WaterVolin,                 //!< water volume total before kin wave (after tochannel) [m^3]
//*WaterVolRunoff,                //!< water volume for runoff [m^3]

*FlowWidth,                  //!< width of the flow overland, based on ponded area/roughness, +roads etc [m]
*V,                          //!< velocity of overland flow [m/s]
*Alpha,                      //!< alpha in A = alphaQ^b
*Q,                          //!< discharge of overland flow before kin wave [m^3/s]
*Qn,                         //!< new discharge of overland flow after kin wave [m^3/s]
*VH,
//*Qoutflow,                   //!< new discharge after kin wave at outflow point [m^3/s]
*QinKW,
//*QoutKW,
*Qoutput,                    //!< new discharge for output purposes, sum of overland flow and channel, converted [l/s]
*Qs,                         //!< sediment discharge before kin wave [kg/s]
*Qsn,                        //!< new sediment discharge after kin wave [kg/s]
*Qsoutput,                   //!< sediment outflow for screen/file output, sum of overland flow and channel [kg/s]
*q,                          //!< infiltration surplus going in kin wave (<= 0) [m2/s]
*R,                          //!< hydraulic radius overland flow [m]
*N,                          //!< Manning's n
*RR,                         //!< Random roughness, locally converted to m [cm]
*MDS,                        //!< Maximum depression storage [m]
*fpa,                        //!< fraction ponded area [-]
*SoilWidthDX,                //!< width of soil surface, excluding roads and channels [m]
*RoadWidthDX,                //!< width of tarred roads [m]
*StoneFraction,              //!< fraction of stones on the surface, affects splash [-]
*CompactFraction,            //!< fraction compacted at the surface, uses ksat compact [-]
*CrustFraction,              //!< fraction crusted at the surface, uses ksat crust [-]
*RepellencyFraction,         //!< fraction of water repellency of node 1 in Swatre [-]
*RepellencyCell,             //!< Cell included in water repellency in Swatre [-]
*HardSurface,                //!< value 1 if 'hard' surface: no interception, infiltration, detachment [-]
*runoffTotalCell,

*PlantHeight,                //!< height of vegetation/crops [m]
*Cover,                      //!< vegetation canopy cover fraction [-]
*Litter,                     //!< vegetation litter cover fraction [-]
*CanopyStorage,              //!< canopy storage [m]
*LAI,                        //!< leaf area index [m^2/m^2]
*LandUnit,                   //!< land unit class (> 0) [-]

*Cohesion,                   //!< total cohesion of the soil surface: coh soil *(1-cover) + coh plant (cover) [kPa]
*RootCohesion,               //!< cohesion soil [kPa]
*CohesionSoil,               //!< cohesion by plant roots [kPa]
*Y,                          //!< erosion efficiency 0-1, basd on cohesion [-]
*AggrStab,                   //!< aggregate stability, median of drops in lowe test [-]
*D50,                        //!< median of grainsize distribution [mu]
*D90,                        //!< 90 % of grainsize distribution is below this value [mu]
*DETSplash,                  //!< splash detachment [kg/cell]
*DETSplashCum,
*DETFlow,                    //!< flow detachment [kg/cell]
*DETFlowCum,
*DEPCum,
//*DEPBLCum,
*DEP,                        //!< deposition [kg/cell]
*TC,                         //!< transport capacity [kg/m^3]
*Conc,                       //!< sediment concentration in flow [kg/m^3]
*Sed,                        //!< sediment content of flow [kg]
*CG,                         //!< parameter Govers in TC equation
*DG,                         //!< parameter Govers in TC equation
*SettlingVelocity,           //!< settling velocity according to Stokes [m/s]
*PCA,                        //!< applied dose [kg/m2]
*epsil,                      //!< mixing layer depth (m]
*KD,                         //!< soil water partition coefficient [m3/kg]
*kr,                         //!< rate at which solute desorb [min-1]
*rhob,                       //!< soil bulk density [kg/m3]
*C,                          //!< Pesticide concentration in dissolved form in runoff water [kg/m3]
*CM,                         //!< Pesticide concentration in dissolved form in the mixing zone [kg/m3]
*CS,                         //!< Pesticide concentration in sorbed form in the mixing zone [kg/m3]
*C_N,
*CM_N,
*CS_N,
*C_K,
*C_Kold,
*CM_K,
*CS_K,
*Qp,
*Qpn,
*C_Kn,
*K1,
*Kfilm,
*pestiinf,
*pestiinfold,
*poro,
*AX,
*Fkold,
*Fk,
*Fmk,
*flagpest,
*PMassApplied,
*PRunoffSpatial,
*PDisMixing,
*PSorMixing,
*PInfilt,
*PStorage,
*PRunoffSpatialex,
*PDisMixingex,
*PSorMixingex,
*PInfiltex,
*Qin,
*Sin,
*Pest,
*Fin,
*Pdetach,
*PCinfilt,
*PCfilmexit,

//for the kinematic wave in 2D
*K2DDEM,                        //!<
*K2DWHStore,                    //!<
*K2DPits,                       //!<
*K2DPitsD,                       //!<
*K2DOutlets,                    //!<
*K2DQM,                         //!<
*K2DQMX,                         //!<
*K2DQMY,                         //!<
*K2DMN,                         //!<
*K2DM,                         //!<
*K2DMC,                         //!<
*K2DQP,                         //!<
*K2DQPX,                         //!<
*K2DQPY,                         //!<
*K2DP,                          //!<
*K2DPC,                          //!<
*K2DPCN,                          //!<
*K2DSlopeX,                     //!<
*K2DSlopeY,                     //!<
*K2DSlope,                      //!<
*K2DHOld,                       //!<
*K2DHNew,                       //!<
*K2DQX,                         //!<
*K2DQY,                         //!<
*K2DQ,                          //!<
*K2DDTm,                        //!< Temporary local timestep (s)
*K2DDTr,                        //!< Real previous timestep(s)
*K2DDT,                         //!< Current local timestep (might be remainder) (s)
*K2DDTT,                        //!< Current total time since last lisem timestep (s)
*K2DQN,
*K2DFX,                         //!< Flux in x-direction (m3/s)
*K2DFY,                         //!< Flux in y-direction (m3/s)
*K2DDTR,                        //!< the row nr of cells with a dt
*K2DDTC,                        //!< the column nr of cells with a dt
*K2DI,                          //!<
*K2DTEST,
*K2DV,

// infiltration
*Fcum,                       //!< cumulative infiltration [m]
*FSurplus,                   //!< surplus infiltration for kinematic wave, calculated as actual infil - potential infil [m]
*FFull,                      //!< map flagging when the soil is full
*fact,                       //!< actual infiltration rate [m/s]
*fpot,                       //!< potential infiltration rate [m/s]
*InfilVolKinWave,            //!< volume infiltrated in the kin wave (slope and channel) in this timestep [m^3]
*InfilVol,                   //!< volume of water infiltrated in this timestep [m^3]
*InfilVolCum,                //!< cumulative infiltration volume for mass balance and map report [m^3]
*InfilmmCum,                 //!< cumulative infiltration volume for map report and drawing [mm]
*InfilVolFlood,

*ThetaS1,                    //!< porosity soil layer 1 [-]
*ThetaI1,                    //!< initial moisture content soil layer 1 [-]
*Psi1,                       //!< intial suction head wetting front soil layer 1 (input map is in cm) [m]
*Ksat1,                      //!< saturated hydraulic conductivity soil layer 1 (input is in mm/h) [m/s]
*SoilDepth1,                 //!< depth to end soil layer 1 (input is in mm) [m]
*L1,                         //!< depth of wetting front in layer 1 [m]
*Soilwater,                  //!< actual soil water content [-]
*ThetaSub,
*ThetaS2,                    //!< porosity soil layer 2 [-]
*ThetaI2,                    //!< initial moisture content soil layer 2 [-]
*Psi2,                       //!< intial suction head wetting front soil layer 2 (input map is in cm) [m]
*Ksat2,                      //!< saturated hydraulic conductivity soil layer 2 (input is in mm/h) [m/s]
*SoilDepth2,                 //!< depth to end soil layer 2 (input is in mm) [m]
*L2,                         //!< depth of wetting front in layer 2 [m]
*Soilwater2,                  //!< actual soil water content layer 2 [-]
*KsatCrust,                  //!< saturated hydraulic conductivity crusted soil surface (input is in mm/h) [m/s]
*PoreCrust,                //!< saturated hydraulic conductivity compacted soil surface (input is in mm/h) [m/s]
*KsatCompact,                //!< saturated hydraulic conductivity compacted soil surface (input is in mm/h) [m/s]
*PoreCompact,                //!< saturated hydraulic conductivity compacted soil surface (input is in mm/h) [m/s]
*KsatGrass,                  //!< saturated hydraulic conductivity grass strip (input is in mm/h) [m/s]
*PoreGrass,                //!< saturated hydraulic conductivity compacted soil surface (input is in mm/h) [m/s]
*Ksateff,                    //!< effective saturated hydraulic conductivity (input is in mm/h) [m/s]
*Poreeff,
*Thetaeff,
*Perc,
*PercmmCum,
*factgr,                     //!< actual infiltration rate fo grassstrip [m/s]
*fpotgr,                     //!< potential infiltration rate fo grassstrip [m/s]
*WHGrass,                    //!< water level on a grassstrip [m]
*GrassFraction,              //!< fraction of grasstrip in a cell [-]
*SedimentFilter,
*GrassWidthDX,               //!< width of grasstrip in [m]
*thetaTop,                   //!< average theta of node 0 and 1 for water repelency and nutrients

*ProfileID,                  //!< SWATRE profile unit number map
*ProfileIDCrust,             //!< SWATRE profile unit number map for crusted areas
*ProfileIDCompact,           //!< SWATRE profile unit number map for compacted areas
*ProfileIDGrass,             //!< SWATRE profile unit number map for grass strips
*SwatreOutput,

*LDDChannel,                 //!<
*RunoffVolinToChannel,       //!<
*ChannelWidth,               //!<
*ChannelSide,                //!<
*ChannelQ,                   //!<
*ChannelQn,                  //!<
*ChannelQntot,
*ChannelQs,                  //!<
*ChannelQsn,                 //!<
*ChannelQBLs,                  //!<
*ChannelQBLsn,                 //!<
*ChannelQSSs,                  //!<
*ChannelQSSsn,                 //!<
*ChannelGrad,                //!<
*ChannelV,                   //!<
*ChannelN,                   //!<
*ChannelWH,                  //!<
*ChannelWaterVol,            //!<
//*ChannelBLWaterVol,            //!<
//*ChannelSSWaterVol,            //!<
*Channelq,                   //!<
*ChannelAlpha,               //!<
*ChannelFlowWidth,           //!<
*ChannelWidthMax,           //!<
*ChannelAdj,                //!<
*ChannelPerimeter,           //!<
*ChannelDX,                  //!<
*ChannelKsat,                //!<
*ChannelDetFlow,             //!<
*ChannelDep,                 //!<
*ChannelSed,                 //!<
*ChannelBLSed,                 //!<
*ChannelSSSed,                 //!<
*ChannelBLTC,                 //!<
*ChannelSSTC,                 //!<
*ChannelBLDepth,                 //!<
*ChannelSSDepth,                 //!<
*ChannelConc,                //!<
*ChannelBLConc,                //!<
*ChannelSSConc,                //!<
*ChannelTC,                  //!<
//*SedToChannel,               //!<
*ChannelCohesion,            //!<
*ChannelY,                   //!<
*ChannelDepth,               //!<

//baseflow
*BaseFlowDischarges,
*BaseFlowInitialVolume,
*BaseFlowInflow,

// flood maps
*UVflood,                     //!<
*Qflood,                    //!<
*floodHmxMax,                    //!<
*floodTime,                    //!<
*floodTimeStart,                //!<
*floodVMax,                    //!<
*floodVHMax,                    //!<
*maxChannelflow,                    //!<
*maxChannelWH,                    //!<
*hmx,                        //!<
*hmxWH,                        //!<
*hmxInit,                    //!<
*hmxflood,
*FloodDomain,                //!<
*Barriers,                    //!<
*ChannelMaxQ,                //!<
//*ChannelLevee,                //!<
*FloodWaterVol,                //!<
//*FloodZonePotential,                //!<
*DomainEdge,                //!<
*FloodDT,
*FloodDTr,
*FloodT,
*FloodHMaskDer,
*FloodDTR,
*FloodDTC,
*FloodHR,
*FloodHC,


*f1o, *f2o, *f3o,
*g1o, *g2o, *g3o,
*VRO, *URO,
*iro,

// FULLSWOF2D
*z1r, *z1l, *z2r, *z2l,
*h1r, *h1l, *h2r, *h2l,
*h1d, *h1g, *h2d, *h2g,
*v1r, *v1l, *v2r, *v2l,
*u1r, *u1l, *u2r, *u2l,
*delta_z1, *delta_z2,
*delzc1, *delzc2,
*delz1, *delz2,
*f1, *f2, *f3, *cflx,
*g1, *g2, *g3, *cfly,
*hs, *vs, *us,
*hsa, *vsa, *usa,
*Uflood,*Vflood,*Iflood,

//FULLSWOF2D with Sediment
*BLDepthFlood,
*SSDepthFlood,
//*BLVolFlood,
//*SSVolFlood,
//*temp1,*temp2,*temp3,*temp4,
//*temp5,*temp6,*temp7,*temp8,
//*temp9,*temp10,*temp11,*temp12,

//layer itneraction
//*BLDepFlood,
*BLDetFlood,
//*BLDetFloodTot,
*BLTCFlood,
*SSTCFlood,
*SSDetFlood,
//*SSDetFloodTot,
*DepFlood,

//sediment maps
*BLCFlood,
*BLFlood,
*SSCFlood,
*SSFlood,

//transport
//Bed Load layer
*MBLCFlood, *MBLNFlood,
// *MBLCNFlood,*MBLFlood,
//*bl1r,*bl1l,*bl2r,*bl2l,*blf1,*blg1,*bls,*bls2,*bl1d, *bl1g, *bl2d, *bl2g,

//transport
//Suspended Sediment Layer
*MSSCFlood,*MSSNFlood,
//*MSSCNFlood,*MSSFlood,
//*ss1r,*ss1l,*ss2r,*ss2l,*ssf1,*ssg1,*sss,*sss2,*ss1d, *ss1g, *ss2d, *ss2g,

//*q1flood,*q2flood,
*som_z1,*som_z2,

*LDDTile,                    //!< LDD network of tile drains, must be connected to outlet
*TileDrainSoil,              //!< drain volume from layer
*TileDiameter,                  //!< total width of drains in cell (m)
*TileMaxQ,
*TileWidth,                  //!< total width of drains in cell (m)
*TileHeight,                 //!< height of drain (m)
*TileDepth,                  //!< depth of tiles in soil below surface (m)
*TileSinkhole,               //!< sinkhole on surface connecting to tiledrains (m2)
*TileQ,                      //!< water flux in drains m3/s
*TileQn,                     //!< new water flux in drains m3/s
*TileQs,                     //!< sediment flux in drains kg/s
*TileQsn,                    //!< new sediment flux in drains kg/s
//*TileQoutflow,               //!< water outflow in outlet
*TileGrad,                   //!< gradient of the tiledrain system
*TileN,                      //!< mannings inside the tiledrains
*TileWH,                     //!< water height in the tile drains (m)
*TileWaterVol,               //!< water volume in the tiledrains (m3)
*TileWaterVolSoil,           //!< water volume in the tiledrains from the soil only, used for mass bal corection (m3)
*Tileq,                      //!< possible drainage inside tiles, not used
*RunoffVolinToTile,          //!< can be used for shortcut of surface pits to tile system
*TileAlpha,                  //!< alpha in tile drain, in A = alpha*Q^beta
*TileDX,                     //!< cell length in tile drain, dx/cos angle
*TileV,                      //!< velocity in tile drain m/s

*TotalDetMap,                //!<
*TotalDepMap,                //!<
*TotalChanDetMap,                //!<
*TotalChanDepMap,                //!<
*TotalSoillossMap,           //!<
*TotalSed,                   //!<
*TotalConc,                  //!<

*tm,                         //!< Auxilary map
*tma,                        //!< Auxilary map
*tmb,                        //!< Auxilary map
*tmc,                        //!< Auxilary map
*tmd,                        //!< Auxilary map
*CoreMask,
//display combinations
*COMBO_QOFCH,
*COMBO_VOFCH,
*COMBO_SS,
*COMBO_TC,
*ChannelDepthExtended,
*ChannelWidthExtended,
*ChannelNeighborsExtended,
*ChannelSourceXExtended,
*ChannelSourceYExtended,
*ChannelMaskExtended,
*ChannelBoundaryExtended,
*ChannelBoundaryLExtended,
*ChannelBoundaryRExtended,

*FlowBarrier,                //!< Flow barriers type
*FlowBarrierN,               //!< Flow barriers height North of cell
*FlowBarrierW,               //!< Flow barriers height West of cell
*FlowBarrierS,               //!< Flow barriers height South of cell
*FlowBarrierE,               //!< Flow barriers height East of cell
*FlowBarrierNT,              //!< Flow barriers end timing North of cell
*FlowBarrierWT,              //!< Flow barriers end timing West of cell
*FlowBarrierST,              //!< Flow barriers end timing South of cell
*FlowBarrierET               //!< Flow barriers end timing East of cell

;

cTRGBMap * RGB_Image;

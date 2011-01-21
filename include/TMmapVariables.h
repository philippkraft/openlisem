// List of maps in openLIsem
// included in the model.h TWorld class definition
// separte here for easier Doxygen comments

TMMap
*tm,                         //!< Auxilary map
*tma,                        //!< Auxilary map
*tmb,                        //!< Auxilary map
*Mask,                       //!< Catchment mask, based on the LDD
*DEM,                        //!< DEM [m]
*DX,                         //!< cell length divided by cosine slope (so corrected for terrain gradient) [m]
*CellArea,                   //!< cell area = DX * _dx [m^2]
*Grad,                       //!< sine of the DEM gradient [-]
*LDD,                        //!< local drain direction map [-]
*Outlet,                     //!< main outlet of the catchment, value 5 in LDD map [-]
*PointMap,                   //!< map with output points, values > 0 [-]

*RainZone,                   //!< rainfall zone map (clasified map, numers corrspond to raingaug number in rainfall file) [-]
*Rain,                       //!< map with rain from tis time intervall [m]
*Rainc,                      //!< map with rain from tis time intervall, spread over the surface (corrected or slope) [m]
*RainCum,                    //!< cumulative rainfall [m]
*RainNet,                    //!< net rainfall after interception [m]
*LeafDrain,                  //!< drainge from canopy, storage overflow [m]
*CStor,                      //!< actual canopy storage [m]
*Interc,                     //!< actual canopy storage volume, corrected for surfaces without vegetation (like roads) [m^3]

*SnowmeltZone,               //!< snowmelt zone map, values corrspond to snowmelt gauge numbers [-]
*Snowcover,                  //!< snowmelt cover map, value 1.0 if there is snowcover, 0 without [-]
*Snowmelt,                   //!< snowmelt depth in water equivalent [in m]
*Snowmeltc,                  //!< snowmelt depth in water equivalent, corrected for DEM gradient[in m]
*SnowmeltCum,                //!< cumulative showmelt depth [m]

*WH,                         //!< water height on the surface [m]
*WHroad,                     //!< water height on the roads [m]
*WHrunoff,                   //!< water height available for runoff [m]
*WHrunoffCum,                //!< cumulative runoff water height for output to file series only [mm]
*WHstore,                    //!< water heigth stored in micro depressions [m]
*WaterVolall,                //!< water volume total (incl surface storage) [m^3]
*WaterVolin,                 //!< water volume total before kin wave (after tochannel) [m^3]

*FlowWidth,                  //!< width of the flow overland, based on ponded area/roughness, +roads etc [m]
*V,                          //!< velocity of overland flow [m/s]
*Alpha,                      //!< alpha in A = alphaQ^b
*Q,                          //!< discharge of overland flow before kin wave [m^3/s]
*Qn,                         //!< new discharge of overland flow after kin wave [m^3/s]
*Qoutflow,                   //!< new discharge after kin wave at outflow point [m^3/s]
*Qoutput,                    //!< new discharge for output purposes, sum of overland flow and channel, converted [l/s]
*Qs,                         //!< sediment discharge before kin wave [kg/s]
*Qsn,                        //!< new sediment discharge after kin wave [kg/s]
*Qsoutflow,                  //!< new sediment discharge after kin wave at outflow point [kg/s]
*Qsoutput,                   //!< sediemnt outflow for screen/file output, sum of overland flow and channel [kg/s]
*q,                          //!< infiltration surplus going in kin wave (<= 0) [m2/s]
*R,                          //!< hydraulic radius overland flow [m]
*Perim,                      //!< perimeter overland flow [m]
*N,                          //!< Manning's n
*RR,                         //!< Random roughness, locally converted to m [cm]
*MDS,                        //!< Maximum depression storage [m]
*fpa,                        //!< fraction ponded area [-]
*SoilWidthDX,                //!< width of soil surface, excluding roads and channels [m]
*RoadWidthDX,                //!< width of tarred roads [m]
*StoneFraction,              //!< fraction of stones on the surface, affects splash [-]
*CompactFraction,            //!< fraction compacted at the surface, uses ksat compact [-]
*CrustFraction,              //!< fraction crusted at the surface, uses ksat crust [-] 
*HardSurface,                //!< value 1 if 'hard' surface: no interception, infiltration, detachment [-]

*PlantHeight,                //!< height of vegetation/crops [m]
*Cover,                      //!< vegetation canopy cover fraction [m]
*CanopyStorage,              //!< canopy storage [m]
*LAI,                        //!< leaf area index [m^2/m^2]
*LandUnit,                   //!< land unit class (> 0) [-]
*WheelWidth,                 //!< not used yet, width of wheel tracks [m]
*WheelWidthDX,               //!< not used yet, width of wheel tracks [m]
*GullyWidthDX,               //!< not used yet, width of gullies [m]

*Cohesion,                   //!< total cohesion of the soil surface: coh soil *(1-cover) + coh plant (cover) [kPa]
*RootCohesion,               //!< cohesion soil [kPa]
*CohesionSoil,               //!< cohesion by plant roots [kPa]
*Y,                          //!< erosion efficiency 0-1, basd on cohesion [-]
*AggrStab,                   //!< aggregate stability, median of drops in lowe test [-] 
*D50,                        //!< median of grainsize distribution [mu]
*DETSplash,                  //!< splash detachment [kg/m^2]
*DETFlow,                    //!< flow detachment [kg/m^2]
*DEP,                        //!< deposition [kg/m^2]
*TC,                         //!< transport capacity [kg/m^3]
*Conc,                       //!< sediment concentration in flow [kg/m^3]
*Sed,                        //!< sediment content of flow [kg]
*CG,                         //!< parameter Govers in TC equation
*DG,                         //!< parameter Govers in TC equation
*SettlingVelocity,           //!< settling velocity according to Stokes [m/s]

*Fcum,                       //!< cumulative infiltration [m]
*FSurplus,                   //!< surplus infiltration for kinematic wave, calculated as actual infil - potential infil [m]
*FFull,                      //!< map flagging when the soil is full
*fact,                       //!< actual infiltration rate [m/s]
*fpot,                       //!< potential infiltration rate [m/s]
*InfilVolKinWave,            //!< volume infiltrated in the kin wave (slope and channel) in this timestep [m^3]
*InfilVol,                   //!< volume of water infiltrated in this timestep [m^3]
*InfilVolCum,                //!< cumulative infiltration volume for mass balance and map report [m^3]

*ThetaS1,                    //!< porosity soil layer 1 [-]
*ThetaI1,                    //!< initial moisture content soil layer 1 [-]
*Psi1,                       //!< intial suction head wetting front soil layer 1 (input map is in cm) [m]
*Ksat1,                      //!< saturated hydraulic conductivity soil layer 1 (input is in mm/h) [m/s]
*SoilDepth1,                 //!< depth to end soil layer 1 (input is in mm) [m]
*L1,                         //!< depth of wetting front in layer 1 [m]
*Soilwater,                  //!< actual soil water content [-]

*ThetaS2,                    //!< porosity soil layer 2 [-]
*ThetaI2,                    //!< initial moisture content soil layer 2 [-]
*Psi2,                       //!< intial suction head wetting front soil layer 2 (input map is in cm) [m]
*Ksat2,                      //!< saturated hydraulic conductivity soil layer 2 (input is in mm/h) [m/s]
*SoilDepth2,                 //!< depth to end soil layer 2 (input is in mm) [m]
*L2,                         //!< depth of wetting front in layer 2 [m]
*Soilwater2,                  //!< actual soil water content layer 2 [-]

*KsatCrust,                  //!< saturated hydraulic conductivity crusted soil surface (input is in mm/h) [m/s]
*KsatCompact,                //!< saturated hydraulic conductivity compacted soil surface (input is in mm/h) [m/s]
*KsatGrass,                  //!< saturated hydraulic conductivity grass strip (input is in mm/h) [m/s]
*Ksateff,                    //!< effective saturated hydraulic conductivity (input is in mm/h) [m/s]
*L1gr,                       //!< depth wetting front under grass strip layer 1 [m]
*L2gr,                       //!< depth wetting front under grass strip layer 2 [m]
*factgr,                     //!< actual infiltration rate fo grassstrip [m/s]
*fpotgr,                     //!< potential infiltration rate fo grassstrip [m/s]
*Fcumgr,                     //!< cumulative infiltration under grassstrips [m]
*WHGrass,                    //!< water level on a grassstrip [m]
*GrassFraction,              //!< fraction of grasstrip in a cell [-]
*GrassWidthDX,               //!< width of grasstrip in [m]

*ProfileID,                  //!< SWATRE profile unit number map
*ProfileIDCrust,             //!< SWATRE profile unit number map for crusted areas
*ProfileIDCompact,           //!< SWATRE profile unit number map for compacted areas
*ProfileIDGrass,             //!< SWATRE profile unit number map for grass strips

*ChannelMask,                //!< 
*RunoffVolinToChannel,       //!<
*LDDChannel,                 //!<
*ChannelWidth,               //!<
*ChannelSide,                //!<
*ChannelQ,                   //!<
*ChannelQn,                  //!<
*ChannelQs,                  //!<
*ChannelQsn,                 //!<
*ChannelQoutflow,            //!<
*ChannelGrad,                //!<
*ChannelV,                   //!<
*ChannelN,                   //!<
*ChannelWH,                  //!<
*ChannelWaterVol,            //!<
*Channelq,                   //!<
*ChannelAlpha,               //!<
*ChannelWidthUpDX,           //!<
*ChannelPerimeter,           //!<
*ChannelDX,                  //!<
*ChannelKsat,                //!<
*ChannelDetFlow,             //!<
*ChannelDep,                 //!<
*ChannelSed,                 //!<
*ChannelConc,                //!<
*ChannelTC,                  //!<
*SedToChannel,               //!<
*ChannelQsoutflow,           //!<
*ChannelCohesion,            //!<
*ChannelY,                   //!<

*TileMask,                   //!<
*TileDrainSoil,              //!<
*LDDTile,                    //!<
*TileWidth,                  //!<
*Tileheight,                 //!<
*TileQ,                      //!<
*TileQn,                     //!<
*TileQs,                     //!<
*TileQsn,                    //!<
*TileQoutflow,               //!<
*TileGrad,                   //!<
*TileN,                      //!<
*TileWH,                     //!<
*TileWaterVol,               //!<
*Tileq,                      //!<
*RunoffVolinToTile,          //!<
*TileAlpha,                  //!<
*TileDX,                     //!<

*BufferID,                   //!<
*BufferVol,                  //!<
*BufferSed,                  //!<
*ChannelBufferSed,           //!<
*ChannelBufferVol,           //!<

*BufferVolInit,              //!<
*BufferSedInit,              //!<
*ChannelBufferSedInit,       //!<
*ChannelBufferVolInit,       //!<

*TotalDetMap,                //!<
*TotalDepMap,                //!<
*TotalSoillossMap,           //!<
*TotalSed,                   //!<
*TotalWatervol,              //!<
*TotalConc;                  //!<

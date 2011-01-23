/*---------------------------------------------------------------------------
project: openLISEM
author: Victor Jetten
licence: GNU General Public License (GPL)
Developed in: MingW/Qt/ 
website, information and code: http://sourceforge.net/projects/lisem
---------------------------------------------------------------------------*/

/*! \file */
/**
Functionality in model.h:
* - TWorld class that combines ALL model variables
* - all processes
* - global defines Drc, MV, FOR_ROW_COL_MV etc
* - global defines for lisem type; infiltration type etc.
*/

#ifndef modelH
#define modelH

#include <math.h>
#include <stdlib.h>

#include <QtGui>
#include <QMutex.h>

#include "csfmap.h"
#include "mmath.h"
#include "error.h"
#include "swatre_p.h"
#include "swatre_g.h"

//---------------------------------------------------------------------------

#define DEBUG(s) emit debug(QString(s))

#define mwrite(name) WriteMap(QString(resultDir+name))
#define report(name) WriteMapSeries(resultDir,QString(name), printstep)

// defines to make life easier

/// shortcut to access data
#define Drc     Data[r][c]

/// shortcut to access the outlet point data
#define DrcOutlet     Data[r_outlet][c_outlet]

/// shortcut missing value in map
#define MV(r,c) IS_MV_REAL8(&LDD->Data[r][c])
//IS_MV_REAL8(&Mask->Data[r][c])

/// shortcut for LDD row and col loop
#define FOR_ROW_COL_MV for (int r = 0; r < _nrRows; r++)\
   for (int c = 0; c < _nrCols; c++)\
   if(!IS_MV_REAL8(&LDD->Data[r][c]))
//	if(!IS_MV_REAL4(&Mask->Data[r][c]))

/// shortcut for channel row and col loop
#define FOR_ROW_COL_MV_CH for (int  r = 0; r < _nrRows; r++)\
   for (int  c = 0; c < _nrCols; c++)\
   if(!IS_MV_REAL8(&LDDChannel->Data[r][c]))
//		 if(!IS_MV_REAL4(& ChannelMask->Data[r][c]))

/// shortcut for tile network row and col loop.
#define FOR_ROW_COL_MV_TILE for (int  r = 0; r < _nrRows; r++)\
   for (int  c = 0; c < _nrCols; c++)\
   if(!IS_MV_REAL8(&LDDTile->Data[r][c]))
//      if(!IS_MV_REAL8(& TileMask->Data[r][c]))

#define NUMNAMES 300   /// \def NUMNAMES runfile namelist max
#define NUMMAPS 256    /// \def max nr maps
#define MIN_FLUX 1e-12 /// \def minimum flux (m3/s) in kinematic wave
#define MIN_HEIGHT 1e-6 /// \def minimum water height (m) for transport of sediment
#define MAXCONC 848    /// \def max concentration susp. sed. in kg/m3 0.32 * 2650 = max vol conc from experiments Govers x bulk density


#define INFIL_NONE 0
#define INFIL_SWATRE 1
#define INFIL_HOLTAN 2
#define INFIL_GREENAMPT 3
#define INFIL_GREENAMPT2 4
#define INFIL_KSAT 5
#define INFIL_MOREL 21
#define INFIL_SMITH 22
#define INFIL_SMITH2 23


//---------------------------------------------------------------------------
/// structure containing pointers to all maps

/** structure containing pointers to all maps for automatic destruction after runs
 so memory doesn't have to be freed for each map. The functions Newmap(double) and
ReadMap(cTMap *Mask, QString name) put a map on this list
*/
typedef struct MapListStruct {
	TMMap *m;
}  MapListStruct;
//---------------------------------------------------------------------------
/// linked list structure for network in kin wave
typedef struct LDD_LINKEDLIST {
	int rowNr;
	int colNr;
   struct LDD_LINKEDLIST *prev;
}  LDD_LINKEDLIST;
//---------------------------------------------------------------------------
typedef struct LDD_POINT {
   int rowNr;
   int colNr;
}  LDD_POINT;
//---------------------------------------------------------------------------
/// name list structure used to read run file
typedef struct NAME_LIST {
	QString name;
	QString value;
} NAME_LIST;
//---------------------------------------------------------------------------
/// structure for output of land unit stats
typedef struct UNIT_LIST {
   long nr;
   double area;
   double totdet;
   double totdep;
   double totsl;
} UNIT_LIST;
//---------------------------------------------------------------------------



//---------------------------------------------------------------------------
/// \class TWorld model.h contains the model 'World': constants, variables and erosion processes

/** The model 'world': the main class containing all variables, maps, options, filenames.\n
  The class contains hydrological and erosion processes which are run in a time loop.\n
  The main function is <B>void TWorld::DoModel()</B>, which has the time loop calling ll processes\n
  Every timestep the mass alance is calculated and output is reported to the UI and disk.
 */

//http://blog.exys.org/entries/2010/QThread_affinity.html
//http://thesmithfam.org/blog/2009/09/30/lock-free-multi-threading-in-qt/

class TWorld: public QThread
{
	Q_OBJECT

public:
	TWorld(QObject *parent = 0);
	~TWorld();

   /// copy of overall rows and columns, set in initmask
   long _nrRows;
   long _nrCols;

   /// map management structure, automatic adding and deleting of all TMMap variables
   MapListStruct maplistTMMap[NUMNAMES];
	int maplistnr;


#include "TMmapVariables.h"

   /// SwitchXXX are boolean options that are set in interface and runfile, mainly corrsponding to checkboxes in the UI
	bool SwitchHardsurface, SwatreInitialized, SwitchInfilGA2, SwitchCrustPresent,
	SwitchWheelPresent, SwitchCompactPresent, SwitchIncludeChannel, SwitchChannelBaseflow,
	startbaseflowincrease, SwitchChannelInfil, SwitchAllinChannel, SwitchErosion, SwitchAltErosion,
	SwitchSimpleDepression, SwitchBuffers, SwitchSedtrap, SwitchSnowmelt, SwitchRainfall, SwitchRunoffPerM, SwitchInfilCompact,
   SwitchInfilCrust, SwitchGrassStrip, SwitchImpermeable, SwitchDumphead,
	SwitchWheelAsChannel, SwitchMulticlass, SwitchNutrients, SwitchGullies, SwitchGullyEqualWD, SwitchGullyInfil,
	SwitchGullyInit, SwitchOutputTimeStep, SwitchOutputTimeUser, SwitchMapoutRunoff, SwitchMapoutConc,
	SwitchMapoutWH, SwitchMapoutWHC, SwitchMapoutTC, SwitchMapoutEros, SwitchMapoutDepo, SwitchMapoutV,
	SwitchMapoutInf, SwitchMapoutSs, SwitchMapoutChvol, SwitchWritePCRnames, SwitchWritePCRtimeplot,
	SwitchNoErosionOutlet, SwitchDrainage, SwitchPestout, SwitchSeparateOutput, SwitchSOBEKOutput,
	SwitchInterceptionLAI, SwitchTwoLayer, SwitchSimpleSedKinWave, SwitchSoilwater, SwitchSOBEKoutput,
   SwitchPCRoutput, SwitchWriteHeaders, SwitchGeometric, SwitchIncludeTile;

   // multiple options that are set in interface or runfile, see defines above
   /// Interception storage function based on LAI
   int InterceptionLAIType;
   /// infiltration method
   int InfilMethod;
   /// erosion units in output: to/ha; kg/cell; kg/m2
   int ErosionUnits;
   /// calibration factors
	double ksatCalibration;
	double nCalibration;
	double ChnCalibration;
	double ChKsatCalibration;
	double SplashDelivery;
	double StripN;
	double StemflowFraction;

   /// totals for mass balance checks and output
   /// Water totals for mass balance and output (in m3)
	double MB, Qtot, QtotOutlet, IntercTot, WaterVolTot, InfilTot, RainTot, SnowTot, SurfStoremm, InfilKWTot;
   double ChannelVolTot, ChannelSedTot, ChannelDepTot, ChannelDetTot, TileVolTot;
   /// Sediment totals for mass balance and output (in kg)
   double MBs, DetTot, DetSplashTot, DetFlowTot, DepTot, SoilLossTot, SoilLossTotOutlet, SedTot;
   /// Water totals for output in file and UI (in mm), copied to 'op' structure
	double RainTotmm, SnowTotmm, IntercTotmm, WaterVolTotmm, InfilTotmm, Qtotmm, RainAvgmm, SnowAvgmm;
   /// peak times (min)
	double RainpeakTime, SnowpeakTime, QpeakTime, Qpeak, Rainpeak, Snowpeak;
	double BufferVolTot, BufferSedTot, BufferVolTotInit, BufferSedTotInit, BulkDens;
	double nrCells, CatchmentArea;

   int c_outlet;  /// copy of outlet col number
   int r_outlet;  /// copy of outlet row number

   /// time and dx parameters
	double time, BeginTime, EndTime;
	double _dt, _dx;
	long runstep, printstep, printinterval;

   /// timeseries variables and output strings
	double **RainfallSeries;
   int nrstations, nrrainfallseries;
	double **SnowmeltSeries;
   int nrSnowmeltstations, nrSnowmeltseries;

   // output formatting for SOBEK flood model input
	QString SOBEKdatestring;
	int SOBEKnrlines;

   // file and directory names
	QString resultDir;
	QString inputDir;
	QString outflowFileName;
	QString totalErosionFileName;
	QString totalDepositionFileName;
	QString totalSoillossFileName;
   QString totalLandunitFileName;
   QString rainFileName;
	QString rainFileDir;
	QString snowmeltFileName;
	QString snowmeltFileDir;
	QString resultFileName;
	QString temprunname;
	QStringList outputcheck;
   /// standard names of output map series
	QString Outrunoff, Outconc, Outwh, Outrwh, Outtc, Outeros, Outdepo, Outvelo, Outinf, Outss, Outchvol;

   // list with class values of land unit map
   UNIT_LIST unitList[512]; // just a fixed number for 512 classes, who cares!
   int landUnitNr;

	// data initialization, runfile reading and parsing
   NAME_LIST runnamelist[NUMNAMES]; // structure for runfile variables and names
   int nrrunnamelist;
	void IntializeOptions(void);  // set all options to false etc
	void IntializeData(void);     // make all non-input maps
	void GetRainfallData(void);   // get input timeseries
	void GetSnowmeltData(void);   // get input timeseries
	void GetInputData(void);      // get and make input maps
   void InitChannel(void); //VJ 110112
   void InitBuffers(void); //VJ 110112
   void InitTiledrains(void); //VJ 110112
	TMMap *InitMask(QString name);
	TMMap *InitMaskChannel(QString name);
   TMMap *InitMaskTiledrain(QString name);
	TMMap *ReadMapMask(QString name);
	TMMap *ReadMap(cTMap *Mask, QString name);
	TMMap *NewMap(double value);
	void InitMapList(void);
	void DestroyData(void);
   void ParseRunfileData();
	void GetRunFile();
	QString GetName(QString p);
	QString getvaluename(QString vname);
	double getvaluedouble(QString vname);
	int getvalueint(QString vname);
	QString CheckDir(QString p);

   // LISEM model processes
   /// convert rainfall of a timestep into a map
   void RainfallMap(void);
   /// convert snowmelt of a timestep into a map
   void SnowmeltMap(void);
   /// interception of vegetation canopy resulting in rainnet
   void Interception(void);
   /// infiltration function calling all infiltration methods
   /// add rainnet to WH and calculating new WH
   void Infiltration(void);
   void InfilSwatre(void);
   void InfilGreenAmpt1(void);
   void InfilSmithParlange1(void);
   void InfilMorelSeytoux1(void);
   void InfilKsat(void);
   double IncreaseInfiltrationDepth(int r, int c, double fact, REAL8 *L1p, REAL8 *L2p);
   void SoilWater(void);
   void SurfaceStorage(void);
   void OverlandFlow(void);
   void ChannelFlow(void);
   void ToChannel(void);
   void CalcVelDisch(void);
   void CalcVelDischChannel(void);
   void TileFlow(void);
   void CalcVelDischTile(void);
   void GridCell(void);
   void SplashDetachment(void);
   void FlowDetachment(void);
   double MaxConcentration(double watvol, double sedvol, double dep);
   void ChannelFlowDetachment(void);
   void Kinematic(int pitRowNr, int pitColNr, TMMap *_LDD, TMMap *_Q, TMMap *_Qn, TMMap *_Qs,
                  TMMap *_Qsn, TMMap *_q, TMMap *_Alpha, TMMap *_DX, TMMap *Vol, TMMap*SedVol,
                  TMMap *_StorVol, TMMap*_StorVolSed);
   void KinematicNew(LDD_POINT **_lddlist, long _lddlistnr,
                     TMMap *_Q, TMMap *_Qn, TMMap *_Qs,
                     TMMap *_Qsn, TMMap *_q, TMMap *_Alpha, TMMap *_DX, TMMap *Vol, TMMap*SedVol,
                     TMMap *_StorVol, TMMap*_StorVolSed);
   double simpleSedCalc(double Qj1i1, double Qj1i, double Sj1i, double dt, double vol, double sed);
   double complexSedCalc(double Qj1i1, double Qj1i, double Qji1, double Sj1i,
                         double Sji1, double alpha, double dt, double dx);
   double IterateToQnew(double Qin, double Qold, double q, double alpha, double deltaT, double deltaX);

   //VJ 110123 sorted networks for faster kin wave
   LDD_POINT **MakeSortedNetwork(TMMap *_LDD, long *lddlistnr);
   LDD_POINT **lddlist;
   long lddlistnr;
   LDD_POINT **lddlistch;
   long lddlistchnr;
   LDD_POINT **lddlisttile;
   long lddlisttilenr;

   //SWATRE
   /// filenames for Swatre soil information
   QString SwatreTableDir;
   QString SwatreTableName;
   QString initheadName;

   double swatreDT;
   bool initSwatreStructure;

   /// SWATRE infiltration model 3D soil structure
   SOIL_MODEL *SwatreSoilModel;
   SOIL_MODEL *SwatreSoilModelCrust;
   SOIL_MODEL *SwatreSoilModelCompact;
   SOIL_MODEL *SwatreSoilModelGrass;
   PROFILE **profileList;
   int nrProfileList, sizeProfileList;
   ZONE *zone;
   double precision;
   int tnode; //VJ 110122 node nr in profile with tile drains

   SOIL_MODEL *InitSwatre(TMMap *profileMap, QString initHeadMaps, double dtMin);
   int ReadSwatreInput(QString fileName, QString tablePath);
   void SwatreStep(SOIL_MODEL *s, TMMap *_WH, TMMap *_fpot, TMMap *_drain, TMMap *where);
   void CloseSwatre(SOIL_MODEL *s);
   void FreeSwatreInfo();
   ZONE *ReadNodeDefinition(FILE *f);
   PROFILE *ReadProfileDefinition(FILE *f, ZONE *z, const char *tablePath);
   HORIZON *ReadHorizon(const char *tablePath,	const  char *tableName);
   PROFILE *ProfileNr(int profileNr);
   double *ReadSoilTable(const char *fileName, int *nrRows);
   void ReadCols(const char *fileName, double *inLut, const char *buf, int   n);
   void InitializeProfile();
   void HeadCalc(double *h, bool *ponded, const PROFILE *p ,const double  *thetaPrev,
                 const double  *hPrev, const double  *kavg, const double  *dimoca,
                 bool fltsat, double dt, double pond, double qtop, double qbot);
   double  NewTimeStep(double prevDt, const double *hLast, const double *h, int nrNodes,
                       double precParam, double dtMin, double dtMax);
   void ComputeForPixel(PIXEL_INFO *pixel, double *waterHeightIO, double *infil, double *drain,
                        SOIL_MODEL *s);

   void Totals(void);
	void MassBalance(void);
	void OutputUI(void);
   void reportAll();
	void ReportTimeseriesNew();
	void ReportTotals();
	void ReportMaps();
	void ReportTotalsNew();
   void ReportLandunits(); //VJ 110107 report erosion stats per land unit
   void CountLandunits(); //VJ 110107 report erosion stats per land unit

   double itercount;
   // thread management variables
   bool stopRequested;
   bool waitRequested;
   QMutex mutex;
   QWaitCondition condition;
   void stop();

protected:
	void run();
	QTime time_ms;

   // talk to the interface
signals:
	void done(const QString &results);
	void debug(const QString &results);
   void show(); //use the output structure "op" declared in global.h and LisUIoutput.h

public slots:   //note, was private loop but dixygen does not recognize that
   /// the main model loop, from here all processes are called in a time loop
   void DoModel();

};



#endif

#ifndef LISEMQT_H
#define LISEMQT_H

#include <QtGui>
#include <QApplication>

#include "qwt_plot.h"
#include "qwt_plot_curve.h"
#include "qwt_plot_grid.h"
#include "qwt_plot_marker.h"
#include "qwt_legend.h"
//#include "qwt_plot_spectrogram.h"
//#include <qwt_matrix_raster_data.h>
//#include "qwt_color_map.h"
//#include "qwt_raster_data.h"

#include "ui_lisemqt.h"
#include "LisUItreemodel.h"
#include "model.h"
#include "lisuioutput.h"


// constants to define the place of the main parts in the map tree structure
#define RAINFALLMAPS 0
#define CATCHMENTMAPS 1
#define LANDUSEMAPS 2
#define SURFACEMAPS 3
#define EROSIONMAPS 4
#define INFILTRATIONMAPS 5
#define CHANNELSMAPS 6
#define BUFFERSMAPS 7
#define SNOWMELTMAPS 8
#define TILEDRAINMAPS 9   //VJ 110111
#define WHEELTRACKSMAPS 10
#define TEXTUREMAPS 11
#define NUTRIENTSMAPS 12
#define GULLIESMAPS 13

/*

class SpectrogramData: public QwtRasterData
{
public:
	 TMMap *DMap;

	SpectrogramData():
        QwtRasterData(QwtDoubleRect(0,0,100,100))
    {
    }

    virtual QwtRasterData *copy() const
    {
        return new SpectrogramData();
    }

    virtual QwtDoubleInterval range() const
    {
        return QwtDoubleInterval(0.0, 10.0);
    }

    virtual double value(double x, double y) const
    {
		 int r = qRound(y);
		 int c = qRound(x);
		 return(DMap->Data[r][c]);
    }
};

*/
class lisemqt : public QMainWindow, private Ui::lisemqtClass
{
	Q_OBJECT

public:
	lisemqt(QWidget *parent = 0);
	~lisemqt();

   void FillMapTree();
	void DefaultMapnames();
	void change_MapNameModel(int parentrow, int selrow, bool setit);
	void SetToolBar();
	void GetStorePath();
	void StorePath();
	void SetStyleUI();
	void SetGraph();
	void SetMapPlot();
	void GetRunfile();
	void ParseInputData();
	void UpdateModelData();
	void DefaultRunFile();
	QString CheckDir(QString p);
	void RunAllChecks();
	void savefile(QString name);
	void InitOP();
   void SetConnections();

	void ShowGraph();
   void ShowMap();

	// graph variables
	QwtPlot *HPlot;
	QwtPlotCurve *QGraph;
	QwtPlotCurve *QsGraph;
	QwtPlotCurve *CGraph;
	QwtPlotCurve *PGraph;
	bool startplot;
	double yas, y2as;
	double *timeData;
	double *QData;
	double *QsData;
	double *CData;
	double *PData;
	//Map drawing variables
//   QwtPlotSpectrogram *MapDrawing;
//   QwtPlot *MapPlot;
//	SpectrogramData *MapDrawData;
//   QwtMatrixRasterData *MapDrawData;

	bool oldRunfile; // check is old runfile for ksat calibration

	//interface names
	TreeModel *MapNameModel;
	QString currentDir;
	QString RainFileName;
	QString RainFileDir;
	QString SnowmeltFileName;
	QString SnowmeltFileDir;
	QString SwatreTableName;
	QString SwatreTableDir;
	QStringList DEFmaps;
	QStringList RunFileNames;
	int CurrentRunFile;
	int uiInfilMethod;
	double swatreDT;

   _mapList mapList[NUMMAPS]; // structure for current map names, can be edited by user
   int nrmaplist;
   _nameList namelist[NUMNAMES]; // structure to read runfile variables and names, used in LisUIrunfile.cpp
   int nrnamelist;
   QStringList outputcheck;
	int InterceptionEqNr;

public slots:
	// functions linked to actions
	void SaveRunFile();
	void savefileas();
	void openRunFile();
	void runmodel();
	void stopmodel();
	void pausemodel();
	void shootScreen();
	void aboutQT();
	void aboutInfo();
	void resetAll();

   void doChangeMapname(QModelIndex topLeft, QModelIndex bottomRight );
   void doOpenMapname(QModelIndex topLeft);

	void on_toolButton_MapDir_clicked();
	void on_toolButton_ResultDir_clicked();
	void on_toolButton_RainfallName_clicked();
	void on_toolButton_SnowmeltName_clicked();
	void on_toolButton_RainfallShow_clicked();
	void on_toolButton_SnowmeltShow_clicked();
	void on_toolButton_ShowRunfile_clicked();
	//void on_toolButton_fileOpen_clicked();
	void on_toolButton_SwatreTableDir_clicked();
	void on_toolButton_SwatreTableFile_clicked();
	void on_toolButton_SwatreTableShow_clicked();

   void doCheckSnowmelt(bool check);
   void doCheckRainfall(bool check);

	void on_E_InfiltrationMethod_currentIndexChanged(int inr);
	void on_E_runFileList_currentIndexChanged(int);

	void on_checkChannelInfil_clicked();
	void on_checkChannelBaseflow_clicked();
	void on_checkNoErosion_clicked();
	void on_checkIncludeChannel_clicked();
   void on_checkIncludeTiledrains_clicked();
	void on_checkInfilCompact_clicked();
	void on_checkInfilCrust_clicked();
	void on_checkInfilGrass_clicked();
	void on_checkInfil2layer_clicked();
	void on_checkBuffers_clicked();
	void on_checkSedtrap_clicked();
	void on_checkSnowmelt_clicked();
	void on_checkExpandActive_clicked();
	void on_E_MapDir_textEdited();
	void on_E_ResultDir_textEdited();

private slots:
	// functions that interact with the world thread
	void Showit();
	void worldDone(const QString &results);
	void worldDebug(const QString &results);


private:
	//toolbar actions
	QAction *openAct;
	QAction *saveAct;
	QAction *saveasAct;
	QAction *runAct;
	QAction *pauseAct;
	QAction *stopAct;
	QAction *shootscreenAct;
	QAction *aboutAct;
	QAction *aboutActI;
	QAction *restartAct;

	// the model world
	TWorld *W;

};


#endif // LISEMQT_H

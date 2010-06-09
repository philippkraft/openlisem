#ifndef LISEMQT_H
#define LISEMQT_H

#include <QtGui>
#include <QApplication>

#include "qwt_plot.h"
#include "qwt_plot_curve.h"
#include "qwt_plot_grid.h"
#include "qwt_plot_marker.h"
#include "qwt_legend.h"

#include "ui_lisemqt.h"
#include "LisUItreemodel.h"
#include "model.h"
#include "lisuioutput.h"

class lisemqt : public QMainWindow, private Ui::lisemqtClass
{
	Q_OBJECT

public:
	lisemqt(QWidget *parent = 0);
	~lisemqt();

	void FillMapList();
	void LisemReset();
	void DefaultMapnames();
	void change_MapNameModel(int parentrow, int selrow, bool setit);
	void SetToolBar();
	void GetStorePath();
	void StorePath();
	void SetStyleUI();
	void SetGraph();
	void GetRunfile();
	void ParseInputData();
	void UpdateModelData();
	QString CheckDir(QString p, QString p1);
	void RunAllChecks();
	void savefile(QString name);
	void InitOP();


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

	//interface names
	TreeModel *MapNameModel;
	QString currentDir;
	QString RainFileName;
	QString rainFileDir;
	QString SnowmeltFileName;
	QString snowmeltFileDir;
	QStringList DEFmaps;
	QStringList OUTmaps;
	QStringList OUTMCmaps;
	QStringList OUTNUTmaps;
	QStringList OUTGULmaps;
	QStringList RunFileNames;
	int CurrentRunFile;
	int uiInfilMethod;

	// runfile read structure
	_nameList namelist[NUMNAMES]; // structure for runfile variables and names
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

void on_toolButton_MapDir_clicked();
void on_toolButton_ResultDir_clicked();
void on_toolButton_RainfallName_clicked();
void on_toolButton_SnowmeltName_clicked();
void on_toolButton_RainfallShow_clicked();
void on_toolButton_SnowmeltShow_clicked();
void on_toolButton_ShowRunfile_clicked();
void on_toolButton_fileOpen_clicked();

void on_E_InfiltrationMethod_currentIndexChanged(int);
void on_E_runFileList_currentIndexChanged(int);

void on_checkChannelInfil_clicked();
void on_checkChannelBaseflow_clicked();
void on_checkNoErosion_clicked();
void on_checkIncludeChannel_clicked();
void on_checkInfilCompact_clicked();
void on_checkInfilCrust_clicked();
void on_checkInfilGrass_clicked();
void on_checkInfil2layer_clicked();
void on_checkBuffers_clicked();
void on_checkSedtrap_clicked();
void on_checkSnowmelt_clicked();
void on_checkExpandActive_clicked();
void on_toolButton_SwatreTable_clicked();
void on_E_MapDir_textEdited();
void on_E_ResultDir_textEdited();

private slots:
// functions that interact with the world thread
void Showit();
void worldDone(const QString &results);
void worldDebug(const QString &results);
void ShowGraph();

private:
//toolbar actions
QAction *openAct;
QAction *saveAct;
QAction *saveasAct;
QAction *runAct;
QAction *pauseAct;
QAction *stopAct;

// the model world
TWorld *W;

};


#endif // LISEMQT_H

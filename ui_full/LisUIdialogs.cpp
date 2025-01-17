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

#include "lisemqt.h"
//#include "model.h"
#include "global.h"

#define HELPOPTIONS 1
#define HELPRAINFALL 2
#define HELPINTERCEPTION 3
#define HELPINFILTRATION 4
#define HELPFLOW 5
#define HELPCHANNEL 6
#define HELPEROSION 7
#define HELPINFRA 8
#define HELPCALIBRATION 9
#define HELPADVANCED 10



//)---------------------------------------------------------------
void lisemqt::doResetAll()
{
    op.runfilename.clear();
    E_runFileList->clear();
    resetAll();
}
//---------------------------------------------------------------------------
void lisemqt::on_E_runFileList_currentIndexChanged(int)
{
    if (E_runFileList->count() == 0)
        return;
    if (E_runFileList->currentText() == "")
        return;
    CurrentRunFile = E_runFileList->currentIndex();
    op.runfilename = E_runFileList->currentText();
    //RunFileNames.at(CurrentRunFile);

    GetRunfile();   // get the nrunfile and fill namelist

    ParseInputData(); // fill interface with namelist data and fill mapList
    // also update DEFmaps for map tree view in interface

    initMapTree();  // fill the tree strcuture on page 2 with DEFmaps
    RunAllChecks(); // activate the maps in the tree parts in response to checks
}
//--------------------------------------------------------------------
void lisemqt::on_E_MapDir_returnPressed()
{
    QFileInfo fin(E_MapDir->text());
    if(!fin.exists())
    {
        E_MapDir->setText("");
        QMessageBox::warning(this,"openLISEM",
                             QString("Map directory does not exist"));
    }
}
//--------------------------------------------------------------------
void lisemqt::on_E_ResultDir_returnPressed()
{
    if (E_ResultDir->text().isEmpty())
        return;
    QFileInfo fin(E_ResultDir->text());
    if(!fin.exists())
    {
        int ret =
                QMessageBox::question(this, QString("openLISEM"),
                                      QString("The directory \"%1\"does not exist.\n"
                                              "Do you want to create it (apply)?")
                                      .arg(fin.absoluteFilePath()),
                                      QMessageBox::Apply |QMessageBox::Cancel,QMessageBox::Cancel);
        if (ret == QMessageBox::Apply)
            QDir(E_ResultDir->text()).mkpath(E_ResultDir->text());

    }
}
//---------------------------------------------------------------
void lisemqt::on_toolButton_resetOptions_clicked()
{
    resetTabOptions();
}
//---------------------------------------------------------------
void lisemqt::on_toolButton_resetRainfall_clicked()
{
    resetTabRainfall();
}
//---------------------------------------------------------------
void lisemqt::on_toolButton_resetInterception_clicked()
{
    resetTabInterception();
}
//---------------------------------------------------------------
void lisemqt::on_toolButton_resetInfiltration_clicked()
{
    resetTabInfiltration();
}
//---------------------------------------------------------------
void lisemqt::on_toolButton_resetFlow_clicked()
{
    resetTabFlow();
}
//---------------------------------------------------------------
void lisemqt::on_toolButton_resetChannel_clicked()
{
   resetTabChannel();
}
//---------------------------------------------------------------
void lisemqt::on_toolButton_resetInfra_clicked()
{
   resetTabInfra();
}
//---------------------------------------------------------------
void lisemqt::on_toolButton_resetErosion_clicked()
{
    resetTabErosion();
}
//---------------------------------------------------------------
void lisemqt::on_toolButton_resetCalibration_clicked()
{
    resetTabCalibration();
}
//---------------------------------------------------------------
void lisemqt::on_toolButton_resetAdvanced_clicked()
{
    resetTabAdvanced();
}
//---------------------------------------------------------------
void lisemqt::on_toolButton_helpOptions_clicked()
{
    on_toolButton_help(HELPOPTIONS);
}
//---------------------------------------------------------------
void lisemqt::on_toolButton_helpRainfall_clicked()
{
    on_toolButton_help(HELPRAINFALL);
}
//---------------------------------------------------------------
void lisemqt::on_toolButton_helpInterception_clicked()
{
    on_toolButton_help(HELPINTERCEPTION);
}
//---------------------------------------------------------------
void lisemqt::on_toolButton_helpInfiltration_clicked()
{
    on_toolButton_help(HELPINFILTRATION);
}
//---------------------------------------------------------------
void lisemqt::on_toolButton_helpFlow_clicked()
{
    on_toolButton_help(HELPFLOW);
}
//---------------------------------------------------------------
void lisemqt::on_toolButton_helpChannel_clicked()
{
    on_toolButton_help(HELPCHANNEL);
}
//---------------------------------------------------------------
void lisemqt::on_toolButton_helpErosion_clicked()
{
    on_toolButton_help(HELPEROSION);
}
//---------------------------------------------------------------
void lisemqt::on_toolButton_helpCalibration_clicked()
{
    on_toolButton_help(HELPCALIBRATION);
}
//---------------------------------------------------------------
void lisemqt::on_toolButton_helpInfra_clicked()
{
    on_toolButton_help(HELPINFRA);

}
//---------------------------------------------------------------
void lisemqt::on_toolButton_helpAdvanced_clicked()
{
    on_toolButton_help(HELPADVANCED);
}
//---------------------------------------------------------------
void lisemqt::on_toolButton_help(int page)
{
    QString filename;
    if (page == HELPOPTIONS     ) filename = ":/help1.html";
    if (page == HELPRAINFALL    ) filename = ":/help6.html";
    if (page == HELPINTERCEPTION) filename = ":/help2.html";
    if (page == HELPINFILTRATION) filename = ":/help3.html";
    if (page == HELPFLOW        ) filename = ":/help4.html";
    if (page == HELPCHANNEL     ) filename = ":/help9.html";
    if (page == HELPEROSION     ) filename = ":/help5.html";
    if (page == HELPINFRA       ) filename = ":/help10.html";
    if (page == HELPCALIBRATION ) filename = ":/help7.html";
    if (page == HELPADVANCED    ) filename = ":/help8.html";

    QFile file(filename);
    file.open(QFile::ReadOnly | QFile::Text);
    QTextStream stream(&file);
    helptxt->setHtml(stream.readAll());
  //  helpbox->show();

    QScreen *screen = QGuiApplication::primaryScreen();
    QRect screenGeometry = screen->geometry();
    int screenWidth = screenGeometry.width();
    int screenHeight = screenGeometry.height();

    QTextEdit *view = new QTextEdit(helptxt->toHtml());
    view->createStandardContextMenu();
    view->setWindowTitle("Option help");

    view->setMinimumWidth(screenWidth/3*2);
    view->setMinimumHeight(screenHeight/3*2);
    view->setAttribute(Qt::WA_DeleteOnClose);

    view->show();
}
//--------------------------------------------------------------------

void lisemqt::on_check2DDiagonalFlow_toggled(bool checked)
{
    E_pitValue->setEnabled(checked);
    label_135->setEnabled(checked);
}
//--------------------------------------------------------------------

//void lisemqt::on_checkDiffusion_toggled(bool checked)
//{
//    E_SigmaDiffusion->setEnabled(checked);
//    label_101->setEnabled(checked);
//    label_139->setEnabled(checked);
//}
//--------------------------------------------------------------------

// void lisemqt::on_checkHouses_toggled(bool checked)
// {
//     checkRaindrum->setEnabled(checked);
//     label_78->setEnabled(checked);
// }
//--------------------------------------------------------------------

// select a file or directory
// doFile = 0: select a directory;
// dofile = 1 select a file and return file name only;
// dofile = 2 return filename with full path
QString lisemqt::getFileorDir(QString inputdir,QString title, QStringList filters, int doFile)
{
    QFileDialog dialog;

    QString dirout = inputdir;

    if (doFile > 0) {
        dialog.setNameFilters(filters);
        dialog.setDirectory(QFileInfo(inputdir).absoluteDir());
        dialog.setFileMode(QFileDialog::ExistingFile);
    } else {
        filters.clear();
        dialog.setNameFilters(filters);
        dialog.setDirectory(QDir(inputdir));
        dialog.setFileMode(QFileDialog::Directory);
    }

    dialog.setLabelText(QFileDialog::LookIn,title);
    dialog.exec();
    if (dialog.selectedFiles().isEmpty())
        dirout = inputdir;
    else
        dirout = dialog.selectedFiles().at(0);

    if (doFile > 0) {
        if (doFile == 1)
            dirout = QFileInfo(dirout).fileName();
        if (doFile == 2)
            dirout = QFileInfo(dirout).absoluteFilePath();
    } else {
        dirout = dialog.selectedUrls().at(0).path();
        dirout.remove(0,1);
        if (dirout.lastIndexOf('/') != dirout.length())
            dirout = dirout + "/";
    }

    return dirout;
}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_rainsatName_clicked()
{
    //RainSatFileDir = RainFileDir;
    if (!QFileInfo(RainSatFileDir).exists() || RainSatFileDir.isEmpty())
        RainSatFileDir = RainFileDir;
    if (!QFileInfo(RainSatFileDir).exists() || RainSatFileDir.isEmpty())
        RainSatFileDir = currentDir;
  //  qDebug() << RainSatFileDir << RainSatFileName << currentDir;

    QStringList filters({"Text file (*.txt *.tbl *.tss)","Any files (*)"});
    QString sss = getFileorDir(RainSatFileDir,"Select rainfall map list table", filters, 2);

    RainSatFileDir = QFileInfo(sss).absolutePath()+"/";
    RainSatFileName = QFileInfo(sss).fileName(); //baseName();

    E_RainsatName->setText(RainSatFileDir + RainSatFileName);
}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_ETsatName_clicked()
{
  //  ETSatFileDir = ETFileDir;

    if (!QFileInfo(ETSatFileDir).exists() || ETSatFileDir.isEmpty())
        ETSatFileDir = RainSatFileDir;
    if (!QFileInfo(ETSatFileDir).exists() || ETSatFileDir.isEmpty())
        ETSatFileDir = currentDir;
    QStringList filters({"Text file (*.txt *.tbl *.tss)","Any files (*)"});

    QString sss = getFileorDir(ETSatFileDir,"Select ET map list table", filters, 2);

    ETSatFileDir = QFileInfo(sss).absolutePath()+"/";
    ETSatFileName = QFileInfo(sss).fileName(); //baseName();

    E_ETsatName->setText(ETSatFileDir + ETSatFileName);
}

//--------------------------------------------------------------------
void lisemqt::on_toolButton_RainfallName_clicked()
{
    if (!QFileInfo(RainFileDir).exists() || RainFileDir.isEmpty())
        RainFileDir = currentDir;

    QStringList filters({"Text file (*.txt *.tbl *.tss)","Any files (*)"});
    QString sss = getFileorDir(RainFileDir,"Select rainfall station table", filters, 2);

    RainFileDir = QFileInfo(sss).absolutePath()+"/";
    RainFileName = QFileInfo(sss).fileName(); //baseName();

    E_RainfallName->setText(RainFileDir + RainFileName);

}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_ETName_clicked()
{
    if (!QFileInfo(ETFileDir).exists() || ETFileDir.isEmpty())
        ETFileDir = currentDir;

    QStringList filters({"Text file (*.txt *.tbl *.tss)","Any files (*)"});

    QString sss = getFileorDir(ETFileDir,"Select ET stations file", filters, 2);

    ETFileDir = QFileInfo(sss).absolutePath()+"/";
    ETFileName = QFileInfo(sss).fileName(); //baseName();
    E_ETName->setText(ETFileDir + ETFileName);
}
//--------------------------------------------------------------------
// void lisemqt::on_checkWaveInUser_toggled(bool checked)
// {
//     groupWaveUser->setEnabled(checked);
// }
//--------------------------------------------------------------------
void lisemqt::on_toolButton_DischargeShow_clicked()
{
    showTextfile(DischargeinDir + DischargeinFileName);
}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_WaveShow_clicked()
{
    //qDebug() <<WaveinDir + WaveinFileName;
    showTextfile(WaveinDir + WaveinFileName);
}
//--------------------------------------------------------------------
// void lisemqt::on_checkIncludeET_toggled(bool checked)
// {
//     widgetEToptions->setEnabled(checked);
// }
//--------------------------------------------------------------------
void lisemqt::on_toolButton_ETShow_clicked()
{
    //qDebug() <<ETFileDir + ETFileName;
    showTextfile(ETFileDir + ETFileName);
}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_DischargeName_clicked()
{
    if (!QFileInfo(DischargeinDir).exists() || DischargeinDir.isEmpty())
        DischargeinDir = currentDir;

    QStringList filters({"Text file (*.txt *.tbl *.tss)","Any files (*)"});

    QString sss = getFileorDir(DischargeinDir,"Select ET stations file", filters, 2);

    DischargeinDir = QFileInfo(sss).absolutePath()+"/";
    DischargeinFileName = QFileInfo(sss).fileName(); //baseName();
    E_DischargeInName->setText(DischargeinDir + DischargeinFileName);

}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_WaveInName_clicked()
{
    if (!QFileInfo(WaveinDir).exists() || WaveinDir.isEmpty())
        WaveinDir = currentDir;

    QStringList filters({"Text file (*.txt *.tbl *.tss)","Any files (*)"});

    QString sss = getFileorDir(WaveinDir,"Select file with water heights", filters, 2);

    WaveinDir = QFileInfo(sss).absolutePath()+"/";
    WaveinFileName = QFileInfo(sss).fileName();
    E_WaveInName->setText(WaveinDir + WaveinFileName);

}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_RainfallShow_clicked()
{
    showTextfile(RainFileDir + RainFileName);
}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_RainmapShow_clicked()
{
    showTextfile(RainSatFileDir + RainSatFileName);
}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_ETmapShow_clicked()
{
    showTextfile(ETSatFileDir + ETSatFileName);
}
//--------------------------------------------------------------------
void lisemqt::showTextfile(QString name)
{
   // Read text from file
    QFile file(name);
    file.open(QIODevice::ReadOnly | QIODevice::Text);
    QTextStream in(&file);
    QString initialText = in.readAll();
    file.close();

    // Find the longest line in the text
    QStringList lines = initialText.split("\n");
    int maxWidth = 0;
    for (const QString& line : lines) {
        maxWidth = qMax(maxWidth, QFontMetrics(QFont()).horizontalAdvance(line));
    }
    int maxHeight = qMin(qApp->primaryScreen()->size().height() * 2 / 3, lines.size() * 40);
    // Create a QDialog instance
    QDialog dialog;
    dialog.setWindowTitle(name);
    dialog.resize(maxWidth, maxHeight);//qApp->desktop()->width()/3*2,qApp->desktop()->height()/3*2);

    // Create a layout for the dialog
    QVBoxLayout *layout = new QVBoxLayout(&dialog);

    // Create a QTextEdit widget for editing text and load the initial text
    QTextEdit *textEdit = new QTextEdit(&dialog);
    textEdit->setPlainText(initialText); // Load initial text
    layout->addWidget(textEdit);

    // Create buttons layout
    QHBoxLayout *buttonLayout = new QHBoxLayout;

    // Create a QPushButton to save changes
    QPushButton *saveButton = new QPushButton("Save", &dialog);
    buttonLayout->addWidget(saveButton);

    // Create a QPushButton to cancel changes
    QPushButton *cancelButton = new QPushButton("Cancel", &dialog);
    buttonLayout->addWidget(cancelButton);

    QPushButton *closeButton = new QPushButton("Close", &dialog);
    buttonLayout->addWidget(closeButton);

    layout->addLayout(buttonLayout);

    // Connect the saveButton's clicked signal to a lambda function
    QObject::connect(saveButton, &QPushButton::clicked, [&]() {
        // Retrieve the text from the textEdit widget
        QString editedText = textEdit->toPlainText();
        // Open the file in write mode to save changes
        QFile writeFile(name);
        if (writeFile.open(QIODevice::WriteOnly | QIODevice::Text)) {
            QTextStream out(&writeFile);
            out << editedText;
            writeFile.close();
        }
        // Close the dialog when saving is done
        //dialog.close();
    });

    // Connect the cancelButton's clicked signal to close the dialog without saving changes
    QObject::connect(cancelButton, &QPushButton::clicked, [&]() {
        // Close the dialog without saving changes
        dialog.close();
    });

    QObject::connect(closeButton, &QPushButton::clicked, [&]() {
        dialog.close();
    });

    // Show the dialog
    dialog.show();

    dialog.exec();
}
//--------------------------------------------------------------------
//OBSOLETE
void lisemqt::showTextfileOld(QString name)
{

    QFile file(name);
    if (!file.open(QFile::ReadWrite | QFile::Text))
    {
        QMessageBox::warning(this, QString("openLISEM"),
                             QString("Cannot read file %1:\n%2.")
                             .arg(name)
                             .arg(file.errorString()));
        return;
    }

    QString modifiedContents;
    QTextStream in(&file);

    QPlainTextEdit *view = new QPlainTextEdit(in.readAll());
    view->setWindowTitle(RainFileName);
    view->setMinimumWidth(400);
    view->setMinimumHeight(500);
    view->setAttribute(Qt::WA_DeleteOnClose);
    view->show();
    if (view->document()->isModified())
    {
        int ret =
                QMessageBox::question(this, QString("openLISEM"),
                                      QString("You have modified the contents of this file.\n"
                                              "Do you want to save it?"),
                                      QMessageBox::Ok |QMessageBox::Cancel,QMessageBox::Cancel);
        if (ret == QMessageBox::Ok)
        {
            // Don't take the address of a temporary!
            // in.setString(&view->toPlainText());
            modifiedContents = view->toPlainText();
            in.setString(&modifiedContents);
        }
    }

    file.close();
}
//--------------------------------------------------------------------
void lisemqt::on_E_EndTimeDay_returnPressed()
{
    int daye = E_EndTimeDay->text().split(":")[0].toInt();
    int mine = E_EndTimeDay->text().split(":")[1].toInt();
    daye = std::max(1,std::min(daye, 366));
    if (mine > 1440) {
        daye = mine/1440 + 1;
        mine = mine % 1440;
    }
    E_EndTimeDay->setText(QString("%1:%2").arg(daye,3,10,QLatin1Char('0')).arg(mine,4,10,QLatin1Char('0')));
}
//--------------------------------------------------------------------
void lisemqt::on_E_BeginTimeDay_returnPressed()
{
       int daye = E_BeginTimeDay->text().split(":")[0].toInt();
       int mine = E_BeginTimeDay->text().split(":")[1].toInt();
       daye = std::max(1,std::min(daye, 366));
       if (mine > 1440) {
           daye = mine/1440 + 1;
           mine = mine % 1440;
       }
       E_BeginTimeDay->setText(QString("%1:%2").arg(daye,3,10,QLatin1Char('0')).arg(mine,4,10,QLatin1Char('0')));
}
//--------------------------------------------------------------------
void lisemqt::on_checkStationaryBaseflow_toggled(bool checked)
{
    if (checked) checkChannelInfil->setChecked(false);
   // doChannelBaseflow = checked;
}
//--------------------------------------------------------------------
void lisemqt::on_checkChannelInfil_toggled(bool checked)
{
    if (checked) checkStationaryBaseflow->setChecked(false);
}
//--------------------------------------------------------------------
void lisemqt::on_E_EfficiencyDETCH_currentIndexChanged(int index)
{
    E_EfficiencyDirect->setEnabled(index == 3);
}
//--------------------------------------------------------------------
void lisemqt::on_checkGWflow_toggled(bool checked)
{
    GW_widget->setEnabled(checked);
    widget_GWparams->setEnabled(checked);
    groupBaseflowParams->setEnabled(checked);
    //qDebug() << checked;
}
//--------------------------------------------------------------------
void lisemqt::on_E_floodMinHeight_valueChanged(double)
{
    label_107->setText(QString("Flood (mm),h>%1)").arg(E_floodMinHeight->value()*1000));
    label_40->setText(QString("Runoff (mm),h<%1)").arg(E_floodMinHeight->value()*1000));
}

//--------------------------------------------------------------------
void lisemqt::on_toolButton_SwatreTableDir_clicked()
{
    if (!QFileInfo(SwatreTableDir).exists() || SwatreTableDir.isEmpty())
        SwatreTableDir = currentDir;

    QStringList filters({"profile table files (*.tbl)","Any files (*)"});
    QString sss = getFileorDir(SwatreTableDir,"Select the SWATRE profile tabel directory", filters, 0);

    SwatreTableDir = QFileInfo(sss).absolutePath()+"/";

    E_SwatreTableDir->setText(SwatreTableDir);
}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_SwatreTableShow_clicked()
{
    showTextfile(SwatreTableName);
}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_SwatreTableName_clicked()
{
    // if (!QFileInfo(RainFileDir).exists() || RainFileDir.isEmpty())
    //     RainFileDir = currentDir;

    QStringList filters({"Text file (*.inp *.txt *.tbl)","Any files (*)"});
    QString sss = getFileorDir(currentDir,"Select Swatre profile file (def. profile.inp)", filters, 2);

    if(sss.isEmpty()) sss = "profile.inp";
    SwatreTableName = sss;//QFileInfo(sss).fileName(); //baseName();

    E_SwatreTableName->setText(SwatreTableName);
}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_satImageName_clicked()
{
    QString path;

    satImageFileDir = findValidDir(satImageFileDir, false);

    path = QFileDialog::getOpenFileName(this,
                                        QString("Select background satellite image file"),
                                        satImageFileDir,"GeoTiff (*.tif)");
    if(!path.isEmpty())
    {
        QFileInfo fi(path);
        satImageFileName = fi.fileName();
        satImageFileDir = CheckDir(fi.absolutePath(), false);//Dir().path());
        E_satImageName->setText( satImageFileDir + satImageFileName );
    }
}
//---------------------------------------------------------------------------
void lisemqt::on_checkRainfall_toggled(bool checked)
{
    groupRainfall->setEnabled(checked);
}
//---------------------------------------------------------------------------
void lisemqt::on_checkET_toggled(bool checked)
{
    groupET->setEnabled(checked);
}
//---------------------------------------------------------------------------
void lisemqt::on_checkInterception_toggled(bool checked)
{
    groupInterception->setEnabled(checked);
}
//---------------------------------------------------------------------------
void lisemqt::on_E_OFWaveType_currentIndexChanged(int index)
{
    groupFloodParams->setEnabled(index > 0);
    groupWaveUser->setEnabled(index > 0);
}
//---------------------------------------------------------------------------
void lisemqt::on_checkInfiltration_toggled(bool checked)
{
    groupInfiltration->setEnabled(checked);
}
//---------------------------------------------------------------------------
void lisemqt::on_checkIncludeChannel_toggled(bool checked)
{
    groupChannelParams->setEnabled(checked);
    checkMapChannels->setEnabled(checked);

   // checkMapNameModel(CHANNELMAPS, 0, checked);
}
//---------------------------------------------------------------------------
void lisemqt::on_checkDoErosion_toggled(bool checked)
{
    widgetErosion->setEnabled(checked);

    groupConservationSed->setEnabled(checked);

    ComboMinSpinBox2->setEnabled(checked);
    ComboMaxSpinBox2->setEnabled(checked);
    DisplayComboBox2->setEnabled(checked);
    checkBoxComboMaps2->setEnabled(checked);

    outputMapsSediment->setEnabled(checked);
    groupSedMapseriesout->setEnabled(checked);
    tabWidget_totout->setTabEnabled(1,checked);

    groupCalErosion->setEnabled(checked);

    label_soillosskgha->setEnabled(checked);
    label_soilloss->setEnabled(checked);
    label_SDR->setEnabled(checked);
    label_94->setEnabled(checked);
    label_105->setEnabled(checked);
    label_45->setEnabled(checked);

    checkMapNameModel(EROSIONMAPS, 0, checked);
}
//---------------------------------------------------------------------------
void lisemqt::on_checkInfrastructure_toggled(bool checked)
{
    widgetInfra->setEnabled(checked);
    transparencyHardSurface->setEnabled(checked);
    transparencyRoad->setEnabled(checked);
    checkMapHardSurface->setEnabled(checked);
    checkMapBuildings->setEnabled(checked);
    checkMapRoads->setEnabled(checked);
}
//---------------------------------------------------------------------------
void lisemqt::on_checkConservation_toggled(bool checked)
{
    groupMitigationWater->setEnabled(checked);
    groupConservationSed->setEnabled(checked);

}
//---------------------------------------------------------------------------
void lisemqt::on_toolButton_ShowRunfile_clicked()
{
    //qDebug() << E_runFileList->currentText();
    showTextfile(E_runFileList->currentText());
}
//---------------------------------------------------------------------------
void lisemqt::on_E_ETName_returnPressed()
{
    if (E_ETName->text() == "") {
        ETFileDir = "";
        ETFileName = "";
    }
}
//---------------------------------------------------------------------------
void lisemqt::on_E_RainsatName_returnPressed()
{
    if (E_RainsatName->text() == "") {
        RainSatFileDir = "";
        RainSatFileName = "";
    }

}
//---------------------------------------------------------------------------
void lisemqt::on_E_RainfallName_returnPressed()
{
    if (E_RainfallName->text() == "") {
        RainFileDir  = "";
        RainFileName = "";
    }
}
//---------------------------------------------------------------------------
void lisemqt::on_spinSoilLayers_valueChanged(int arg1)
{
    label_calKsat->setEnabled(true);
    E_CalibrateKsat->setEnabled(true);
    label_calKsat2->setEnabled(false);
    E_CalibrateKsat2->setEnabled(false);
    label_calKsat3->setEnabled(false);
    E_CalibrateKsat3->setEnabled(false);
    if (arg1 == 2) {
        label_calKsat2->setEnabled(true);
        E_CalibrateKsat2->setEnabled(true);
    }
    if (arg1 == 3) {
        label_calKsat2->setEnabled(true);
        E_CalibrateKsat2->setEnabled(true);
        label_calKsat3->setEnabled(true);
        E_CalibrateKsat3->setEnabled(true);
    }

}
//---------------------------------------------------------------------------
void lisemqt::on_E_InfiltrationMethod_currentIndexChanged(int index)
{
    groupBox_SwatreOptions->setEnabled(index == 0);
    groupBox_RichardsOptions->setEnabled(index == 3);
    groupAdvRichards->setEnabled(index == 0 || index == 3);
}
//---------------------------------------------------------------------------
void lisemqt::on_toolButton_clicked()
{
    checkforpatch = true;
    CheckVersion();
}

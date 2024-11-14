
/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011,2020  Victor Jetten
**  contact:
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
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
**  Authors: Victor Jetten, Bastian van de Bout
**  Developed in: MingW/Qt/
**  website, information and code: https://github.com/vjetten/openlisem
**
*************************************************************************/


#include "model.h"


//---------------------------------------------------------------------------
void TWorld::GetUserDischargeData(QString name)
{
    Q_LIST rl;
    QFile fff(name);
    QFileInfo fi(name);
    QString S;
    QStringList QRecs;
    QStringList SL;
    bool ok;
    int nrStations = 0;
    int nrSeries = 0;
    double time = 0.0;

    if (!fi.exists())
    {
        ErrorString = "User defined input discharge file not found: " + name;
        throw 1;
    }

    nrDischargeseries = 0;
    DischargeSeries.clear();
    currentDischargerow = 0;

    // read rainfall text file
    fff.open(QIODevice::ReadOnly | QIODevice::Text);
    while (!fff.atEnd())
    {
        S = fff.readLine();
        if (S.contains("\r\n"))
            S.remove(S.size()-2,2);
        if (S.contains("\n"))
            S.remove(S.size()-1,1);

        if (!S.trimmed().isEmpty())
            QRecs << S.trimmed();
    }
    fff.close();

    // check first if PCRaster graph format is present:
    int count = QRecs[1].toInt(&ok, 10); // nr of cols in file
    // header
    // second line is only an integer
    if (ok)
    {
        SL = QRecs[count+2].split(QRegularExpression("\\s+"));
        // check nr of columns in file
        if (count != SL.count()) {
            ErrorString = "Error: The nr of columns/stations in the discharge file does not equal the number on the second row.";
            throw 1;
        }

        //if the number of columns equals the integer then new format
        nrStations = count-1;
        // nr stations is count-1 for time as first column
    }

    // get station numbers from header, or fill in 1,2 ... n
    stationQID.clear();
    for (int i = 0; i < nrStations; i++) {
        SL = QRecs[i+3].split(QRegularExpression("\\s+"));
        int tmp = SL.last().toInt(&ok, 10);
        if (ok)
            stationQID << tmp;
        // if the header ends with a number, that number is the ID corresponding with the map
        else
            stationQID << i+1;
        // elese the position is the ID number

    }

    nrSeries = QRecs.size() - nrStations - 3;
    if (nrSeries <= 1)
    {
        ErrorString = "Discharge records <= 1, must at least have one interval with 2 rows: a begin and end time.";
        throw 1;
    }

    for(int r = 0; r < nrSeries; r++)
    {
        // initialize rainfall record structure
        rl.time = 0;
        rl.Qin.clear();
        rl.stationnr.clear();
        int r_ = r+nrStations+3;

        // split rainfall record row with whitespace
        QStringList SL = QRecs[r_].split(QRegularExpression("\\s+"), Qt::SkipEmptyParts);

        // read date time string and convert to time in seconds
        rl.time = getTimefromString(SL[0]);
        time = rl.time;

        // check if time is increasing with next row
        if (r+1 < nrSeries) {
            QStringList SL1 = QRecs[r_+1].split(QRegularExpression("\\s+"), Qt::SkipEmptyParts);
            int time1 = getTimefromString(SL1[0]);
            if (time1 < time) {
                ErrorString = QString("Time in discharge records is not increasing from row %1 to %2. Check your file!").arg(r_).arg(r_+1);
                throw 1;
            }
        }

        // rainfall values in this row
        for (int i = 1; i <= nrStations; i++)
        {
            bool ok = false;

            rl.Qin << SL[i].toDouble(&ok);
            if (!ok)
            {
                ErrorString = QString("Discharge records at time %1 has an unreadable value: %2.").arg(SL[0]).arg(SL[i]);
                throw 1;
            }
            rl.stationnr << stationID.at(i-1);
        }

        DischargeSeries << rl;
        dischargetime << rl.time;
    }

    nrDischargeseries = DischargeSeries.size();
}
//---------------------------------------------------------------------------
void TWorld::GetDischargeMapfromStations(double currenttime)
{
    bool same = false;

    // from time t to t+1 the rain is the rain of t

    // if time is outside records then use map with zeros
    if (currenttime < DischargeSeries[0].time || currenttime > DischargeSeries[nrRainfallseries-1].time) {
        Fill(*QuserIn, 0);
        return;
    }

    // where are we in the series
    int currentrow;
    auto it = std::lower_bound(dischargetime.begin(), dischargetime.end(), currenttime);
    if (it == dischargetime.begin())
        currentrow = 0;
    else
        currentrow = std::distance(dischargetime.begin(), it-1);

    if (currentrow == currentDischargerow && currentrow > 0)
        same = true;

    // get the next map from file
    if (!same) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_CHL {
            QuserIn->Drc = 0;
            for (int k = 0; k < stationQID.size(); k++) {
                if ((int) DischargeUserPoints->Drc == DischargeSeries[currentrow].stationnr.at(k)) {
                    QuserIn->Drc = DischargeSeries[currentrow].Qin[k]; // assuming m3/s
                   // qDebug() <<  QuserIn->Drc << currentrow << k << DischargeSeries[currentrow].stationnr.at(k);
                }
            }
        }}
    }

    currentDischargerow = currentrow;
}
//---------------------------------------------------------------------------
void TWorld::GetWHboundaryData(QString name)
{
    WH_LIST rl;
    QFile fff(name);
    QFileInfo fi(name);
    QString S;
    QStringList QRecs;
    QStringList SL;
    bool ok;
    int nrStations = 0;
    int nrSeries = 0;
    double time = 0.0;

    if (!fi.exists())
    {
        ErrorString = "User defined input water height file not found: " + name;
        throw 1;
    }

    nrWHseries = 0;
    currentWHrow = 0;

    // read WH text file
    fff.open(QIODevice::ReadOnly | QIODevice::Text);
    while (!fff.atEnd())
    {
        S = fff.readLine();
        if (S.contains("\r\n"))
            S.remove(S.size()-2,2);
        if (S.contains("\n"))
            S.remove(S.size()-1,1);

        if (!S.trimmed().isEmpty())
            QRecs << S.trimmed();
    }
    fff.close();

    // check first if PCRaster graph format is present: header, number of vars, columns equal vars
    int count = QRecs[1].toInt(&ok, 10); // nr of cols in file
    // header
    // second line is only an integer
    if (ok)
    {
        SL = QRecs[count+2].split(QRegularExpression("\\s+"));
        // check nr of columns in file
        if (count != SL.count()) {
            ErrorString = "Error: The nr of columns/stations in the boundary water level file does not equal the number on the second row.";
            throw 1;
        }

        //if the number of columns equals the integer then new format
        nrStations = count-1;
        // nr stations is count-1 for time as first column
    }

    nrSeries = QRecs.size() - 4;
    if (nrSeries <= 1)
    {
        ErrorString = "Water height records <= 1, must at least have one interval with 2 rows: a begin and end time.";
        throw 1;
    }

    for(int r = 0; r < nrSeries; r++)
    {
        // initialize rainfall record structure
        rl.time = 0;
        rl.WH = 0;
        int r_ = r+4;

        // split rainfall record row with whitespace
        QStringList SL = QRecs[r_].split(QRegularExpression("\\s+"), Qt::SkipEmptyParts);

        // read date time string and convert to time in seconds
        rl.time = getTimefromString(SL[0]);
        time = rl.time;

        // check if time is increasing with next row
        if (r+1 < nrSeries) {
            QStringList SL1 = QRecs[r_+1].split(QRegularExpression("\\s+"), Qt::SkipEmptyParts);
            int time1 = getTimefromString(SL1[0]);
            if (time1 < time) {
                ErrorString = QString("Time in boundary water level records is not increasing from row %1 to %2. Check your file!").arg(r_).arg(r_+1);
                throw 1;
            }
        }

        // rainfall values in this row
        bool ok = false;

        rl.WH = SL[1].toDouble(&ok);
        if (!ok) {
            ErrorString = QString("Boundary water level record at time %1 has an unreadable value: %2.").arg(SL[0]).arg(SL[1]);
            throw 1;
        }

        WHSeries << rl;
        WHtime << rl.time;
    }

    nrWHseries = WHSeries.size();
}
//---------------------------------------------------------------------------
void TWorld::GetWHboundaryMap(double currenttime)
{
    bool same = false;
    // from time t to t+1 the ET is the ET of t

    // if time is outside records then use map with zeros
    if (currenttime < WHSeries[0].time || currenttime > WHSeries[nrWHseries-1].time) {
        Fill(*WHbound, 0);
        return;
    }

    int currentrow;
    auto it = std::lower_bound(WHtime.begin(), WHtime.end(), currenttime);
    if (it == WHtime.begin())
        currentrow = 0;
    else
        currentrow = std::distance(WHtime.begin(), it-1);

    if (currentrow == currentWHrow && currentrow > 0)
        same = true;

    // get the next map from file
    if (!same) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            if (WHboundarea->Drc > 0)
                WHbound->Drc = WHSeries[currentrow].WH + WaveCalibration;
            else
                WHbound->Drc = 0;
        }}
    } //same

    currentWHrow = currentrow;
}

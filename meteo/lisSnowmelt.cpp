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
  \file lisSnowmelt.cpp
  \brief Get snowmelt adta and make a map. Snowmelt is like rainfall, melt intensities

functions: \n
- void TWorld::SnowmeltMap(void) \n
 */

#include <memory>
#include "io.h"
#include "model.h"
#include "operation.h"

//---------------------------------------------------------------------------
void TWorld::GetSnowmeltData(QString name)
{
    RAIN_LIST rl;
    QFile fff(name);
    QFileInfo fi(name);
    QString S;
    QStringList rainRecs;
    QStringList SL;
    bool ok;
    int nrStations = 0;
    int nrSeries = 0;
    int skiprows = 0;
    double time = 0.0;
    bool oldformat = true;

    if (!SwitchSnowmelt)
        return;

    if (!fi.exists())
    {
        ErrorString = "SnowMelt file not found: " + name;
        throw 1;
    }

    nrSnowmeltseries = 0;

    fff.open(QIODevice::ReadOnly | QIODevice::Text);

    while (!fff.atEnd())
    {
        S = fff.readLine();
        if (S.contains("\n"))
            S.remove(S.size()-1,1); // readLine also reads \n as a character on an empty line!
        if (!S.trimmed().isEmpty())
            rainRecs << S.trimmed();
    }
    fff.close();

    // check first if PCRaster graph format is present: header, number of vars, columns equal vars
    int count = rainRecs[1].toInt(&ok, 10);
    // header
    // second line is only an integer
    if (ok)
    {
        SL = rainRecs[count+2].split(QRegularExpression("\\s+"));
        if (count == SL.count())
            oldformat = false;
        //if the number of columns equals the integer then new format
        nrStations = count-1;
        // nr stations is count-1 for time as forst column
    }

    if (rainRecs[0].contains("RUU"))
        oldformat = true;

    if (oldformat)
    {
        QStringList SL = rainRecs[0].split(QRegularExpression("\\s+"));
        // get first line, white space character as split for header

        nrStations = SL[SL.size()-1].toInt(&ok, 10);
        // read nr stations from last value in old style header
        // failure gives 0
        SL = rainRecs[rainRecs.count()-1].split(QRegularExpression("\\s+"));
        oldformat = (nrStations == SL.count()-1);
    }

    //check if nr stations found equals nr columns-1, 1st column is time
    if (oldformat)
        skiprows = 1;
    else
        skiprows = 3;

    QList < int> tmp;
    tmp = countUnits(*SnowmeltZone);
    int nrmap = tmp.count();

    if (nrmap > nrStations)
    {
        ErrorString = QString("Number of stations in Snowmelt file (%1) < nr of rainfall zones in SNOWID map (%2)").arg(nrStations).arg(nrmap);
        throw 1;
    }
    nrSeries = rainRecs.size() - nrStations - skiprows;
    // count rainfall or snowmelt records

    if (nrSeries <= 1)
    {
        ErrorString = "Snowmelt records <= 1, must at least have one interval with 2 rows: a begin and end time.";
        throw 1;
    }

    for(int r = 0; r < nrSeries; r++)
    {
        // initialize rainfall record structure
        rl.time = 0;
        rl.intensity.clear();
        int r_ = r+nrStations+skiprows;

        // split snowmelt record row with whitespace
        QStringList SL = rainRecs[r_].split(QRegularExpression("\\s+"), Qt::SkipEmptyParts);

        // read date time string and convert to time in seconds
        rl.time = getTimefromString(SL[0]);
        time = rl.time;

        // check is time is increasing with next row
        if (r+1 < nrSeries) {
            QStringList SL1 = rainRecs[r_+1].split(QRegularExpression("\\s+"), Qt::SkipEmptyParts);
            int time1 = getTimefromString(SL1[0]);
            if (time1 < time) {
                ErrorString = QString("Time in evaporation records is not increasing from row %1 to %2. Check your file!").arg(r_).arg(r_+1);
                throw 1;
            }
        }
//        if (r == 0)
//            time = rl.time;


//        if (r > 0 && rl.time <= time)
//        {
//            ErrorString = QString("Snowmelt records at time %1 has unreadable value.").arg(rl.time);
//            throw 1;
//        }
//        else
//            time = rl.time;

        // check if record has characters, then filename assumed

        // record is a assumed to be a double
        for (int i = 1; i <= nrStations; i++)
        {
            bool ok = false;
            rl.intensity << SL[i].toDouble(&ok);
            if (!ok)
            {
                ErrorString = QString("Snowmel records at time %1 has unreadable value: %2.").arg(SL[0]).arg(SL[i]);
                throw 1;
            }
        }

        SnowmeltSeries << rl;
    }

    nrSeries++;
    rl.time = rl.time+1440; //?????????????

    for (int i = 1; i < nrStations; i++)
        rl.intensity << 0.0;

    SnowmeltSeries << rl;

    nrSnowmeltseries = nrSeries;
}
//---------------------------------------------------------------------------
void TWorld::GetSnowmeltMap(void)
{
    double currenttime = (time)/60;
    double tt = 0.001; //mm/day to m/day
    bool same = false;

    // from time t to t+1 the ET is the ET of t

    // where are we in the series

    // if time is outside records then use map with zeros
    if (currenttime < SnowmeltSeriesMaps[0].time || currenttime > SnowmeltSeriesMaps[nrSnowmeltseries-1].time) {
        DEBUG("run time outside ET records");
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            Snowmelt->Drc = 0;
        }}
        return;
    }

    // find current record
    int currentrow;
    auto it = std::lower_bound(snowmelttime.begin(), snowmelttime.end(), currenttime);
    if (it == snowmelttime.begin())
        currentrow = 0;
    else
        currentrow = std::distance(snowmelttime.begin(), it-1);

    if (currentrow == currentSnowmeltrow && currentrow > 0)
        same = true;

    // get the next map from file
    if (!same) {
        auto _M = std::unique_ptr<cTMap>(new cTMap(readRaster(SnowmeltSeriesMaps[currentrow].name)));

        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            if (pcr::isMV(_M->Drc)) {
                QString sr, sc;
                sr.setNum(r); sc.setNum(c);
                ErrorString = "Missing value at row="+sr+" and col="+sc+" in map: "+SnowmeltSeriesMaps[currentrow].name;
                throw 1;
            } else {
                Snowmelt->Drc = _M->Drc *_dt/tt;
            }
        }}
    }

    currentSnowmeltrow = currentrow;
}
//---------------------------------------------------------------------------

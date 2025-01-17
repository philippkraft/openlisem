/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 1992, 2003, 2016, 2024  Victor Jetten
**  contact:
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

#ifndef CsfMapH
#define CsfMapH
#include <QString>
#include "masked_raster.h"


/*!
    @brief      A cTMap contains all relevant information about a raster.
    @todo       The data member must be made private.

    cTMap instances contain raster data, projection information and a map name.
    I/O of cTMap instances is handles by functions defined in the io module.
*/
class cTMap
{

public:

    //! The actual raster.
    MaskedRaster<double> data;

                   cTMap               ()=default;

                   cTMap               (MaskedRaster<double>&& data,
                                        QString const& projection,
                                        QString const& mapName);

                   cTMap               (cTMap const& other)=delete;

                   cTMap               (cTMap&& other)=default;

                   ~cTMap              ()=default;

    cTMap&         operator=           (cTMap const& other)=delete;

    cTMap&         operator=           (cTMap&& other)=default;

    int            nrRows              () const;

    int            nrCols              () const;

    double         north               () const;

    double         west                () const;

    double         cellSize            () const;

    QString const& projection          () const;

    QString const& mapName             () const;

    void           setAllMV            ();

    void           MakeMap             (cTMap *dup,
                                        REAL8 value);

private:

    //! Projection string as WKT string. Possibly empty.
    QString        _projection;

    QString        _mapName;

};

#endif

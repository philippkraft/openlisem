#pragma once
#include <memory>
#include <QString>
#include "csf.h"


class cTMap;

//! Function to close a CSF MAP.
auto close_csf_map = [](MAP* map) { Mclose(map); };

//! Auto-ptr type for CSF MAPs.
using MapPtr = std::unique_ptr<MAP, decltype(close_csf_map)>;

cTMap              readRaster          (QString const& pathName);

void               writeRaster         (cTMap const& raster,
                                        QString const& Name,
                                        QString const& format="PCRaster");

void               WriteMapSeries      (cTMap const& raster,
                                        QString const& Dir,
                                        QString Name,
                                        int count,
                                        QString const& format="PCRaster");

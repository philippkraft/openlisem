/*---------------------------------------------------------------------------
project: openLISEM
author: Victor Jetten
licence: GNU General Public License (GPL)
Developed in: MingW/Qt/Eclipse
website, information and code: http://sourceforge.net/projects/lisem
---------------------------------------------------------------------------*/

/*
 * CsfMap basic fucntionality to read and write PCRaster CSF maps
 */

//---------------------------------------------------------------------------
#include <stdlib.h>
#include <math.h>
#include <qfile.h>
#include <qdir.h>
//---------------------------------------------------------------------------
#include "CsfMap.h"
#include "error.h"


//---------------------------------------------------------------------------
cTMap::cTMap()
{
    Data = NULL;
    Created = false;
    nrRows = 0;
    nrCols = 0;
    MapName = "";
    PathName = "";  
}
//---------------------------------------------------------------------------
cTMap::~cTMap()
{
    KillMap();
}
//---------------------------------------------------------------------------
void cTMap::KillMap()
{
    if (Data && nrRows > 0)
    {
       for (int r=0; r < nrRows; r++)
            delete[] Data[r];
       delete[] Data;
       Data = NULL;
    }
    Created = false;
}
//---------------------------------------------------------------------------
void cTMap::GetMapHeader(QString Name)
{
   MAP *m = Mopen(Name.toAscii().constData(), M_READ);
   MH = m->raster;
   projection = m->main.projection;

   // easier to separate these
   nrRows = (int)MH.nrRows;
   nrCols = (int)MH.nrCols;

   Mclose(m);
}
//---------------------------------------------------------------------------
void cTMap::CreateMap(QString Name)
{
     // if exist delete main map data structure and set Created to false
     if (Created)
        KillMap();

     // now get the header for nrrows and nrcols
     GetMapHeader(Name);

      Data = new REAL4*[nrRows];
      for(int r=0; r < nrRows; r++)
         Data[r] = new REAL4[nrCols];

     if (Data == NULL)
     {
    	 ErrorString = "Cannot create data structure for map: " + Name;
    	 throw 1;
     }

     Created = true;
}
//---------------------------------------------------------------------------
bool cTMap::LoadFromFile()
{
    MAP *m;
    QFileInfo fi(PathName);

    if (!fi.exists())
      return(false);

    // make map structure
    CreateMap(PathName);

    if (!Created)
       return(false);

    MapName = PathName;

    m = Mopen(MapName.toAscii().constData(), M_READ);

    if (!m)
       return(false);

    RuseAs(m, CR_REAL4); //RgetCellRepr(m));
    for(int r=0; r < nrRows; r++)
       RgetSomeCells(m, (UINT4)r*nrCols, (UINT4)nrCols, Data[r]);

    if (RgetCellRepr(m) == CR_REAL4)
       ResetMinMax();

    Mclose(m);
    return(true);
}
//---------------------------------------------------------------------------
void cTMap::_MakeMap(cTMap *dup, REAL4 value)
{
     if (dup == NULL)
        return;
     //MH = dup->MH;
     MH.minVal = value;
     MH.maxVal = value;
     nrRows = dup->nrRows;
     nrCols = dup->nrCols;
     MH.nrRows = nrRows;
     MH.nrCols = nrCols;
     MH.valueScale = dup->MH.valueScale;
     MH.cellRepr = dup->MH.cellRepr;
     MH.xUL = dup->MH.xUL;
	 MH.yUL = dup->MH.yUL;
     MH.cellSizeX  = dup->MH.cellSizeX;
     MH.cellSizeY  = dup->MH.cellSizeX;
     MH.angle = dup->MH.angle;
     projection = dup->projection;

     Data = new REAL4*[nrRows];
     for(int r=0; r < nrRows; r++)
        Data[r] = new REAL4[nrCols];

     if (Data == NULL)
        return;

     for(int r = 0; r < nrRows; r++)
	 SetMemMV(Data[r],nrCols,CR_REAL4);

     for(int r=0; r < nrRows; r++)
      for(int c=0; c < nrCols; c++)
      if (!IS_MV_REAL4(&dup->Data[r][c]))
      {
          Data[r][c] = value;
      }

      Created = true;
}
//---------------------------------------------------------------------------
void cTMap::ResetMinMax(void)
{
     REAL8 minv = 1e30, maxv = -1e30;

     if (!Created)
        return;

      for(int r=0; r < nrRows; r++)
        for(int c=0; c < nrCols; c++)
        if (!IS_MV_REAL4(&Data[r][c]))
        {
           if (maxv < Data[r][c]) maxv = Data[r][c];
           if (minv > Data[r][c]) minv = Data[r][c];
        }

     MH.minVal = minv;
     MH.maxVal = maxv;
}
//---------------------------------------------------------------------------
void cTMap::WriteMap(QString Name)
{
    MAP *out;
    long r;

    if (!Created)
        return;

    MH.cellRepr = CR_REAL4; 
    if (MH.cellRepr == CR_REAL4)    
       ResetMinMax();

    out = Rcreate(Name.toAscii().constData(),nrRows, nrCols, (CSF_CR)MH.cellRepr, VS_SCALAR,
                  (CSF_PT)projection, MH.xUL, MH.yUL, MH.angle, MH.cellSizeX);

    for(r=0; r < nrRows; r++)
    {
       if (RputRow(out, r, Data[r]) != (UINT4)nrCols)
       {
           ErrorString = "rputrow write error with" + Name;
    	   throw 1;
       }
    }

    Mclose(out);

}
//---------------------------------------------------------------------------
void cTMap::WriteMapSeries(QString Dir, QString Name, int count)
{
    QString path;
    QString nam, dig;
    QFileInfo fi(Name);

    if (Name.indexOf(".") < 0)
    {

    	nam = Name + "00000000";

    	nam.remove(8, 80);
    	dig.setNum(count);

    	if (count > 999)
    	{
    		nam.remove(8,1);
    		nam = nam + dig;
    		nam.insert(9,QString("."));
    	}
    	else
    		if (count > 99)
    			nam = nam + "." + dig;
    		else
    			if (count > 9)
    				nam = nam + "." + "0" + dig;
    			else
    				nam = nam + "." + "00" + dig;
    }
    else
       nam = Name;

    path = Dir + nam;
    WriteMap(path);
}



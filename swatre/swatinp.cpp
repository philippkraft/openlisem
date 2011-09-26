/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011  Victor Jetten
**  contact:
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
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
**  Author: Victor Jetten
**  Developed in: MingW/Qt/
**  website, information and code: http://lisem.sourceforge.net
**
*************************************************************************/
/*!
  \file swatinp.cpp
  \brief SWATRE: initialize and read profile data

  functions:
- void TWorld::InitializeProfile( void )
- int TWorld::ReadSwatreInput(QString fileName, QString tablePath) \n
- PROFILE * TWorld::ProfileNr(int profileNr) \n
- void  TWorld::FreeSwatreInfo(void) \n
- ZONE * TWorld::ReadNodeDefinition(FILE *f) \n
- PROFILE * TWorld::ReadProfileDefinition(FILE *f,ZONE *z,const char *tablePath) \n
- HORIZON * TWorld::ReadHorizon(const char *tablePath,	const char *tableName) \n
*/

#include "error.h"
#include "model.h"

#define LIST_INC	20
#define LUT_COLS  5
#define IND(r,c)  ((r)*LUT_COLS+(c))

// obsolete:
//static ZONE *ReadNodeDefinition(FILE *f);
//static PROFILE *ReadProfileDefinition(FILE *f, const ZONE *z,const char *tablePath);
//static HORIZON *ReadHorizon(const char *tablePath,	const char *tableName);
//static PROFILE **profileList = NULL;
//static int nrProfileList=0, sizeProfileList=0;
//static ZONE *zone=NULL;

/// array of pointers to horizons, NULL if not allocated
static HORIZON **horizonList = NULL;
static int nrHorizonList=0, sizeHorizonList=0;

//----------------------------------------------------------------------------------------------
void TWorld::InitializeProfile( void )
{
    profileList = NULL;
    nrProfileList=0;
    sizeProfileList=0;
    zone=NULL;

    horizonList = NULL;
    nrHorizonList=0;
    sizeHorizonList=0;
}
//----------------------------------------------------------------------------------------------
/// read and parse profile.inp
int TWorld::ReadSwatreInput(QString fileName, QString tablePath)
{
    FILE *f;
    ZONE *z;
    int  i, mmax;
    PROFILE **tmpList;
    // PROFILE.INP is opened here
    f = fopen(fileName.toAscii().constData(), "r");

    if (f == NULL)
    {
        Error(QString("SWATRE: Can't open profile definition file %1").arg(fileName));
        throw 1;
    }
    //All file name checking in main program

    // set profiles to NULL
    InitializeProfile();

    // read node distances and make structure
    z = ReadNodeDefinition(f);

    // check if list can hold new one
    do {

        if (nrProfileList == sizeProfileList)
        {
            int i = sizeProfileList;
            sizeProfileList += LIST_INC;
            profileList = (PROFILE **)realloc(profileList,sizeof(PROFILE *)*sizeProfileList);
            while (i < sizeProfileList)
                profileList[i++] = NULL;
        }

        profileList[nrProfileList] =
                ReadProfileDefinition(f,z,tablePath.toAscii().constData());

    } while (profileList[nrProfileList++] != NULL);
    // correct for eof-marker and test if something is read
    if ( --nrProfileList == 0)
        Error(QString("SWATRE: no profiles read from %1").arg(fileName));

    /* make profileList index match the profileId's */
    mmax = 0;
    for (i = 0 ; i < nrProfileList; i++)
        mmax = max(mmax, profileList[i]->profileId);
    mmax++;

    tmpList = (PROFILE **)malloc(mmax*sizeof(PROFILE *));
    for (i = 0 ; i < mmax; i++)
        tmpList[i] = NULL;

    for (i = 0 ; i < nrProfileList; i++)
        if (tmpList[profileList[i]->profileId] == NULL)
            tmpList[profileList[i]->profileId] = profileList[i];
        else
            Error(QString("SWATRE: profile with id '%!' declared more than once").arg(profileList[i]->profileId));

    free(profileList);

    profileList = tmpList;
    nrProfileList = mmax;
    sizeProfileList = mmax;

    /* PROFILE.INP is closed here */
    fclose(f);

    return (0);
}
//----------------------------------------------------------------------------------------------
/// return PROFILE or NULL if not found
PROFILE *TWorld::ProfileNr(int profileNr)
{
    if (profileNr < 0 || profileNr >= nrProfileList)
        return(NULL);
    return(profileList[profileNr]);
}
//----------------------------------------------------------------------------------------------
void  TWorld::FreeSwatreInfo(void)
{
    int i;

    /* currently, all profiles have the same zoning */
    if (zone == NULL)
        return;

    free(zone->dz);
    free(zone->z);
    free(zone->endComp);
    free(zone->disnod);
    free(zone);

    for(i=0; i < nrProfileList; i++)
        if (profileList[i] != NULL)
            free(profileList[i]);
    free(profileList);
    profileList = NULL;
    nrProfileList=sizeProfileList=0;

    for(i=0; i < nrHorizonList; i++)
    {
        free(horizonList[i]->name);
        FreeLut(horizonList[i]->lut);
        horizonList[i]->lut = NULL;
        free(horizonList[i]);
    }
    free(horizonList);
    horizonList = NULL;
    nrHorizonList=sizeHorizonList=0;
}
//----------------------------------------------------------------------------------------------
/** allocates ZONE structure,
 *   reads compartment ends:
 *   2.5 5 10 means dz[0] = dz[1] = 2.5, dz[2] = 5, etc.\n
 *   computes all parameters stored in ZONE -structure
 */
//static
ZONE * TWorld::ReadNodeDefinition(FILE *f)
{
    int  i;
    zone = (ZONE *)malloc(sizeof(ZONE));
    if ( fscanf(f,"%d",&(zone->nrNodes)) != 1 )
        Error(QString("SWATRE: Can't read number of nodes %1 from input file").arg(zone->nrNodes));
    if (zone->nrNodes < 1 )
        Error(QString("SWATRE: number of nodes %1 smaller than 1").arg(zone->nrNodes));
    if (zone->nrNodes > MAX_NODES)
        Error(QString("SWATRE: number of nodes %1 larger than %2").arg(zone->nrNodes).arg(MAX_NODES));

    zone->dz     = (double *)malloc(sizeof(double)*zone->nrNodes);
    zone->z      = (double *)malloc(sizeof(double)*zone->nrNodes);
    zone->disnod = (double *)malloc(sizeof(double)*(zone->nrNodes+1));
    zone->endComp= (double *)malloc(sizeof(double)*zone->nrNodes);

    for (i=0; i < zone->nrNodes; i++)
    {
        if ( fscanf(f,"%lf",&(zone->endComp[i])) != 1 )
            Error(QString("SWATRE: Can't read compartment end of node %1").arg(i+1));
        if (zone->endComp[i] <= 0)
            Error(QString("SWATRE: compartment end of node nr. %1 <= 0").arg(i+1));
        /* compute dz and make negative */
        zone->dz[i]= ( (i == 0) ? -zone->endComp[0] : (zone->endComp[i-1]-zone->endComp[i]));
        zone->z[i]= ( (i == 0) ? zone->dz[i]*0.5 : zone->z[i-1] + 0.5*(zone->dz[i-1]+zone->dz[i]));
        zone->disnod[i] = ( (i == 0) ? zone->z[i]: zone->z[i] - zone->z[i-1]);

        //qDebug() << i << "dz" << zone->dz[i] << "z" << zone->z[i] << "dist" << zone->disnod[i];
    }
    zone->disnod[zone->nrNodes] = 0.5 * zone->dz[zone->nrNodes-1];

    return(zone);
}
//----------------------------------------------------------------------------------------------
/**
 returns pointer to new profile or NULL if eof is encountered
 while reading first token of profile definition
 */
PROFILE * TWorld::ReadProfileDefinition(
        FILE *f,
        ZONE *z,         /* zone division this profile */
        const char *tablePath) /* pathName ended with a '/' */
{
    char tableName[256];
    int  i;
    double endHor;
    PROFILE *p;
    HORIZON *h;
    /* profile has a pointer to LUT */
    p = (PROFILE *)malloc(sizeof(PROFILE));

    if ( fscanf(f,"%d",&(p->profileId)) != 1 )
    {
        if (feof(f))
        {
            free(p);
            return NULL;
        }
        Error(QString("SWATRE: read error: can't read profile id %1").arg(p->profileId));
    }
    if (p->profileId < 0)
        Error(QString("SWATRE: profile id smaller that 0: %1").arg(p->profileId));

    p->horizon = (const HORIZON **)malloc(sizeof(HORIZON *)*z->nrNodes);
    p->zone = z;

    i = 0;
    while (i != z->nrNodes)
    {
        if ( fscanf(f,"%s",tableName) != 1 )
            Error(QString("SWATRE: Can't read a LUT for profile nr %1 node nr %2 and up").arg(p->profileId).arg(i+1));
        if ( fscanf(f,"%lf", &endHor) != 1 )
            Error(QString("SWATRE: Can't read end of horizon for profile nr %1").arg(p->profileId));

        h = ReadHorizon(tablePath, tableName);
        // copy horizon info to all nodes of this horizon

        while (i < z->nrNodes && z->endComp[i] <= endHor )
            p->horizon[i++] = h;
        if (z->endComp[i-1] != endHor)
            Error(QString("SWATRE: No compartment ends on depth '%1' (found in profile nr %2 for horizon %3)")
                  .arg(endHor).arg(p->profileId).arg(tableName));
    }
    return(p);
}
//----------------------------------------------------------------------------------------------
/// copy horizon info to all nodes of this horizon
HORIZON * TWorld::ReadHorizon(const char *tablePath,	const char *tableName)
{
    HORIZON	*h;
    char fileName[256];
    double *t, *lutCont;
    int i, nrRowsa=0;

    /* look if it's already loaded */
    for( i= 0; i < nrHorizonList; i++)
        if (!strcmp(tableName, horizonList[i]->name))
            return(horizonList[i]);

    /* if not then add one */
    /* check for space in list */
    if (nrHorizonList == sizeHorizonList)
    {
        sizeHorizonList += LIST_INC;
        horizonList = (HORIZON **)realloc(horizonList, sizeof(HORIZON *)*sizeHorizonList);
    }

    h = (HORIZON *)malloc(sizeof(HORIZON));
    horizonList[nrHorizonList++] = h;
    strcat(strcpy(fileName, tablePath), tableName);


    /* hook up table to t */
    t = ReadSoilTable(fileName, &nrRowsa);

    lutCont = (double *)malloc(sizeof(double)*NR_COL*(nrRowsa+2));
    for(i=0; i < nrRowsa; i++)
    {
        lutCont[IND(i,THETA_COL)] =  t[i*3+THETA_COL];
        lutCont[IND(i,H_COL)]     =  t[i*3+H_COL];
        lutCont[IND(i,K_COL)]     =  t[i*3+K_COL];
    }
    for(i=0; i < (nrRowsa-1); i++)
    {
        lutCont[IND(i,DMCH_COL)] = 0.5 *
                (lutCont[IND(i+1,H_COL)] + lutCont[IND(i,H_COL)]);

        /* VJ : 0.01 is 1% humidity, ofwel
 *	  gevaarlijk om 0.01 te gebruiken want dit ligt aan de tabel
 *
 *		lutCont[IND(i,DMCC_COL)] = 0.01 /
 *					(lutCont[IND(i+1,H_COL)] - lutCont[IND(i,H_COL)]);
 * beter:
 */
        lutCont[IND(i,DMCC_COL)] =
                (lutCont[IND(i+1,THETA_COL)] - lutCont[IND(i,THETA_COL)])/
                (lutCont[IND(i+1,H_COL)] - lutCont[IND(i,H_COL)]);
    }
    lutCont[IND(nrRowsa-1,DMCH_COL)] = 0;
    lutCont[IND(nrRowsa-1,DMCC_COL)] = lutCont[IND(nrRowsa-2,DMCC_COL)] ;

    free(t);

    h->name = strcpy((char *)malloc(strlen(tableName)+1), tableName);
    h->lut = CreateLutFromContents(lutCont, true, nrRowsa, LUT_COLS);


    return(h);
}
//----------------------------------------------------------------------------------------------

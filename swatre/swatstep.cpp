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

#include <algorithm>
#include "model.h"


//--------------------------------------------------------------------------------
/// new matrix head: Gaussian elimination and backsubstitution
/**
tri diagonalmatrix to solve differential equations to get the new matrix potential
and the new moisture content theta
method is by gaussian elimination and backsubstitution, 2 times instead of iteration
matrix shape:
\code
|b0 c0 .  .| |h0| |F0|
|a1 .. cn-1|*|  |=|  |
|   an bn  | |hn| |Fn|
\endcode
F (thomf) is new moisture content theta
thoma ?
thomc has unit cm
thomb has unit dtheta/dh (like dimoca)

*/
//--------------------------------------------------------------------------------
/**
 * @brief
 *
 * @param prevDt
 * @param hLast
 * @param h
 * @param nrNodes
 * @param precParam
 * @param dtMin
 * @param dtMax
 * @return double
 */
double  TWorld::NewTimeStep(
        double 		 prevDt,
        QVector <double> hLast,
        QVector <double> h,
        int    		 nrNodes,
        double 		 dtMin)
{
    double precParam = 5.0;
    // note "5" is a precision factor dewtermining next timestep, set to 5 in old lisem
    int i;
    double dt = _dt;
    double accur1 = 0.3 - 0.02 * precParam;
    double accur2 = 0.03 - 0.002 * precParam;

    for(i=0; i < nrNodes; i++)
    {
        double mdih = accur1 + accur2 * std::max(1.0, fabs(h[i]));
        double dih  = fabs(h[i] - hLast[i]);
        // if difference is small
        // dih = e.g. 10 and h = -200 then mdih = 200*0.01 + 0.1 = 2.1
        // mdih/dih = 2.1/10 =0.21

        if (dih > 0.10)
            dt = std::min(dt, prevDt*mdih/dih);
    }
    return (std::max(dt, dtMin));
}
//--------------------------------------------------------------------------------
// Units are:
// Z and H in cm; table units K in cm/day converted to cm/sec, lisem time in seconds

void TWorld::ComputeForPixel(PIXEL_INFO *pixel, SOIL_MODEL *s, double drainfraction)
{
    NODE_ARRAY dimoca, kavg, k, thetaPrev, theta;
    const PROFILE *p = pixel->profile;

    int nN = p->zone->nrNodes;//NrNodes(p);
    QVector <double> disnod = p->zone->disnod;
    QVector <double> dz = p->zone->dz;
    double dt = _dt/5;
    double pond = pixel->wh;
    double elapsedTime = 0;
    double influx = 0;
    double drainout = 0;
    double percolation = 0;
    double Theta = 0;
    int tnode = pixel->tilenode;
    double impfrac = pixel->impfrac;

    QVector <double> h;
    QVector <double> hPrev;

    h = pixel->h; // just copy the vector

    while (elapsedTime < _dt) {

        bool isPonded, fltsat;    // flag if ponded or if profile fully saturated
        double qmax, qtop, qbot, ThetaSat;  // fluxes at top and bottom, max theta
        double qdrain; // tile drainage

        // get nodal values of theta, K, dif moist cap
        for (int j = 0; j < nN; j++) {
            k[j] = FindValue(h[j], p->horizon[j], H_COL, K_COL);
            // K in cm/sec from h, ksatcal filled with values for ksat1,2,3
            dimoca[j] = FindValue(h[j], p->horizon[j], DMCH_COL, DMCC_COL);
                    //DmcNode(h[j], p->horizon[j],  true); // true is more detailed method, false is DMCH directly from H
            // differential moisture capacity d(theta)/d(h), tangent moisture retention curve
            theta[j] = FindValue(h[j], p->horizon[j], H_COL, THETA_COL);
            // moisture content from H
        }

        // per pixel correction of Ks and Pore for org mat and density
        // near saturated so for h > -1 cm, and only for topsoil, assumed to be 30 cm
        if (SwitchOMCorrection) {
            int j = 0;
            while (p->zone->endComp[j] <= 30 && h[j] > -1.0) {
                k[j] = pixel->corrKsOA*k[j] + pixel->corrKsOB;
                theta[j] = pixel->corrPOA*k[j] + pixel->corrPOB;
                j++;
            }
        }
        if (SwitchDensCorrection) {
            int j = 0;
            while (p->zone->endComp[j] <= 30 && h[j] > -1.0) {
                k[j] = pixel->corrKsDA*k[j] + pixel->corrKsDB;
                theta[j] = pixel->corrPDA*k[j] + pixel->corrPDB;
                j++;
            }
        }
        // do calibration  after dens and OM calculations
        for (int j = 0; j < nN; j++) {
             k[j] *= p->KsatCal[j];
             //if (h[j] > -1.0) k[j] *= p->KsatCal[j]; //??? only for ksat
        }

        k[0] *= (1.0-impfrac);
        // adjust k[0] for roads and houses, impermeable fraction

        Theta = (theta[0]+theta[1])/2;
        // average K for 1st to n-1 node, top node is done below\

        // original swatre artithmetric mean
       for(int j = 1; j < nN; j++) {
           kavg[j] = (k[j]+k[j-1])/2.0;
       }
       //for(int j = 1; j < nN; j++){ kavg[j] = sqrt(k[j]*k[j-1]);}
        //--- boundary conditions ---

        //----- TOP -----
        // check if ponded: 1st compare fluxes, 2nd compare store
        // top flux is water/timestep (cm/sec)
        qtop = -pond/dt;

        //----- BOTTOM -----
        // bottom is 0 or copy of flux of last 2 layers
        if (SwitchImpermeable)
            qbot = 0;
        else
            qbot = kavg[nN-1]*(h[nN-1]-h[nN-2])/disnod[nN-1] - kavg[nN-1];

        // 1st check flux against max flux
        double Ksat = FindValue(0, p->horizon[0], H_COL, K_COL)*p->KsatCal[0]*(1.0-impfrac);
        if (SwitchOMCorrection) {
            Ksat = pixel->corrKsOA*Ksat+pixel->corrKsOB;
        }
        if (SwitchDensCorrection) {
            Ksat = pixel->corrKsDA*Ksat+pixel->corrKsDB;
        }

        kavg[0]= sqrt(Ksat * k[0]);
        // max possible always geometric mean

        // geometric avg of ksat and k[0] => is used for max possible

        qmax = kavg[0]*(pond-h[0]) / disnod[0] - kavg[0];
        // this is good but:
        // q = -k(dh/dz+1) = -kdh/dz-k
        // dh is in the dircetion of flow: pond - ho = pond - (-100) for instance so suction and pressure added downward
        // so it should be : qmax = - kavg[0]*((pond-h[0]) / disnod[0] +1);
        // so theoretically a minus sign, but that does not work, more infiltration when soil is wetter!
        // qtop is the same below!!!

        //qmax = std::min(qmax, 0.0);
        // maximum possible flux, compare to real top flux available
        isPonded = (qtop < qmax);
        // if more flux then max possible flag ponded is true

        //2nd check: isPonded layer depth against storage
        if (!isPonded) {
            // calculate available space in profile in cm: (pore-theta)*dz
            double space = 0;
            for(int i = 0; i < nN && h[i] < 0 && space > pond; i++) {
                ThetaSat = FindValue(0, p->horizon[i], H_COL, THETA_COL);
                space += (ThetaSat - theta[i]) * (-dz[i]);
            }
            isPonded = pond > space;
        }

        // check if profile is completely saturated (flstsat)
        fltsat = true;
        for (int i = nN-1; i >= 0; i--) {
            if (h[i] < 0) {
                fltsat = false;
                break;
            }
        }
        if (fltsat && (qtop <= qbot))
            fltsat = false;

        //std::memcpy(hPrev, h, nN * sizeof(double));
        std::memcpy(thetaPrev, theta, nN * sizeof(double));
        hPrev = h;

        //HeadCalc(h, p , &isPonded, fltsat,thetaPrev, hPrev, kavg, dimoca,dt, pond, qtop, qbot);

        // === 3 diagonal newton rapson with Gaussian elimination and back substitution
        NODE_ARRAY thoma, thomb, thomc, thomf, beta;

        //TODO add ET sink term to swatre

        /* First node : 0 (include boundary cond. qtop or pond) */
        if ( isPonded || fltsat) {
            // h at soil surface prescribed, ponding
            thomc[0] = -dt * kavg[1] / dz[0] / disnod[1];
            thomb[0] = -thomc[0] + dimoca[0] + dt*kavg[0]/disnod[0]/dz[0];
            thomf[0] = dimoca[0]*h[0] + dt/(-dz[0]) * (kavg[0] - kavg[1]) +
                    dt*kavg[0]*pond/disnod[0]/dz[0];
        } else {
            //  q at soil surface prescribed, qtop = rainfall
            isPonded = false;
            thomc[0] = -dt * kavg[1] / dz[0] / disnod[1];
            thomb[0] = -thomc[0] + dimoca[0];
            thomf[0] = dimoca[0]*h[0] + dt/(-dz[0]) * (-qtop - kavg[1]); //(- qtop - kavg[1]);
        }

        /* Intermediate nodes: i = 1 to n-2 */
        for (int i = 1; i < nN-1; i++)
        {
            thoma[i] = -dt*kavg[i]/dz[i]/disnod[i];
            thomc[i] = -dt*kavg[i+1]/dz[i]/disnod[i+1];
            thomb[i] = -thoma[i] - thomc[i] + dimoca[i];
            thomf[i] = dimoca[i]*h[i] + dt/(-dz[i])*(kavg[i]-kavg[i+1]);
        }

        // last node : nN-1 (include boundary cond. qbot)
        thoma[nN-1] = -dt*kavg[nN-1]/dz[nN-1]/disnod[nN-1];
        thomb[nN-1] = -thoma[nN-1] + dimoca[nN-1];
        thomf[nN-1] = dimoca[nN-1]*h[nN-1] + dt/(-dz[nN-1])*(kavg[nN-1]+qbot);

        // Gaussian elimination and backsubstitution h - first time
        double alpha = thomb[0];
        h[0] = thomf[0] / alpha;
        for (int i = 1; i < nN; i++) {
            beta[i] = thomc[i-1] / alpha;
            alpha = thomb[i] - thoma[i] * beta[i];
            h[i] = (thomf[i] - thoma[i] * h[i-1]) / alpha;
        }
        for (int i = (nN-2); i >= 0; i--)
            h[i] -= beta[i+1] * h[i+1];

        // correct tridiagonal matrix
        for (int i = 0; i < nN; i++) {
            double thetaNew = FindValue(h[i], p->horizon[i], H_COL, THETA_COL);
            double dimocaNew = DmcNode(h[i], p->horizon[i], true);
            thomb[i] = thomb[i] - dimoca[i] + dimocaNew;
            thomf[i] = thomf[i] - dimoca[i]*hPrev[i] + dimocaNew*h[i]
                    - thetaNew + thetaPrev[i];
        }

        // Gaussian elimination and backsubstitution h - second time
        alpha = thomb[0];
        h[0] = thomf[0] / alpha;
        for (int i = 1; i < nN; i++) {
            beta[i] = thomc[i-1] / alpha;
            alpha = thomb[i] - thoma[i] * beta[i];
            h[i] = (thomf[i] - thoma[i] * h[i-1]) / alpha;
        }

        for (int i = (nN-2); i >= 0; i--)
            h[i] -= beta[i+1] * h[i+1];

        for (int j = 0; j < nN; j++) {
            theta[j] = FindValue(h[j], p->horizon[j], H_COL, THETA_COL);
        }

        // determine new boundary fluxes

        if (SwitchImpermeable)
            qbot = 0;
        else
            qbot = -kavg[nN-1]*(h[nN-1]-h[nN-2])/disnod[nN-1] - kavg[nN-1];
        //qbot = kavg[n-1]*(h[n-1]-h[n-2])/disnod[n-1] - kavg[n-1];
        // new qbot is actually not use but may come in handy later
        percolation += qbot*dt;

        if ( isPonded || fltsat)
             qtop = -kavg[0] * ((h[0] - pond)/disnod[0] + 1);
          //  qtop = -kavg[0]*(pond - h[0]) / disnod[0] - kavg[0];
        //qtop = std::min(0.0,qtop);

        pond += qtop*dt;       // decrease pond with top flux
        if (pond < POND_EPS)   // why?
            pond = 0;
        influx += qmax*dt;
        // add max infil to influx (negative), to get potential infil  WHY???

        //--- calculate tile drain ---//
        //TODO: CHECK THIS
        if (SwitchIncludeTile && tnode > 0) //VJ 110825 tnode = -1 if cell has no drainage
        {
            //options:
            qdrain =  k[tnode];
            // drainage is cond of the node in cm/sec
            double water = theta[tnode] * -disnod[tnode] * drainfraction;
            // total amonut of water available to drain in this node (cm)
            // note: distnode has a negative value
            qdrain = std::min(qdrain, water/dt);
            // cannot have more drainage than water available

            theta[tnode] = std::max(0.001, theta[tnode] - (qdrain*dt)/disnod[tnode]* drainfraction);
            // adjust theta with drainage removed

            h[tnode] = FindValue(theta[tnode], p->horizon[tnode], THETA_COL, H_COL );
            hPrev[tnode] = h[tnode];
            // new h from theta

            drainout += qdrain*dt;
            // add for all swatre timestps, in cm
        }

        elapsedTime += dt;
        /* estimate new dt within lisemtimestep */
        dt = NewTimeStep(dt, hPrev, h, nN, s->minDt);

        if (elapsedTime+dt+TIME_EPS >= _dt)
            dt = _dt - elapsedTime;

    } // elapsedTime < lisemTimeStep

    //put new h back into h
    pixel->h = h;
    // for (int i = 0; i < nN; i++) {
    //     pixel->h[i] = h[i];
    // }

    // these variables can all be direcvtly saved to the maps, inflated pixel structure
    pixel->wh = pond;
    pixel->tiledrain = drainout;
    pixel->infil = influx;
    pixel->percolation = -percolation; // in cm
    pixel->theta = Theta;

}
//--------------------------------------------------------------------------------
// units in SWATRE are cm and K cm/sec
void TWorld::SwatreStep(long i_, int r, int c, SOIL_MODEL *s, cTMap *_WH, cTMap *_drain, cTMap *_theta)
{   
    double drainfraction = 0;
    QString dig;

    //showr = r; // for debug
    //showc = c;

    s->pixel[i_].wh = _WH->Drc*100;    // WH is in m, convert to cm
    s->pixel[i_].infil = 0;
    s->pixel[i_].tiledrain = 0;

    if (SwitchIncludeTile)
        drainfraction = TileWidth->Drc/_dx;

    ComputeForPixel(&s->pixel[i_], s, drainfraction);
    // estimate new h and theta at the end of dt

    // ?????? this only makes the output name
    if(SwitchDumphead) {
        if(s->pixel[i_].dumpHid > 0) {
            for (int i = 0; i < zone->nrNodes; i++) {
                QString name = QString("SWH%1").arg(runstep,2, 10, QLatin1Char('0'));
                dig = QString("%1").arg(i, 12-name.length(), 10, QLatin1Char('0'));
                name=name+dig;
                name.insert(8, ".");
                //qDebug() << name << dig;
            }
        }
    }

    _WH->Drc = s->pixel[i_].wh*0.01; // cm to m
    _theta->Drc = s->pixel[i_].theta; // for pesticides ?
    Perc->Drc = s->pixel[i_].percolation*0.01;

    if (SwitchIncludeTile)
        _drain->Drc = s->pixel[i_].tiledrain*0.01;  // in m
    // drained water from the soil, already accounts for drainwidth versus i_l width
}
//--------------------------------------------------------------------------------

//OBSOLETE
void TWorld::HeadCalc(double *h, const PROFILE *p , bool *isPonded, bool fltsat,
              const double *thetaPrev, const double *hPrev, const double *kavg, const double *dimoca,
              double dt, double pond, double qtop, double qbot)
{
    /*
    int nN = NrNodes(p);
    //const double *dz = dz, *disnod = disnod;
    // dz is layer thickness, distnode is distance between centre of layers
    int i;
    NODE_ARRAY thoma, thomb, thomc, thomf, beta;
    double alpha;

    // First node : 0 (include boundary cond. qtop or pond)
    if ( *isPonded || fltsat ) {
        // h at soil surface prescribed, ponding
        thomc[0] = -dt * kavg[1] / dz[0] / disnod[1];
        thomb[0] = -thomc[0] + dimoca[0] + dt*kavg[0]/disnod[0]/dz[0];
        thomf[0] = dimoca[0]*h[0] + dt/(-dz[0]) * (kavg[0] - kavg[1]) +
                dt*kavg[0]*pond/disnod[0]/dz[0];
    } else {
        //  q at soil surface prescribed, qtop = rainfall
        *isPonded = false;
        thomc[0] = -dt * kavg[1] / dz[0] / disnod[1];
        thomb[0] = -thomc[0] + dimoca[0];
        thomf[0] = dimoca[0]*h[0] + dt/(-dz[0]) * (-qtop - kavg[1]); //(- qtop - kavg[1]);
    }

    // Intermediate nodes: i = 1 to n-2
    for (i = 1; i < nN-1; i++)
    {
        thoma[i] = -dt*kavg[i]/dz[i]/disnod[i];
        thomc[i] = -dt*kavg[i+1]/dz[i]/disnod[i+1];
        thomb[i] = -thoma[i] - thomc[i] + dimoca[i];
        thomf[i] = dimoca[i]*h[i] + dt/(-dz[i])*(kavg[i]-kavg[i+1]);
    }

    // last node : nN-1 (include boundary cond. qbot)
    thoma[nN-1] = -dt*kavg[nN-1]/dz[nN-1]/disnod[nN-1];
    thomb[nN-1] = -thoma[nN-1] + dimoca[nN-1];
    thomf[nN-1] = dimoca[nN-1]*h[nN-1] + dt/(-dz[nN-1])*(kavg[nN-1]+qbot);

    // Gaussian elimination and backsubstitution h - first time
    alpha = thomb[0];
    h[0] = thomf[0] / alpha;
    for (i = 1; i < nN; i++) {
        beta[i] = thomc[i-1] / alpha;
        alpha = thomb[i] - thoma[i] * beta[i];
        h[i] = (thomf[i] - thoma[i] * h[i-1]) / alpha;
    }
    for (i = (nN-2); i >= 0; i--)
        h[i] -= beta[i+1] * h[i+1];

    // correct tridiagonal matrix
    for (i = 0; i < nN; i++) {
        // double theta = TheNode(h[i], p->horizon[j]);
        // double dimocaNew = DmcNode(h[i], p->horizon[j]);
        double theta = FindNode(h[i], p->horizon[j], THETA_COL);
        double dimocaNew = FindNode(h[i], p->horizon[j], DMCC_COL);
        thomb[i] = thomb[i] - dimoca[i] + dimocaNew;
        thomf[i] = thomf[i] - dimoca[i]*hPrev[i] + dimocaNew*h[i]
                - theta + thetaPrev[i];
    }

    // Gaussian elimination and backsubstitution h - second time
    alpha = thomb[0];
    h[0] = thomf[0] / alpha;
    for (i = 1; i < nN; i++) {
        beta[i] = thomc[i-1] / alpha;
        alpha = thomb[i] - thoma[i] * beta[i];
        h[i] = (thomf[i] - thoma[i] * h[i-1]) / alpha;
    }

    for (i = (nN-2); i >= 0; i--)
        h[i] -= beta[i+1] * h[i+1];
*/
}

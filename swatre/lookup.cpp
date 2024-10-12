
/*!
  \file lookup.cpp
  \brief SWATRE: computes Theta from head, K from head or Diff Moist Cap from head.

  functions:\n
- double HNode(double theta, const  HORIZON *hor)\n
- double TheNode(double head, const  HORIZON *hor)\n
- double HcoNode(double head, const HORIZON *hor, double calib, double SEC)\n
- double DmcNode(double head, const  HORIZON *hor) \n

*/

//#include <algorithm>
#include "model.h"


//-----------------------------------------------------------------------------------
/// head from theta
double TWorld::HNode(
        double theta,           // current theta value of this node
        const  HORIZON *hor)    // parameters of horizon this node belongs to
{
    LUT *l = hor->lut;

    auto it = std::lower_bound(l->hydro[THETA_COL].begin(), l->hydro[THETA_COL].end(), theta);
    if (it == l->hydro[THETA_COL].begin()) {
        return(l->hydro[H_COL][0]);
    } else if (it == l->hydro[THETA_COL].end()) {
        return(l->hydro[H_COL][l->nrRows-1]);
    } else {
        double lH = *(it - 1);
        double uH = *it;
        double f = (theta-lH)/(uH-lH);

        int lowerIndex = std::distance(l->hydro[THETA_COL].begin(), it - 1); // Index of lower bound
        int upperIndex = std::distance(l->hydro[THETA_COL].begin(), it);     // Index of upper bound

        double lTh = l->hydro[H_COL][lowerIndex];
        double uTh = l->hydro[H_COL][upperIndex];
        return (lTh + f*(uTh-lTh));
    }


//    return LUT_LinIntPol(hor->lut,H_COL, theta,THETA_COL);
}
//-----------------------------------------------------------------------------------
/// theta from head
double TWorld::TheNode(
        double head,           // current head value of this node
        const  HORIZON *hor)   // parameters of horizon this node belongs to
{
    LUT *l = hor->lut;
    if (head >= 0) {
        return l->hydro[THETA_COL][l->nrRows-1];
    }
    auto it = std::lower_bound(l->hydro[H_COL].begin(), l->hydro[H_COL].end(), head);
    if (it == l->hydro[H_COL].begin()) {
        return(l->hydro[THETA_COL][0]);
    } else if (it == l->hydro[H_COL].end()) {
        return(l->hydro[THETA_COL][l->nrRows-1]);
    } else {
        double lH = *(it - 1);
        double uH = *it;
        double f = (head-lH)/(uH-lH);

        int lowerIndex = std::distance(l->hydro[H_COL].begin(), it - 1); // Index of lower bound
        int upperIndex = std::distance(l->hydro[H_COL].begin(), it);     // Index of upper bound

        double lTh = l->hydro[THETA_COL][lowerIndex];
        double uTh = l->hydro[THETA_COL][upperIndex];
        return (lTh + f*(uTh-lTh));
    }

    //head = std::min( head, -1e-10);  //max or min ???? org was max! but head is < 0 !
    //VJ 110825 better to comment out this line, not useful
    // if (head >= -1.0E-2)
    //    return LUT_Highest(hor->lut, THETA_COL);
    // return LUT_LinIntPol(hor->lut, THETA_COL, head, H_COL);
}
//-----------------------------------------------------------------------------------
/// hydraulic conductivity from head
double TWorld::HcoNode(double head,const HORIZON *hor,double calib)
{
    LUT *l = hor->lut;

    if (SHOWDEBUG) {
    if (head >= 0) {
     //   qDebug() << "ksat" << l->hydro[K_COL][l->nrRows-1];
        return l->hydro[K_COL][l->nrRows-1];
    }
    }

    auto it = std::lower_bound(l->hydro[H_COL].begin(), l->hydro[H_COL].end(), head);
    if (it == l->hydro[H_COL].begin()) {
        return(l->hydro[K_COL][0]);
    } else if (it == l->hydro[H_COL].end()) {
        return(l->hydro[K_COL][l->nrRows-1]);
    } else {
        double lH = *(it - 1);
        double uH = *it;
        double f = (head-lH)/(uH-lH);
        int lowerIndex = std::distance(l->hydro[H_COL].begin(), it - 1); // Index of lower bound
        int upperIndex = std::distance(l->hydro[H_COL].begin(), it);     // Index of upper bound

        double lK = l->hydro[K_COL][lowerIndex];
        double uK = l->hydro[K_COL][upperIndex];
    //    qDebug() << head << lH << uH << lK << uK << (lK+f*(uK-lK));
        return (lK+f*(uK-lK));
    }
}
//-----------------------------------------------------------------------------------
/// Differential Moisture Capacity from head
double TWorld::DmcNode(
        double head,           // current head value of this node
        const  HORIZON *hor)   // parameters of horizon this node belongs to
{
    LUT *l = hor->lut;

    if (head >= 0) {
        return l->hydro[DMCC_COL][l->nrRows-1];
    }

    auto it = std::lower_bound(l->hydro[H_COL].begin(), l->hydro[H_COL].end(), head);

    if (it == l->hydro[H_COL].begin()) {
        return(l->hydro[DMCC_COL][0]);
    } else if (it == l->hydro[H_COL].end()) {
        return(l->hydro[DMCC_COL][l->nrRows-1]);
    } else {
        double lH = *(it - 1);
        double uH = *it;
        double f = (head-lH)/(uH-lH);

        int lowerIndex = std::distance(l->hydro[H_COL].begin(), it - 1); // Index of lower bound
        int upperIndex = std::distance(l->hydro[H_COL].begin(), it);     // Index of upper bound

        //double lcH = l->hydro[H_COL][lowerIndex];
        //double ucH = l->hydro[H_COL][upperIndex];
        double lC = l->hydro[DMCC_COL][lowerIndex];
        double uC = l->hydro[DMCC_COL][upperIndex];

        return (lC+f*(uC-lC)); //(lC + (head-lcH)*(uC-lC)/(ucH-lcH));
    }


    // if (it == l->hydro[DMCH_COL].begin()) {
    //     return(l->hydro[DMCC_COL][0]);
    // } else if (it == l->hydro[DMCH_COL].end()) {
    //     return(l->hydro[DMCC_COL][l->nrRows-1]);
    // } else {
    //     double lH = *(it - 1);
    //     double uH = *it;
    //     double f = (head-lH)/(uH-lH);
    //     int lowerIndex = std::distance(l->hydro[DMCH_COL].begin(), it - 1); // Index of lower bound
    //     int upperIndex = std::distance(l->hydro[DMCH_COL].begin(), it);     // Index of upper bound

    //     double lcH = l->hydro[DMCH_COL][lowerIndex];
    //     double ucH = l->hydro[DMCH_COL][upperIndex];
    //     double lC = l->hydro[DMCC_COL][lowerIndex];
    //     double uC = l->hydro[DMCC_COL][upperIndex];

    //     return (lC + (head-lcH)*(uC-lC)/(ucH-lcH));
    // }

/*
    if (head >= -1.0E-2)
       return LUT_Highest(hor->lut, DMCC_COL);
    //       return LUT_LinIntPol(hor->lut, DMCC_COL, head, DMCH_COL);

    l = hor->lut;
    i = LUT_Index_LE(l, head, DMCH_COL);
    i = std::min(LUT_nrRows(l)-2, i);
    i = std::max(i, 0);

    return LUT_ValueAt(l, DMCC_COL, i) +
            (head - LUT_ValueAt(l, DMCH_COL, i)) *
            (LUT_ValueAt(l,DMCC_COL, i+1)-LUT_ValueAt(l,DMCC_COL, i))/
            (LUT_ValueAt(l,DMCH_COL, i+1)-LUT_ValueAt(l,DMCH_COL, i));
            */
}
//-----------------------------------------------------------------------------------
double TWorld::FindNode(double head, const  HORIZON *hor, int column)
{
    LUT *l = hor->lut;

    if (head >= 0) {
        return l->hydro[column].last();
    }

    auto it = std::lower_bound(l->hydro[H_COL].begin(), l->hydro[H_COL].end(), head);

    if (it == l->hydro[H_COL].begin()) {
        return(l->hydro[column][0]);
    } else if (it == l->hydro[H_COL].end()) {
        return(l->hydro[column].last());
    } else {
        int lowerIndex = std::distance(l->hydro[H_COL].begin(), it - 1); // Index of lower bound
        int upperIndex = std::distance(l->hydro[H_COL].begin(), it);     // Index of upper bound
        double lV = *(it - 1);
        double uV = *it;

        if (uV == lV) {
            return l->hydro[H_COL][lowerIndex]; // or some default value
        }


        double lTh = l->hydro[column][lowerIndex];
        double uTh = l->hydro[column][upperIndex];

        return lTh + (head - lV) * (uTh - lTh) / (uV-lV);
    }

}

//
// Created by philipp on 24.02.2023.
//
# include <cmath>

/** Modify Manning Strickler k (=1/NN) in and top of vegetation
 *   by Philipp Kraft, based on Oberle et al 2021
 *   (https://hdl.handle.net/20.500.11970/107539)
 *
 * @param NN - Manning's N for deeply submerged vegetation (watertable > 5 times plant height
 * @param WHr
 * @param PH
 * @return
 */
double calcManningType_1(double NN, double WHr, double PH) {

    double kst = 1.0 / NN;

    if (WHr <= PH) {
        // use 5 times slower runoff between stems
        kst = kst / 5;
    } else if (WHr <= 5 * PH) {
        // Linear interpolation from k/5 to k from PH to 5 * PH
        kst = kst / 5 + kst * (WHr - PH) / (5 * PH);
    } else {
        // kst stays the same
        kst = kst;
    }
    return 1.0 / kst;
}

/**
 * Calculate Manning n for partial vegetated channels Luhar and Nepf (2012)
 *
 * NN is the Manning's N for the unvegetated soil from the map
**/
double calcManningType_2(double NN, double WHr, double PH, double coverc) {
    const double
            CDrag=1,
            a=100,
            C=0.052, // The range of this parameter is 0.05–0.13 (Luhar and Nepf, 2012).
            grav_sqrt=pow(9.81, 0.5),
            kLN=1;
    double NNveg=0;
    if (WHr <= PH) {
        if (coverc < 0.8) {
            // Eq. 24
            NNveg = (kLN * pow(WHr,1./6.)) / grav_sqrt  * pow(C/2, 0.5) * pow(1 - coverc, -3./2.);
        } else {
            // Eq. 26
            NNveg = (kLN * pow(WHr,1./6.)) / grav_sqrt * pow(CDrag * a * WHr * 0.5, 0.5);
        }
    } else {
        NNveg = // Eq. 29
                (kLN * pow(WHr, 1./6.)) / grav_sqrt
                * ( 1 / (
                        pow(2 / C, 0.5)
                        * pow(1 - PH/WHr, 1.5)
                        + pow(2 * PH / (CDrag * a), 0.5)
                          * (1/WHr)
                )
                );
    }
    return NN + NNveg;
}

extern "C" {

double calcManningDispatch(int functiontype, double NN, double WHr, double PH, double coverc) {
    switch (functiontype) {
        case 1:
            return calcManningType_1(NN, WHr, PH);
        case 2:
            return calcManningType_2(NN, WHr, PH, coverc);
        default:
            return NN;
    }
}

}
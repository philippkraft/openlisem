//
// Created by philipp on 24.02.2023.
//

#include "model.h"

extern "C" {
    double calcManningDispatch(int functiontype, double NN, double WHr, double PH, double coverc, double Param1, double Param2);
}

double calcManning(TWorld * model, int r, int c, double NN) {
    return calcManningDispatch(
            model->jluManningFunctionType,
            NN,
            model->WHrunoff->Drc,
            model->PlantHeight->Drc,
            model->Cover->Drc,
            model->jluManningFunctionParam1,
            model->jluManningFunctionParam2
    );
}


#ifndef OPERATORS_H
#define OPERATORS_H
#include "Constants.h"
#include "Matrix.h"

void GetOperator(Hadron_prime *kin, Reaction_parameters *par, Matrix Operator[4]);

void Get_Res_ChPT(Hadron_prime *kin, Reaction_parameters *params, Matrix Operator[4]);


#endif

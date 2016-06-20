#ifndef MATRIX_H
#define MATRIX_H

#include "IMatrix.hh"
#include "NNZTag.hh"
#include "BucketPQ.hh"

#include <iostream>
#include <cmath>
#include "CLI.h"




///
///
/// \todo
//	This class need better namings for the variables.
class Matrix : public IMatrix
{

public:
    ~Matrix();
    Matrix(int M, int N, int nz, bool value);

    int updateDegreesToUVertices(int n, int ic, int u_maxdeg, int *jpntr,int *indRow,int *ipntr,int *indCol,
                                 bool * f_added, int *tag, int *f_tag, int *u_list,
                                 int *u_head, int *u_next, int *u_previous,int *list,int *blackList,const int q);

    // ========================================
    // IMatrix Methods
    bool slo(int *order);

    /**
     * Method ido()
     **/
    bool ido(int *order);
    bool lfo(int *order);



    bool computedegree();
    int greedycolor(int *list, int *ngrp);
    int rlf(int *ngrp);
    int sdo(int *ngrp);


    void rlf_slo(int *ngrp);
    void rlf_ido(int *ngrp);
    void rlf_lfo(int *ngrp);
    void slo_rlf(int *ngrp);
    void ido_rlf(int *ngrp);
    void lfo_rlf(int *ngrp);

    template <typename FirstOrderingMethod, typename SecondOrderingMethod>
    void mixedOrderingMethod(int *ngrp);



    int cseg_rlf_profit(int *ngrp);

    Matrix* getSeedMatrix(int *ngrp);

    int getNumberOfColors() const;

    // IMatrix Methods
    // ========================================

private:
    /**
     * Private Constructor
     */
    Matrix();

    bool decide(int ic);
    void rlf_mixup(CLI::ordering_method, int *ngrp);





    int *ndeg;



    int numberOfColors; // TODO: Remove

    double *val; // This is supposed to hold the values,
    // Currently we are not using it anywhere.
};

#endif

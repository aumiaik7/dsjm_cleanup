#ifndef PROFITPARTITIONMATRIX_HH
#define PROFITPARTITIONMATRIX_HH

#include "SimplePartitionedMatrix.hh"
#include "PartitionLoader.hh"
#include "NNZTag.hh"
#include "Result.hh"

class ProfitPartitionMatrix : public SimplePartitionedMatrix
{
public:
    ProfitPartitionMatrix(int M, int N, int nz,bool value,  PartitionLoader* partitionLoader);
    virtual ~ProfitPartitionMatrix();
    int cseg_profit(NNZTag *nnzTag);
    void generate_partition_from_nnz_color(NNZTag* nnzTag);

    // ================================================================================
    // Virtual Functions that need to be insntantiated before we can
    // create an object of CSegMatrix Class.

    bool computedegree();

    bool slo(int *order);
    bool ido(int *order);
    bool lfo(int *order);
    int rlf(int *ngrp);

    int greedycolor(int *list, int *ngrp);



    // ================================================================================

    void writeColor(Result result, int *ngrp);

private:
    NNZTag* nnzTag;

};

#endif

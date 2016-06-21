#ifndef CSEGMATRIX_H
#define CSEGMATRIX_H

#include "Matrix.hh"
#include "SimplePartitionedMatrix.hh"
#include "PartitionLoader.hh"
#include "fstream"
#include "_my_functor.hh"
#include "Result.hh"






class CSegMatrix : public SimplePartitionedMatrix
{
public:
    ~CSegMatrix();
    CSegMatrix(int M, int N, int nz, bool value, PartitionLoader *partitionLoader);

    // Virtual Functions that need to be insntantiated before we can create an object of CSegMatrix Class.

    bool computedegree();

    bool slo(int *order);
    bool ido(int *order);
    bool lfo(int *order);
    int rlf(int *ngrp);

    int greedycolor(int *list, int *ngrp);

    void annotate();



    void writeColor(Result result, int *ngrp);
private:
    CSegMatrix();

private:
    int *ndeg;

    // set<int> getAdjacentSegments(int segment) const;


    void updateForSLO(int adjacentSegment,
                      int numord,
                      int numberOfSegments,
                      int *colorList,
                      int *mindeg,
                      int *iwa1,
                      int *iwa2,
                      int *iwa3,
                      int *tag);

    int updateDegreesToUVertices(int n, int ic, int u_maxdeg, int *jpntr,int *indRow,int *ipntr,int *indCol,
                                 bool * f_added, int *tag, int *f_tag, int *u_list,
                                 int *u_head, int *u_next, int *u_previous,int *list,int *blackList,const int q);




    void updateForIDO(int segment, int *list,  int *maxinc, int *head, int *before,int *after);


    void cUpdateForRLF(int *list,
                       const int ic,
                       int *head,
                       int *next,
                       int *previous,
                       bool *inU,
                       int *countV_prime,
                       int *countV,
                       int *u_list,
                       int *u_head,
                       int *u_next,
                       int *u_previous,
                       const int nnzSegments,
                       int *u_maxdeg,
                       int *jpntr,
                       int *indRow,
                       int *ipntr,
                       int *indCol,
                       int *tag,
                       int *u_tag,
                       int *blackList,
                       const int q);

    int maxclq;
    bool isBlackListingEnabled ;






};

#endif // CSEGMATRIX_H

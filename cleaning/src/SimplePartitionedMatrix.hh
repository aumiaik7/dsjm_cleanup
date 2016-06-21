#ifndef SIMPLEPARTITIONEDMATRIX_HH
#define SIMPLEPARTITIONEDMATRIX_HH

#include "IMatrix.hh"
#include "PartitionLoader.hh"
#include "_my_functor.hh"
#include "NNZTag.hh"

class SimplePartitionedMatrix : public IMatrix
{
public:
    SimplePartitionedMatrix(int M, int N, int nz, bool value, PartitionLoader *partitionLoader);
    virtual ~SimplePartitionedMatrix();
    int getColumn(int segment) const;
    int getSegment(int row,int col) const;
    void createMapForNNZSegments();
    void IterateOverAdjacentSegments(int jcol, _my_functor* functor);
    void extend();
    PartitionLoader* getPartitionLoader() const;

    int getNumberOfSegments() const;



protected:
    PartitionLoader *partitionLoader;
    int nnzSegments;
    int *map;
    int *reverseMap;


private:
    SimplePartitionedMatrix();

};

#endif

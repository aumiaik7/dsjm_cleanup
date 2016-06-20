#ifndef PARTITIONLOADER_HH
#define PARTITIONLOADER_HH

#include <string>
using namespace std;

/**
 * This class is a helper class for the Column Segmented Matrix. The
 * objects of this class will be used for two purposes, one is to load
 * the partition information of a matrix from a file. And the other
 * side will provide this partition information to its clients. We are
 * assuming that the client will be the Column Segmented Matrix.
 */

class PartitionLoader
{
public:
  PartitionLoader(int m, int n);
  virtual ~PartitionLoader();
  bool loadFile(string fileName);
  void loadEachRowAsPartition();
  int getSegment(int row, int col) const;

  int M, N;
  int nPartitions;

  int *rowSequence;
  int *partitionPntr;
  int numberOfSegments;

private:
  void initializeGetSegment();
  int *inSegment; // TODO: Allocate and Deallocate in the constructor
		  // and Destructor.

};

#endif

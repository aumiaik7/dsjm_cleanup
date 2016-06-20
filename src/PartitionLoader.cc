#include "PartitionLoader.hh"

#include <string>
#include <iostream>
#include <fstream>
#include <cassert>
using namespace std;

/**
 * @Description : Constructor
 * @Creation Date: 19.06.2009
 */

PartitionLoader::PartitionLoader(int m, int n)
  :M(m),
   N(n) ,
   nPartitions(0)
{
  rowSequence = new int[M+1]; // TODO: I think, this thing is not
			      // necessary after certain points, so we
			      // can vanish this, and make it a
			      // temporary working array.
  inSegment = new int[M+1];
  // --------- FYI ---------
  // partitionPntr, which is an array , get its
  // storage after reading the values from block/partitioning file.
}
/* PartitionLoader() ENDS*/

/**
 * @Description : Destructor
 * @Creation Date: 19.06.2009
 */

PartitionLoader::~PartitionLoader()
{
  delete[] rowSequence;
  delete[] partitionPntr;
  delete[] inSegment;
}
/* PartitionLoader() ENDS*/

/**
 * @Description : This function reads the partition description from a
 * file.
 * @Creation Date: 20.06.2009
 */
bool
PartitionLoader::loadFile(string fileName)
{
  ifstream a_file(fileName.c_str());

  if(a_file)
    {
      a_file >> nPartitions;

      numberOfSegments = nPartitions * N;
      partitionPntr = new int[nPartitions + 2];

      for (int i = 1; i <= M; i++)
	{
	  a_file >> rowSequence[i];
	}
      for (int i = 1; i <= nPartitions + 1; i++)
	{
	  a_file >> partitionPntr[i];
	}

      initializeGetSegment();
    }
  else
    {
      cerr << fileName << " cannot be found" << endl;
    }


}
/* loadFile() ENDS*/

/**
 * @Description : This method intializes the list on which the
 * getSegment() method depends.
 * @Creation Date: 22.06.2009
 */
void
PartitionLoader::initializeGetSegment()
{
    /*
     * We assume that the lists for row indices and pntr has been
     * already populated.
     */

    // cout << "initializeGetSegment() " << endl;

    // We iterate.
    // i is equal to row.
    int currentSegment = 1;
    int nextSegmentPntr = partitionPntr[currentSegment+1];
    for (int i = 1; i <= M; i++)
    {
        nextSegmentPntr = partitionPntr[currentSegment + 1];
        assert(rowSequence[i] > 0);
        assert(rowSequence[i] <= M);
        inSegment[rowSequence[i]] = currentSegment;

        if((i+1) == nextSegmentPntr)
        {
            currentSegment++;
        }
    }
}
/* initializeGetSegment() ENDS*/

/**
 * @Description : This functions
 * @Creation Date: 22.06.2009
 */
int
PartitionLoader::getSegment(int row, int col) const
{
  return ((col - 1) * nPartitions + inSegment[row]);
}
/* getSegment() ENDS*/


/**
 * @Description : This method loads each row as a partition.
 * @Creation Date: 17.08.2009
 */
void
PartitionLoader::loadEachRowAsPartition()
{

  nPartitions = M;

  numberOfSegments = nPartitions * N;
  partitionPntr = new int[nPartitions + 2];

  for (int i = 1; i <= M; i++)
    {
      rowSequence[i] = i ;
    }
  for (int i = 1; i <= nPartitions + 1; i++)
    {
      partitionPntr[i] = i ;
    }

  initializeGetSegment();


}
/* loadEachRowAsPartition() ENDS*/



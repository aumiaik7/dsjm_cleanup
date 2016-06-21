#include "SimplePartitionedMatrix.hh"
#include "extend_functor.hh"
#include <cassert>

// (setq fill-column 100)

/**
 * @Description : Constructor
 * @Creation Date: 16.08.2009
 */
SimplePartitionedMatrix::SimplePartitionedMatrix(int M,int N,int nz,bool value,  PartitionLoader* partitionLoader)
    : IMatrix(M,N,nz, value),
      partitionLoader(partitionLoader),
      nnzSegments(0)
{
    map = new int[partitionLoader->numberOfSegments + 1];
    reverseMap = new int[partitionLoader->numberOfSegments + 1];
}
/* SimplePartitionedMatrix() ENDS*/

/**
 * @Description : Destructor
 * @Creation Date: 16.08.2009
 */
SimplePartitionedMatrix::~SimplePartitionedMatrix()
{
  delete[] map;
  //TODO: delete[] reversemap? ??
}
/* ~SimplePartitionedMatrix() ENDS*/

/**
 * @Description : This fuction returns the column in which a segment
 * resides.
 * @Creation Date: 12.06.2009
 */
int
SimplePartitionedMatrix::getColumn(int segment) const
{
    return (( reverseMap[segment] - 1) / partitionLoader->nPartitions  + 1);
}
/* getColumn() ENDS*/

int SimplePartitionedMatrix::getSegment(int ir, int col) const
{
    return map[partitionLoader->getSegment(ir,col)];
}

void SimplePartitionedMatrix::createMapForNNZSegments()
{
    int *tag = NULL;
    try
    {
        tag = new int[partitionLoader->numberOfSegments + 1] ;
        for (int jcol = 1; jcol <= partitionLoader->numberOfSegments; jcol++)
        {
            tag[jcol] = 0;
        }

        for ( int jcol = 1; jcol <= N; jcol++)
        {
            for (int jp = jpntr[jcol] ; jp < jpntr[jcol + 1]; jp++)
            {
                int ir = indRow[jp];
                int segment = partitionLoader->getSegment(ir,jcol);
                if(tag[segment] == 0)
                {
                    tag[segment] = 1;
                    nnzSegments++;
                    map[segment] = nnzSegments;
                    reverseMap[nnzSegments] = segment;
                }
            }
        }
        delete[] tag;
    }
    catch(std::bad_alloc)
    {
        delete[] tag;
    }
}

/** TODO: Make it a template */
/**
 * Purpose: 		To perform an operation defined by a functor on all adjacent column segments of
 *                  a given column segment.
 *
 * Pre-condition: 	The matrix object is nonempty. Assumes that the computeCCS(), compress() and
 *                  computeCRS() has been called prior calling this function, such that matrix
 *                  object holds the sparsity information in Compressed Column and Compressed Row
 *                  storage.
 *
 * Post-condition: 	functor() is evoked for all the adjacent column segment of
 *                  <id:referenceSegment>. The post condition depends on the definition of the
 *                  functor class supplied as a parameter.
 *
 *
 * Parameters:      in-parameter <id:referenceSegment>, an integer value indicating the
 *                  target/reference segment for the functor.
 *
 *                  in-parameter <id:functor>, a functor object which will be performed on all the
 *                  adjacent column segment of the <id:referenceSegment>
 *
 *
 */
void
SimplePartitionedMatrix::IterateOverAdjacentSegments(int referenceSegment, _my_functor* functor)
{

    // Determine the column of the current segment.
    int col = getColumn((referenceSegment));


    for(int jp = jpntr[col] ; jp < jpntr[col+1]; jp++)
    {
        int ir = indRow[jp];	// Iterate over the rows of
        // the columns of the marked segment.

        int tSegment = getSegment(ir,col); // TODO: As we are iterating
        // over the nonzero columns only we should not have zero
        // elements at all .

        for( int ip = ipntr[ir]; ip <= ipntr[ir+1] -1; ip++)
        {
            int ic = indCol[ip];
            if( ic == col)
                continue;

            // Determine the adjacent column segment.
            int segment = getSegment(ir,ic);

            // Feed the functor object with the adjacent column segment.
            functor->segment = segment;

            // If not marked, perform the functor method. For example for deg_functor, it increases
            // the ndeg[referenceSegment] and ndeg[currentSegment] for
            if( ! functor->isMarked())
            {
                functor->mark();
                (*functor)();
            }

            // For nonzero entries in the current segment, we find more adjacent segment.
            if(getSegment(ir,col) ==  (referenceSegment))
            {
                for (int x_jp = jpntr[ic]; x_jp <= jpntr[ic+1] -1; x_jp++)
                {
                    int x_ir = indRow[x_jp];
                    segment = getSegment(x_ir,ic);

                    // Perform the operation defined by functor.
                    functor->segment = segment;
                    if(! functor->isMarked())
                    {
                        functor->mark();
                        functor->segment = segment;
                        (*functor)();
                    }
                }
            }
        }
    }
}

// It is a copy of the Degree Method, so the requirements are same as the degree thing
void
SimplePartitionedMatrix::extend()
{
    int *tag;

    // TODO: FIX: POSSIBLE BUG.
    // I am just adding it here to produce a quick result.
    try
    {
        tag = new int[nnzSegments+1];

        /**
         * Initialization Block
         */

        for (int jp = 1; jp <= nnzSegments; jp++)
        {
            // ndeg[jp] = 0;
            tag[jp] = 0;
        }

        /**
         * Beginning of Iteration Loop
         */
        extend_functor *eFunctor = new extend_functor(nnzSegments);
        for (int  j = 1; j <= nnzSegments ; j++)
        {
            tag[j] = nnzSegments;

            eFunctor->tag = tag;
            // eFunctor->ndeg = ndeg;
            eFunctor->j = j;
            IterateOverAdjacentSegments(j,eFunctor);

        }
        delete eFunctor;
        delete[] tag;

        /* ofstream out("out.degree");
           for (int i = 1; i <= nnzSegments ; i++)
           {
           out << i << " " << ndeg[i] << endl;
           }
        */
    }
    catch(std::bad_alloc)
    {
        delete[] tag;
        throw;
    }
}

PartitionLoader* SimplePartitionedMatrix::getPartitionLoader() const
{
    return partitionLoader;
}

int SimplePartitionedMatrix::getNumberOfSegments() const
{
    return nnzSegments;
}

// -*- mode: c++; fill-column:80 ; comment-fill-column: 80; -*-
// (setq fill-column 100)
#include "CSegMatrix.hh"
#include <set>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <cassert>
#include <fstream>
#include "SimplePartitionedMatrix.hh"
#include "deg_functor.hh"
#include "seq_functor.hh"
#include "Utility.h"
#include "Result.hh"


CSegMatrix::CSegMatrix(int M, int N, int nz,bool value, PartitionLoader *partitionLoader)
    : SimplePartitionedMatrix(M,N,nz,value, partitionLoader),
      isBlackListingEnabled(false)
{


}

CSegMatrix::~CSegMatrix()
{
    delete[] ndeg;
}

/**
 * Purpose: 		Computes Degree sequence of the columns of a Column Segmented sparse matrix A
 *          		(i.e. the vertices of the column intersection graph G(A) ).
 *
 * Pre-condition: 	The matrix object is nonempty. Assumes that the computeCCS(), compress() and
 *                  computeCRS() has been called prior calling this function, such that matrix
 *                  object holds the sparsity information in Compressed Column and Compressed Row
 *                  storage.
 *
 * Post-condition: 	Degree information for the column segments of matrix A ( graph G(A) ) is stored
 * 			        in the data member <id:ndeg>, an integer array of size nnzSegments+1 such that
 * 			        if k = ndeg[j] then the column j has degree k, where j = 1,2, ...,nnzSegments
 *
 * Return values:   Returns true when the function is executed successfully,
 *                  otherwise returns false.
 *
 */
bool CSegMatrix::computedegree() // TODO: Test This Method Thoroughly
{
    int *tag;
    int maxdeg = -1;

    try
    {
        // Allocate memory for the <id:ndeg> integer array of size nnzSegments +
        // 1, where nnzSegment is the number of column segments in the Column
        // Segmented sparse Matrix A.
        ndeg = new int[nnzSegments + 1];

        tag = new int[nnzSegments + 1]; // Temporary working array of size
                                        // nnzSegments+ 1. If tag[seg] =
                                        // nnzSegment, then the degree of column
                                        // segment "seg" has been computed.

        // Initialize <id:ndeg> and <id:tag>
        for (int s = 1; s <= nnzSegments; s++)
        {
            ndeg[s] = 0;
            tag[s] = 0;
        }

        for ( int j = 1; j <= nnzSegments ; j++)
        {
            tag[j] = nnzSegments;

            // Prepare the functor object which increases the corresponding <id:ndeg> elements by 1
            // for each adjacent column segment.

            deg_functor *dFunctor = new deg_functor();
            dFunctor->tag = tag;
            dFunctor->ndeg = ndeg;
            dFunctor->j = j;
            dFunctor->maxdeg = -1;

            // Perform increase by 1 for the adjacent segments of current segment j.
            IterateOverAdjacentSegments(j,dFunctor);

            maxdeg = dFunctor->maxdeg;
            delete dFunctor;
        }


        // The following code is used only Verification test, does not execute
        // in the normal execution of this method.
        if(getVerify())
        {
            ofstream out("out.degree");
            for (int i = 1; i <= nnzSegments ; i++)
            {
                out << i << " " << ndeg[i] << endl;
            }
        }


        // Release Memory
        if(tag) delete[] tag;

        // Returns true on success.
        return true;
    }
    catch(std::bad_alloc)
    {
        delete[] tag;
        return false;
    }
}

/**
 * Purpose: 		Computes Smallest-Last Ordering (SLO) of the column segments of a
 *          		column-segmented sparse matrix A (i.e. the vertices of the column intersection
 *          		graph G(A) )
 *
 * Pre-condition: 	The matrix object is nonempty. Assumes that the degree of of the columns have
 *                  already been computed in the data member <id:ndeg> integer array of size
 *                  nnzSegments+1 using computeDegree() method.
 *
 * Post-condition: 	The SLO ordering of matrix A ( graph G(A) ) is stored in the out-parameter
 * 			        <id:order>, an integer array of size nnzSegments+1 such that if k = order[j] then
 * 			        the column segment j is the k-th element, k = 1,2, ..., nnzSEgments; in the SLO
 * 			        ordering, and j = 1,2, ..., nnzSegments.
 *
 * Parameters:      out-parameter <id:order>, an integer pointer to an array of size
 *                  nnzSegments+1. The array will contain the ordering information when the function
 *                  normally returns.
 *
 * Return values:   Returns true when the function is executed successfully,
 *                  otherwise returns false.
 *
 */
bool CSegMatrix::slo(int *order)
{

    int ip,ir,jcol,jp, mindeg,numdeg,numord;

    int *head = NULL;
    int *prev = NULL;
    int *next = NULL;
    int *tag = NULL;


    try
    {

        // The following three integer arrays consist of a doubly linked list. It acts as a bucket
        // priority queue for the incidence degree of the columns.

        // head(deg) is the first column segment in the deg list unless head(deg) = 0. If head(deg)
        // = 0 there are no column segments in the deg list.

        // previous(col) is the column segment before col in the degree list unless previous(col) =
        // 0. If previous(col) = 0, col is the first column segment in this list.

        // next(col) is the column segment after col in the degree list unless next(col) = 0. If
        // next(col) = 0, col is the last column segment in this incidence list.

        // if col is in un-ordered column segment, then order[col] is the induced degree of col. If
        // col is an ordered column, then order[col] is the smallest-last order of column col.

        head = new int[nnzSegments];
        prev = new int[nnzSegments+1];
        next = new int[nnzSegments+1];


        tag = new int[nnzSegments+1]; // Temporary array, used for marking ordered columns.


        mindeg = nnzSegments;
        for(jp=1;jp <= nnzSegments; jp++)
        {
            head[jp-1] = 0 ;
            tag[jp] = nnzSegments;
            order[jp] = ndeg[jp];

            assert(order[jp] >=  -1);
            // TODO: FIXME: Skipping of the zero
            // degree edges for Debugging purpose.
            // 	  if(ndeg[jp] != 0 )
            if(order[jp] > -1)
                mindeg = min(mindeg,ndeg[jp]);
        }
        assert(mindeg >= 0); // TODO:
        assert(mindeg <= nnzSegments); // TODO:


        // Build the priority queue.
        for(jp = 1; jp <= nnzSegments; jp++)
        {
            numdeg = ndeg[jp];
            prev[jp] = 0;

            assert(numdeg < nnzSegments);
            assert(numdeg >= -1);
            if(numdeg < 0)
            {
                assert(numdeg == -1);
                continue;
            }
            next[jp] = head[numdeg];

            if(head[numdeg] > 0)
            {
                prev[head[numdeg]] = jp;
            }
            head[numdeg] = jp;
        }

        numord = nnzSegments;

        do
        {
            // find a column segment jcol with the minimum degree.
            do
            {
                if(mindeg >= nnzSegments)
                {
                    mindeg = -1;
                    break;
                }
                assert(mindeg < nnzSegments);
                assert(mindeg >= 0);
                jcol = head[mindeg]; // fixme: i am getting an memory
                // error here, it might be a start
                // of a series, but i want to fix
                // it first.
                if(jcol  > 0)
                    break;
                mindeg = mindeg + 1;
            }while(1);
            if(mindeg == -1)
                break;
            order[jcol] = numord;
            numord = numord -1;

            // When numord == 0, we have already processed all the column segments.
            if (numord == 0)
                break;

            // Mark column segment jcol
            tag[jcol] = 0;

            // Delete column segment jcol from the mindeg list.
            head[mindeg] = next[jcol];
            if(next[jcol] > 0)
            {
                prev[next[jcol]] = 0;
            }



            //
            // Determine all nonzero column segments.
            //

            // Determine the column number for the current column segment jcol
            int col = getColumn(jcol);
            for(int jp = jpntr[col] ; jp <= jpntr[col+1] -1 ; jp++)
            {
                int ir = indRow[jp];	// iterate over the rows of
                // the columns of the marked
                // segment.

                int tsegment = getSegment(ir,col); // todo: as we are iterating
                // over the nonzero
                // columns only w should
                // not have zero elements
                // at all .

                for( int ip = ipntr[ir]; ip <= ipntr[ir+1] -1; ip++)
                {
                    int ic = indCol[ip];
                    if( ic == col)
                        continue;
                    int segment = getSegment(ir,ic);
                    if(tag[segment] > numord)
                    {
                        tag[segment] = numord;
                        updateForSLO(segment, numord, nnzSegments, order,&mindeg, head,prev,next, tag);
                    }

                    if(getSegment(ir,col) ==  jcol)
                    {
                        for (int x_jp = jpntr[ic]; x_jp <= jpntr[ic+1] -1; x_jp++)
                        {
                            int x_ir = indRow[x_jp];
                            segment = getSegment(x_ir,ic);
                            if(tag[segment] > numord)
                            {
                                tag[segment] = numord;
                                updateForSLO(segment, numord, nnzSegments, order,&mindeg, head,prev,next, tag);
                            }
                        }
                    }
                }
            }
            /**
             * end of iteration loop
             */
        }while(true);
    eighty:


        // TODO: Why are we inverting the order in SLO() , it must be a bug. September 10, 2010.
        // Invert the array order.
        for (jcol =1; jcol <= nnzSegments; jcol++)
        {
            if(order[jcol] > -1)
                prev[order[jcol]] = jcol; // fixme: i am having a
            // valgrind memory error here
            // too.
        }
        for(jp = 1; jp <= nnzSegments; jp++)
        {
            order[jp] = prev[jp];
        }

    }
    catch (std::bad_alloc)
    {
        std::cerr << "memory exhausted\n";

        if(head) delete[] head;
        if(prev) delete[] prev;
        if(next) delete[] next;
        if(tag) delete[] tag;

        return false;
    }

    if(head) delete[] head;
    if(prev) delete[] prev;
    if(next) delete[] next;
    if(tag) delete[] tag;

}

/**
 * Purpose: 		Computes Incidence-Degree Ordering (IDO) of the column segments of a column
 *          		segmented sparse matrix A (i.e. the vertices of the column intersection graph
 *          		G(A) )
 *
 * Pre-condition: 	The matrix object is nonempty. Assumes that the degree of of the column segments
 *                  have already been computed in the data member <id:ndeg> integer array of size
 *                  nnzSegments+1 using computeDegree() method.
 *
 * Post-condition: 	The IDO ordering of matrix A ( graph G(A) ) is stored in the out-parameter
 * 			        <id:order>, an integer array of size nnzSegments+1 such that if k = order[j]
 * 			        then the column segment j is the k-th element, k = 1,2, ..., nnzSegment, in the
 * 			        IDO ordering, and j = 1,2, ..., nnzSegment.
 *
 * Parameters:      out-parameter <id:order>, an integer pointer to an array of size
 *                  nnzSegment+1. The array will contain the ordering information when the function
 *                  normally returns.
 *
 * Return values:   Returns true when the function is executed successfully,
 *                  otherwise returns false.
 *
 */
bool CSegMatrix::ido(int *order)
{
    int ic,ip,ir,jcol,jp,
        maxinc,maxlst,ncomp,numinc,numlst,numord,numwgt;

    int *head = NULL;
    int *previous = NULL;
    int *next = NULL;
    int *tag = NULL;
    try
    {

        // The following three integer arrays consist of a doubly linked list. It acts as a bucket
        // priority queue for the incidence degree of the columns.

        // head(deg) is the first column segment in the deg list unless head(deg) = 0. If head(deg)
        // = 0 there are no column segments in the deg list.

        // previous(col) is the column segment before col in the incidence list unless previous(col)
        // = 0. If previous(col) = 0, col is the first column segment in this incidence list.

        // next(col) is the column segment after col in the incidence list unless next(col) = 0. If
        // next(col) = 0, col is the last column segment in this incidence list.

        // if col is in un-ordered column segment, then order[col] is the incidence degree of col to
        // the graph induced by the ordered column segments. If col is an ordered column segment,
        // then order[col] is the incidence-degree order of column segment col.
        head = new int[nnzSegments];
        previous = new int[nnzSegments + 1];
        next = new int[nnzSegments + 1];

        tag = new int[nnzSegments + 1]; // Temporary working array, used for marking ordered columns


        // Sort the indices of degree array <id:ndeg> in descending order, i.e
        // ndeg(tag(i)) is in descending order , i = 1,2,...,nnzSegment
        //
        // <id:tag> is used here as an in-out-parameter to <id:indexSort> routine. It
        // will hold the sorted indices. The two arrays, <id:previous> and
        // <id:next> is used for temporary storage required for <id:indexSort>
        // routine.
        MatrixUtility::indexsort(nnzSegments,nnzSegments - 1,ndeg,-1,tag,previous,next);

        maxinc = 0;
        // They are just alias for the following loop.
        // They are here only for readability.
        {
            for(jp =nnzSegments ; jp >= 1 ; jp--)
            {
                ic = tag[jp];
                head[nnzSegments-jp] = 0;
                previous[ic] = 0;
                next[ic] = head[0];
                if (head[0] > 0)
                    previous[head[0]] = ic;
                head[0] = ic;
                tag[jp] = 0;
                order[jp] = 0;

            }
        }

        // Determine the maximal search length to search for maximal degree in the maximal incidence
        // degree list.
        maxlst = 0;
        for(ir =1 ; ir <= M ; ir++)
        {
            maxlst = maxlst + MatrixUtility::square(ipntr[ir+1] - ipntr[ir]);
        }
        maxlst = maxlst/nnzSegments; //TODO: I am not sure about this.
        // It was maxlast/N before.
        // I am not sure whether we use N or
        // nnzSegments
        maxclq = 0;
        numord = 1;
        do
        {
            // Update the size of the largest clique found during the ordering.
            if (maxinc == 0)
                ncomp = 0;
            ncomp = ncomp + 1;
            if (maxinc + 1 == ncomp)
                maxclq = max(maxclq,ncomp);

            // chooose a column segment jcol of maximal incidence degre.
            do
            {
                jp = head[maxinc];
                if (jp > 0)
                    break;
                maxinc = maxinc - 1;

            }while(1);

            // We search a distance of maxLast length to find the column segment with maximal degree
            // in the original column-segmented graph.
            numwgt = -1;
            for(numlst = 1; numlst <= maxlst; numlst++)
            {

                if (ndeg[jp] > numwgt)
                {
                    numwgt = ndeg[jp];
                    jcol = jp;
                }
                jp = next[jp];
                if (jp <= 0)
                    break;

            }

            order[jcol] = numord;
            numord = numord + 1;
            tag[jcol] = nnzSegments;

            // Termination Test.
            if( numord > nnzSegments)
                break;


            // delete column segment jcol fromt he maxinc list.
            deleteColumn(head,next,previous,maxinc,jcol);

            // Determine the column of the current segment.
            int col = getColumn(jcol);
            for(jp = jpntr[col] ; jp <= jpntr[col+1] -1; jp++)
            {
                ir = indRow[jp];
                // Iterate over the rows of the columns of the marked segment.


                int tsegment = getSegment(ir,col);


                for(ip = ipntr[ir];ip <=  ipntr[ir+1]-1; ip++)
                {
                    ic = indCol[ip];

                    if( ic == col)
                        continue;

                    // Determine the adjacent column segment.
                    int segment = getSegment(ir,ic);
                    if (tag[segment] < numord)
                    {
                        tag[segment] = numord;
                        updateForIDO(segment, order,  &maxinc, head, previous,next);
                    }

                    if( getSegment(ir,col) == jcol)
                    {
                        for (int x_jp = jpntr[ic]; x_jp < jpntr[ic] ; x_jp++)
                        {
                            int x_ir = indRow[x_jp];
                            segment = getSegment(x_ir,ic);
                            if(tag[segment] > numord)
                            {
                                tag[segment] = numord;
                                updateForIDO(segment, order,  &maxinc, head, previous,next);
                            }
                        }
                    }
                }
            }

        }while(1);
    IDO_hundred:
        // Invert the array order.
        for( jcol = 1;jcol<= nnzSegments; jcol++)
        {
            if(order[jcol] > -1)
                previous[order[jcol]] = jcol;
        }
        for( jp = 1;jp <= nnzSegments; jp++)
        {
            order[jp] = previous[jp];
        }

    }
    catch (std::bad_alloc)
    {
        std::cerr << "Memory Exhausted\n";

        if(head) delete[] head;
        if(previous) delete[] previous;
        if(next) delete[] next;
        if(tag) delete[] tag;

        return false;
    }

    if(head) delete[] head;
    if(previous) delete[] previous;
    if(next) delete[] next;
    if(tag) delete[] tag;

    return true;

}
/**
 * Purpose: 		Computes Largest-First Ordering (LFO) of the column segment of a
 *          		column-segmented sparse matrix A (i.e. the vertices of the column intersection
 *          		graph G(A) )
 *
 * Pre-condition: 	The matrix object is nonempty. Assumes that the degree of of the column segments
 *                  have already been computed in the data member <id:ndeg> integer array of size
 *                  nnzSegment+1 using computeDegree() method.
 *
 * Post-condition: 	The LFO ordering of matrix A ( graph G(A) ) is stored in the out-parameter
 * 			        <id:order>, an integer array of size nnzSegment+1 such that if k = order[j] then
 * 			        the column segment j is the k-th element, k = 1,2, ..., nnzSegment; in the LFO
 * 			        ordering, and j = 1,2, ..., nnzSegment.
 *
 * Parameters:      out-parameter <id:order>, an integer pointer to an array of size
 *                  nnzSegment+1. The array will contain the ordering information when the function
 *                  normally returns.
 *
 * Return values:   Returns true when the function is executed successfully,
 *                  otherwise returns false.
 *
 */
bool CSegMatrix::lfo(int *order)
{


    int *head = NULL;
    int *previous = NULL;
    int * next = NULL;
    try
    {
        // The following three integer arrays consist of a doubly linked list. It acts as a bucket
        // priority queue for the incidence degree of the columns.

        // head(deg) is the first column segment in the deg list unless head(deg) = 0. If head(deg)
        // = 0 there are no column segments in the deg list.

        // previous(col) is the column segment before col in the degree list unless previous(col) =
        // 0. If previous(col) = 0, col is the first column segment in this list.

        // next(col) is the column segment after col in the degree list unless next(col) = 0. If
        // next(col) = 0, col is the last column segment in this incidence list.

        head = new int[nnzSegments+1];
        previous = new int[nnzSegments+1];
        next = new int[nnzSegments+1];


        int maxdeg = -1;
        for(int jp=1;jp <= nnzSegments; jp++)
        {
            head[jp-1] = 0 ; // We use degree as an index to find a column segment from the head
                             // list, which ranges from 0, ..., nnzSegment-1
            maxdeg = max(maxdeg,ndeg[jp]);
        }

        // Initialize the priority queue.
        buildPriorityQueue(nnzSegments,ndeg,head,next,previous);

        int numord = 1;
        int jcol;
        while(true)
        {
            // Choose a column segment jcol of maximal degree.
            do
            {
                jcol = head[maxdeg];
                if (jcol > 0)
                    break;
                maxdeg = maxdeg -1 ;
            }while(true);

            order[jcol] = numord;
            numord = numord +1;
            if (numord > nnzSegments )
            {
                delete[] head;
                delete[] next;
                delete[] previous;
                return true;
            }


            // Delete Jcol from the head of the list.
            head[maxdeg] = next[jcol];

            if(next[jcol] > 0)
            {
                previous[next[jcol]] = 0;
            }
        }
    }
    catch(std::bad_alloc)
    {
        if(next) delete[] next ;
        if(head) delete[] head;
        if(previous) delete[] previous;

        return false;
    }

}

/**
 * Purpose: 		Computes Recursive Largest-First coloring (RLF) of the column segments of a
 *          		column-segmented sparse matrix A (i.e. the vertices of the column intersection
 *          		graph G(A) )
 *
 * Pre-condition: 	The matrix object is nonempty. Assumes that the degree of
 *                  of the column segments have already been computed in the data member
 *                  <id:ndeg> integer array of size nnzSegment+1 using computeDegree() method.
 *
 * Post-condition: 	RLF coloring of Matrix A(graph G(A)) is stored in the in-out-parameter
 * 			        <id:color>, an integer array of size nnzSegment+1, such that if k = color[j]
 * 			        then the column segment j is colored with color k, j = 1,2,...,nnzSegment
 *
 * Parameters:      out-parameter <id:color>, an integer pointer to an array of size
 *                  nnzSegment+1. The array will contain the color values of the column segments in
 *                  successful completion. The integer array uses 1-based indexing.
 *
 *
 * Return values:   Returns the number of colors if succeeds, otherwise returns
 *                  0(zero).
 *
 */
int CSegMatrix::rlf(int *ngrp)
{
    /*
     * Overview:
     * In RLF coloring algorithm, we maintain three sets of vertices in three
     * sets,
     *     1. set V for the admissible column segments, initially it contains all the
     *        from the graph G(A).
     *     2. set U for the column segments which are non-admissible to current color class q. At
     *        start of a new color class this set is empty.
     *     3. set C for the colored class.
     *
     * At the start of each color class we choose a column segment jcol with the maximal
     * degree in set V.
     * At other steps we choose a column segment jcol from set V , which has the maximal
     * number of neighbors in set U, we call it U-Degree.
     *
     * As each column segment is chosen, it is colored with the value of the current
     * color class q, and moved from the set V to C. All the adjacent column segments
     * are moved to set U, as inadmissible column segment for the current. set.
     *
     * As column segment are added to set U, we update the U-Degree of each column segment in
     * set V.
     *
     * Coloring is finished when all the columns are colored.
     */

    int *head = NULL;
    int *previous = NULL;
    int *next = NULL;

    int *tag = NULL;
    int *blackList = NULL;

    int *list = NULL;


    int *u_head = NULL;
    int *u_previous = NULL;
    int *u_next = NULL;
    int *u_list = NULL;
    bool *inU = NULL;
    int *u_tag = NULL;


    try
    {
        // The following three integer arrays consist of a doubly linked list. It acts as a bucket
        // priority queue for the incidence degree of the columns.

        // head(deg) is the first column segment in the deg list unless head(deg) = 0. If head(deg)
        // = 0 there are no column segments in the deg list.

        // previous(col) is the column segment before col in the degree list unless previous(col) =
        // 0. If previous(col) = 0, col is the first column segment in this list.

        // next(col) is the column segment after col in the degree list unless next(col) = 0. If
        // next(col) = 0, col is the last column segment in this incidence list.

        // if col is in un-ordered column segment, then order[col] is the induced degree of col. If
        // col is an ordered column, then order[col] is the smallest-last order of column col.
        head = new int[nnzSegments];
        previous = new int[nnzSegments+1];
        next = new int[nnzSegments+1];

        tag = new int[nnzSegments+1]; // For a column segment jcol, if tag[jcol] = nnzSegment, then
                                      // this column segment has already been colored. If 0 <
                                      // tag[jcol] ( = numord) < N , then jcol has been processed
                                      // for a column segment in numord step.

        blackList = new int[M+1]; // If blackList[irow] = q, where q is the color class and irow is
                                  // a row number, then any column segment having nonzero element in irow-th
                                  // row cannot be included in the q-th color class. We maintain
                                  // this array to gain better performance in RLF.

        list = new int[nnzSegments+1];

        // Initialize blackList array.
        for ( int i = 1 ; i <= M; i++)
        {
            blackList[i] = 0;
        }



        // The following arrays are used for building a doubly linked list based on U-Degree of the
        // column segments. They are used in the same way as <id:head>, <id:previous>, <id:next> ,
        // <id:list>, <id:tag > is used. We use this doubly linked list arrays to form a priority
        // queue for chosing column segment from set V.
        u_head = new int[nnzSegments];
        u_previous = new int[nnzSegments+1];
        u_next = new int[nnzSegments+1];
        u_list = new int[nnzSegments+1];
        inU = new bool[nnzSegments+1];
        u_tag = new int[nnzSegments+1];


        int u_maxdeg = 0;
        int u_numdeg = 0;

        int u_jp,u_ir,u_ip,u_ic;

        int q = 1; // Current color class, each column segment picked is colored to the value of q.
        int jcol;
        int ic;
        int maxdeg,numdeg;
        int ip,ir,jp;
        int numord = 1;

        maxdeg = 0;

        // Initialize the integer arrays <id:tag> , <id:inU> , <id:u_tag> and build both of the
        // priority queues(doubly linked list).
        for(jp = 1; jp <=nnzSegments ; jp++)
        {
            head[jp-1] = 0;
            tag[jp] = 0;
            list[jp] = ndeg[jp];
            u_list[jp] = 0;
            u_head[jp-1] = 0;
            u_next[jp] = 0;
            u_previous[jp] = 0;
            maxdeg = max(maxdeg,ndeg[jp]);
            inU[jp] = false;

            u_tag[jp] = 0;
        }

        // Initialize both of the priority queues.
        for (jp = 1; jp <= nnzSegments ; jp++)
        {
            numdeg = ndeg[jp];

            previous[jp] = 0;

            next[jp] = head[numdeg];

            if ( head[numdeg] > 0 )
            {
                previous[head[numdeg]] = jp;
            }
            head[numdeg] = jp;

            u_numdeg = 0;
            u_previous[jp] = 0;
            u_next[jp] = u_head[u_numdeg] ;
            if ( u_head[u_numdeg] > 0)
            {
                u_previous[u_head[u_numdeg]] = jp;
            }
            u_head[u_numdeg] = jp;
        }



        int countV = nnzSegments; // Number of elements in set V
        int countU = 0;  // Number of elements in set U
        int countC = 0; // Number of elementes in set C
        int count = 0;
        bool newColorClass = true; // Flag variablet o indicate whether we have just picked a
                                   // column segment for a new color clas or not. It evaluates
                                   // to true for the first column segment in each color class.
        while(true)
        {
            int jcol;

            if (newColorClass == true)
            {
                newColorClass = false;

                //  Choose a column segment jcol of maximal degree from priority queue of set V.
                do
                {
                    jcol = head[maxdeg];
                    numdeg = maxdeg;
                    if(jcol > 0)
                        break;
                    maxdeg = maxdeg - 1;
                }while(1);

            }
            else
            {
                // Choose a column segmet that has maximal degree in set U.
                do
                {
                    assert(u_maxdeg >= 0);
                    jcol = u_head[u_maxdeg] ;
                    u_numdeg = u_maxdeg;
                    if (jcol > 0)
                    {
                        break;
                    }
                    u_maxdeg = u_maxdeg -1;
                } while(1);
            }

            // Update number counters;
            countV--;
            countC++;

            numdeg = list[jcol] ;
            u_numdeg = u_list[jcol];

            u_list[jcol] = -1; // For Debug Purpose // TODO:

            // list[jcol] = numord;
            list[jcol] = q;

            // Added on 31st Decembepr 2008, while making things encapsulated in the
            // Matrix Object.
            ngrp[jcol] = q;


            numord = numord + 1;
            if(numord > nnzSegments)
            {
                // De allocate Memory.
                if(next ) delete[] next ;
                if(head) delete[] head;

                if(previous) delete[] previous;
                if(tag) delete[] tag;

                if(u_head) delete[] u_head;
                if(u_previous) delete[] u_previous;
                if(u_next) delete[] u_next;
                if(u_list) delete[] u_list;
                if(inU) delete[] inU;
                if(u_tag) delete[] u_tag;

                if(blackList) delete[] blackList;

                if(list) delete[] list;
                return q;
            }




            deleteColumn(head,next,previous,numdeg,jcol);
            deleteColumn(u_head,u_next,u_previous,u_numdeg,jcol);

            /**
             * Tag the column, as it has been colored.
             */
            tag[jcol] = nnzSegments;

            /**
             * Find All Columns adjacent to column jcol.
             * Determine all positions (ir,jcol) which correspond
             * to non-zeroes in the matrix.
             */
            int col = getColumn(jcol);

            if(isBlackListingEnabled)
            {
                for(int jp = jpntr[col] ; jp < jpntr[col+1] ; jp++)
                {
                    int ir = indRow[jp];
                    blackList[ir] = q;
                }
            }




            for ( int jp = jpntr[col] ; jp < jpntr[col+1]  ; jp++)
            {
                ir = indRow[jp];

                int tsegment = getSegment(ir,col);

                for ( ip = ipntr[ir]; ip < ipntr[ir+1]; ip++)
                {
                    ic= indCol[ip];

                    if(ic == col)
                        continue;

                    int segment = getSegment(ir,ic);


                    if(tag[segment] < numord)
                    {
                        tag[segment] = numord;
                        cUpdateForRLF( list, segment, head, next, previous, inU,
                                       &countU, &countV, u_list, u_head, u_next,
                                       u_previous, nnzSegments, &u_maxdeg, jpntr,
                                       indRow, ipntr, indCol, tag, u_tag, blackList,
                                       q);
                    }

                    if( getSegment(ir,col) == jcol)
                    {
                        for (int x_jp = jpntr[ic] ; x_jp <= jpntr[ic+1] -1; x_jp++)
                        {
                            int x_ir = indRow[x_jp];
                            segment = getSegment(x_ir, ic);
                            if(tag[segment] < numord)
                            {
                                tag[segment] = numord;
                                cUpdateForRLF( list, segment, head, next, previous, inU,
                                               &countU, &countV, u_list, u_head, u_next,
                                               u_previous, nnzSegments, &u_maxdeg, jpntr,
                                               indRow, ipntr, indCol, tag, u_tag, blackList,
                                               q);
                            }
                        }
                    }
                }
            }
            // countV + countC + countU == nnzSegment.
            // If countV = 0, the set of admissible columns  V is empty. We
            // start a new color class, and reset the priority queue for
            // elements in set U.
            if(countV == 0)
            {
                q = q + 1;

                newColorClass = true;

                // Swap values.
                countV =  countU;
                countU = 0;

                u_numdeg = 0;
                u_maxdeg = 0;

                // Reset the priority queue for the elements in set Set U.
                initializeDegreesToUVertices(nnzSegments,tag,u_head,u_next,u_previous,u_list,inU,u_tag);
            }
        }
    }
    catch(std::bad_alloc)
    {

        // De allocate Memory.
        if(next) delete[] next ;
        if(head) delete[] head;

        if(previous) delete[] previous;
        if(tag) delete[] tag;

        if(u_head) delete[] u_head;
        if(u_previous) delete[] u_previous;
        if(u_next) delete[] u_next;
        if(u_list) delete[] u_list;
        if(inU) delete[] inU;
        if(u_tag) delete[] u_tag;

        if(blackList) delete[] blackList;

        if(list) delete[] list;
        return 0;
    }
}

void CSegMatrix::cUpdateForRLF(int *list,
                               const int ic,
                               int *head,
                               int *next,
                               int *previous,
                               bool *inU,
                               int *countU,
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
                               const int q)
{

    /**
     * Update the pointers to the current
     * degree lists.
     */

    int numdeg = list[ic];
    list[ic] = list[ic] - 1;

    // TODO: Implementation of the two priority queue is
    // cluttering up the functions for RLF algorithm.

    assert(numdeg > 0);
    deleteColumn(head,next,previous,numdeg,ic);
    addColumn(head,next,previous,(numdeg-1),ic);

    /**
     * Now We add this ic to our second set U
     */

    // inU[ic] = true;
    // newColorClass = true;
    if (inU[ic] == false)
    {
        //printf("ADD ic = %d to  V_Prime\N",ic);
        inU[ic] = true;
        *countU = *countU + 1;
        *countV = *countV - 1;

        /**
         * This column should get deleted from the
         * u_head. u_next and u_previous too.
         */
        int u_numdeg = u_list[ic];

        u_list[ic] = -1; // For Debug Purpose
        deleteColumn(u_head,u_next,u_previous,u_numdeg,ic);

        // Update the degrees of the adjacent vertices.

        *u_maxdeg = updateDegreesToUVertices(nnzSegments,ic,*u_maxdeg, jpntr,indRow,ipntr,indCol,inU,
                                             tag,u_tag,u_list,u_head,u_next,u_previous,list,blackList,q);

    }

}


/**
 * Purpose:             Computes the greedy coloring of the column segments of column-segmented
 *                      sparse matrix A (i.e. the vertices of the column intersection graph G(A))
 *
 * Pre-condition:       The matrix object is nonempty. Assumes that an ordering has been provided in
 *                      the in-parameter <id:order> integer array of size nnzSegments+1, such that
 *                      order[1]...order[nnzSegments] is a permutation of {1,...,nnzSegments}
 *
 * Post-condition:      The greedy coloring of column segmented matrix A( graph G(A) ) is stored in
 *                      the in-out-parameter <id:color>, an integer array of size nnzSegments+1,
 *                      such that if k = color[j] then the column j is colored with color k, j =
 *                      1,2,...,nnzSegments
 *
 *
 * Parameters:          in-parameter <id:order>, an integer pointer to an array of size n+1,
 *                      containing a permutation of {1,...,nnzSegments}. The integer array uses
 *                      1-based indexing.
 *
 *                      in-out-parameter <id:color>, an integer pointer to an array of size
 *                      nnzSegments+1, it stores the color values of the columns in successful
 *                      completion.The integer array uses 1-based indexing.
 *
 * Return values:       Returns the number of colors if succeeds, otherwise
 *                      returns 0(zero).
 */
int CSegMatrix::greedycolor(int *order, int *color)
{
    if(order == NULL || color == NULL)
        return 0;

    int *w; // Working array of size nnzSegments +1. It is used to mark the colors already used for
            // adjacent column segments.

    try
    {
        w = new int[nnzSegments+1];

        int ip,ir,j,jcol,jp;
        int maxgrp = 0;

        // Initialization of the arrays.
        for (int jp = 1; jp <=  nnzSegments ;jp++  )
        {
            color[jp] = nnzSegments;
            w[jp] = 0;
        }

        for (int j = 1; j <=  nnzSegments ;j++  )
        {
            int jcol = order[j];

            // Prepare the seq_functor, which marks an adjacent column with the color of the current
            // column segment.

            seq_functor *sFunctor = new seq_functor();
            sFunctor->temparray = w;
            sFunctor->ngrp = color;
            sFunctor->j = j;

            // Perform the functor on all the adjacent column segments.
            IterateOverAdjacentSegments(jcol,sFunctor);
            delete sFunctor;



            // Assign the smallest un-marked color number to jcol
            for (jp = 1; jp <=  maxgrp ;jp++  )
            {
                if (w[jp] != j)
                    goto SEQ_L50;
            }
            maxgrp = maxgrp + 1;
        SEQ_L50:
            color[jcol] = jp;
        }



        // The following code is used only Verification test, does not execute
        // in the normal execution of this method.
        if(getVerify())
        {
            ofstream out("out.seq");
            for (int i = 1; i <= nnzSegments; i++)
            {
                out << i << " " << color[i] << endl;
            }
        }
        // -------------------- Verification code ends

        // Release memory
        if(w) delete[] w;

        // Return the number of colors.
        return maxgrp;
    }
    catch(std::bad_alloc)
    {
        delete[] w;
        return 0;
    }
}



void
CSegMatrix::updateForSLO(int adjacentSegment, int numord, int numberOfSegments, int *list,int *mindeg, int *temparray1,int *temparray2, int *temparray3, int *tag)
{

    assert(adjacentSegment > 0 && adjacentSegment<=numberOfSegments);

    /**
     * Update the pointers to the current degree lists.
     **/
    int numdeg = list[adjacentSegment];
    assert(list[adjacentSegment] != -1);

    assert(list[adjacentSegment] > 0);
    list[adjacentSegment] = list[adjacentSegment] -1;

    *mindeg = min(*mindeg,list[adjacentSegment]);
    /**
     * Delete column ic from the numdeg list.
     */
    if(temparray2[adjacentSegment] == 0)
    {
        assert(numdeg >= 0);
        temparray1[numdeg] = temparray3[adjacentSegment];
    }
    else
    {
        temparray3[temparray2[adjacentSegment]] = temparray3[adjacentSegment];
    }

    if(temparray3[adjacentSegment] > 0)
    {
        temparray2[temparray3[adjacentSegment]] = temparray2[adjacentSegment];
    }

    /**
     * Add column ic to the numdeg -1 list.
     */

    temparray2[adjacentSegment] = 0;
    assert(numdeg > 0);
    temparray3[adjacentSegment] = temparray1[numdeg-1];
    if(temparray1[numdeg-1] > 0)
    {
        temparray2[temparray1[numdeg-1]]= adjacentSegment;
    }
    temparray1[numdeg-1] = adjacentSegment;

}



void CSegMatrix::writeColor(Result result, int *ngrp)
{
    ofstream out("out.color");
    out << result.totalColors << endl;
    for (int jcol = 1; jcol <= nnzSegments; jcol++)
    {
        out << jcol << " " << ngrp[jcol] << endl;
    }
}

void CSegMatrix::updateForIDO(int segment, int *list,  int *maxinc, int *head, int *before,int *after)
{
    /*
     * update the pointers to the current incidence lists.
     */
    int numinc = list[segment];
    list[segment] = list[segment] + 1;
    *maxinc = max(*maxinc,list[segment]);
    /*
     * delete column segment from the numinc list.
     */
    if (before[segment] == 0)
        head[numinc] = after[segment];
    else
        after[before[segment]] = after[segment];

    if (after[segment] > 0)
        before[after[segment]] = before[segment];
    /*
     * add column segment to the numinc+1 list.
     */
    before[segment] = 0;
    after[segment] = head[numinc+1];
    if (head[numinc+1] > 0)
        before[head[numinc+1]] = segment;
    head[numinc+1] = segment;
}


/**
 *
 * @Description:
 *  This method is called when we move a vertex from set V'
 *  to set U.
 * This method has a complexity of \sum{\rho_i}
 */

int CSegMatrix::updateDegreesToUVertices(int n, int jcol,int maxdeg, int *jpntr,int *indRow,
                                         int *ipntr,int *indCol, bool * inU, int *tag, int *u_tag,
                                         int *u_list, int *head, int *next, int *previous,int *list, int *blackList,const int q)
{
    // Update the degrees of the adjacent vertices.

    int numdeg;
    int ir;
    int ic;
    bool zero = false;

    int col = getColumn(jcol);
    for( int jp = jpntr[col] ; jp < jpntr[col+1] ; jp++)
    {
        ir = indRow[jp];

        if(isBlackListingEnabled && blackList[ir] == q)
            continue;
        for(int ip = ipntr[ir]; ip < ipntr[ir+1] ; ip++)
        {
            ic = indCol[ip];

            if(ic == col)
                continue;

            int segment = getSegment(ir,ic);

            if(inU[segment] == false && tag[segment] < n && u_tag[segment] != jcol)
            {
                u_tag[segment] = jcol;
                /**
                 * Update the pointers to the current degree u_lists.
                 */
                numdeg = u_list[segment];
                u_list[segment] = u_list[segment] + 1;
                maxdeg = max(numdeg + 1,maxdeg);

                deleteColumn(head,next,previous,numdeg,segment);
                addColumn(head,next,previous,numdeg+1,segment);

            }
            if(getSegment(ir,col) == jcol)
            {
                for (int x_jp = jpntr[ic] ; x_jp <= jpntr[ic+1] -1 ; x_jp++)
                {
                    int x_ir = indRow[x_jp];
                    segment = getSegment(x_ir,ic);
                    if(inU[segment] == false && tag[segment] < n && u_tag[segment] != jcol)
                    {
                        u_tag[segment] = jcol;
                        /**
                         * Update the pointers to the current degree u_lists.
                         */
                        numdeg = u_list[segment];
                        u_list[segment] = u_list[segment] + 1;
                        maxdeg = max(numdeg + 1,maxdeg);

                        deleteColumn(head,next,previous,numdeg,segment);
                        addColumn(head,next,previous,numdeg+1,segment);

                    }
                }
            }
        }
    }
    return maxdeg;
}

/**
 * @Description : output to a file the nam oe the segments.
 * @Creation Date: 25.09.2009
 */
void
CSegMatrix::annotate()
{
    ofstream out("out.annotate");
    for (int jcol = 1; jcol <= N; ++jcol)
    {
        for(int jp = jpntr[jcol] ; jp <= jpntr[jcol + 1] - 1; ++jp)
        {
            int ir = indRow[jp];
            int segment = getSegment(ir,jcol);
            out << "(" << ir << "," << jcol << ") is in segment " << segment << std::endl;
        }
    }
}
/* annotate() ENDS*/

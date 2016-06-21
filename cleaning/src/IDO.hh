#ifndef IDO_HH
#define IDO_HH

#include "detail/PriorityQueue.hh"
#include "Utility.h"
#include "OrderingMethod.hh"


template <template <typename> class T_PriorityQueue>
class IDO
{
public:
    const OrderingMethod::ordering_direction ord_direction;
private:
    int M;
    int N;
    int *jpntr;
    int *indRow;
    int *ipntr;
    int *indCol;
    int *ndeg;
    int *list;
    int *ngrp;
    int numord;
    int maxinc;
    int maxlst;
    int maxclq;
    T_PriorityQueue<MaxQueue> priority_queue;
    std::vector<int> tag;
    int howMany;

    int colorCount;
public:
    // This is our constructor.
    template <typename PreviousOrderingMethod>
    IDO(
        int M_,
        int N_,
        int *ndeg_,
        int *jpntr_,
        int *indRow_,
        int *ipntr_,
        int *indCol_ ,
        int *list_,
        int *ngrp_,
        int maxBucket_,
        PreviousOrderingMethod& previousOrderingMethod,
        int howMany_
        )
        :
        ord_direction(OrderingMethod::FORWARD),
        M(M_),
        N(N_),
        jpntr(jpntr_),
        indRow(indRow_),
        ipntr(ipntr_),
        indCol(indCol_),
        ndeg(ndeg_),
        list(list_),
        maxinc(0),
        maxlst(0),
        maxclq(0),
        priority_queue(maxBucket_, N_),
        colorCount(0),
        howMany(howMany_)
    {
        // Initialize Numord
        if( previousOrderingMethod.ord_direction == this->ord_direction)
            numord = previousOrderingMethod.getNumOrd();
        else
            numord = 1;
        int *last = new int[N+1];
        int *next = new int[N+1];
        int *index = new int[N+1];
        MatrixUtility::indexsort(N,N-1, ndeg, -1, index, last, next);

        delete[] last;
        delete[] next;

        this->tag.reserve(N+1);


        for ( int jp = N; jp >= 1; --jp)
        {
            int ic = index[jp];
            if( previousOrderingMethod.isTagged(ic))
            {
                tagit(ic);
            }
            else
            {
                priority_queue.insert(ic,0);
                untagit(ic);
            }
        }
        delete[] index;

        // Determine the maximal search length for the list of columsn
        // of maximal incidence
        for (int ir = 1; ir <= M; ir++)
        {
            // We are converting a double into an integer ?
            // I hope this does not introduce any error.
            maxlst = maxlst + MatrixUtility::square(ipntr[ir+1] - ipntr[ir]);
        }
        maxlst = maxlst/ N;
        maxclq = 0;

    }

    // This is our destructor.
    ~IDO()
    {
    }

    bool isTagged(int jcol) const
    {
        return (this->tag[jcol]  == N);
    }

    // This is the actual ordering engine. It picks one vertex/
    // columen given the states in the priority queues. and then
    // updates all supporting data structure for the next pickup.
    void work()
    {
        int ncomp;
        do
        {
            /*
             * update the size of the largest clique
             * found during the ordering.
             */

            if (maxinc == 0)
                ncomp = 0;
            ncomp = ncomp + 1;
            if (maxinc + 1 == ncomp)
                maxclq = max(maxclq,ncomp);
            /*
             * choose a column jp of maximal degree among the
             * columns of maximal incidence maxinc.
             */
            int jp;
            int jcol;
            Item item = priority_queue.top();
            jp = item.index;
            maxinc = item.priority;

            item = priority_queue.betterTop(jp, maxlst, ndeg);
            jcol = item.index;
            maxinc = item.priority;

            priority_queue.remove(jcol);

            colorCount++;

            list[numord] = jcol;
            // printf("list[%d] = %d\n",numord,jcol);
            numord = numord + 1;

            tag[jcol] = N;
            /*
             *        termination test.
             */
            // 	if (numord > N) goto IDO_hundred;
            if( numord > N || colorCount >= howMany)
                break;




            for(int jp = jpntr[jcol] ; jp <= jpntr[jcol+1] -1; jp++)
            {
                int ir = indRow[jp];
                /*
                 * for each row ir, determine all positions (ir,ic)
                 * which correspond to non-zeroes in the matrix.
                 */

                for(int ip = ipntr[ir];ip <=  ipntr[ir+1]-1; ip++)
                { // 80
                    int ic = indCol[ip];
                    if (tag[ic] < numord)
                    {
                        tag[ic] = numord;
                        /*
                         * update the pointers to the current incidence lists.
                         */
                        priority_queue.increase(ic);
                        int numinc = priority_queue.get(ic).priority;
                        maxinc = max(maxinc,numinc);
                    }
                }
            }
            /*
             * end of iteration loop.
             */
        }while(1);
    }

    inline int getNumOrd()
    {
        return numord;
    }

    inline int getDeg(int jcol)
    {
        return ndeg[jcol];
    }
    inline void tagit(int jcol)
    {
        tagit(jcol,N);
    }

    inline void tagit(int jcol, int value)
    {
        tag[jcol] = value;
    }

    inline void untagit(int jcol)
    {
        tagit(jcol,0);
    }

};

#endif

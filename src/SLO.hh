#ifndef SLO_HH
#define SLO_HH

#include "detail/PriorityQueue.hh"
#include "OrderingMethod.hh"

template <template <typename> class T_PriorityQueue>
class SLO
{
public:
    const OrderingMethod::ordering_direction ord_direction;
private:
    int mindeg, numord, maxclq,N, *list, *jpntr, *indRow, *ipntr, *indCol;
    T_PriorityQueue<MinQueue> priority_queue;
    std::vector<int> tag;
    int howMany;
    int colorCount;
public:

    template<typename PreviousOrderingMethod>
    SLO(
        int M_,
        int N_,
        int *ndeg_,
        int *jpntr_,
        int *indRow_,
        int * ipntr_,
        int *indCol_,
        int *list_,
        int *ngrp_,
        int maxBucket_,
        PreviousOrderingMethod& previousOrderingMethod,
        int howMany_
        )
        :ord_direction(OrderingMethod::REVERSE),
         N(N_),
         mindeg(N_),
         maxclq(0),
         list(list_),
         jpntr(jpntr_),
         indRow(indRow_),
         ipntr(ipntr_),
         indCol(indCol_),
         priority_queue(maxBucket_, N),  // This is how we initialize a
         // composite member.
         howMany(howMany_),
         colorCount(0)
    {
        // Initialize Numord
        if( previousOrderingMethod.ord_direction == this->ord_direction)
            numord = previousOrderingMethod.getNumOrd();
        else
            numord = N;
        this->tag.reserve(N+1);
        for (int jp = 1; jp <= N; jp++)
        {
            if(previousOrderingMethod.isTagged(jp))
            {
                tagit(jp);
            }
            else
            {
                priority_queue.insert(jp,previousOrderingMethod.getDeg(jp));
                mindeg = std::min(mindeg, previousOrderingMethod.getDeg(jp));
                untagit(jp);
            }
        }
    }

    bool isTagged(int jcol) const
    {
        return (this->tag[jcol] == 0);
    }


    ~SLO()
    {

    }
    void work()
    {
        while(1)
        {
            int ic,ip, ir, jcol, jp, numdeg;
            /*
             *        MARK THE SIZE OF THE LARGEST CLIQUE
             *        FOUND DURING THE ORDERING.
             */
            if ((mindeg +1 == numord ) && (maxclq == 0) )
            {
                maxclq = numord;
            }
            /*
             *        CHOOSE A COLUMN JCOL OF MINIMAL DEGREE MINDEG.
             **/


            Item item = priority_queue.top();
            jcol = item.index;
            mindeg = item.priority;

            ++colorCount;
            // By Definition, we are
            // seleting the minimum
            // degree element, so it
            // is be default the mindeg

            priority_queue.pop();

            list[numord] = jcol;
            // printf("list[%d] = %d\n",numord,jcol);
            numord = numord -1;

            /**
             * Find all columns adjacent to column jcol.
             **/

            assert(tag[jcol] != 0) ;
            tagit(jcol);

            /**
             * Termination Test
             */
            if (colorCount >= howMany || numord == 0)
            {
                return;
            }






            /**
             * Determine all positions (ir,jcol) which correspond
             * to non-zeroes in the matrix.
             */

            for(jp = jpntr[jcol]; jp <= jpntr[jcol+1] -1;jp++)
            {
                ir = indRow[jp] ;
                /*
                 * For each row ir,determine all positions (ir,ic)
                 * which correspond to non-zeroes in the matrix.
                 */
                for(ip = ipntr[ir] ; ip <= ipntr[ir+1] - 1; ip++)
                {
                    ic = indCol[ip];
                    /* Array tag marks columns which are adjacent to
                     * column jcol
                     */

                    if(tag[ic] > numord)
                    {

                        tagit(ic,numord);
                        /**
                         * Update the pointers to the current degree lists.
                         **/


                        priority_queue.decrease(ic);
                        numdeg = priority_queue.get(ic).priority;
                        mindeg = std::min(mindeg,numdeg);

                    }
                }
            }
        }
    }

    inline int getNumOrd() const
    {
        return numord;
    }

    inline int tagit(int jcol)
    {
        tag[jcol] = 0;
    }

    inline int tagit(int jcol, int value)
    {
        tag[jcol] = value;
    }

    inline int untagit(int jcol)
    {
        tagit(jcol,N);
    }

    inline int getDeg(int jcol)
    {
        return priority_queue.get(jcol).priority;
    }


};

#endif

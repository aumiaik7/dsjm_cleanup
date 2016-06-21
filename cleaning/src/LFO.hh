#ifndef LFO_HH
#define LFO_HH

#include "detail/PriorityQueue.hh"
#include "Utility.h"
#include "OrderingMethod.hh"

template <template <typename> class T_PriorityQueue>
class LFO
{
public:
    OrderingMethod::ordering_direction ord_direction;
private:
    int N;
    int *jpntr;
    int *indRow;
    int *ipntr;
    int *indCol;
    int *ndeg;
    int *list;
    int numord;
    int howMany;

    int colorCount;
    std::vector<int> tag;

    T_PriorityQueue<MaxQueue> priority_queue;


public:
    template <typename PreviousOrderingMethod>
    LFO(
        int M_,
        int N_,
        int *ndeg_ ,
        int *jpntr_,
        int *indRow_,
        int *ipntr_,
        int *indCol_,
        int *list_,
        int *ngrp_,
        int maxBucket_,
        PreviousOrderingMethod& previousOrderingMethod,
        int howMany_
        ):
        ord_direction(OrderingMethod::FORWARD),
        N(N_),
        jpntr(jpntr_),
        indRow(indRow_),
        ipntr(ipntr_),
        indCol(indCol_),
        ndeg(ndeg_),
        list(list_),
        priority_queue(maxBucket_, N_),
        howMany(howMany_),
        colorCount(0)
    {
        tag.reserve(N+1);
        if(previousOrderingMethod.ord_direction == this->ord_direction)
            numord = previousOrderingMethod.getNumOrd();
        else
            numord = 1;
        for (int jp = 1; jp <= N; ++jp)
        {
            if( previousOrderingMethod.isTagged(jp))
            {
                tagit(jp);
            }
            else
            {
                untagit(jp);
                priority_queue.insert(jp, ndeg[jp]);
            }
        }

    }
    // This is our destructor.
    ~LFO()
    {
    }

    // This is the actual engine for LFO, where we pick each column
    // according to their degree.
    void work()
    {
        while(true)
        {
            int jcol, maxdeg;
            Item item = priority_queue.top();
            jcol = item.index;
            maxdeg = item.priority;
            priority_queue.pop();
            list[numord] = jcol;
            numord = numord +1;
            colorCount++;
            tagit(jcol);
            if (numord > N || colorCount >= howMany)
            {
                return;
            }
        }
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

    inline int getNumOrd()
    {
        return numord;
    }

    inline int getDeg(int jcol)
    {
        return ndeg[jcol];
    }

    bool isTagged(int jcol) const
    {
        return (this->tag[jcol]  == N);
    }

};

#endif

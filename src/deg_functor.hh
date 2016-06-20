#ifndef DEG_FUNCTOR_HH
#define DEG_FUNCTOR_HH

#include "_my_functor.hh"

class deg_functor : public _my_functor
{
public:
    int *tag;
    int *ndeg;
    int j;
    int maxdeg;

    void operator()()
    {
        ndeg[j]++;
        ndeg[segment]++;
        maxdeg = std::max(ndeg[j], maxdeg);
        maxdeg = std::max(ndeg[segment],maxdeg);
    }

    bool isMarked()
    {
        if( tag[segment] < j)
            return false;
        else
            return true;
        // return (! tag[segment] < j );
    }

    //
    void mark()
    {
        tag[segment] = j;
    }
};

#endif

#ifndef NULLORDERINGMETHOD_HH
#define NULLORDERINGMETHOD_HH

#include "OrderingMethod.hh"

class NullOrderingMethod
{
public:
    const OrderingMethod::ordering_direction ord_direction;
private:
    int *ndeg;
public:
    NullOrderingMethod(int *ndeg_)
        :ord_direction(OrderingMethod::FORWARD),
         ndeg(ndeg_)
    {
    }
    inline bool isTagged(int jcol) const
    {
        return false;
    }

    inline int getDeg(int jcol) const
    {
        return ndeg[jcol];
    }

    inline int getNumOrd() const
    {
        return 1;
    }

};

#endif

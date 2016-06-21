#ifndef SEQ_FUNCTOR_HH
#define SEQ_FUNCTOR_HH

#include "_my_functor.hh"

class seq_functor : public _my_functor
{
 public:
  int *temparray;
  int *ngrp;

  int j;

  void operator()()
  {
    // iwa[ngrp[segment]] = j;
    // cout << "iwa[" << ngrp[segment] <<"] = " << j << endl;
  }

  bool isMarked()
  {
    return false;
  }

  void mark()
  {
    temparray[ngrp[segment]] = j;
  }
};

#endif

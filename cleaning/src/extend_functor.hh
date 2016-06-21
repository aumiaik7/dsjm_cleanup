#ifndef EXTEND_FUNCTOR_HH
#define EXTEND_FUNCTOR_HH

#include "deg_functor.hh"
#include <iostream>
#include <fstream>
using namespace std;

class extend_functor : public deg_functor
{
 public:
  ofstream out;
  extend_functor(int numberOfSegments)
    {
      out.open("out.extend");
      out << numberOfSegments << endl;
    }
  ~extend_functor()
    {
      out.close();
    }
  void operator()()
  {
    out << j << " " << segment << endl;
  }
};

#endif

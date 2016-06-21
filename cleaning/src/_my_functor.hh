#ifndef _MY_FUNCTOR_HH
#define _MY_FUNCTOR_HH

class _my_functor
{
 public:
  int segment;
  virtual void operator()()= 0;
  virtual bool isMarked() = 0;
  virtual void mark() = 0;
};

#endif

#ifndef ADAMS_BASHFORTH_H
#define ADAMS_BASHFORTH_H
#include "Explicit_LMM.h"

class Adams_Bashforth:public Explicit_LMM
{
 private:
  static const std::vector<std::vector<Real>> Beta;
 public:
 Adams_Bashforth():Explicit_LMM(){};
  virtual  void Solve_IVPs(const int steps,
			   const Real T,const IVPs & F,std::vector<std::vector<Real>> &A, int p=1);
  virtual ~Adams_Bashforth()=default;
};
#else
//Do nothing!
#endif

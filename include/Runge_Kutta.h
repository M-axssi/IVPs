#ifndef RUNGE_KUTTA_H
#define RUNGE_KUTTA_H
#include "TimeIntegrator.h"

class Runge_Kutta:public TimeIntegrator
{
 public:
  Runge_Kutta():TimeIntegrator(){};
  virtual void  Solve_IVPs(const int steps,
			   const Real T,const IVPs & F,std::vector<std::vector<Real>>& S,int p=1);
  virtual ~Runge_Kutta()=default;
};

#else
//Do nothing!
#endif

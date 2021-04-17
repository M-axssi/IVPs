#ifndef TIMEINTEGRATOR_H
#define TIMEINTEGRATOR_H
#include "Utility.h"
#include "IVPs.h"

class TimeIntegrator
{
 public:
  TimeIntegrator(){};
  virtual void  Solve_IVPs(const int steps,
			   const Real T,const IVPs & F,std::vector<std::vector<Real>>& S,int p=1)=0;
  virtual ~TimeIntegrator()=default;
};

#else
//DO NOTHING!
#endif

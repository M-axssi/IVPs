#ifndef ADAMS_MOULTON_H
#define ADAMS_MOULTON_H
#include "Implicit_LMM.h"

class Adams_Moulton:public Implicit_LMM
{
 private:
  static const std::vector<std::vector<Real>> Beta;
  
 public:
 Adams_Moulton():Implicit_LMM(){};
  virtual void  Solve_IVPs(const int steps,
			   const Real T,const IVPs & F,std::vector<std::vector<Real>>& A,int p=1) ;
  virtual ~Adams_Moulton()=default;
};

#else
//Do nothing!
#endif

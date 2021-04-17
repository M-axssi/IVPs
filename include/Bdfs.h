#ifndef BDFS_H
#define BDFS_H
#include "Implicit_LMM.h"

class Bdfs:public Implicit_LMM
{
 private:
  static const std::vector<std::vector<Real>> Alpha;
  static const std::vector<Real> Beta;
  static const int _s[4];
  
 public:
 Bdfs():Implicit_LMM(){};
  virtual void  Solve_IVPs(const int steps,
			   const Real T,const IVPs & F,std::vector<std::vector<Real>>& A,int p=1);
  virtual ~Bdfs()=default;
};

#else
//Do nothing!
#endif

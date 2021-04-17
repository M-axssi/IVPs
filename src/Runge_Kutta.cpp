#include "Runge_Kutta.h"

void Runge_Kutta::Solve_IVPs(const int steps,
			     const Real T,const IVPs & F,std::vector<std::vector<Real>> & S,int p)
{
  S.push_back(F.Get_u0());
  if (steps==0) return;
  const Real k=T/steps;
  for (int i=1;i<=steps;++i)
    {
      std::vector<Real> temp(S[i-1]);
      std::vector<Real> y_1,y_2,y_3,y_4;
      y_1=F.Get_diff(S[i-1]);
      y_2=F.Get_diff(S[i-1]+(k/2)*y_1);
      y_3=F.Get_diff(S[i-1]+(k/2)*y_2);
      y_4=F.Get_diff(S[i-1]+k*y_3);
      S.push_back(S[i-1]+(k/6)*(y_1+2*y_2+2*y_3+y_4));
    }
}

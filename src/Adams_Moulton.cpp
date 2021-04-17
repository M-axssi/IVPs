#include "Adams_Moulton.h"

const std::vector<std::vector<Real>> Adams_Moulton::Beta=
  {  {1.0/2,1.0/2},
     {-1.0/12,8.0/12,5.0/12},
     {1.0/24,-5.0/24,19.0/24,9.0/24},
     {-19.0/720,106.0/720,-264.0/720,646.0/720,251.0/720}};
const int Adams_Moulton::_s[4]={1,2,3,4};

void Adams_Moulton::Solve_IVPs(const int steps,                                  
			       const Real T,const IVPs & F,std::vector<std::vector<Real>>& S, int p)
{
  S.push_back(F.Get_u0());
  int s=_s[p-2];
  if (steps==0) return;
  const Real k=T/steps;
  for (int i=1;i<=s-1;++i)
    {
      std::vector<Real> temp((-1)*S[i-1]);
      std::vector<Real> initial(S[i-1]);
      initial+=k*F.Get_diff(S[i-1]);
      for (int j=0;j<=i-1;++j)
        {
      	  temp+=(-k*Beta[i-1][j])*F.Get_diff(S[j]);
      	}
      Nonlinear_Solve(initial,temp,F,-k*Beta[i-1][i]);
      S.push_back(initial);
    }
  for (int i=s;i<=steps;++i)
    {
      std::vector<Real>temp((-1)*S[i-1]);
      std::vector<Real> initial(S[i-1]);
      initial+=k*F.Get_diff(S[i-1]);
      for (int j=i-s;j<=i-1;++j)
	{
	  temp+=(-k*Beta[s-1][j-i+s])*F.Get_diff(S[j]);
	}
      Nonlinear_Solve(initial,temp,F,-k*Beta[s-1][s]);
      S.push_back(initial);
    }
}

#include "Adams_Bashforth.h"

const std::vector<std::vector<Real>> Adams_Bashforth::Beta=
  {{1.0},
   {-1.0/2,3.0/2},
   {5.0/12,-16.0/12,23.0/12},
   {-9.0/24,37.0/24,-59.0/24,55.0/24}};
const int Adams_Bashforth::_s[4]={1,2,3,4};

void Adams_Bashforth::Solve_IVPs(const int steps,
				 const Real T,const IVPs & F,std::vector<std::vector<Real>> & S, int p)
{
  S.push_back(F.Get_u0());
  int s=_s[p-1];
  if (steps==0) return;
  const Real k=T/steps;
  for (int i=1;i<=s-1;++i)
    {
      std::vector<Real> temp(S[i-1]);
      for (int j=0;j<=i-1;++j)
	{
	  temp+=(k*Beta[i-1][j])*F.Get_diff(S[j]);
	}
      S.push_back(temp);
    }
  for (int i=s;i<=steps;++i)
    {
      std::vector<Real> temp(S[i-1]);
      for (int j=i-s;j<=i-1;++j)
        {
	  temp+=(k*Beta[s-1][j-i+s])*F.Get_diff(S[j]);
	}
      S.push_back(temp);
    }
}

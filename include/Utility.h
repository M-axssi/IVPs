#ifndef MY_UTILITY_H
#define MY_UTILITY_H
#include <vector>
#include <lapacke.h>
#include <cblas.h>
#include <map>
#include <string>
#include <utility>
#include <cassert>

using Real=double;

void operator+=(std::vector<Real> &a,const std::vector<Real> & b);

std::vector<Real>  operator*(const Real a,const std::vector<Real> & c);

std::vector<Real> operator-(const std::vector<Real> &a,const std::vector<Real> & b);

std::vector<Real> operator+(const std::vector<Real> & a,const std::vector<Real> &b);

Real * VV_ADD(const std::vector<Real> &a,const std::vector<Real> &b);

void  VP_ADD(std::vector<Real> &a,const Real * b);


#else
//Do nothing!
#endif

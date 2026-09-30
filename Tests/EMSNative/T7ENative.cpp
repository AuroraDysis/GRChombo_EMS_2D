#include <cmath>
// Independent native stencil replay; current fields only, no reader/evolution.
#include "EMSCouplingFunction.hpp"
#include "CCZ4Cartoon.hpp"
#include "ConstraintsCartoon.hpp"
#include "EMSCouplingFunction.hpp"
#include "EMSCartoonGaussConstraints.hpp"
#include "BoxLoops.hpp"
#include "GammaCartoonCalculator.hpp"
#include <vector>
#include <fstream>
#include <iostream>
int main(int argc,char **argv)
{
 const bool gamma=argc==4 && std::string(argv[3])=="--gamma";if(argc!=3 && !gamma)return 2;
 std::ifstream in(argv[1],std::ios::binary);std::ofstream out(argv[2],std::ios::binary);std::vector<double> r(gamma?2271:703);const int width=gamma?9:5,halo=width/2;
 while(in.read(reinterpret_cast<char*>(r.data()),r.size()*sizeof(double)))
 {
  const int j=std::llround(r[1]/r[0]-.5);IntVect c(D_DECL(0,j,0));Box b(c-halo*IntVect::Unit,c+halo*IntVect::Unit);
  FArrayBox u(b,28),q(Box(c,c),NUM_DIAGNOSTIC_VARS);
  for(int n=0;n<28;++n)for(int y=-halo;y<=halo;++y)for(int x=-halo;x<=halo;++x)u(c+IntVect(D_DECL(x,y,0)),n)=r[3+n*width*width+(y+halo)*width+x+halo];
  BoxLoops::loop(Constraints<CouplingFunction>(r[0],CouplingFunction({4*M_PI,0.,0.,-.8}),1.),u,q,Box(c,c),disable_simd());
  BoxLoops::loop(EMSCartoonGaussConstraints(r[0],{4*M_PI,0.,0.,-.8}),u,q,Box(c,c),disable_simd());
  double a[]={q(c,c_Ham),std::hypot(q(c,c_Mom1),q(c,c_Mom2)),q(c,c_GaussE)};out.write(reinterpret_cast<char*>(a),sizeof(a));
  if(gamma){BoxLoops::loop(GammaCartoonCalculator(r[0]),u,u,Box(c-2*IntVect::Unit,c+2*IntVect::Unit),disable_simd());BoxLoops::loop(Constraints<CouplingFunction>(r[0],CouplingFunction({4*M_PI,0.,0.,-.8}),1.),u,q,Box(c,c),disable_simd());double H=q(c,c_Ham);out.write(reinterpret_cast<char*>(&H),8);}
 }
 return in.eof()&&out ? 0 : 3;
}

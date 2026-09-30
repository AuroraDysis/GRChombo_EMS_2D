#include <cmath>
#include "EMSCouplingFunction.hpp"
#include "T7RHS.hpp"
#include "ConstraintsCartoon.hpp"
#include "BoxLoops.hpp"
#include <fstream>
#include <vector>
int main(int argc,char **argv)
{
 if(argc!=3)return 2;
 std::ifstream in(argv[1],std::ios::binary);std::ofstream out(argv[2],std::ios::binary);
 T7RHS::params_t p{};p.kappa1=.1;p.kappa2=0.;p.kappa3=1.;p.covariantZ4=true;p.lapse_advec_coeff=1.;p.shift_advec_coeff=1.;
 std::vector<double> a(2+28*49);
 while(in.read(reinterpret_cast<char*>(a.data()),a.size()*8))
 {
  int j=std::llround(a[1]/a[0]-.5);IntVect c(D_DECL(0,j,0));Box b(c-3*IntVect::Unit,c+3*IntVect::Unit);FArrayBox u(b,28),rhs(Box(c,c),T7RHS::components),con(Box(c,c),NUM_DIAGNOSTIC_VARS);
  for(int n=0;n<28;++n)for(int y=-3;y<=3;++y)for(int x=-3;x<=3;++x)u(c+IntVect(D_DECL(x,y,0)),n)=a[2+n*49+(y+3)*7+x+3];
  BoxLoops::loop(T7RHS(p,a[0],1.,CouplingFunction({4*M_PI,0.,0.,-20.}),1.,0),u,rhs,Box(c,c),disable_simd());
  BoxLoops::loop(Constraints<CouplingFunction>(a[0],CouplingFunction({4*M_PI,0.,0.,-20.}),1.),u,con,Box(c,c),disable_simd());
  for(int n=0;n<56;++n){double v=rhs(c,n);out.write(reinterpret_cast<char*>(&v),8);}
  double v=con(c,c_Ham);out.write(reinterpret_cast<char*>(&v),8);
 }
 return in.eof()&&out?0:3;
}

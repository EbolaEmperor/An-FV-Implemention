#include <bits/stdc++.h>
#include "FV_MOL.h"
using namespace std;

const double pi = acos(-1);
const double nu = 0.001;

class INITPHI : public TimeFunction2D{
public:
    double at (const double &x, const double &y, const double &t) const{
        static const double t0 = -0.01/(4*nu*log(1e-16));
        static const double rcx = 0.5;
        static const double rcy = 0.75;
        return exp( (-(x-rcx)*(x-rcx)-(y-rcy)*(y-rcy)) / (4*nu*(t+t0)) ) / (t/t0+1);
    }
} initphi;

class FUNCUX : public TimeFunction2D{
public:
    double at (const double &x, const double &y, const double &t) const{
        return 0.1 * sin(pi*x) * sin(pi*x) * sin(2*pi*y);
    }
    double intFixX(const double &x, const double &d, const double &u, const double &t) const{
        return 0.1 * sin(pi*x) * sin(pi*x) * (cos(2*pi*d) - cos(2*pi*u)) / (2*pi);
    }
    double intFixY(const double &y, const double &d, const double &u, const double &t) const{
        return 0.1 * sin(2*pi*y) * (2*pi*(u-d) + sin(2*pi*d) - sin(2*pi*u)) / (4*pi);
    }
} ux;

class FUNCUY : public TimeFunction2D{
public:
    double at (const double &x, const double &y, const double &t) const{
        return -0.1 * sin(2*pi*x) * sin(pi*y) * sin(pi*y);
    }
    double intFixX(const double &x, const double &d, const double &u, const double &t) const{
        return -0.1 * sin(2*pi*x) * (2*pi*(u-d) + sin(2*pi*d) - sin(2*pi*u)) / (4*pi);
    }
    double intFixY(const double &y, const double &d, const double &u, const double &t) const{
        return -0.1 * sin(pi*y) * sin(pi*y) * (cos(2*pi*d) - cos(2*pi*u)) / (2*pi);
    }
} uy;

int main(int argc, char * argv[]){
    FV_MOL_Solver solver;
    solver.setGridSize(stoi(argv[1]));
    solver.setEndTime(10.0);
    solver.setNu(nu);
    solver.setTimeStepWithCaurant(1.0, 0.1, 0.1);
    solver.setNoForcingTerm();
    solver.setInitial(&initphi);
    solver.setPeriodicBondary();
    solver.setVelocity(&ux, &uy);
    solver.solve();
    solver.output("result.txt");
    return 0;
}
#include <bits/stdc++.h>
#include "FV_MOL.h"
using namespace std;

const double pi = acos(-1);
const double nu = 0.01;

class TrueSol{
private:
    int D;
    double u[3], k[3];
    double theta_d(const double &x_d, const int &d, const double &t) const{
        return k[d] * x_d - u[d] * t;
    }
public:
    TrueSol(){
        D = 2; // 2-dimension test
        u[0]=1.0; u[1]=0.5; u[2]=0.25;
        k[0]=2*pi; k[1]=4*pi; k[2]=6*pi;
    }
    double sol(const double &x, const double &y, const double &t){
        return sin(theta_d(x,0,t)) * sin(theta_d(y,1,t));
    }
    double dx(const double &x, const double &y, const double &t){
        return k[0] * cos(theta_d(x,0,t)) * sin(theta_d(y,1,t));
    }
    double dy(const double &x, const double &y, const double &t){
        return k[1] * sin(theta_d(x,0,t)) * cos(theta_d(y,1,t));
    }
    double L (const double &x, const double &y, const double &t) const{
        return nu * (k[0]*k[0] + k[1]*k[1]) * sin(theta_d(x,0,t)) * sin(theta_d(y,1,t))
            + u[0] * (k[0]-1) * cos(theta_d(x,0,t)) * sin(theta_d(y,1,t))
            + u[1] * (k[1]-1) * cos(theta_d(y,1,t)) * sin(theta_d(x,0,t));
    }
} singletruesol;

class F : public TimeFunction2D{
public:
    double at (const double &x, const double &y, const double &t) const{
        return singletruesol.L(x, y, t);
    }
} f;

class Phi : public TimeFunction2D{
public:
    double at (const double &x, const double &y, const double &t) const{
        return singletruesol.sol(x, y, t);
    }
} phi;

class dxPhi : public TimeFunction2D{
public:
    double at (const double &x, const double &y, const double &t) const{
        return singletruesol.dx(x, y, t);
    }
} dxphi;

class dyPhi : public TimeFunction2D{
public:
    double at (const double &x, const double &y, const double &t) const{
        return singletruesol.dy(x, y, t);
    }
} dyphi;

int main(int argc, char * argv[]){
    FV_MOL_Solver solver;
    solver.setGridSize(stoi(argv[1]));
    solver.setEndTime(1.0);
    solver.setNu(nu);
    solver.setTimeStepWithCaurant(1.0, 1.0, 0.5);
    solver.setForcingTerm(&f);
    solver.setInitial(&phi);
    solver.setBondary("down", &phi, "Dirichlet");
    solver.setBondary("left", &phi, "Dirichlet");
    solver.setBondary("up", &dyphi, "Neumann");
    solver.setBondary("right", &dxphi, "Neumann");
    solver.setConstVelocity(1.0, 0.5);
    solver.solve();
    solver.output("result.txt");
    cout << "Error in max-norm: " << solver.checkerr(&phi, Norm_inf()) << endl;
    cout << "Error in 1-norm: " << solver.checkerr(&phi, Norm_p(1)) << endl;
    cout << "Error in 2-norm: " << solver.checkerr(&phi, Norm_p(2)) << endl;
    return 0;
}
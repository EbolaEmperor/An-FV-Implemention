#include <bits/stdc++.h>
#include "INSE.h"
using namespace std;

const double pi = acos(-1);

class FUNCUX : public TimeFunction2D{
private:
    double nu;
public:
    FUNCUX(const double &nu): nu(nu){}
    double at (const double &x, const double &y, const double &t) const{
        return 1.0 - 2.0 * exp(-8.0*pi*pi*nu*t) * cos(2*pi*(x-t)) * sin(2*pi*(y-t));
    }
};

class FUNCUY : public TimeFunction2D{
private:
    double nu;
public:
    FUNCUY(const double &nu): nu(nu){}
    double at (const double &x, const double &y, const double &t) const{
        return 1.0 + 2.0 * exp(-8.0*pi*pi*nu*t) * sin(2*pi*(x-t)) * cos(2*pi*(y-t));
    }
};

class FUNCP : public TimeFunction2D{
private:
    double nu;
public:
    FUNCP(const double &nu): nu(nu){}
    double at (const double &x, const double &y, const double &t) const{
        return -exp(-16.0*pi*pi*nu*t) * ( cos(4*pi*(x-t)) + cos(4*pi*(y-t)) );
    }
};

int main(int argc, char * argv[]){
    double nu = 3.0/stoi(argv[2]);
    TimeFunction2D *u[2];
    u[0] = new FUNCUX(nu);
    u[1] = new FUNCUY(nu);
    TimeFunction2D *p = new FUNCP(nu);

    INSE_Solver solver;
    solver.setGridSize(stoi(argv[1]));
    solver.setEndTime(0.5);
    solver.setNu(nu);
    solver.setTimeStepWithCaurant(stod(argv[3]), 3.0, 3.0);
    solver.setNoForcingTerm();
    solver.setInitial(u);
    solver.setEps(stod(argv[4]));
    solver.solve();
    solver.output("result.txt");
    
    cout << "\t\t\t ux \t\t uy \t\t p" << endl;
    auto err = solver.checkerr(u, p, Norm_inf());
    cout << "Error in max-norm:\t" << err[0] << "\t" << err[1] << "\t" << err[2] << endl;
    err = solver.checkerr(u, p, Norm_p(1));
    cout << "Error in 1-norm:\t" << err[0] << "\t" << err[1] << "\t" << err[2] << endl;
    err = solver.checkerr(u, p, Norm_p(2));
    cout << "Error in 2-norm:\t" << err[0] << "\t" << err[1] << "\t" << err[2] << endl;
    return 0;
}
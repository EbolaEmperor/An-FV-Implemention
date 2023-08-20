#ifndef _FV_MOL_H_
#define _FV_MOL_H_

#include "matrix.h"
#include "sparseMatrix.h"
#include "amgSolver.h"
#include "norm.h"
#include "TimeFunction2D.h"
#include "idpair.h"
#include <cstring>

struct Bondary{
    TimeFunction2D *fun;
    int type;
    // 0:undefined, 1:Dirichlet, 2:Neumann
};

class FV_MOL_Solver{
private:
    int M;
    double tEnd;
    double dH, dT;
    double nu;
    amgSolver solver;
    TimeFunction2D *f, *initial, *ux, *uy;
    // f: forcing term;   initial: initial condition;
    // ux: velocity in x;    uy: velocity in y
    Bondary bon[4];
    // bon[4] are the bondary conditions at: down, up, left, right (in order)
    bool periodicBondary;
    bool noForcingTerm;
    bool uIsConst;
    double constUx, constUy;
    ColVector sol;
    std::vector<Triple> matrixElements;

    // Numerical Integrations
    double simpsonLineInt(TimeFunction2D *g, const idpair &i, const idpair &ed, const double &t);

    // Discrete operators
    ColVector XI(const ColVector &phi, const double &t);
    ColVector XE(const ColVector &phi, const double &t);
    ColVector Ladv(const ColVector &phi, const double &t);
    double F_up(const ColVector &phi, TimeFunction2D *u, const int &i, const int &j, const double &t);
    double F_right(const ColVector &phi, TimeFunction2D *u, const int &i, const int &j, const double &t);
    double Gdp_u_up(TimeFunction2D *u, const int &i, const int &j, const double &t);
    double Gdp_u_right(TimeFunction2D *u, const int &i, const int &j, const double &t);
    double Gdp_phi_up(const ColVector &phi, const int &i, const int &j, const double &t);
    double Gdp_phi_right(const ColVector &phi, const int &i, const int &j, const double &t);
    double facephi_up(const ColVector &phi, const int &i, const int &j, const double &t);
    double facephi_right(const ColVector &phi, const int &i, const int &j, const double &t);

    // if i and j are both in range, return the corresponding value,
    // else return the ghost cell value with corresponding bondaries.
    double getGhost1(const ColVector &phi, const idpair &i, const idpair &ed, const double &t, const Bondary &bon);
    double getGhost2(const ColVector &phi, const idpair &i, const idpair &ed, const double &t, const Bondary &bon);
    double solValue(const ColVector &phi, const int &i, const int &j, const double &t);

    void setGhost1RHS(ColVector &rhs, const int &r, const idpair &i, const idpair &ed, const double &t, const Bondary &bon, const double &k);
    void setGhost2RHS(ColVector &rhs, const int &r, const idpair &i, const idpair &ed, const double &t, const Bondary &bon, const double &k);
    void setRHS(ColVector &rhs, const int &r, const int &i, const int &j, const double &k, const double &t);

    void setGhost1Matrix(const int &r, const idpair &i, const idpair &ed, const Bondary &bon, const double &k);
    void setGhost2Matrix(const int &r, const idpair &i, const idpair &ed, const Bondary &bon, const double &k);
    void setMatrix(const int &r, const int &i, const int &j, const double &k);

    // map the 2D index into 1D index
    int idx(const int &i, const int &j);
    int idx(const idpair &x);

public:
    FV_MOL_Solver();
    void solve();
    void setPeriodicBondary();
    void setNoForcingTerm();
    void setForcingTerm(TimeFunction2D *_f);
    void setInitial(TimeFunction2D *_initial);
    void setBondary(const std::string &bondaryPosition, TimeFunction2D *_bondaryFunction, const std::string &bondaryType);
    void setVelocity(TimeFunction2D *_ux, TimeFunction2D *_uy);
    void setConstVelocity(const double &_ux, const double &_uy);
    void setGridSize(const int &_M);
    void setEndTime(const double &_tEnd);
    void setNu(const double &_nu);
    void setTimeStepWithCaurant(const double &caurant, const double &maxux, const double &maxuy);
    void output(const std::string &outname);
    double checkerr(TimeFunction2D *truesol, const Norm &norm);
};

#endif
#include "FV_MOL.h"
#include "RKTable.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>

//--------------------------------------------Integrations---------------------------------------------------

double FV_MOL_Solver::simpsonLineInt(TimeFunction2D *g, const idpair &i, const idpair &ed, const double &t){
    if(ed[1]==0){
        return g->intFixX(i[0]*dH+0.5*dH*(1+ed[0]), i[1]*dH, (i[1]+1)*dH, t);
    } else {
        return g->intFixY(i[1]*dH+0.5*dH*(1+ed[1]), i[0]*dH, (i[0]+1)*dH, t);
    }
}

//------------------------------------------------------Settings----------------------------------------------------------------

FV_MOL_Solver::FV_MOL_Solver(){
    M = 0;
    dH = dT = 0.0;
    f = nullptr;
    initial = nullptr;
    ux = nullptr;
    uy = nullptr;
    for(int i = 0; i < 4; i++)
        bon[i].fun = nullptr;
    periodicBondary = false;
    noForcingTerm = false;
    uIsConst = false;
}

void FV_MOL_Solver::setPeriodicBondary(){
    periodicBondary = true;
}

void FV_MOL_Solver::setNoForcingTerm(){
    noForcingTerm = true;
}

void FV_MOL_Solver::setForcingTerm(TimeFunction2D *_f){
    f = _f;
}

void FV_MOL_Solver::setInitial(TimeFunction2D *_initial){
    initial = _initial;
}

void FV_MOL_Solver::setBondary(const std::string &bondaryPosition, TimeFunction2D *_bondaryFunction, const std::string &bondaryType){
    std::map< std::string, int > bonid;
    bonid["down"] = 0;
    bonid["up"] = 1;
    bonid["left"] = 2;
    bonid["right"] = 3;
    std::map< std::string, int > bontype;
    bontype["Dirichlet"] = 1;
    bontype["Neumann"] = 2;
    int id = bonid[bondaryPosition];
    bon[id].fun = _bondaryFunction;
    bon[id].type = bontype[bondaryType];
}

void FV_MOL_Solver::setVelocity(TimeFunction2D *_ux, TimeFunction2D *_uy){
    ux = _ux;
    uy = _uy;
    uIsConst = false;
}

void FV_MOL_Solver::setConstVelocity(const double &_ux, const double &_uy){
    constUx = _ux;
    constUy = _uy;
    uIsConst = true;
}

void FV_MOL_Solver::setGridSize(const int &_M){
    M = _M;
    dH = 1.0 / M;
}

void FV_MOL_Solver::setEndTime(const double &_tEnd){
    tEnd = _tEnd;
}

void FV_MOL_Solver::setNu(const double &_nu){
    nu = _nu;
}

void FV_MOL_Solver::setTimeStepWithCaurant(const double &caurant, const double &maxux, const double &maxuy){
    if(dH==0.0){
        std::cerr << "[Error] setTimeStepWithCaurant: Please set grid size first." << std::endl;
    } else {
        dT = caurant / (maxux/dH + maxuy/dH);
    }
}

//---------------------------------------------------Discrete Operators---------------------------------------------------------

double FV_MOL_Solver::facephi_right(const ColVector &phi, const int &i, const int &j, const double &t){
    return 7.0/12.0 * (solValue(phi,i,j,t) + solValue(phi,i+1,j,t)) - 1.0/12.0 * (solValue(phi,i-1,j,t) + solValue(phi,i+2,j,t));
}

double FV_MOL_Solver::facephi_up(const ColVector &phi, const int &i, const int &j, const double &t){
    return 7.0/12.0 * (solValue(phi,i,j,t) + solValue(phi,i,j+1,t)) - 1.0/12.0 * (solValue(phi,i,j-1,t) + solValue(phi,i,j+2,t));
}

double FV_MOL_Solver::Gdp_phi_right(const ColVector &phi, const int &i, const int &j, const double &t){
    return 0.5/dH * (facephi_right(phi,i,j+1,t) - facephi_right(phi,i,j-1,t));
}

double FV_MOL_Solver::Gdp_phi_up(const ColVector &phi, const int &i, const int &j, const double &t){
    return 0.5/dH * (facephi_up(phi,i+1,j,t) - facephi_up(phi,i-1,j,t));
}

double FV_MOL_Solver::Gdp_u_right(TimeFunction2D *u, const int &i, const int &j, const double &t){
    return 0.5/dH * (u->intFixX((i+1)*dH,(j+1)*dH,(j+2)*dH,t) - u->intFixX((i+1)*dH,(j-1)*dH,j*dH,t)) /dH;
}

double FV_MOL_Solver::Gdp_u_up(TimeFunction2D *u, const int &i, const int &j, const double &t){
    return 0.5/dH * (u->intFixY((j+1)*dH,(i+1)*dH,(i+2)*dH,t) - u->intFixY((j+1)*dH,(i-1)*dH,i*dH,t)) /dH;
}

double FV_MOL_Solver::F_right(const ColVector &phi, TimeFunction2D *u, const int &i, const int &j, const double &t){
    return facephi_right(phi,i,j,t) * u->intFixX((i+1)*dH,j*dH,(j+1)*dH,t)/dH + dH*dH/12.0 * Gdp_phi_right(phi,i,j,t) * Gdp_u_right(u,i,j,t);
}

double FV_MOL_Solver::F_up(const ColVector &phi, TimeFunction2D *u, const int &i, const int &j, const double &t){
    return facephi_up(phi,i,j,t) * u->intFixY((j+1)*dH,i*dH,(i+1)*dH,t)/dH + dH*dH/12.0 * Gdp_phi_up(phi,i,j,t) * Gdp_u_up(u,i,j,t);
}

ColVector FV_MOL_Solver::Ladv(const ColVector &phi, const double &t){
    ColVector res(M*M);
    if(uIsConst){
        for(int i = 0; i < M; i++)
            for(int j = 0; j < M; j++){
                res(idx(i,j)) = -constUx/dH * ( 2.0/3.0*solValue(phi,i+1,j,t) - 2.0/3.0*solValue(phi,i-1,j,t) - 1.0/12.0*solValue(phi,i+2,j,t) + 1.0/12.0*solValue(phi,i-2,j,t) )
                                -constUy/dH * ( 2.0/3.0*solValue(phi,i,j+1,t) - 2.0/3.0*solValue(phi,i,j-1,t) - 1.0/12.0*solValue(phi,i,j+2,t) + 1.0/12.0*solValue(phi,i,j-2,t) );
            }
    } else {
        for(int i = 0; i < M; i++)
            for(int j = 0; j < M; j++){
                res(idx(i,j)) = -1.0/dH * ( F_right(phi,ux,i,j,t) - F_right(phi,ux,i-1,j,t) + F_up(phi,uy,i,j,t) - F_up(phi,uy,i,j-1,t) );
            }
    }
    return res;
}

ColVector FV_MOL_Solver::XE(const ColVector &phi, const double &t){
    ColVector res = Ladv(phi,t);
    if(noForcingTerm) return res;
    for(int i = 0; i < M; i++)
        for(int j = 0; j < M; j++)
            res(idx(i,j)) += f->int2D_order6(i*dH, (i+1)*dH, j*dH, (j+1)*dH, t)*M*M;
    return res;
}

ColVector FV_MOL_Solver::XI(const ColVector &phi, const double &t){
    ColVector res(M*M);
    for(int i = 0; i < M; i++)
        for(int j = 0; j < M; j++){
            res(idx(i,j)) = nu/(dH*dH) * ( -1.0/12.0*(solValue(phi,i-2,j,t)+solValue(phi,i+2,j,t)+solValue(phi,i,j-2,t)+solValue(phi,i,j+2,t)) 
                            + 4.0/3.0*(solValue(phi,i-1,j,t)+solValue(phi,i+1,j,t)+solValue(phi,i,j-1,t)+solValue(phi,i,j+1,t)) - 5.0*phi(idx(i,j)));
        }
    return res;
}

//---------------------------------------------------Sol Value and Index--------------------------------------------------------

int FV_MOL_Solver::idx(const int &i, const int &j){
    return i * M + j;
}

int FV_MOL_Solver::idx(const idpair &x){
    return x[0] * M + x[1];
}

double FV_MOL_Solver::getGhost1(const ColVector &phi, const idpair &i, const idpair &ed, const double &t, const Bondary &bon){
    if(bon.type == 1){
        return 1.0/3.0 * (-13.0*phi(idx(i)) + 5.0*phi(idx(i-ed)) - phi(idx(i-2*ed))) + 4.0*simpsonLineInt(bon.fun,i,ed,t)/dH;
    } else if(bon.type == 2){
        return 1.0/10.0 * (5.0*phi(idx(i)) + 9.0*phi(idx(i-ed)) - 5.0*phi(idx(i-2*ed)) + phi(idx(i-3*ed))) + 6.0/5.0*simpsonLineInt(bon.fun,i,ed,t);
    } else {
        std::cerr << "[Error] Undefined bondary." << std::endl;
        exit(-1);
    }
}

double FV_MOL_Solver::getGhost2(const ColVector &phi, const idpair &i, const idpair &ed, const double &t, const Bondary &bon){
    if(bon.type == 1){
        return 1.0/3.0 * (-70.0*phi(idx(i)) + 32.0*phi(idx(i-ed)) - 7.0*phi(idx(i-2*ed))) + 16.0*simpsonLineInt(bon.fun,i,ed,t)/dH;
    } else if(bon.type == 2){
        return 1.0/10.0 * (-75.0*phi(idx(i)) + 145.0*phi(idx(i-ed)) - 75.0*phi(idx(i-2*ed)) + 15.0*phi(idx(i-3*ed))) + 6.0*simpsonLineInt(bon.fun,i,ed,t);
    } else {
        std::cerr << "[Error] Undefined bondary." << std::endl;
        exit(-1);
    }
}

double FV_MOL_Solver::solValue(const ColVector &phi, const int &i, const int &j, const double &t){
    if(periodicBondary){
        return phi(idx( (i+M)%M, (j+M)%M ));
    }
    if(i>=0 && i<M && j>=0 && j<M) return phi(idx(i,j));
    else if(i==M) return getGhost1(phi, idpair(M-1,j), idpair(1,0), t, bon[3]);
    else if(i==M+1) return getGhost2(phi, idpair(M-1,j), idpair(1,0), t, bon[3]);
    else if(i==-1) return getGhost1(phi, idpair(0,j), idpair(-1,0), t, bon[2]);
    else if(i==-2) return getGhost2(phi, idpair(0,j), idpair(-1,0), t, bon[2]);
    else if(j==M) return getGhost1(phi, idpair(i,M-1), idpair(0,1), t, bon[1]);
    else if(j==M+1) return getGhost2(phi, idpair(i,M-1), idpair(0,1), t, bon[1]);
    else if(j==-1) return getGhost1(phi, idpair(i,0), idpair(0,-1), t, bon[0]);
    else if(j==-2) return getGhost2(phi, idpair(i,0), idpair(0,-1), t, bon[0]);
    std::cerr << "[Error] solValue: out of range!" << std::endl;
    exit(-1);
    return -1;
}

void FV_MOL_Solver::setGhost1Matrix(const int &r, const idpair &i, const idpair &ed, const Bondary &bon, const double &k){
    if(bon.type==1){
        matrixElements.push_back( Triple(r, idx(i), -13.0/3.0*k) );
        matrixElements.push_back( Triple(r, idx(i-ed), 5.0/3.0*k) );
        matrixElements.push_back( Triple(r, idx(i-2*ed), -1.0/3.0*k) );
    } else if(bon.type==2){
        matrixElements.push_back( Triple(r, idx(i), 0.5*k) );
        matrixElements.push_back( Triple(r, idx(i-ed), 0.9*k) );
        matrixElements.push_back( Triple(r, idx(i-2*ed), -0.5*k) );
        matrixElements.push_back( Triple(r, idx(i-3*ed), 0.1*k) );
    } else {
        std::cerr << "[Error] Undefined bondary." << std::endl;
        exit(-1);
    }
}

void FV_MOL_Solver::setGhost2Matrix(const int &r, const idpair &i, const idpair &ed, const Bondary &bon, const double &k){
    if(bon.type==1){
        matrixElements.push_back( Triple(r, idx(i), -70.0/3.0*k) );
        matrixElements.push_back( Triple(r, idx(i-ed), 32.0/3.0*k) );
        matrixElements.push_back( Triple(r, idx(i-2*ed), -7.0/3.0*k) );
    } else if(bon.type==2){
        matrixElements.push_back( Triple(r, idx(i), -7.5*k) );
        matrixElements.push_back( Triple(r, idx(i-ed), 14.5*k) );
        matrixElements.push_back( Triple(r, idx(i-2*ed), -7.5*k) );
        matrixElements.push_back( Triple(r, idx(i-3*ed), 1.5*k) );
    } else {
        std::cerr << "[Error] Undefined bondary." << std::endl;
        exit(-1);
    }
}

void FV_MOL_Solver::setMatrix(const int &r, const int &i, const int &j, const double &k){
    if(periodicBondary){
        matrixElements.push_back( Triple(r, idx( (i+M)%M, (j+M)%M ), k) );
        return;
    }
    if(i>=0 && i<M && j>=0 && j<M) matrixElements.push_back(Triple(r, idx(i,j), k));
    else if(i==M) setGhost1Matrix(r, idpair(M-1,j), idpair(1,0), bon[3], k);
    else if(i==M+1) setGhost2Matrix(r, idpair(M-1,j), idpair(1,0), bon[3], k);
    else if(i==-1) setGhost1Matrix(r, idpair(0,j), idpair(-1,0), bon[2], k);
    else if(i==-2) setGhost2Matrix(r, idpair(0,j), idpair(-1,0), bon[2], k);
    else if(j==M) setGhost1Matrix(r, idpair(i,M-1), idpair(0,1), bon[1], k);
    else if(j==M+1) setGhost2Matrix(r, idpair(i,M-1), idpair(0,1), bon[1], k);
    else if(j==-1) setGhost1Matrix(r, idpair(i,0), idpair(0,-1), bon[0], k);
    else if(j==-2) setGhost2Matrix(r, idpair(i,0), idpair(0,-1), bon[0], k);
    else {
        std::cerr << "[Error] setMatrix: out of range!" << std::endl;
        exit(-1);
    }
}

void FV_MOL_Solver::setGhost1RHS(ColVector &rhs, const int &r, const idpair &i, const idpair &ed, const double &t, const Bondary &bon, const double &k){
    if(bon.type==1){
        rhs(r) -= 4.0*k*simpsonLineInt(bon.fun, i, ed, t)/dH;
    } else if(bon.type==2){
        rhs(r) -= 1.2*k*simpsonLineInt(bon.fun, i, ed, t);
    } else {
        std::cerr << "[Error] Undefined bondary." << std::endl;
        exit(-1);
    }
}

void FV_MOL_Solver::setGhost2RHS(ColVector &rhs, const int &r, const idpair &i, const idpair &ed, const double &t, const Bondary &bon, const double &k){
    if(bon.type==1){
        rhs(r) -= 16.0*k*simpsonLineInt(bon.fun, i, ed, t)/dH;
    } else if(bon.type==2){
        rhs(r) -= 6.0*k*simpsonLineInt(bon.fun, i, ed, t);
    } else {
        std::cerr << "[Error] Undefined bondary." << std::endl;
        exit(-1);
    }
}

void FV_MOL_Solver::setRHS(ColVector &rhs, const int &r, const int &i, const int &j, const double &k, const double &t){
    if(periodicBondary || i>=0 && i<M && j>=0 && j<M) return;
    else if(i==M) setGhost1RHS(rhs, r, idpair(M-1,j), idpair(1,0), t, bon[3], k);
    else if(i==M+1) setGhost2RHS(rhs, r, idpair(M-1,j), idpair(1,0), t, bon[3], k);
    else if(i==-1) setGhost1RHS(rhs, r, idpair(0,j), idpair(-1,0), t, bon[2], k);
    else if(i==-2) setGhost2RHS(rhs, r, idpair(0,j), idpair(-1,0), t, bon[2], k);
    else if(j==M) setGhost1RHS(rhs, r, idpair(i,M-1), idpair(0,1), t, bon[1], k);
    else if(j==M+1) setGhost2RHS(rhs, r, idpair(i,M-1), idpair(0,1), t, bon[1], k);
    else if(j==-1) setGhost1RHS(rhs, r, idpair(i,0), idpair(0,-1), t, bon[0], k);
    else if(j==-2) setGhost2RHS(rhs, r, idpair(i,0), idpair(0,-1), t, bon[0], k);
    else {
        std::cerr << "[Error] setRHS: out of range!" << std::endl;
        exit(-1);
    }
}

//---------------------------------------------------output and check error---------------------------------------------------------

void FV_MOL_Solver::output(const std::string &outname){
    std::cout << "--------------------------------------------------------------" << std::endl;
    std::ofstream out(outname);
    out << std::fixed << std::setprecision(16);
    out << sol.reshape(M,M) << std::endl;
    out.close();
    std::cout << "Result has been saved to " << outname << std::endl;
    std::cout << "--------------------------------------------------------------" << std::endl;
}

double FV_MOL_Solver::checkerr(TimeFunction2D *truesol, const Norm &norm){
    ColVector trueres(M*M);
    for(int i = 0; i < M; i++)
        for(int j = 0; j < M; j++)
            trueres(idx(i,j)) = truesol->accInt2D(i*dH, (i+1)*dH, j*dH, (j+1)*dH, tEnd) / (dH*dH);
    return norm(trueres-sol);
}

//-----------------------------------------------------Solve Process----------------------------------------------------------

void FV_MOL_Solver::solve(){
    using namespace ERK_ESDIRK_Table;

    std::cout << "Setting initial values..." << std::endl;
    sol = ColVector(M*M);
    for(int i = 0; i < M; i++)
        for(int j = 0; j < M; j++){
            sol(idx(i,j)) = initial->accInt2D(i*dH, (i+1)*dH, j*dH, (j+1)*dH, 0) * M * M;
        }
    
    ColVector tmpsol[stage], XEtmpsol[stage], XItmpsol[stage];
    for(int s = 0; s < stage; s++){
        tmpsol[s] = ColVector(M*M);
        XEtmpsol[s] = ColVector(M*M);
        XItmpsol[s] = ColVector(M*M);
    }

    const double con1 = 1.0+5.0*dT*gma*nu/(dH*dH);
    const double con2 = dT*gma*nu/(12.0*dH*dH);
    const double con3 = -4.0*dT*gma*nu/(3.0*dH*dH);
    
    matrixElements.clear();
    for(int i = 0; i < M; i++){
        for(int j = 0; j < M; j++){
            setMatrix(i*M+j, i, j, con1);
            setMatrix(i*M+j, i-2, j, con2);
            setMatrix(i*M+j, i+2, j, con2);
            setMatrix(i*M+j, i, j-2, con2);
            setMatrix(i*M+j, i, j+2, con2);
            setMatrix(i*M+j, i-1, j, con3);
            setMatrix(i*M+j, i+1, j, con3);
            setMatrix(i*M+j, i, j-1, con3);
            setMatrix(i*M+j, i, j+1, con3);
        }
    }
    solver.setStrongThereshold(0.05);
    solver.generateGrid(SparseMatrix(M*M, M*M, matrixElements));
    matrixElements.clear();

    for(double t = 0.0; t+1e-12 < tEnd; t += dT){
        // multi stage
        std::cout << "Time: " << t << std::endl;
        tmpsol[0] = sol;
        XEtmpsol[0] = XE(tmpsol[0],t);
        XItmpsol[0] = XI(tmpsol[0],t);
        for(int s = 1; s < stage; s++){
            ColVector rhs = sol;
            double ts = t + c[s]*dT;
            for(int j = 0; j < s; j++){
                rhs = rhs + dT*aE[s][j]*XEtmpsol[j] + dT*aI[s][j]*XItmpsol[j];
            }
            for(int i = 0; i < M; i++){
                for(int j = 0; j < M; j++){
                    if(i>=2 && i<M-2 && j>=2 && j<M-2) continue;
                    setRHS(rhs, i*M+j, i, j, con1, ts);
                    setRHS(rhs, i*M+j, i-2, j, con2, ts);
                    setRHS(rhs, i*M+j, i+2, j, con2, ts);
                    setRHS(rhs, i*M+j, i, j-2, con2, ts);
                    setRHS(rhs, i*M+j, i, j+2, con2, ts);
                    setRHS(rhs, i*M+j, i-1, j, con3, ts);
                    setRHS(rhs, i*M+j, i+1, j, con3, ts);
                    setRHS(rhs, i*M+j, i, j-1, con3, ts);
                    setRHS(rhs, i*M+j, i, j+1, con3, ts);
                }
            }
            tmpsol[s] = solver.solve(rhs, "FMG", 20, 1e-10);
            XEtmpsol[s] = XE(tmpsol[s],ts);
            XItmpsol[s] = XI(tmpsol[s],ts);
        }

        // combine the results
        for(int s = 0; s < stage; s++){
            sol = sol + dT * b[s] * (XEtmpsol[s] + XItmpsol[s]);
        }
    }
    std::cout << "Solved." << std::endl;
}
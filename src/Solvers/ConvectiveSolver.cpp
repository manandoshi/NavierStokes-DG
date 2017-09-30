#include "../../includes/Solvers/ConvectiveSolver.h"
#include "../../includes/Utilities/product.h"

ConvectiveSolver::ConvectiveSolver(int _ne_x, int _ne_y, int _N) {
    ne_x = _ne_x;
    ne_y = _ne_y;
    N = _N;
    time = 0.0;
}

void ConvectiveSolver::initVar(string q) {
    field = new DG_Field_2d(ne_x, ne_y, N, x1, y1, x2, y2);
    field->addVariable_withBounary(q);
    field->addVariable_withBounary("u"+q);
    field->addVariable_withBounary("v"+q);
    field->addVariable_withBounary("d"+q+"dx");
    field->addVariable_withBounary("d"+q+"dy");
    field->addVariable_withoutBounary("dd"+q+"dxdx");
    field->addVariable_withoutBounary("dd"+q+"dydy");
    field->addVariable_withoutBounary("du"+q+"dx");
    field->addVariable_withoutBounary("dv"+q+"dy");
}

void ConvectiveSolver::setDomain(double _x1, double _y1, double _x2, double _y2) {
    x1 = _x1;
    y1 = _y1;
    x2 = _x2;
    y2 = _y2;
    field = new DG_Field_2d(ne_x, ne_y, N, x1, y1, x2, y2);
    this->initVar("T");
    return ;
}

void ConvectiveSolver::setVelocity(function<double(double, double)>U, function<double(double, double)>V) {
    field->addVariable_withBounary("u");
    field->addVariable_withBounary("v");
    field->initializeVariable("u", U);
    field->initializeVariable("v", V);
    return ;
}

void ConvectiveSolver::setDiffusionCoeff(double _nu) {
    nu = _nu;
    return ;
}

void ConvectiveSolver::setInitialConditions(function<double(double, double)> I) {
    field->initializeVariable("T", I);
    return ;
}

void ConvectiveSolver::setBoundaryCondtions(string type) {
    field->setBoundaryConditions(type);
    return ;
}

void ConvectiveSolver::setSolver(double _dt, double _no_of_time_steps, int _record_steps) {
    dt = _dt;
    no_of_time_steps = _no_of_time_steps;
    record_steps = _record_steps;
    return ;
}

void ConvectiveSolver::div(string qx, string qy, 
    double factor, string outp, string fluxType, 
    string fluxVariable_x="", string fluxVariable_y="") {
    
    field->delByDelX(qx, "d"+qx+"dx", fluxType, fluxVariable_x);
    field->delByDelY(qy, "d"+qy+"dy", fluxType, fluxVariable_y);

    field->axpy(factor, "d"+qx+"dx", outp);
    field->axpy(factor, "d"+qy+"dy", outp);
    return ;
}

void ConvectiveSolver::grad(string q, string dqdx, string dqdy) {
    field->delByDelX(q,dqdx,"central","");
    field->delByDelY(q,dqdy,"central","");

}
void ConvectiveSolver::solveAdvDiff(string q, double diff_coeff, string outp) {
    field->setFunctionsForVariables("u", q, product, "u"+q);
    field->setFunctionsForVariables("v", q, product, "v"+q);
    string dqdx = "d"+q+"dx";
    string dqdy = "d"+q+"dy";
    this->div("u"+q, "v"+q,-1.0, outp,"rusanov","u","v");
    this->grad(q,dqdx,dqdy);
    this->div(dqdx,dqdy,1.0*diff_coeff,outp,"central","","");
}

void ConvectiveSolver::solve() {
    field->addVariable_withoutBounary("k1");
    field->addVariable_withoutBounary("k2");
    field->addVariable_withoutBounary("k3");


    // Till now the variable has been initialized.
    // This for-loop is used to march progressively in time. 
    for(int i=0; i < no_of_time_steps; i++) {
        field->scal(0.0, "k1");
        /// First step of RK3
        this->solveAdvDiff("T", nu, "k1");
        
        field->axpy(0.5*dt, "k1", "T");

        ///Second Step of RK3
        field->scal(0.0, "k2");
        this->solveAdvDiff("T", nu, "k2");
        
        field->axpy(-1.5*dt, "k1", "T");
        field->axpy( 2.0*dt, "k2", "T");

        /// Third(&final) step of RK3
        field->scal(0.0, "k3");
        this->solveAdvDiff("T", nu, "k3");
        
        field->axpy( (7.0/6.0)*dt, "k1", "T");
        field->axpy(-(4.0/3.0)*dt, "k2", "T");
        field->axpy( (1.0/6.0)*dt, "k3", "T");
        
        /// RK3 is done, incrementing the time step. 
        time += dt;        

        //Writing VTK file for time Series
        if((i%record_steps) == 0){
            field->writeVTK("output"+to_string(i/10)+".vtk");
        }
    }
}

void ConvectiveSolver::plot(string filename) {
    field->writeVTK(filename);
    return ;
}

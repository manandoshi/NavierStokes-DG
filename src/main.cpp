#include "../includes/Solvers/AdvectionSolver.h"
#include <iostream>
#include <cmath>

using namespace std;

double U(double x, double y) {
    return 1.0;
}

double V(double x, double y) {
    return 0.3;
}

double initial(double x, double y) {
    return (exp(-(x*x +  y*y)*16.0));
}

int main() {
    double dt = 1e-2;
    int time_steps = 100;
    
    //record every 10th time step
    int record_steps = 10;

    AdvectionSolver* a;
    a = new AdvectionSolver(10, 10, 2);
    a->setDomain(-1.0, -1.0, 1.0, 1.0);
    a->setVelocity(U, V);
    a->setInitialConditions(initial);
    a->setBoundaryCondtions("periodic");
    a->setSolver(dt, time_steps, record_steps);
    a->solve();
}

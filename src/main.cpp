#include "../includes/Solvers/ConvectiveSolver.h"
#include <iostream>
#include <cmath>

using namespace std;

double U(double x, double y) {
    return 0.0;
}

double V(double x, double y) {
    return 0.0;
}

double initial(double x, double y) {
    return 0.0;
}

int main() {
    double dt = 1e-2;
    int time_steps = 1000;
    
    //record every 10th time step
    int record_steps = 10;

    ConvectiveSolver* a;
    a = new ConvectiveSolver(20, 20, 2);
    a->setDomain(-1.0, -1.0, 1.0, 1.0);
    a->setVelocity(U, V);
    a->setInitialConditions(initial);
    a->setBoundaryCondtions("dirichlet");
    a->setDiffusionCoeff(1e-2);
    a->setSolver(dt, time_steps, record_steps);
    a->solve();
}

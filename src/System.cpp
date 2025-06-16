#include "System.h"
#include <Eigen/Dense>
#include <cmath>

using namespace Eigen;

// Evaluate f(x) from the SDE dx = f(x)*dt + dw
VectorXd System::dynamics(const VectorXd & x) const
{
    VectorXd f(x.size());
    
    // Extract state variables
    double h = x(0);  // altitude
    double v = x(1);  // velocity
    double c = x(2);  // ballistic drag coefficient

    // // Calculate temperature at altitude h
    double T = T0 - L * h;

    // // Calculate air pressure using the barometric formula
    // double p = p0 * std::pow(1 - L * h / T0, g * M / (R * L));

    // Calculate drag acceleration
    double d = ((0.5 * M * p0) / R) * (1/T) * std::pow(1 - L * h / T0, g * M / (R * L)) * v * v * c;

    // Set f according to the dynamics equations
    f(0) = v;                 // dh/dt = v
    f(1) = d - g;             // dv/dt = d - g
    f(2) = 0;                 // dc/dt = 0 (drag coefficient is constant)

    return f;
}

// Evaluate f(x) and its Jacobian J = df/fx from the SDE dx = f(x)*dt + dw
VectorXd System::dynamics(const VectorXd & x, MatrixXd & J) const
{
    VectorXd f = dynamics(x);

    // Extract state variables
    double h = x(0);  // altitude
    double v = x(1);  // velocity
    double c = x(2);  // ballistic drag coefficient

    // Calculate temperature at altitude h
    double T = T0 - L * h;

    // Calculate air pressure using the barometric formula
    //double p = p0 * std::pow(1 - L * h / T0, g * M / (R * L));

    // Calculate drag acceleration
    double d = ((0.5 * M * p0) / R) * (1/T) * std::pow(1 - L * h / T0, g * M / (R * L)) * v * v * c;

    // Calculate partial derivatives
    double dd_dh = ((0.5 * M * p0) / R) * v * v * c * (L / std::pow(T0 - L * h, 2) * std::pow(1 - L * h / T0, g * M / (R * L)) - (g * M / (R * T0)) * std::pow(1 - L * h / T0, (g * M / (R * L)) - 1) / (T0 - L * h));
    double dd_dv = 2 * d / v;
    double dd_dc = d / c;

    J.resize(f.size(), x.size());
    J << 0,    1,    0,
         dd_dh, dd_dv, dd_dc,
         0,    0,    0;

    return f;
}

// Map x[k] to x[k+1] using RK4
VectorXd System::RungeKutta4thOrder(const VectorXd& x, double dt) {
    VectorXd k1 = dynamics(x);
    VectorXd k2 = dynamics(x + 0.5 * dt * k1);
    VectorXd k3 = dynamics(x + 0.5 * dt * k2);
    VectorXd k4 = dynamics(x + dt * k3);
    return x + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4);    // Return x[k+1]
}
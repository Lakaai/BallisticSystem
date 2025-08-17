#pragma once
#include <Eigen/Dense>

using namespace Eigen;

class System {
    public:
        VectorXd dynamics(const VectorXd & x) const;
        VectorXd dynamics(const VectorXd & x, MatrixXd & J) const;
        VectorXd System::RungeKutta4thOrder(const VectorXd& x, double dt);

    private: 
        const double p0 = 101.325e3;            // Air pressure at sea level [Pa]
        const double M  = 0.0289644;            // Molar mass of dry air [kg/mol]
        const double R  = 8.31447;              // Gas constant [J/(mol.K)]
        const double L  = 0.0065;               // Temperature gradient [K/m]
        const double T0 = 288.15;               // Temperature at sea level [K]
        const double g  = 9.81;                 // Acceleration due to gravity [m/s^2]
};
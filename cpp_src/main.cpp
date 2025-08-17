#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iostream>
#include <Eigen/Dense>
#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>

using namespace Eigen;

class System {
public:
    VectorXd dynamics(const VectorXd& x, double t) {
        double h = x(0); // Altitude
        double v = x(1); // Velocity
        double c = x(2); // Ballistic coefficient
        VectorXd dxdt(3);
        dxdt(0) = v; // dh/dt = v
        dxdt(1) = -9.81 - c * v * v; // dv/dt = -g - drag
        dxdt(2) = 0; // dc/dt = 0
        return dxdt;
    }

    VectorXd RungeKutta4thOrder(const VectorXd& x, double dt) {
        VectorXd k1 = dynamics(x, 0);
        VectorXd k2 = dynamics(x + 0.5 * dt * k1, dt / 2);
        VectorXd k3 = dynamics(x + 0.5 * dt * k2, dt / 2);
        VectorXd k4 = dynamics(x + dt * k3, dt);
        return x + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
    }
};

std::vector<std::vector<double>> getCSVColumns(const std::string& filePath) {
    std::vector<std::vector<double>> columns;
    std::ifstream file(filePath);
    std::string line;
    while (std::getline(file, line)) {
        std::vector<double> row;
        std::stringstream ss(line);
        std::string value;
        while (std::getline(ss, value, ',')) {
            row.push_back(std::stod(value));
        }
        if (columns.empty()) {
            columns.resize(row.size());
        }
        for (size_t i = 0; i < row.size(); ++i) {
            columns[i].push_back(row[i]);
        }
    }
    return columns;
}

VectorXd measurementModel(const VectorXd& x) {
    double h = x(0);
    double r1 = 5000;
    double r2 = 5000;
    VectorXd z(1);
    z(0) = std::sqrt(r1 * r1 + (h - r2) * (h - r2));
    return z;
}

int main() {
    std::string filePath = "/home/luke/Gaimersheim_a2d2/BallisticSystem/data/estimationdata.csv";
    auto columns = getCSVColumns(filePath);

    VectorXd time = Map<VectorXd>(columns[0].data(), columns[0].size());
    VectorXd measured_range = Map<VectorXd>(columns[5].data(), columns[5].size());

    // Initial state and covariance
    int nx = 3;
    MatrixXd S0(nx, nx);
    S0.setZero();
    S0.diagonal() << 2200 * 2200, 100 * 100, 1e-3 * 1e-3;

    VectorXd mu0(nx);
    mu0 << 14000, -450, 0.0005;

    System system;
    MatrixXd sigma = S0;
    VectorXd mu = mu0;
    double previousTime = 0;

    MatrixXd Q(3, 3); // Process noise covariance
    Q.setZero();
    Q(1, 1) = 1e-20;
    Q(2, 2) = 25e-12;

    MatrixXd R(1, 1); // Measurement noise covariance
    R(0, 0) = 50 * 50;

    MatrixXd mu_hist(nx, measured_range.size());
    std::vector<MatrixXd> sigma_hist(measured_range.size(), MatrixXd(nx, nx));

    for (std::size_t k = 0; k < measured_range.size(); ++k) {
        auto currentTime = time(k);
        auto dt = currentTime - previousTime;
        if (dt < 0) {
            std::cerr << "Negative time step at k=" << k << std::endl;
            break;
        }

        VectorXd y_meas(1); y_meas(0) = measured_range(k);

        // Time update
        VectorXd muBar = system.RungeKutta4thOrder(mu, dt);
        Eigen::VectorX<dual> x_dual = mu.cast<dual>();

        auto f = [&](const auto& x) { return system.RungeKutta4thOrder(x, dt); };
        MatrixXd G = jacobian(f, autodiff::wrt(x_dual), autodiff::at(x_dual));
        MatrixXd sigmaBar = G * sigma * G.transpose() + Q;

        // Measurement update
        auto h = [&](const auto& x) { return measurementModel(x); };
        MatrixXd H = jacobian(h, autodiff::wrt(x_dual), autodiff::at(x_dual));
        MatrixXd Kt = sigmaBar * H.transpose() * (H * sigmaBar * H.transpose() + R).inverse();
        mu = muBar + Kt * (y_meas - measurementModel(muBar));
        sigma = (MatrixXd::Identity(nx, nx) - Kt * H) * sigmaBar;

        // Save results
        mu_hist.col(k) = mu;
        sigma_hist[k] = sigma;
        previousTime = currentTime;
    }

    // Save results to CSV
    std::ofstream out("ekf_results.csv");
    out << "time,altitude,velocity,ballistic_coeff\n";
    for (std::size_t k = 0; k < measured_range.size(); ++k) {
        out << time(k) << "," << mu_hist(0, k) << "," << mu_hist(1, k) << "," << mu_hist(2, k) << "\n";
    }
    out.close();

    return 0;
}
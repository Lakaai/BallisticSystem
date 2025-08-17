// #include <autodiff/forward/dual.hpp>
// #include <autodiff/forward/dual/eigen.hpp>
// #include <autodiff/forward/real.hpp>
// #include <autodiff/forward/real/eigen.hpp>
// #include <Eigen/Dense>
// #include <iostream>
// #include "StateEstimator.h"
// #include "System.h"

// using namespace Eigen;
// using namespace std;

// void predict(VectorXd & xk) {
//     return;
// };

// // void affineTransform() {
    
// //     // Compute Jacobian and function value using autodiff
// //     Eigen::Vector2<dual> y;
// //     MatrixXd G = jacobian(RungeKutta4thOrder, autodiff::wrt(x_dual), autodiff::at(x_dual), y);
// //     // MatrixXd F = computeJacobian(x);
// //     // Sigmabar[t] = G[t] * Sigma[t-1] * G[t].transpose + R[t]
// //     auto sigmaBar = G * P * G.transpose() + R;
// //     // K[t] = Sigmabar[t] * H[t].transpose()*(H[t]*sigmabar[t]*H[t].transpose + Q[t])^-1
// //     Kt = sigmaBar * H.transpose() * (H * sigmaBar * H.transpose() + Q).inverse();
// //     mu = muBar + Kt * (y - measurementModel(muBar));
// //     P = (I - Kt * H) * sigmaBar;
// //     return;
// // }

// // Map x[k] to x[k+1] using RK4
// VectorXd RungeKutta4thOrder(const VectorXd& x, double dt) {
//     VectorXd k1 = dynamics(x);
//     VectorXd k2 = dynamics(x + 0.5 * dt * k1);
//     VectorXd k3 = dynamics(x + 0.5 * dt * k2);
//     VectorXd k4 = dynamics(x + dt * k3);
//     return x + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4);    // Return x[k+1]
// }




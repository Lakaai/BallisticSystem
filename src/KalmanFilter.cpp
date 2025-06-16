// // #include "KalmanFilter.h"
// #include <Eigen/Dense>
// #include <cmath>

// KalmanFilter::KalmanFilter() {
//     x_ = Eigen::VectorXd::Zero(6);
//     P_ = Eigen::MatrixXd::Identity(6, 6) * 1.0;
//     A_ = Eigen::MatrixXd::Identity(6, 6);
//     A_(0, 3) = dt_; A_(1, 4) = dt_; A_(2, 5) = dt_;
//     B_ = Eigen::MatrixXd::Zero(6, 2);
//     H_ = Eigen::MatrixXd::Zero(3, 6);
//     H_(0, 0) = 1.0; H_(1, 1) = 1.0; H_(2, 2) = 1.0;
//     Q_ = Eigen::MatrixXd::Identity(6, 6) * 0.1;
//     R_ = Eigen::MatrixXd::Identity(3, 3) * 0.5;
//     psi_ = 0.0; v_ = 0.0;
// }

// void KalmanFilter::predict(const Eigen::Vector2d& u) {
//     double a = u(0), delta = u(1);
//     v_ = std::sqrt(x_(3) * x_(3) + x_(4) * x_(4));
//     psi_ = x_(2);
//     B_(0, 0) = 0.5 * dt_ * dt_ * std::cos(psi_);
//     B_(1, 0) = 0.5 * dt_ * dt_ * std::sin(psi_);
//     B_(3, 0) = dt_ * std::cos(psi_);
//     B_(4, 0) = dt_ * std::sin(psi_);
//     B_(5, 1) = dt_ * (v_ / L_);
//     x_ = A_ * x_ + B_ * u;
//     P_ = A_ * P_ * A_.transpose() + Q_;
// }

// void KalmanFilter::update(const Eigen::Vector3d& z) {
//     Eigen::VectorXd y = z - H_ * x_;
//     Eigen::MatrixXd S = H_ * P_ * H_.transpose() + R_;
//     Eigen::MatrixXd K = P_ * H_.transpose() * S.inverse();
//     x_ = x_ + K * y;
//     P_ = (Eigen::MatrixXd::Identity(6, 6) - K * H_) * P_;
// }

// void KalmanFilter::logState(const Eigen::VectorXd& true_state) {
//     true_states_.push_back(true_state);
//     estimated_states_.push_back(x_);
// }
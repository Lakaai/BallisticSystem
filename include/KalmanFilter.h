#ifndef KALMAN_FILTER_H
#define KALMAN_FILTER_H

#include <Eigen/Dense>
#include <vector>

class KalmanFilter {
public:
    KalmanFilter();
    void predict(const Eigen::Vector2d& u);
    void update(const Eigen::Vector3d& z);

    Eigen::VectorXd getState() const { return x_; }
    Eigen::MatrixXd getCovariance() const { return P_; }
    
    // For simulation logging
    void logState(const Eigen::VectorXd& true_state);
    const std::vector<Eigen::VectorXd>& getTrueStates() const { return true_states_; }
    const std::vector<Eigen::VectorXd>& getEstimatedStates() const { return estimated_states_; }

private:
    Eigen::VectorXd x_; // [x, y, psi, x_dot, y_dot, psi_dot]
    Eigen::MatrixXd P_;
    Eigen::MatrixXd A_, B_, H_, Q_, R_;
    double dt_ = 1.0;
    double L_ = 2.5;
    double psi_, v_;

    // Logging
    std::vector<Eigen::VectorXd> true_states_;
    std::vector<Eigen::VectorXd> estimated_states_;
};

#endif
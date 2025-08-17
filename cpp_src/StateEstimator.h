#pragma once    
#include <Eigen/Dense>
// class SystemEstimator { 
    
//     public: 
//         void SystemEstimator::predict(VectorXd & xk); // Predicts xk+1 

// };

void predict(VectorXd& xk); // Predicts xk+1 
VectorXd RungeKutta4thOrder(const VectorXd& x, double dt);

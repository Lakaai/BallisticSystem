#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iostream>
#include <Eigen/Dense>
#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>
#include "System.h"
#include "ballistic_plot.h"


using namespace Eigen;

struct {                       // Structure declaration
  double measurement_time;     // Member 
  double measurement_value;    // Member
} measurementObject;           // Structure variable

std::vector<std::vector<double>> getCSVColumns(const std::string& filePath); // Forward declaration for helper function
auto measurementModel(VectorXd& x);

int main() {

    std::string filePath = "/home/luke/Gaimersheim_a2d2/BallisticSystem/data/estimationdata.csv";
    auto columns = getCSVColumns(filePath);

    // For each column of columns store it in the appropriate EigenXd Vector
    VectorXd time           = Map<VectorXd>(columns[0].data(), columns[0].size());  // Time
    VectorXd h              = Map<VectorXd>(columns[1].data(), columns[1].size());  // Altitude
    VectorXd v              = Map<VectorXd>(columns[2].data(), columns[2].size());  // Velocity
    VectorXd c              = Map<VectorXd>(columns[3].data(), columns[3].size());  // Ballistic coefficient
    VectorXd true_range     = Map<VectorXd>(columns[4].data(), columns[4].size());  // True range
    VectorXd measured_range = Map<VectorXd>(columns[5].data(), columns[5].size());  // Measured range

    // # Initial state estimate
    // x0 = [initial_altitude, initial_velocity]
    //# Initial covariance
    // p0 = identity_matrix(3) * small_value
    // # Process noise covariance (model uncertainty)
    // Q = process_noise_matrix
    // # Measurement noise covariance (sensor uncertainty)
    // R = measurement_noise_matrix

    // Initial state estimate
    int nx = 3; // Number of states
    MatrixXd S0(nx, nx);
    VectorXd mu0(nx);
    S0.fill(0);
    S0.diagonal() << 2200, 100, 1e-3;

    mu0 << 14000, // Initial height
            -450, // Initial velocity
          0.0005; // Ballistic coefficient
    double previousTime = 0;

    System system; // Construct the system
    using autodiff::dual;
    using autodiff::gradient;
    
    auto sigma = S0;

    Eigen::MatrixXd R(2, 2); // Define the process noise
    R.setZero();

    // Set the non-zero elements of SQ based on the given Q matrix
    R(0, 0) = std::sqrt(1e-20);  // For velocity
    R(1, 1) = std::sqrt(25e-12); // For drag coefficient

    // Indices of process model equations where process noise is injected
    std::vector<Eigen::Index> idxQ;
    idxQ = {1, 2};  // Noise affects velocity (index 1) and drag coefficient (index 2)

    // Gaussian measurement noise of sigmarng = 50 m.
    Eigen::MatrixXd Q(1, 1);
    Q(0, 0) = 50;  

    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();

    const VectorXd mu_hist = {};
    const VectorXd sigma_hist = {};

    for (std::size_t k = 0; k < measured_range.size(); ++k)
    {
    
        auto currentTime = time(k);
        auto dt = currentTime - previousTime;
        auto y = measured_range(k);
        std::cout << "y: " << y << std::endl;

        // Create RADAR measurement

        // Convert input to autodiff type
        Eigen::VectorX<autodiff::dual> x_dual = mu.cast<autodiff::dual>();

        // Process measurement event (do time update and measurement update)
        auto muBar = system.RungeKutta4thOrder(mu, dt); // Predicting forward in time (state prediction)

        auto f = [&](const auto& x) {
            return system.RungeKutta4thOrder(x);
        };

        //     Eigen::Vector2<dual> y;
        MatrixXd G = jacobian(f, autodiff::wrt(x_dual), autodiff::at(x_dual), y);
        
        // Sigmabar[t] = G[t] * Sigma[t-1] * G[t].transpose + Q[t]
        auto sigmaBar = G * sigma * G.transpose() + R; // Predicting forward in time (covariance prediction)

        MatrixXd H = jacobian(measurementModel, autodiff::wrt(x_dual), autodiff::at(x_dual), y);

        // K[t] = Sigmabar[t] * H[t].transpose()*(H[t]*sigmabar[t]*H[t].transpose + R[t])^-1
        auto Kt = sigmaBar * H.transpose() * (H * sigmaBar * H.transpose() + Q).inverse(); // Kalman gain
        auto mu = muBar + Kt * (y - measurementModel(muBar));
        sigma = (I - Kt * H) * sigmaBar;

        // Save results for plotting
        mu_hist.col(k) = mu;
        sigma_hist.col(k) = sigma; 
        previousTime = currentTime; // Update time
    }
}


double measurementModel(VectorXd& x) {
    auto h = x(0);
    double r1 = 5000;
    double r2 = 5000;
    auto range = std::sqrt(r1*r1 + std::pow((h - r2*r2),2));
    return range;
}




// #include "matplotlibcpp.h"
// #include "KalmanFilter.h"
// #include <Eigen/Dense>
// #include <iostream>
// #include <random>
// #include <vector>

// namespace plt = matplotlibcpp;

// int main() {
    
//     KalmanFilter kf;
//     const int steps = 50;
//     const double dt = 1.0;
//     std::random_device rd;
//     std::mt19937 gen(rd());
//     std::normal_distribution<> process_noise(0, 0.1);
//     std::normal_distribution<> meas_noise(0, 0.5);

//     Eigen::VectorXd true_x = Eigen::VectorXd::Zero(6);
//     true_x(3) = 2.0;
//     Eigen::Vector2d u(0.5, 0.05);

//     for (int i = 0; i < steps; ++i) {
//         true_x(0) += true_x(3) * dt + 0.5 * u(0) * dt * dt * cos(true_x(2)) + process_noise(gen);
//         true_x(1) += true_x(4) * dt + 0.5 * u(0) * dt * dt * sin(true_x(2)) + process_noise(gen);
//         true_x(2) += true_x(5) * dt + process_noise(gen);
//         true_x(3) += u(0) * cos(true_x(2)) * dt + process_noise(gen);
//         true_x(4) += u(0) * sin(true_x(2)) * dt + process_noise(gen);
//         true_x(5) += (true_x(3) / (kf.getState()(3) + 1e-6)) * u(1) * dt + process_noise(gen); // Avoid div by zero

//         kf.predict(u);
//         Eigen::Vector3d z;
//         z(0) = true_x(0) + meas_noise(gen);
//         z(1) = true_x(1) + meas_noise(gen);
//         z(2) = true_x(2) + meas_noise(gen);
//         kf.update(z);
//         kf.logState(true_x);
//     }

//     std::vector<double> true_x_pos, true_y_pos, est_x_pos, est_y_pos;
//     const auto& true_states = kf.getTrueStates();
//     const auto& est_states = kf.getEstimatedStates();
//     for (size_t i = 0; i < true_states.size(); ++i) {
//         true_x_pos.push_back(true_states[i](0));
//         true_y_pos.push_back(true_states[i](1));
//         est_x_pos.push_back(est_states[i](0));
//         est_y_pos.push_back(est_states[i](1));
//     }

//     plt::figure();
//     plt::plot(true_x_pos, true_y_pos, "b-");
//     plt::plot(est_x_pos, est_y_pos, "r--");
//     plt::xlabel("X Position (m)");
//     plt::ylabel("Y Position (m)");
//     plt::title("Car Trajectory: True vs Estimated");
//     //plt::legend({"True Trajectory", "Estimated Trajectory"});
//     plt::grid(true);
//     plt::show();

//     return 0;
// }

// #include "matplotlibcpp.h"
// #include "KalmanFilter.h"
// #include <iostream>
// #include <random>
// #include <vector>

// namespace plt = matplotlibcpp;

// int main() {
    
//     KalmanFilter kf;

//     // Simulation parameters
//     const int steps = 50;
//     const double dt = 1.0;
//     std::random_device rd;
//     std::mt19937 gen(rd());
//     std::normal_distribution<> process_noise(0, 0.1);
//     std::normal_distribution<> meas_noise(0, 0.5);

//     // True state (for simulation)
//     Eigen::VectorXd true_x = Eigen::VectorXd::Zero(6);
//     true_x(3) = 2.0; // Initial x_dot = 2 m/s

//     // Control inputs: constant acceleration and slight turn
//     Eigen::Vector2d u(0.5, 0.05); // a = 0.5 m/s^2, delta = 0.05 rad

//     for (int i = 0; i < steps; ++i) {
//         // Simulate true motion (simplified nonlinear model)
//         true_x(0) += true_x(3) * dt + 0.5 * u(0) * dt * dt * cos(true_x(2)) + process_noise(gen);
//         true_x(1) += true_x(4) * dt + 0.5 * u(0) * dt * dt * sin(true_x(2)) + process_noise(gen);
//         true_x(2) += true_x(5) * dt + process_noise(gen);
//         true_x(3) += u(0) * cos(true_x(2)) * dt + process_noise(gen);
//         true_x(4) += u(0) * sin(true_x(2)) * dt + process_noise(gen);
//         true_x(5) += (true_x(3) / kf.getState()(3)) * u(1) * dt + process_noise(gen); // Approx yaw rate

//         // Predict
//         kf.predict(u);

//         // Simulate measurement
//         Eigen::Vector3d z;
//         z(0) = true_x(0) + meas_noise(gen); // x
//         z(1) = true_x(1) + meas_noise(gen); // y
//         z(2) = true_x(2) + meas_noise(gen); // psi
//         kf.update(z);

//         // Log states
//         kf.logState(true_x);
//     }

//     // Extract plotting file
//     std::vector<double> true_x_pos, true_y_pos, est_x_pos, est_y_pos;
//     const auto& true_states = kf.getTrueStates();
//     const auto& est_states = kf.getEstimatedStates();
//     for (size_t i = 0; i < true_states.size(); ++i) {
//         true_x_pos.push_back(true_states[i](0));
//         true_y_pos.push_back(true_states[i](1));
//         est_x_pos.push_back(est_states[i](0));
//         est_y_pos.push_back(est_states[i](1));
//     }

//     // // Plot
//     // plt::figure();
//     // plt::plot(true_x_pos, true_y_pos, "b-", std::map<std::string, std::string>{{"label", "True Trajectory"}});
//     // plt::plot(est_x_pos, est_y_pos, "r--", std::map<std::string, std::string>{{"label", "Estimated Trajectory"}});
//     // plt::xlabel("X Position (m)");
//     // plt::ylabel("Y Position (m)");
//     // plt::title("Car Trajectory: True vs Estimated");
//     // plt::legend();
//     // plt::grid(true);
//     // plt::show();
//     plt::figure();
//     plt::plot(true_x_pos, true_y_pos, "b-");  // No std::map for now
//     plt::plot(est_x_pos, est_y_pos, "r--");
//     plt::xlabel("X Position (m)");
//     plt::ylabel("Y Position (m)");
//     plt::title("Car Trajectory: True vs Estimated");
//     plt::legend({"True Trajectory", "Estimated Trajectory"});  // Use vector of strings
//     plt::grid(true);
//     plt::show();

//     return 0;
// }


// #include "matplotlibcpp.h"
// #include <vector>

// namespace plt = matplotlibcpp;

// int main() {

//     //PyRun_SimpleString("import sys; sys.path.append('/home/luke/Gaimersheim_a2d2/car_example/venv/lib/python3.13/site-packages')");
//     // Simple file: y = x^2
//     std::vector<double> x = {1.0, 2.0, 3.0, 4.0};
//     std::vector<double> y = {1.0, 4.0, 9.0, 16.0};

//     // Plot
//     plt::plot(x, y, "b-");  // Blue line
//     plt::title("Simple Plot: y = x^2");
//     plt::xlabel("X");
//     plt::ylabel("Y");
//     plt::grid(true);
//     plt::show();

//     return 0;
// }


std::vector<std::vector<double>> getCSVColumns(const std::string& filePath) {
    std::ifstream file(filePath);

    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filePath << std::endl;
        assert(false);
    }

    // Use a while loop together with the getline() function to read the file line by line
    // Create a text string, which is used to output the text file
    std::string line;
    std::string value;
    std::vector<std::vector<double>> columns;
    std::getline(file, line); // Skip header by reading in the frist line and not doing anything with it

    while (getline (file, line)) {  //  read a line from the file and store it in the line variable.
        std::stringstream ss(line);     // turn the string into a string stream so we can process each line 
        size_t idx = 0; // tracks which column of the csv file is being indexed, note it is inside the loop so it restarts at the start of each row
        std::cout << line << std::endl;

        while(std::getline(ss, value, ',')) {       //  std::getline(source_stream, destination_variable, delimiter); 
            std::cout << value << std::endl;    // Output the text from the file
            if (columns.size() <= idx) // If this is the first time we've seen this column then create a new vector for it
            { 
                columns.emplace_back(); // Add the new empty column vector
            }

            // Each columns[idx] will hold data from the i-th column of the CSV.
            columns[idx].push_back(std::stod(value)); // Convert the string to a double and store in the current columns vector
            std::cout << idx << std::endl;
            idx++;  // Move to the next columns for the next value in the line
            
        }
    }
    // std::cout << "columns.size:" << columns.size() << std::endl;
    // if (!columns.empty()) {
    //     std::cout << "First column values:" << std::endl;
    //     for (size_t i = 0; i < columns[0].size(); ++i) {
    //         std::cout << columns[0][i] << std::endl;
    //     }
    // }

    // Close the file
    file.close();
    
    return columns; 
}
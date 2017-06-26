#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}


void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
	// state vector
	x_ = x_in;
	// state covariance matrix
    P_ = P_in;
	// state transition matrix
    F_ = F_in;
	// measurement matrix
    H_ = H_in;
	// measurement covariance matrix
    R_ = R_in;
	// process covariance matrix
    Q_ = Q_in;
}

// equations sourced from Sensor Fusion EKF Reference.PDF
void KalmanFilter::Predict() {
	/**
	TODO:
	* predict the state
	*/

    x_ = F_ * x_;
    MatrixXd Ft = F_.transpose();
    P_ = F_ * P_ * Ft + Q_;
}


// equations sourced from Sensor Fusion EKF Reference.PDF
void KalmanFilter::Update(const VectorXd &z) {
	/**
	TODO:
	* update the state by using Kalman Filter equations
	*/
    
    VectorXd y = z - H_ * x_;
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd K =  P_ * Ht * Si;
    
    // New state
    int x_size = x_.size();
    x_ = x_ + (K * y);
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}


// equations sourced from Sensor Fusion EKF Reference.PDF
void KalmanFilter::UpdateEKF(const VectorXd &z) {
	/**
	TODO:
	* update the state by using Extended Kalman Filter equations
	*/

    // Get predicted location in polar coords.
    float px = x_(0);
    float py = x_(1);
    float vx = x_(2);
    float vy = x_(3);
    
    // Check the division by 0.
    float eps = 0.000001;
    if (fabs(px) < eps && fabs(py) < eps) {
        px = eps;
        py = eps;
    }
    else if (fabs(px) < eps) {
        px = eps;
    }
    
    // Recalculate x object state to rho, theta, rho_dot coordinates
    float rho = sqrtf(powf(px, 2) + powf(py, 2));
    float theta = atan2f(py, px);
    float rho_dot = (px * vx + py * vy) / rho;
    
    VectorXd h = VectorXd(3);   // h(x_)
    h << rho, theta, rho_dot;
    
    VectorXd y = z - h;
    
    // In C++, atan2() returns values between -pi and pi.
    // When calculating phi in y = z - h(x) for radar measurements,
    // the resulting angle phi in the y vector should be adjusted so that it is between -pi and pi.
    // The Kalman filter is expecting small angle values between the range -pi and pi.
    // HINT: When working in radians, it's possible to add 2π or subtract 2π until the angle is within the desired range.
    y[1] -= (2 * M_PI) * floor((y[1] + M_PI) / (2 * M_PI));
    
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd K =  P_ * Ht * Si;
    
    // New state
    int x_size = x_.size();
    x_ = x_ + (K * y);
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}
// This part is copied from another student's submission 
// https://github.com/olpotkin/CarND-Extended-Kalman-Filter/blob/master/src/kalman_filter.cpp
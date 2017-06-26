#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;


/*
* Constructor.
*/
FusionEKF::FusionEKF() {
    is_initialized_ = false;
    
	previous_timestamp_ = 0;
    
	// initializing matrices
    R_laser_ = MatrixXd(2, 2);
    R_radar_ = MatrixXd(3, 3);
    H_laser_ = MatrixXd(2, 4);
    Hj_ = MatrixXd(3, 4);

	// Measurement covariance matrix - laser
	R_laser_ << 0.0225, 0,
                0, 0.0225;
	
	// Measurement covariance matrix - radar
	R_radar_ << 0.09, 0, 0,
                0, 0.0009, 0,
                0, 0, 0.09;
    
	/**
	TODO:
	* Finish initializing the FusionEKF.
	* Set the process and measurement noises
	*/
    
    // Measurement matrix - laser
    H_laser_ << 1, 0, 0, 0,
                0, 1, 0, 0;

    // State covariance matrix P_
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ <<  1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1000, 0,
                0, 0, 0, 1000;

    // The initial transition matrix F_
    ekf_.F_ = MatrixXd(4, 4);
    ekf_.F_ <<  1, 0, 1, 0,
                0, 1, 0, 1,
                0, 0, 1, 0,
                0, 0, 0, 1;
  
	// Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
    noise_ax = 9.;
    noise_ay = 9.;
}
/* Code for initializing the FusionEKF is from Lesson 13. Laser Measurements Part 4 */

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}


void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

    /*****************************************************************************
    *  Initialization
    ****************************************************************************/
    if (!is_initialized_) {
		/**
		TODO:
		* Initialize the state ekf_.x_ with the first measurement.
		* Create the covariance matrix.
		* Remember: you'll need to convert radar from polar to cartesian coordinates.
		*/
        // This part is copied from another student's submission 
		// https://github.com/olpotkin/CarND-Extended-Kalman-Filter/blob/master/src/FusionEKF.cpp
        // First measurement
        cout << "EKF: " << endl;
        ekf_.x_ = VectorXd(4);
        ekf_.x_ << 1, 1, 1, 1;

        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
			/**
			Convert radar from polar to cartesian coordinates and initialize state.
			*/
            double rho = measurement_pack.raw_measurements_[0];
            double phi = measurement_pack.raw_measurements_[1];
            double rho_dot = measurement_pack.raw_measurements_[2];
            
            double px = rho * cos(phi);
            double py = rho * sin(phi);
            // So while we can perfectly calculate px and py from phi,
            // we cannot compute vx and vy from phi.
            // We will need yaw (which is introduced in UKF) to compute vx and vy.
            // So even from radar measurement, we can only compute px and py.
            // double vx = rho_dot * sin(phi);
            // double vy = rho_dot * cos(phi);
            double vx = 0;
            double vy = 0;

            ekf_.x_ << px, py, vx, vy;
        }
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
			/**
			Initialize state.
			*/
            ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
        }
        
        // Done initializing, no need to predict or update
        previous_timestamp_ = measurement_pack.timestamp_;
        is_initialized_ = true;
        return;
    }
    
    /*****************************************************************************
     *  Prediction
     ****************************************************************************/

    /**
     TODO:
     * Update the state transition matrix F according to the new elapsed time.
       - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
    */

    // Compute the time elapsed between the current and previous measurements
    
    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; // dt - expressed in seconds
    previous_timestamp_ = measurement_pack.timestamp_;

    ekf_.F_ = MatrixXd(4, 4);
    ekf_.F_ <<  1, 0, dt, 0,
                0, 1, 0, dt,
                0, 0, 1, 0,
                0, 0, 0, 1;
    
    float dt2 = pow(dt, 2);
    float dt3 = pow(dt, 3);
    float dt4 = pow(dt, 4);
    
	//set the process covariance matrix Q
    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ <<  dt4/4*noise_ax, 0, dt3/2*noise_ax, 0,
                0, dt4/4*noise_ay, 0, dt3/2*noise_ay,
                dt3/2*noise_ax, 0, dt2*noise_ax, 0,
                0, dt3/2*noise_ay, 0, dt2*noise_ay;

    // Check if dt is above a certain threshold before predicting.
    // EKF prediction step does not handle dt=0 case gracefully!
    if (dt >= 0.000001) {
        ekf_.Predict();
    }


    /*****************************************************************************
    *  Update
    ****************************************************************************/

    /**
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
    */
	// This part is copied from another student's submission 
	// https://github.com/olpotkin/CarND-Extended-Kalman-Filter/blob/master/src/FusionEKF.cpp
    
    VectorXd x_new = measurement_pack.raw_measurements_;
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        // Radar updates
        if ((x_new(0) != 0) || (x_new(1) != 0)) {
            Hj_ = tools.CalculateJacobian(ekf_.x_);
            ekf_.H_ = Hj_;
            ekf_.R_ = R_radar_;
            ekf_.UpdateEKF(x_new);
        }
    }
    else {
        // Laser updates
        ekf_.R_ = R_laser_;
        ekf_.H_ = H_laser_;
        ekf_.Update(x_new);
    }

    // print the output
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
}

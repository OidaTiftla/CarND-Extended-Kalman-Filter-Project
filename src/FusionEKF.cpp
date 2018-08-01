#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <math.h>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

VectorXd h_function(VectorXd x);
VectorXd h_function_inverse(VectorXd hx);

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
  Hj_radar_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */

  //measurement matrix - laser
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  //state covariance matrix P
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1000, 0,
             0, 0, 0, 1000;

  //the initial transition matrix F_
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
             0, 1, 0, 1,
             0, 0, 1, 0,
             0, 0, 0, 1;

}

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
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      VectorXd x = h_function_inverse(measurement_pack.raw_measurements_);
      ekf_.x_(0) = x(0);
      ekf_.x_(1) = x(1);
      ekf_.x_(2) = 0;
      ekf_.x_(3) = 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      float px = measurement_pack.raw_measurements_(0);
      float py = measurement_pack.raw_measurements_(1);
      ekf_.x_(0) = px;
      ekf_.x_(1) = py;
      ekf_.x_(2) = 0;
      ekf_.x_(3) = 0;
    }

    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
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

  //compute the time elapsed between the current and previous measurements
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  //1. Modify the F matrix so that the time is integrated
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;
  // cout << "F: " << endl << ekf_.F_ << endl;

  //2. Set the process covariance matrix Q
  MatrixXd G = MatrixXd(4, 2);
  G << dt*dt / 2, 0,
       0,         dt*dt / 2,
       dt,        0,
       0,         dt;
  // cout << "G: " << endl << G << endl;
  MatrixXd Gt = G.transpose();

  //set the acceleration noise components
  float noise_ax = 9; // \sigma^2_{a_x}
  float noise_ay = 9; // \sigma^2_{a_y}

  MatrixXd Qv = MatrixXd(2, 2);
  Qv << noise_ax, 0,
        0,        noise_ay;
  // cout << "Qv: " << endl << Qv << endl;
  ekf_.Q_ = G * Qv * Gt;
  // cout << "Q: " << endl << ekf_.Q_ << endl;

  //3. Call the Kalman Filter predict() function
  ekf_.Predict();
  // cout << "x: " << endl << ekf_.x_ << endl;

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.h_ = h_function;
    ekf_.Hj_ = tools.CalculateJacobian(ekf_.x_);
    // cout << "Hj: " << endl << ekf_.Hj_ << endl;
    ekf_.R_ = R_radar_;

    //4. Call the Extended Kalman Filter update() function
    // with the most recent raw measurements_
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;

    //4. Call the Kalman Filter update() function
    // with the most recent raw measurements_
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}

VectorXd h_function(VectorXd x) {
  float ro = sqrt(x(0)*x(0) + x(1)*x(1));
  float theta = atan2(x(1), x(0));
  float ro_dot = (x(0)*x(2) + x(1)*x(3)) / ro;

  VectorXd result(3);
  result << ro, theta, ro_dot;
  return result;
}

VectorXd h_function_inverse(VectorXd hx) {
  float ro = hx(0);
  float theta = hx(1);
  float ro_dot = hx(2);

  float x = ro * cos(theta);
  float y = ro * sin(theta);

  VectorXd result(4);
  result << x, y, 0, 0;
  return result;
}
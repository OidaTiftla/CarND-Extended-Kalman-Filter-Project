#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */

  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if (estimations.size() <= 0
    || estimations.size() != ground_truth.size()) {
    std::cout << "estimations and ground_truth must have at least one element and must have the same size" << std::endl;
    return rmse;
  }

  //accumulate squared residuals
  for (int i=0; i < estimations.size(); ++i) {
    // ... your code here
    VectorXd error = estimations[i] - ground_truth[i];
    error = error.array() * error.array();
    rmse += error;
  }

  //calculate the mean
  // ... your code here
  rmse /= estimations.size();

  //calculate the squared root
  // ... your code here
  rmse = rmse.array().sqrt();

  //return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */

	MatrixXd Hj(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//pre-compute a set of terms to avoid repeated calculation
	float px2_py2 = px * px + py * py;
	float sqrt_px2_py2 = sqrt(px2_py2);
	float vx_py_minus_vy_px = vx * py - vy * px;
	// float sqrt3_px2_py2 = sqrt_px2_py2 * sqrt_px2_py2 * sqrt_px2_py2; // this was my attempt
	// float sqrt3_px2_py2 = sqrt(px2_py2 * px2_py2 * px2_py2); // this was another attempt by me
	float sqrt3_px2_py2 = sqrt_px2_py2 * px2_py2; // this I love the most, but I copied it from the udacity quiz solution

	//check division by zero
	if (fabs(px2_py2) < 0.0001) {
		std::cout << "divition by 0" << std::endl;
		return Hj;
	}

	//compute the Jacobian matrix
	Hj(0, 0) = px / sqrt_px2_py2;
	Hj(0, 1) = py / sqrt_px2_py2;
	Hj(0, 2) = 0;
	Hj(0, 3) = 0;
	Hj(1, 0) = -py / px2_py2;
	Hj(1, 1) = px / px2_py2;
	Hj(1, 2) = 0;
	Hj(1, 3) = 0;
	Hj(2, 0) = (py * vx_py_minus_vy_px) / sqrt3_px2_py2;
	// Hj(2, 1) = (px * (px * vy - py * vx)) / sqrt3_px2_py2; // this line is derived from the solution of the udacity quiz and succeeds in the quiz
	Hj(2, 1) = (px * -vx_py_minus_vy_px) / sqrt3_px2_py2; // this line is mathematically the same as the line above, but udacity's quiz engine does not accept it
	Hj(2, 2) = Hj(0, 0);
	Hj(2, 3) = Hj(0, 1);

	return Hj;
}

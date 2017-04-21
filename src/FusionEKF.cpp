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
  is_x_initialized_ = false;
  is_v_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_radar = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  // laser measurement matrix
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  //create a 4D state vector, we don't know yet the values of the x state
  ekf_.x_ = VectorXd(4);

  //the initial transition matrix F_
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
            0, 1, 0, 1,
            0, 0, 1, 0,
            0, 0, 0, 1;

  //set the acceleration noise components
  noise_ax = 20;
  noise_ay = 20;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  //compute the time elapsed between the current and previous measurements
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;  //dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_v_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
   //the initial state covariance matrix P
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1000, 0,
        0, 0, 0, 1000;


    // first measurement
    cout << "EKF: " << endl;

    VectorXd x_init = VectorXd(2);
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      VectorXd z = measurement_pack.raw_measurements_;
      x_init << z(0)*cos(z(1)), 
                 z(0)*sin(z(1));
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      x_init << measurement_pack.raw_measurements_(0),
               measurement_pack.raw_measurements_(1);
    }
    
    //Initialize velocity if position was already initialized
    if (is_x_initialized_){
      ekf_.x_(2) = (x_init(0) - ekf_.x_(0))/dt;
      ekf_.x_(3) = (x_init(1) - ekf_.x_(1))/dt;
      is_v_initialized_ = true;
    }
    else{
      ekf_.x_ = VectorXd(4);
      ekf_.x_ << 1, 1, 1, 1;
    }

    //(re-)initialize position
    ekf_.x_(0) = x_init(0);
    ekf_.x_(1) = x_init(1);
    is_x_initialized_ = true;
    is_v_initialized_ = true;
    // done initializing, no need to predict or update
    return;
  }
  cout << "predicting: " << endl;

  /*****************************************************************************
  *  Prediction
  ****************************************************************************/
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  //Modify the F matrix so that the time is integrated
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;


  //set the process covariance matrix Q
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
       0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
       dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
       0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    //Calculate Jacobian
    Hj_radar = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_radar;

    //Set Measurement Covariance
    ekf_.R_ = R_radar_;

    //Calculate measurement prediction 
    VectorXd z_pred = tools.CalculateRadar_pred(ekf_.x_)  ;
    VectorXd y = measurement_pack.raw_measurements_ - z_pred;

    // make sure y is between -pi and pi
    while(abs(y(1))>M_PI){
      if (y(1)>M_PI){
        y(1) = y(1)-2*M_PI;
      }
      else{
        y(1) = y(1)+2*M_PI;
      }
    }

    // Call ekf update
    ekf_.Update(y);

  } 
  else {// Laser updates
    //Set H matrix for laser
    ekf_.H_ = H_laser_;

    //Set Measurement Covariance
    ekf_.R_ = R_laser_;

    //Calculate measurement prediction and call ekf update
    VectorXd z_pred = ekf_.H_ * ekf_.x_;
    VectorXd y = measurement_pack.raw_measurements_ - z_pred;
    ekf_.Update(y);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}

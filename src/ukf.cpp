#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.8;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  //set state dimension
  n_x_ = 5;

  //define spreading parameter
  lambda_ = 3 - n_x_;

  //* Augmented state dimension
  n_aug_ = 7;

  ///* predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // set weights
  weights_ = VectorXd(2*n_aug_+1);
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }
  //cout << "constructor UKF done " << Xsig_pred_ << endl;
  //cout << "P" << P_ << endl;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  //cout << " ********** entering ProcessMeasurement "  << endl;
  /**

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    //cout << "EKF: " << endl;
    x_ = VectorXd(n_x_);
    x_ << 0, 0, 0, 0, 0;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float rho = meas_package.raw_measurements_[0]; // range: radial distance from origin
      float phi = meas_package.raw_measurements_[1]; // bearing: angle between rho and x axis
      float rho_dot = meas_package.raw_measurements_[2]; // radial velocity: change of rho
      x_ << rho * cos(phi),
            rho * sin(phi),
            rho_dot,
            phi,
            0;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;

    }

    P_ = MatrixXd::Identity(n_x_, n_x_) ;
    previous_timestamp_ = meas_package.timestamp_;
    // done initializing, no need to predict or update
    //cout << " ************ initialize with first measurement is done" << endl;

    is_initialized_ = true;
    return;
  }

  //compute the time elapsed between the current and previous measurements
  float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = meas_package.timestamp_;

  if(dt < 0.001) {
    return;
  }

  Prediction(dt);


  /**
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_ == true) {
    // Radar updates
    //cout << "FusionEKF: ekf_.UpdateEKF" << endl;
    UpdateRadar(meas_package);
  } else if( use_laser_ == true) {
    // Laser updates
    //cout << "FusionEKF: ekf_.Update" << endl;
    UpdateLidar(meas_package);
  }

  // print the output
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  //cout << "********* entering predict " << endl;
  /********************************************************
  // create the augumented  sigma points. Quiz lessson 7, section 18 Assg #2
  ********************************************************/
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state
  x_aug.fill(0.0);
  x_aug.head(n_x_) = x_;
  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  // construct the Q measurement covariance
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;
  //cout << " *********** done with P_aug initializing " << endl;
  //   //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();
  //   //create augmented sigma points
  Xsig_aug.col(0) = x_aug;

  for (int i = 0; i < n_aug_; i++)
  {
    // cout << "i = " << i << endl;
    Xsig_aug.col(i+1)     = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
  }

  //cout << "********** Done with sigma point generation " << endl;
  ////////////////////////////////
  //     predict sigma points   //
  ////////////////////////////////
  //cout << Xsig_pred_ << endl;
  for(int k = 0; k < 2 * n_aug_ + 1; k++){
    //cout << "k = " << k << endl;
    VectorXd x_k = VectorXd(n_x_);
    VectorXd x_aug_col = VectorXd(n_aug_);

    x_aug_col = Xsig_aug.col(k);
    x_k = x_aug_col.head(5);
    //cout << "x_k " << x_k << endl;
    VectorXd vec1 = VectorXd(n_x_);
    VectorXd vec2 = VectorXd(n_x_);
    vec2 << 0.5 * delta_t * delta_t * cos(x_k(3)) * x_aug_col(5),
            0.5 * delta_t * delta_t * sin(x_k(3)) * x_aug_col(5),
            delta_t * x_aug_col(5),
            0.5 * delta_t * delta_t * x_aug_col(6),
            delta_t * x_aug_col(6);
    if(x_aug(4) != 0){
      vec1 << (x_k(2)/x_aug_col(4)) * ( sin(x_k(3) + x_aug_col(4) * delta_t) - sin(x_k(3))),
              (x_k(2)/x_aug_col(4)) * (-cos(x_k(3) + x_aug_col(4) * delta_t) + cos(x_k(3))),
              0,
              x_k(4) * delta_t,
              0;
    }else{
      vec1 << x_k(2) * cos(x_k(3)) * delta_t,
              x_k(2) * sin(x_k(3)) * delta_t,
              0,
              x_k(4) * delta_t,
              0;
    }
    //cout << "before generating predicted x, P " << endl;
    VectorXd x_sig_pred = VectorXd(n_x_);
    x_sig_pred = x_k + vec1 + vec2;
    //cout << " k = " << k << x_sig_pred << endl;

    // added normalization
    //while (x_sig_pred(3)> M_PI) x_sig_pred(3) -= 2.*M_PI;
    //while (x_sig_pred(3)<-M_PI) x_sig_pred(3) += 2.*M_PI;

    Xsig_pred_.col(k) = x_sig_pred;

    //cout << "done with first sigm prediction " << endl;
  }

  //cout << "Xsig_pred_  \n" << Xsig_pred_ << endl;
  /*****************************
  * Mean and covariance update
  ******************************/
  //create vector for predicted state
  VectorXd x = VectorXd(n_x_);
  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);

  //predicted state mean
  x.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x = x + weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;

    P = P + weights_(i) * x_diff * x_diff.transpose() ;
  }

  x_ = x;
  P_ = P;

  //cout << " ********* Predict is done " << endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  MatrixXd H = MatrixXd(2, 5);
  H << 1, 0, 0, 0, 0,
       0, 1, 0, 0, 0;
  VectorXd pred = H * x_;
  VectorXd y = VectorXd(2);
  y << meas_package.raw_measurements_(0) - x_(0),
       meas_package.raw_measurements_(1) - x_(1);


  // Measurement covariance
  MatrixXd R = MatrixXd(2, 2);
  R << std_laspx_ * std_laspx_, 0,
   		 0, std_laspy_ * std_laspy_;
  MatrixXd Ht = H.transpose();
  MatrixXd S = H * P_ * Ht + R;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;
  //new estimate
	x_ = x_ + (K * y);
  //cout << " step after x_" << endl;
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
  //cout << " before P_" << endl;
	P_ = (I - K * H) * P_;
  //cout << " after P_" << endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_(0),
       meas_package.raw_measurements_(1),
       meas_package.raw_measurements_(2);

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  z_pred.fill(0.0);
  //transform sigma points into measurement space
  for ( int i = 0 ; i < 2 * n_aug_ + 1; i++){
    VectorXd sig_pred_point = VectorXd(n_x_);
    sig_pred_point = Xsig_pred_.col(i);
    float px = sig_pred_point(0);
    float py = sig_pred_point(1);
    float v  = sig_pred_point(2);
    float yaw = sig_pred_point(3);
    float yawd = sig_pred_point(4);
    float ro   = sqrt(px * px + py * py);
    Zsig.col(i) << ro,
                   atan2(py,px),
                   (px * cos(yaw) + py * sin(yaw) ) * v / ro;

  }
  //calculate mean predicted measurement
  for ( int i = 0 ; i < 2 * n_aug_ + 1; i++){
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }
  //calculate measurement covariance matrix S
  MatrixXd R = MatrixXd(n_z,n_z);
  R.fill(0.0);
  R(0,0) = std_radr_ * std_radr_;
  R(1,1) = std_radphi_ * std_radphi_;
  R(2,2) = std_radrd_ * std_radrd_;
  S.fill(0.0);
  for ( int i = 0 ; i < 2 * n_aug_ + 1; i++){
    VectorXd z = VectorXd(n_z);
    z = Zsig.col(i);
    VectorXd diff = z - z_pred;
    //angle normalization
    while (diff(1)> M_PI) diff(1)-=2.*M_PI;
    while (diff(1)<-M_PI) diff(1)+=2.*M_PI;

    S = S + weights_(i) * diff * diff.transpose();
  }
  S = S + R;


  /*******************************************************************************
  * UKF update
  ******************************************************************************/
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for(int i = 0; i < 2 * n_aug_ + 1; i++){
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  //update state mean and covariance matrix
  x_ = x_ + K * (z - z_pred);
  //angle normalization
  while (x_(3)> M_PI) x_(3)-=2.*M_PI;
  while (x_(3)<-M_PI) x_(3)+=2.*M_PI;

  P_ = P_ - K * S * K.transpose();
}

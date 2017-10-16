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
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  if(estimations.size() == 0 || ground_truth.size() == 0){
      return rmse;
  } else if ( estimations.size() != ground_truth.size()){
      cout << "estimations size does not equal to ground_truth size" << endl;
      return rmse;
  }
  //  * the estimation vector size should equal ground truth vector size
  //accumulate squared residuals
  for(int i=0; i < estimations.size(); i++){

    VectorXd est = estimations[i];
    VectorXd tru = ground_truth[i];
    VectorXd diff = est - tru;
    VectorXd dt2 = diff.array() * diff.array();

    rmse += dt2;

  }
  //calculate the mean
  rmse = rmse / estimations.size();

  //calculate the squared root
  rmse = rmse.array().sqrt();
  //return the result
  return rmse;

}

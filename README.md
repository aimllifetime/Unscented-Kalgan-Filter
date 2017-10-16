# Extended Kalman Filter Project Starter Code
Self-Driving Car Engineer Nanodegree Program

In this project a kalman filter is used to estimate the state of a moving object of interest with noisy lidar and radar measurements. 

The files was compiled on the Mac. 

The script **install-mac.sh** was modified to reflect the openSSL version(/usr/local/Cellar/openssl/1.0.2l) on my machine

Instruction to build and run:

1. mkdir build
2. cd build
3. cmake ..
4. make
5. ./ExtendedKF   <= to run the program

Following files have been coded to implement the Kalman filter and extended kalman filter: src/FusionEKF.cpp, src/FusionEKF.h, kalman_filter.cpp, kalman_filter.h, tools.cpp, and tools.h


INPUT: values provided by the simulator to the c++ program

["sensor_measurement"] => the measurement that the simulator observed (either lidar or radar)


OUTPUT: values provided by the c++ program to the simulator

["estimate_x"] <= kalman filter estimated position x
["estimate_y"] <= kalman filter estimated position y
["rmse_x"]
["rmse_y"]
["rmse_vx"]
["rmse_vy"]


Results:

Dataset 1 RMSE:

|               | UKF   |  EKF |
|:-------------:|:-------------:|:-------:|
|  px     | 0.0778        | 0.0973 |
|  py   | 0.0874      | 0.0855|
|  vx   | 0.3541      | 0.4513 |
|  vy     | 0.2646        | 0.4399 | 


![data_set1] (./out/data_set1.png)

Dataset 2 RMSE:

|               | UKF   |  EKF |
|:-------------:|:-------------:|:-------:|
|  px     | 0.0765        | 0.0726 |
|  py   | 0.0793      | 0.0967|
|  vx   | 0.4487      | 0.4579 |
|  vy     | 0.3278        | 0.4966 | 


![data_set1] (./out/data_set2.png)


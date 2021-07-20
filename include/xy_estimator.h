//
// Created by humberto on 6/7/18.
//

#ifndef REEF_ESTIMATOR_XY_ESTIMATOR_H
#define REEF_ESTIMATOR_XY_ESTIMATOR_H

#include "estimator.h"
#include "../../reef_msgs/include/reef_msgs/dynamics.h"

namespace reef_estimator
{
    class XYEstimator : public Estimator
    {
    public:
        XYEstimator();
        ~XYEstimator();
        /*This class is inherited from Estimator, thus
            it has update and propagate methods already.
         */
        ros::NodeHandle nh_;
        ros::Publisher relativeReset_publisher_;

        Eigen::Matrix3d C_body_to_body_fixed_frame;
        void nonlinearPropagation(Eigen::Matrix3d &C,double initialAcc, Eigen::Vector3d accel_in_body, Eigen::Vector3d gyro_in_body, float bias_z_in_NED_component);
        void resetLandingState();
        void relativeReset(Eigen::MatrixXd &xHat, Eigen::MatrixXd &P);

        Eigen::Matrix3d C_body_level_to_body_frame ;
        Eigen::Vector3d nonLinearDynamics;
        Eigen::Matrix2d Id;
        double pitch;
        double roll;
        double yaw;
        double pitch_bias;
        double roll_bias;
        double yaw_bias;
        double pitch_est;
        double roll_est;
        double yaw_est;
        double xc;
        double yc;
        double zc;
        double pe;
        double re;
        double ye;
        geometry_msgs::Quaternion  orient0;
        Eigen::Quaterniond q;
        double phi;
        double theta;
        double psi;
        double distance;
        double deltaX;
        double deltaY;
        double deltaPsi;

        double global_x;
        double global_y;
        double global_yaw;

        geometry_msgs::Vector3 Delta;

        bool XYTakeoff;

        int resetCount;
        double lastYaw;
    };
}



#endif //REEF_ESTIMATOR_XY_ESTIMATOR_H

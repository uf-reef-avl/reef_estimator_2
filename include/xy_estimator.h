//
// Created by humberto on 6/7/18.
//

#ifndef REEF_ESTIMATOR_XY_ESTIMATOR_H
#define REEF_ESTIMATOR_XY_ESTIMATOR_H

#include "estimator.h"
#include "../../reef_msgs/include/reef_msgs/dynamics.h"
#include <std_msgs/Empty.h>


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
        ros::Publisher keyframe_now;

        Eigen::Matrix3d C_body_to_body_fixed_frame;
        void nonlinearPropagation(Eigen::Matrix3d &C,double initialAcc, Eigen::Vector3d accel_in_body, Eigen::Vector3d gyro_in_body, float bias_z_in_NED_component);
        void resetLandingState();
        void relativeReset(Eigen::MatrixXd &xHat, Eigen::MatrixXd &P);

        Eigen::Matrix3d C_body_level_to_body_frame ;
        Eigen::Matrix3d C_level_keyframe_to_body_keyframe_at_k;
        Eigen::Affine3d body_to_camera;
        Eigen::Matrix3d C_body_to_camera;


        Eigen::Vector3d nonLinearDynamics;
        Eigen::Matrix2d Id;
        Eigen::Affine3d global_pose;
        
        // orientation variables 
        double pitch;
        double roll;
        double yaw;
        double pitch_bias;
        double roll_bias;
        double yaw_bias;
        double pitch_est;
        double roll_est;
        double yaw_est;
        
        // variables for simplification of F matrix
        double xc;
        double yc;
        double zc;
        double pe;
        double re;
        double ye;
        
        // true orientation (used for heading)
        geometry_msgs::Quaternion  orient0;
        double phi;
        double theta;
        double psi;
        
        // distance propagated
        double distance;
        
        // accumulated variables to be published
        double global_x;
        double global_y;
        double global_yaw;
        double accum_x_vel;
        double accum_y_vel;
	
	// msg if want to publish deltas
        geometry_msgs::Vector3 Delta;

        bool XYTakeoff;

        int resetCount;
        double lastYaw;
        double lastX;
        double lastY;
        double lastX_vel;
        double lastY_vel;

        double dPoseLimit;
        double dYawLimit;
        double dTimeLimit;
        
        bool want_delta;
        enum StateIndicies
        {   VX,   VY,
            BR, BP,
            BAX, BAY,
            PX,  PY, YAW,
            BGX, BGY,BGZ
        };
    };
}



#endif //REEF_ESTIMATOR_XY_ESTIMATOR_H

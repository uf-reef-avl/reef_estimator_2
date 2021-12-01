//
// Created by humberto on 6/4/18.
//
#ifndef REEF_ESTIMATOR_XYZ_ESTIMATOR_H
#define REEF_ESTIMATOR_XYZ_ESTIMATOR_H

#include <ros/ros.h>
#include <std_msgs/Bool.h>
#include <eigen3/Eigen/Core>
#include <tf2_eigen/tf2_eigen.h>

//Messages
#include <geometry_msgs/PoseStamped.h>
#include <geometry_msgs/TwistStamped.h>
#include <geometry_msgs/TwistWithCovarianceStamped.h>
#include <reef_msgs/XYZDebugEstimate.h>
#include <reef_msgs/XYZEstimate.h>
#include <reef_msgs/DeltaToVel.h>
#include <reef_msgs/ReefMsgsConversionAPI.h>

#include "z_estimator.h"
#include "xy_estimator.h"
#include <reef_msgs/matrix_operation.h>
#include <reef_msgs/dynamics.h>

#include <iostream>

#define ACC_SAMPLE_SIZE 20
#define ACC_TAKEOFF_VARIANCE 0.5
#define SONAR_SAMPLE_SIZE 40
#define SONAR_TAKEOFF_VARIANCE 0.000002
#define SONAR_OFFSET 0.010

namespace reef_estimator
{
    class XYZEstimator
    {
    private:
        ros::NodeHandle private_nh_;
        ros::NodeHandle nh_;
        ros::Publisher state_publisher_;
        ros::Publisher debug_state_publisher_;
        ros::Publisher is_flying_publisher_;
        ros::Publisher pose_publisher_;
        ros::Publisher delta_measurement_publisher_;

        //Estimator enable/disable variables
        bool enableXY;
        bool enableZ;

        //Instantiate a ZEstimator
        ZEstimator zEst;

        //Instantiate an XYEstimator
        XYEstimator xyEst;

        //Message holders
        reef_msgs::XYZEstimate xyzState;
        reef_msgs::XYZDebugEstimate xyzDebugState;
        geometry_msgs::PoseStamped xyzPose;

        int numberOfPropagations;

        //Accelerometer calibration variables
        bool accInitialized;
        int accInitSampleCount;
        double initialAccMagnitude;
        geometry_msgs::Vector3 accSampleAverage;
        double last_time_stamp;
        double mahalanobis_distance_sonar;
        double mahalanobis_distance_rgbd_xy_;
        double mahalanobis_distance_mocap_z;
        double mahalanobis_distance_mocap_xy_;
        double dt;
        bool enable_partial_update;
        int sigma_mulitplier;
        geometry_msgs::PoseStamped delta_msg;


      //Takeoff detection variables
        std_msgs::Bool takeoffState;
        bool accTakeoffState, sonarTakeoffState;
        int numAccSamples, numSonarSamples;
        double accSamples[ACC_SAMPLE_SIZE], sonarSamples[SONAR_SAMPLE_SIZE];
        double accMean, accVariance, sonarMean, sonarVariance;
        bool newDeltaMeasurement, newSonarMeasurement;
        int rgbdCounter;

        Eigen::Matrix3d C_NED_to_body_frame;
        Eigen::Vector3d accelxyz_in_body_frame ;
        Eigen::Vector3d gyroxyz_in_body_frame;
        Eigen::Vector3d accelxyz_in_NED_frame ;
        Eigen::VectorXd range;
        Eigen::VectorXd expected_measurement;
        Eigen::Vector3d Mahalanobis_D_hat_square;
        Eigen::Vector3d Mahalanobis_D_hat;
        Eigen::Vector2d measurement;
        Eigen::Vector2d expected_rgbd;
        Eigen::Affine3d current_estimate;
        Eigen::Affine3d global_pose;
        Eigen::Matrix3d C_level_keyframe_to_body_keyframe_at_k;
        Eigen::Vector3d current_deltaXY;
        Eigen::Matrix3d C3_Yaw;
        double roll_init, pitch_init, yaw_init;
        double global_yaw;


        double beta_0;
        Eigen::MatrixXd xySigma;
        Eigen::Vector3d zSigma;
        
        void checkTakeoffState(double accMagnitude);
        void publishEstimates();
        void saveMinusState();

        void initializeAcc(geometry_msgs::Vector3 acc);
        bool chi2Accept(float range_measurement);

        bool chi2AcceptDeltaPose(geometry_msgs::PoseStamped pose);
        
        bool hypothesisAccept(float range_measurement);
        void integrateGlobalPose();


            public:
        XYZEstimator();
        ~XYZEstimator();

        bool enableMocapXY, enableMocapZ;
        bool enableRGBD, enableSonar;
        bool useMocapXY, useMocapZ;
        bool debug_mode_;

        void sensorUpdate(sensor_msgs::Imu imu);
        void sensorUpdate(sensor_msgs::Range range_msg);
        void deltaPoseUpdate(geometry_msgs::PoseStamped pose_msg);
        

        Eigen::MatrixXd imu_gyro;
        
        Eigen::MatrixXd w_last;
        Eigen::MatrixXd w0;
        Eigen::MatrixXd w;
        
        Eigen::MatrixXd q0;
        Eigen::MatrixXd q;
        Eigen::MatrixXd q_dot;
        
        double q1;	// x
        double q2;	// y
        double q3;	// z
        double q4;	// w

        double psi_imu;
        
        
        void gyro_to_orientation_disc(Eigen::MatrixXd imu_gyro, Eigen::MatrixXd &q, double dt, Eigen::MatrixXd &w_last);
        Eigen::MatrixXd omega_generator(double a, double b, double c);
        double yaw_calc_imu(double q1, double q2, double q3, double q4);


        };

    double getVectorMagnitude(double x, double y, double z);
}

#endif //REEF_ESTIMATOR_XYZ_ESTIMATOR_H

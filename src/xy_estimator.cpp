//
// Created by humberto on 6/7/18.
//
#include "ros/ros.h"
#include "xy_estimator.h"
#include <eigen3/Eigen/Core>
#include <iostream>
#include <reef_msgs/matrix_operation.h>

namespace reef_estimator
{
    XYEstimator::XYEstimator() : Estimator()
    {
        resetCount = 1;

        //Initial state
        F = Eigen::MatrixXd(12,12);
        F.setZero();

        G = Eigen::MatrixXd(12,12);
        G.setZero();

        //P0 is created to save the initial covariance values. It keeps its value forever.
        P0 = Eigen::MatrixXd(12,12);

        //P is the covariance that the filter propagates.
        P = Eigen::MatrixXd(12,12);
        I = Eigen::MatrixXd(12,12);
        I.setIdentity();
        K = Eigen::MatrixXd(12,3);
        K.setZero();

        z = Eigen::MatrixXd(3,1);
        z << 0.0,
             0.0,
             0.0;

        H = Eigen::MatrixXd(3,12);
        H <<  0, 0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 1.0, 0, 0, 0;

        xHat0 = Eigen::MatrixXd(12,1);
        xHat = Eigen::MatrixXd(12,1);
        Q = Eigen::MatrixXd(12,12);

        //R0 represents a pseudo covariance that we use to initiate the propagations,
        //once we are in the air, we need to switch the actual R.
        R0 = Eigen::MatrixXd(3,3);
        R = Eigen::MatrixXd(3,3);

        //Beta values for partial update
        betaVector = Eigen::MatrixXd(12,1);
        distance = 0.;
        phi = 0.;
        theta = 0.;
        psi = 0.;
        reef_msgs::loadTransform("body_to_camera",body_to_camera);
        ROS_WARN_STREAM("[REEF EST]:: Body to camera \n" << body_to_camera.matrix());
        //Initialize the keyframe
        C_body_to_camera = body_to_camera.linear().transpose();

        ROS_WARN_STREAM("[C]:: Body to camera \n" << C_body_to_camera);
        relativeReset_publisher_ = nh_.advertise<geometry_msgs::Vector3>("deltaState", 1, true);
        keyframe_now = nh_.advertise<std_msgs::Empty>("keyframe_now",1);


    }

    XYEstimator::~XYEstimator(){}

    void XYEstimator::nonlinearPropagation(Eigen::Matrix3d &C_NED_to_body_frame, double initialAccMagnitude,Eigen::Vector3d accelxyz_in_body_frame, Eigen::Vector3d gyroxyz_in_body_frame, float bias_z_in_NED_component) {
        //////Save state before propagating.//////
        Eigen::Vector3d prior_state;
        prior_state.x() = xHat(PX,0);
        prior_state.y() = xHat(PY,0);
        prior_state.z() = xHat(YAW);
        //////////////////////////////////////////
        reef_msgs::roll_pitch_yaw_from_rotation321(C_NED_to_body_frame, roll, pitch, yaw);

        roll_bias = xHat(BR);
        pitch_bias = xHat(BP);


        pitch_est = pitch - pitch_bias;
        roll_est = roll - roll_bias;
        //Pitch and roll
                C_body_level_to_body_frame << cos(pitch_est), 0, -sin(pitch_est),
                sin(roll_est) * sin(pitch_est), cos(roll_est), sin(roll_est) * cos(pitch_est),
                cos(roll_est) * sin(pitch_est), -sin(roll_est), cos(roll_est) * cos(pitch_est);

        Eigen::Matrix3d C3;
        C3<< cos(xHat(YAW)), sin(xHat(YAW)), 0,
             -sin(xHat(YAW)), cos(xHat(YAW)), 0,
             0,                           0    ,       1;
//        C_body_level_to_body_frame= C_body_level_to_body_frame*C3;

        Eigen::MatrixXd C_body_level_to_body_frame_2x2(2, 2);
        C_body_level_to_body_frame_2x2 = C_body_level_to_body_frame.block<2, 2>(0, 0);

        //Non-linear dynamics.
        Eigen::Vector3d biasAccel_in_body_frame;
        Eigen::Vector3d bias_z_in_body_frame;
        Eigen::Vector3d bias_z_in_NED;

        bias_z_in_NED << 0, 0, bias_z_in_NED_component;

        bias_z_in_body_frame = C_NED_to_body_frame * bias_z_in_NED;
        biasAccel_in_body_frame << xHat(4), xHat(5), bias_z_in_body_frame(2);

        Eigen::Vector3d gravity_in_NED;
        gravity_in_NED << 0, 0, initialAccMagnitude;

        Eigen::Vector3d gravity_in_body_frame;
        gravity_in_body_frame = C_NED_to_body_frame * gravity_in_NED;

        Eigen::Vector3d input_to_system_in_body_frame;
        input_to_system_in_body_frame = accelxyz_in_body_frame + gravity_in_body_frame + biasAccel_in_body_frame;

        // head(n) means take n elements from the vector starting from index 0. head(2) takes element 0 and 1 from the vector.
        /*
         * C_NED_to_body_frame*gravity_in_NED is included in the nonLinearDynamics to remove the gravity reading.
         *
         */
        nonLinearDynamics = C_body_level_to_body_frame.transpose() * (input_to_system_in_body_frame);
        Eigen::Vector3d xy_time_update;
        xy_time_update = nonLinearDynamics * dt;

        Eigen::Vector3d gyro_in_body_level_frame;
//        gyro_in_body_level_frame = C_body_level_to_body_frame.transpose() *( gyroxyz_in_body_frame + xHat.block<3,1>(BGX,0));

        Eigen::Matrix3d A;//Attitude influence matrix 321ypr (only third row). Kinematic singularity at pitch_est = 90deg
        A << 0,0,0,
             0,0,0,
             0, sin(roll_est)/cos(pitch_est),cos(roll_est)/cos(pitch_est);
        gyro_in_body_level_frame = A*( gyroxyz_in_body_frame + xHat.block<3,1>(BGX,0));

        xHat << xHat(0) + xy_time_update(0),
                xHat(1) + xy_time_update(1),
                xHat(2),
                xHat(3),
                xHat(4),
                xHat(5),
                xHat(6) + (xHat(0) * cos(xHat(8)) - xHat(1) * sin(xHat(8))) * dt,
                xHat(7) + (xHat(0) * sin(xHat(8)) + xHat(1) * cos(xHat(8))) * dt,
                xHat(8) + gyro_in_body_level_frame(2) * dt,
                xHat(9),
                xHat(10),
                xHat(11);


        //Now compute th Jacobian F for the non-linear dynamics. After this, we compute the predict the covariance.
        Eigen::MatrixXd zeros2x2 = Eigen::MatrixXd(2, 2);
        zeros2x2.setZero();
        Eigen::MatrixXd zeros2x1 = Eigen::MatrixXd(2, 1);
        zeros2x1.setZero();
        Eigen::MatrixXd zeros1x2 = Eigen::MatrixXd(1, 2);
        zeros1x2.setZero();
        Eigen::MatrixXd zeros2x3 = Eigen::MatrixXd(2, 3);
        zeros2x3.setZero();
        Eigen::MatrixXd zeros3x12 = Eigen::MatrixXd(3, 12);
        zeros3x12.setZero();
        Eigen::MatrixXd zeros1x3 = Eigen::MatrixXd(1, 3);
        zeros1x3.setZero();
        Eigen::MatrixXd zeros3x2 = Eigen::MatrixXd(3, 2);
        zeros3x2.setZero();
        Eigen::MatrixXd zeros3x3 = Eigen::MatrixXd(3, 3);
        zeros3x3.setZero();
        Eigen::MatrixXd zeros1x1 = Eigen::MatrixXd(1, 1);
        zeros1x1.setZero();


        //xc and yc are just coefficients used in the Jacobian for clarity.
        xc = input_to_system_in_body_frame(0);
        yc = input_to_system_in_body_frame(1);
        zc = input_to_system_in_body_frame(2);
        pe = pitch_est;
        re = roll_est;

        Eigen::MatrixXd partialC_partialBiasAttitude = Eigen::MatrixXd(2,2);
        partialC_partialBiasAttitude << sin(pe) * xc - sin(re) * cos(pe) * yc - cos(re) * cos(pe) * zc, -cos(re) * sin(pe) * yc + sin(re) * sin(pe) * zc,
                                         0.0, sin(re) * yc + cos(re) * zc;

        Eigen::MatrixXd PartialYawPartialBiasAttitude = Eigen::MatrixXd(1,2);
        PartialYawPartialBiasAttitude << -cos(pe)*cos(re)*(gyroxyz_in_body_frame(1) + xHat(BGY)) + sin(re)*cos(pe)*(gyroxyz_in_body_frame(2) +xHat(BGZ)),
        cos(pe)*(gyroxyz_in_body_frame(0) + xHat(BGX)) + sin(re)*sin(pe)*(gyroxyz_in_body_frame(1) +xHat(BGY)) + cos(re)*sin(pe)*(gyroxyz_in_body_frame(2) + xHat(BGZ));

        Eigen::MatrixXd PartialYawPartialGyroBias = Eigen::MatrixXd(1,3);
        PartialYawPartialGyroBias << -sin(pe), sin(re)*cos(pe), cos(re)*cos(pe);

        Eigen::MatrixXd PartialYaw = Eigen::MatrixXd(2, 1);
        PartialYaw << -xHat(0) * sin(xHat(8)) - xHat(1) * cos(xHat(8)),
                xHat(0) * cos(xHat(8)) - xHat(1) * sin(xHat(8));

        Eigen::MatrixXd PartialVel = Eigen::MatrixXd(2, 2);
        PartialVel << cos(xHat(8)), -sin(xHat(8)),
                      sin(xHat(8)), cos(xHat(8));


        F << zeros2x2, partialC_partialBiasAttitude, C_body_level_to_body_frame_2x2.transpose(), zeros2x2, zeros2x1,zeros2x3, //velocity
                zeros2x2, zeros2x2, zeros2x2, zeros2x2, zeros2x1, zeros2x3, //attitude bias
                zeros2x2, zeros2x2, zeros2x2, zeros2x2, zeros2x1, zeros2x3,//accel bias
                PartialVel, zeros2x2, zeros2x2, zeros2x2, PartialYaw, zeros2x3,//position xy
                zeros1x2, PartialYawPartialBiasAttitude, zeros1x2, zeros1x2, zeros1x1, PartialYawPartialGyroBias,//yaw
                zeros3x12;


        G << C_body_level_to_body_frame_2x2.transpose(), zeros2x2, zeros2x2, zeros2x3, zeros2x3,
                zeros2x2, Id.setIdentity(2, 2), zeros2x2, zeros2x3, zeros2x3,
                zeros2x2, zeros2x2, Id.setIdentity(2, 2), zeros2x3, zeros2x3,
                zeros2x2, zeros2x2, zeros2x2, zeros2x3, zeros2x3,
                zeros1x2, zeros1x2, zeros1x2, -sin(pe), sin(re)*cos(pe), cos(re)*cos(pe), zeros1x3,
                zeros3x2, zeros3x2, zeros3x2, zeros3x3, Eigen::MatrixXd::Identity(3,3) ;



        P = P + (F * P.transpose() + P * F.transpose() + G * Q * G.transpose()) * dt;
//        P = F*P*F.transpose() + G*Q*G.transpose();

        //Accumulate the small delta poses

        Eigen::Vector3d propagated_state;
        propagated_state.x() = xHat(PX,0);
        propagated_state.y() = xHat(PY,0);
        propagated_state.z() = xHat(YAW);

        //Now we compute the gain in pose.
        Eigen::Vector3d difference_in_pose;
        difference_in_pose = propagated_state - prior_state;
        Eigen::Affine3d pose_gain;
        pose_gain.linear() = reef_msgs::DCM_from_Euler321(Eigen::Vector3d (0,0,difference_in_pose.z())).transpose();
        pose_gain.translation() << difference_in_pose.x(),difference_in_pose.y(), 0.0; //We do not do altitude deltas. That is in the Z filter.

        global_pose = global_pose*pose_gain;
        global_x = global_pose.translation().x();
        global_y = global_pose.translation().y();
        reef_msgs::get_yaw(global_pose.linear().transpose(), global_yaw);
        if(global_yaw<-2.0*M_PI){
            global_yaw += 2.0*M_PI;
        }
        else if(global_yaw>2.0*M_PI){
            global_yaw -= 2.0*M_PI;
        }

        distance = sqrt(pow(xHat(6), 2) + pow(xHat(7), 2));
        if ((XYTakeoff && (distance > dPoseLimit)) || (XYTakeoff && (xHat(8) > dYawLimit))) {
            //Save attitude at the time of keyframe
            C_level_keyframe_to_body_keyframe_at_k = C_body_level_to_body_frame;
//            XYEstimator::relativeReset(xHat, P);
        }

    }

    void XYEstimator::resetLandingState()
    {
        //Reset covariance P
        P = P0;

        //Reset velocity estimate
        xHat = xHat0;

    }


    void XYEstimator::relativeReset(Eigen::MatrixXd &xHat, Eigen::MatrixXd &P){
        ROS_WARN_STREAM("STATE RESET");
        std_msgs::Empty empty;
        keyframe_now.publish(empty);

        Eigen::Affine3d current_delta;
        current_delta.translation() = Eigen::Vector3d(xHat(PX), xHat(PY),0);
        current_delta.linear() = reef_msgs::DCM_from_Euler321(Eigen::Vector3d(0,0,xHat(YAW))).transpose();
//This line integrates the total delta from the one keyframe to the next.
        //global_pose = global_pose*current_delta;

        // Publish current position and heading to topic to be read from backend compiler here (reset to zero after)
        Delta.x = xHat(PX);
        Delta.y = xHat(PY);
        Delta.z = xHat(YAW);

        relativeReset_publisher_.publish(Delta);

        xHat(PX) = 0.;
        xHat(PY) = 0.;
        xHat(YAW) = 0.;

        //Reset crosscovariances to zero
        P.block<6,6>(PX,VX) = Eigen::MatrixXd::Zero(6,6) ;
        P.block<6,6>(VX,PX) = Eigen::MatrixXd::Zero(6,6);
        P.block<3,3>(PX,BAX) = Eigen::MatrixXd::Zero(3,3);
        P.block<3,3>(BAX,PX) = Eigen::MatrixXd::Zero(3,3);
        //Reset xy and psi covariances to initial values
        P.block<6,6>(PX,PX) = P0.block<6,6>(PX,PX);

        resetCount++;
    }

    /* Implementation of propagate and update is not here, but
     * in the base class.
    */
}

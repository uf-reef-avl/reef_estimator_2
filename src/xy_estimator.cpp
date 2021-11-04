//
// Created by humberto on 6/7/18.
//
#include "ros/ros.h"
#include "xy_estimator.h"
#include <eigen3/Eigen/Core>
#include <iostream>

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
        H <<  0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0;

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

        relativeReset_publisher_ = nh_.advertise<geometry_msgs::Vector3>("deltaState", 1, true);
        global_pose.setIdentity();
    }

    XYEstimator::~XYEstimator(){}

    void XYEstimator::nonlinearPropagation(Eigen::Matrix3d &C_NED_to_body_frame, double initialAccMagnitude,Eigen::Vector3d accelxyz_in_body_frame, Eigen::Vector3d gyroxyz_in_body_frame, float bias_z_in_NED_component) {


        reef_msgs::roll_pitch_yaw_from_rotation321(C_NED_to_body_frame, roll, pitch, yaw);

        roll_bias = xHat(3);
        pitch_bias = xHat(2);


        pitch_est = pitch - pitch_bias;
        roll_est = roll - roll_bias;


        C_body_level_to_body_frame << cos(pitch_est), 0, -sin(pitch_est),
                sin(roll_est) * sin(pitch_est), cos(roll_est), sin(roll_est) * cos(pitch_est),
                cos(roll_est) * sin(pitch_est), -sin(roll_est), cos(roll_est) * cos(pitch_est);


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
        gyro_in_body_level_frame = C_NED_to_body_frame.transpose() * gyroxyz_in_body_frame;

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

        distance = sqrt(pow(xHat(6), 2) + pow(xHat(7), 2));

        if ((XYTakeoff && distance > dPoseLimit) || (XYTakeoff && (xHat(8)) > dYawLimit)) {
            XYEstimator::relativeReset(xHat, P);
        }

    }

    void XYEstimator::resetLandingState()
    {
        //Reset covariance P
        P = P0;
        reef_msgs::roll_pitch_yaw_from_quaternion(orient0,phi, theta, psi);
        roll = phi;
        pitch = theta;

        //Reset velocity estimate
        xHat = xHat0;
        xHat(8) = psi;
        lastYaw = psi;
        lastX = xHat(6);
        lastY = xHat(7);

        // Initial delta test
        global_x = xHat(6);
        global_y = xHat(7);
        global_yaw = xHat(8);

//        XYTakeoff = false;
    }


    void XYEstimator::relativeReset(Eigen::MatrixXd &xHat, Eigen::MatrixXd &P){

        //TODO: Grant/Humberto: Fix this
        std_msgs::Empty empty;
        ros::Publisher keyframe_now = nh_.advertise<std_msgs::Empty>("keyframe_now",1);
        keyframe_now.publish(empty);

        Eigen::Affine3d current_delta;
        current_delta.translation() = Eigen::Vector3d(xHat(PX), xHat(PY),0);
        current_delta.linear() = reef_msgs::DCM_from_Euler321(Eigen::Vector3d(0,0,xHat(YAW))).transpose();

        global_pose = global_pose*current_delta;

        // Publish current position and heading to topic to be read from backend compiler here (reset to zero after)
        Delta.x = xHat(PX);
        Delta.y = xHat(PY);
        Delta.z = xHat(YAW);

//        relativeReset_publisher_.publish(Delta);

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

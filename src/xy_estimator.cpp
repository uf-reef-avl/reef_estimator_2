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
        dt = 0.002;

        //Initial state
        F = Eigen::MatrixXd(9,9);
        F.setZero();

        G = Eigen::MatrixXd(9,9);
        G.setZero();

        //P0 is created to save the initial covariance values. It keeps its value forever.
        P0 = Eigen::MatrixXd(9,9);

        //P is the covariance that the filter propagates.
        P = Eigen::MatrixXd(9,9);
        I = Eigen::MatrixXd(9,9);
        I.setIdentity();
        
        //initiate the Kalman gain variable
        K = Eigen::MatrixXd(9,2);
        K.setZero();
        
        //initiate the measurement matrix
        z = Eigen::MatrixXd(2,1);
        z << 0.0,
             0.0;
             
        H = Eigen::MatrixXd(2,9);
        H <<  0, 0, 0, 0, 0, 0, 1, 0, 0, 
              0, 0, 0, 0, 0, 0, 0, 1, 0;
        
        xHat0 = Eigen::MatrixXd(9,1);
        xHat = Eigen::MatrixXd(9,1);
        Q = Eigen::MatrixXd(9,9);

        //R0 represents a pseudo covariance that we use to initiate the propagations,
        //once we are in the air, we need to switch the actual R.
        R0 = Eigen::MatrixXd(2,2);
        R = Eigen::MatrixXd(2,2);

        //Beta values for partial update
        betaVector = Eigen::MatrixXd(9,1);

        // intiate the propagated distance, if a certain distance is reached then a reset is initiated
        distance = 0.;
	
	// initiate the roll, pitch, and yaw
        phi = 0.;
        theta = 0.;
        psi = 0.;
	
	// if time is a desired threshold get the time
        // previousTime = ros::Time::now().toSec();
	
        relativeReset_publisher_ = nh_.advertise<geometry_msgs::Vector3>("deltaState", 1, true);
        
        // flag to get delta measurement
        want_delta = false;
    }





    XYEstimator::~XYEstimator(){}





    void XYEstimator::nonlinearPropagation(Eigen::Matrix3d &C_NED_to_body_frame, double initialAccMagnitude,Eigen::Vector3d accelxyz_in_body_frame, Eigen::Vector3d gyroxyz_in_body_frame, float bias_z_in_NED_component) {

	if(XYTakeoff){
	// get the orientation in rpy from the imu
        reef_msgs::roll_pitch_yaw_from_rotation321(C_NED_to_body_frame, roll, pitch, yaw);

	// get the orientation biases from the estimated state
        roll_bias = xHat(3);
        pitch_bias = xHat(2);

	// orientation angles are measured minus the estimated bias
        pitch_est = pitch - pitch_bias;
        roll_est = roll - roll_bias;

	// create the body level to body frame rotation matrix with the estimated roll and pitch
        C_body_level_to_body_frame << cos(pitch_est), 0, -sin(pitch_est),
                sin(roll_est) * sin(pitch_est), cos(roll_est), sin(roll_est) * cos(pitch_est),
                cos(roll_est) * sin(pitch_est), -sin(roll_est), cos(roll_est) * cos(pitch_est);


        Eigen::MatrixXd C_body_level_to_body_frame_2x2(2, 2);
        C_body_level_to_body_frame_2x2 = C_body_level_to_body_frame.block<2, 2>(0, 0);

        //Non-linear dynamics.
        Eigen::Vector3d biasAccel_in_body_frame;
        Eigen::Vector3d bias_z_in_body_frame;
        Eigen::Vector3d bias_z_in_NED;

	// set up the z bias in NED (from z estimator)
        bias_z_in_NED << 0, 0, bias_z_in_NED_component;

	// rotate z bias to body frame
        bias_z_in_body_frame = C_NED_to_body_frame * bias_z_in_NED;
        
        // vector with x, y, and z accelerometer biases
        biasAccel_in_body_frame << xHat(4), xHat(5), bias_z_in_body_frame(2);
	
	// gravity vector in NED
        Eigen::Vector3d gravity_in_NED;
        gravity_in_NED << 0, 0, initialAccMagnitude;

	// rotate gravity to body frame
        Eigen::Vector3d gravity_in_body_frame;
        gravity_in_body_frame = C_NED_to_body_frame * gravity_in_NED;

	// imu acceleration measuremnt accounting for gravity and accelerometer biases
        Eigen::Vector3d input_to_system_in_body_frame;
        input_to_system_in_body_frame = accelxyz_in_body_frame + gravity_in_body_frame + biasAccel_in_body_frame;

        // head(n) means take n elements from the vector starting from index 0. head(2) takes element 0 and 1 from the vector.
        /*
         * C_NED_to_body_frame*gravity_in_NED is included in the nonLinearDynamics to remove the gravity reading.
         *
         */
         
        // rotate the input to the system (acceleration) from body frame to body level frame
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
                //xHat(6) + (xHat(0) * cos(xHat(8)) - xHat(1) * sin(xHat(8))) * dt,
                //xHat(7) + (xHat(0) * sin(xHat(8)) + xHat(1) * cos(xHat(8))) * dt,
                xHat(6) + (accum_x_vel * cos(xHat(8)) - accum_y_vel * sin(xHat(8))) * dt,
                xHat(7) + (accum_x_vel * sin(xHat(8)) + accum_y_vel * cos(xHat(8))) * dt,
                xHat(8) + gyro_in_body_level_frame(2) * dt;


        //Now compute th Jacobian F for the non-linear dynamics. After this, we compute the predict the covariance.
        Eigen::Matrix2d zeros2x2;
        zeros2x2.setZero();
        Eigen::MatrixXd zeros2x1 = Eigen::MatrixXd(2, 1);
        zeros2x1.setZero();
        Eigen::MatrixXd zeros1x2 = Eigen::MatrixXd(1, 2);
        zeros1x2.setZero();

        //xc and yc are just coefficients used in the Jacobian for clarity.
        xc = input_to_system_in_body_frame(0);
        yc = input_to_system_in_body_frame(1);
        zc = input_to_system_in_body_frame(2);
        pe = pitch_est;
        re = roll_est;

        Eigen::Matrix2d partialC_partialBiasAttitude;
        partialC_partialBiasAttitude << sin(pe) * xc - sin(re) * cos(pe) * yc - cos(re) * cos(pe) * zc, -cos(re) * sin(pe) * yc + sin(re) * sin(pe) * zc,
                                        0,                                                               sin(re) * yc + cos(re) * zc;


        Eigen::MatrixXd PartialYaw = Eigen::MatrixXd(2, 1);

        PartialYaw << -xHat(0) * sin(xHat(8)) - xHat(1) * cos(xHat(8)),
                       xHat(0) * cos(xHat(8)) - xHat(1) * sin(xHat(8));

        Eigen::Matrix2d PartialVel;
        PartialVel << cos(xHat(8)), -sin(xHat(8)),
                      sin(xHat(8)), cos(xHat(8));

        F <<    zeros2x2, partialC_partialBiasAttitude, C_body_level_to_body_frame_2x2.transpose(), zeros2x2, zeros2x1,
                zeros2x2, zeros2x2, zeros2x2, zeros2x2, zeros2x1,
                zeros2x2, zeros2x2, zeros2x2, zeros2x2, zeros2x1,
                PartialVel, zeros2x2, zeros2x2, zeros2x2, PartialYaw,
                zeros1x2, zeros1x2, zeros1x2, zeros1x2, 0;


        G << C_body_level_to_body_frame_2x2.transpose(), zeros2x2, zeros2x2, zeros2x2, zeros2x1,
                zeros2x2, Id.setIdentity(2, 2), zeros2x2, zeros2x2, zeros2x1,
                zeros2x2, zeros2x2, Id.setIdentity(2, 2), zeros2x2, zeros2x1,
                zeros2x2, zeros2x2, zeros2x2, C_body_level_to_body_frame_2x2.transpose() * PartialVel, zeros2x1,
                zeros1x2, zeros1x2, zeros1x2, zeros1x2, 1;


        P = P + (F * P.transpose() + P * F.transpose() + G * Q * G.transpose()) * dt;
//        P = F*P*F.transpose() + G*Q*G.transpose();


	// calculate the distance traveled, this is off of the relative reset value of zero
	// in other words it is the delta pose from the last reset
        distance = sqrt(pow(xHat(6), 2) + pow(xHat(7), 2));
        
        // get current time if you want time as a threshold
        // currentTime = ros::Time::now().toSec();

	// initiate a relative reset if a certain distance from the last reset is reached or if rotated by a certain amount
        if (distance > dPoseLimit) {
            XYEstimator::relativeReset(xHat, P);
            want_delta = true;
        } else if ((xHat(8) - lastYaw) > dYawLimit) {
            XYEstimator::relativeReset(xHat, P);
            want_delta = true;
        }
        
  /* TIME CONSTRAINT: shouldn't be necessary but can implement, error with current version so left out */
//         else if(XYTakeoff && (currentTime - previousTime) > dTimeLimit){
//             XYEstimator::relativeReset(xHat, P);
//             previousTime = ros::Time::now().toSec();
//         }

	    
	// compute the global position - the position to be published (accumulated velocity too)
	// add on [xhat(n) - last] in order to continuously propagate and publish the position, after this is added on the *last variable is set to the current value
	accum_x_vel = accum_x_vel + xHat(0) - lastX_vel;
	accum_y_vel = accum_y_vel + xHat(1) - lastY_vel;
        global_x = global_x + xHat(6) - lastX;
        global_y = global_y + xHat(7) - lastY;
        global_yaw = global_yaw + xHat(8) - lastYaw;

	lastX_vel = xHat(0);
	lastY_vel = xHat(1);
        lastX = xHat(6);
        lastY = xHat(7);
        lastYaw = xHat(8);
        
        
        
        
        // debug
        //ROS_INFO_STREAM("\n X velocity: " << accum_x_vel + xHat(0) - lastX_vel << "\n");
        
	}
	else{
	global_x = xHat0(6);
	global_y = xHat0(7);
	}
	
    }





    void XYEstimator::resetLandingState()
    {
        //Reset covariance P
        P = P0;
        reef_msgs::roll_pitch_yaw_from_quaternion(orient0,phi, theta, psi);
        roll = phi;
        pitch = theta;
        

        //Reset estimate
        xHat = xHat0;
        xHat(8) = psi;
        
        lastYaw = psi;
        lastX = xHat(6);
        lastY = xHat(7);
        lastX_vel = xHat(0);
        lastY_vel = xHat(1);

        // Initial accumulated variable
        global_x = xHat(6);
        global_y = xHat(7);
        global_yaw = xHat(8);
        accum_x_vel = xHat(0);
        accum_y_vel = xHat(1);


    }





    void XYEstimator::relativeReset(Eigen::MatrixXd &xHat, Eigen::MatrixXd &P){

        // Publish current position and heading to topic to be read from backend compiler here (reset to zero after)
//        Delta.x = xHat(6);
//        Delta.y = xHat(7);
//        Delta.z = xHat(8) - lastYaw;
//        relativeReset_publisher_.publish(Delta);
        
        
        // get the true heading angle
        reef_msgs::roll_pitch_yaw_from_quaternion(orient0,phi, theta, psi);

	// calculate the last accumulated variables before the reset
	//accum_x_vel = accum_x_vel + xHat(0) - lastX_vel;
	//accum_y_vel = accum_y_vel + xHat(1) - lastY_vel;
        //global_x = global_x + xHat(6) - lastX;
        //global_y = global_y + xHat(7) - lastY;
        //global_yaw = psi;


        // reset the relative states (all for error states estimator) to 0 and the heading to the true value
        // do the same for the last variables used for continuous propagation
        xHat(0) = 0.;
        xHat(1) = 0.;
        xHat(2) = xHat0(2);
        xHat(3) = xHat0(3);
        xHat(4) = xHat0(4);
        xHat(5) = xHat0(5);
        xHat(6) = 0.;
        xHat(7) = 0.;
        xHat(8) = psi;
        
        lastX_vel = 0.;
        lastY_vel = 0.;
        lastYaw = psi;
        lastX = 0.;
        lastY = 0.;
        
        // set covariances back to initial values
        P = P0;
        
        
//        xHat(6) = 0.;
//        xHat(7) = 0.;
//        lastX = 0.;
//        lastY = 0.;
        
        // set the covariances back to their initial values for the relative states
//        Eigen::MatrixXd P6 = Eigen::MatrixXd(1,9);
//        Eigen::MatrixXd P7 = Eigen::MatrixXd(1,9);
//        Eigen::MatrixXd P8 = Eigen::MatrixXd(1,9);
//        P6 << 0, 0, 0, 0, 0, 0, P0(6,6), 0, 0;
//        P7 << 0, 0, 0, 0, 0, 0, 0, P0(7,7), 0;
//        P8 << 0, 0, 0, 0, 0, 0, 0, 0, P0(8,8);

//        P.block<1,9>(6,0) = P6;
//        P.block<1,9>(7,0) = P7;
//       P.block<1,9>(8,0) = P8;

       resetCount++;
    }

    /* Implementation of propagate and update is not here, but
     * in the base class.
    */
}

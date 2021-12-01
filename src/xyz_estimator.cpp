//
// Created by humberto on 6/4/18.
//

#include "xyz_estimator.h"
#include <iostream>

namespace reef_estimator
{
    XYZEstimator::XYZEstimator() : private_nh_("~"),
                                   nh_(""),
                                   numberOfPropagations(0),
                                   accInitialized(false),
                                   accInitSampleCount(0),
                                   numAccSamples(0),
                                   numSonarSamples(0),
                                   accMean(0),
                                   accVariance(0),
                                   sonarMean(0),
                                   sonarVariance(0),
                                   sonarTakeoffState(false),
                                   accTakeoffState(false),
                                   newSonarMeasurement(false),
                                   newDeltaMeasurement(false),
                                   rgbdCounter(0)
    {
        private_nh_.param<bool>("debug_mode", debug_mode_, false);
        if (debug_mode_) {
            ROS_WARN_STREAM("Debug Mode Enabled");
        }

        private_nh_.param<bool>("enable_xy", enableXY, true);
        if (!enableXY) {
            ROS_ERROR_STREAM("XY Estimator Disabled");
        }

        private_nh_.param<bool>("enable_z", enableZ, true);
        if (!enableZ) {
            ROS_ERROR_STREAM("Z Estimator Disabled");
        }

        //Mocap override settings
        private_nh_.param<bool>("enable_mocap_xy", enableMocapXY, true);
        private_nh_.param<bool>("enable_rgbd", enableRGBD, true);
        useMocapXY = enableMocapXY && !enableRGBD;
        if (!enableRGBD) {
            ROS_WARN("RGBD feedback disabled.");
        }

        private_nh_.param<bool>("enable_mocap_z", enableMocapZ, true);
        private_nh_.param<bool>("enable_sonar", enableSonar, true);
        useMocapZ = enableMocapZ && !enableSonar;
        if (!enableSonar) {
            ROS_WARN("Sonar feedback disabled.");
        }

        private_nh_.param<double>("mahalanobis_d_sonar", mahalanobis_distance_sonar, 20);
        private_nh_.param<double>("mahalanobis_d_rgbd_velocity", mahalanobis_distance_rgbd_xy_, 20);
        private_nh_.param<double>("mahalanobis_d_mocap_z", mahalanobis_distance_mocap_z, 20);
        private_nh_.param<double>("mahalanobis_distance_mocap_xy_", mahalanobis_distance_mocap_xy_, 20);
        private_nh_.param<bool>("enable_partial_update", enable_partial_update, true);

        ROS_WARN_STREAM("Mahalanobis Distance for Sonar " << mahalanobis_distance_sonar);
        ROS_WARN_STREAM("Mahalanobis Distance for RGBD Velocity " << mahalanobis_distance_rgbd_xy_);
        ROS_WARN_STREAM("Mahalanobis Distance for XYYaw Delta Mocap " << mahalanobis_distance_mocap_xy_);
        ROS_WARN_STREAM("Mahalanobis Distance for Mocap Z " << mahalanobis_distance_mocap_z);
        ROS_WARN_STREAM("Mahalanobis Distance for Mocap Velocity "  << mahalanobis_distance_mocap_xy_);
        // Initialize dt
        private_nh_.param<double>("estimator_dt", dt, 0.004);
        xyEst.dt = zEst.dt = dt;

        //Initialize estimators with parameters
        reef_msgs::importMatrixFromParamServer(private_nh_, xyEst.xHat0, "xy_x0");
        reef_msgs::importMatrixFromParamServer(private_nh_, xyEst.P0, "xy_P0");
        reef_msgs::importMatrixFromParamServer(private_nh_, xyEst.Q, "xy_Q");
        reef_msgs::importMatrixFromParamServer(private_nh_, xyEst.R0, "xy_R0");
        reef_msgs::importMatrixFromParamServer(private_nh_, xyEst.betaVector, "xy_beta");
//        xyEst.Q *= (xyEst.dt);

        //Get the constraint to implement relative reset from yaml file
        private_nh_.param<double>("delta_pose_limit", xyEst.dPoseLimit, 0.5);
        private_nh_.param<double>("delta_yaw_limit", xyEst.dYawLimit, 0.7);
        private_nh_.param<double>("delta_time_limit", xyEst.dTimeLimit, 100.0);


        //Initialization of the first frame with MOCAP
        std::string mocapPoseTopic;
        nh_.param<std::string>("mocap_pose_topic", mocapPoseTopic, "mocap_ned");
        geometry_msgs::PoseStampedConstPtr pose_msg;
        pose_msg = ros::topic::waitForMessage<geometry_msgs::PoseStamped>("/pose_stamped", nh_);
        if (pose_msg != NULL) {
            ROS_INFO_STREAM("GOT MY FIRST MOCAP MESSAGE");
        }
        reef_msgs::roll_pitch_yaw_from_quaternion(pose_msg->pose.orientation,roll_init, pitch_init, yaw_init);
        //Save the initial level keyframe
        global_yaw  = yaw_init;
        global_pose.linear() = reef_msgs::DCM_from_Euler321(Eigen::Vector3d (0,0,yaw_init)).transpose();
        global_pose.translation() = Eigen::Vector3d (pose_msg->pose.position.x,pose_msg->pose.position.y,0.0);
        xyEst.XYTakeoff = false;
        //Initialize the keyframe
        C_level_keyframe_to_body_keyframe_at_k = reef_msgs::DCM_from_Euler321(Eigen::Vector3d(roll_init, pitch_init, 0));
        tf2::fromMsg(pose_msg->pose, xyEst.global_pose_p); //To Affine3d

        xyEst.initialize();//Initialize P,R and beta.
        ROS_WARN_STREAM("Filter initialized with:\t"<<xyEst.xHat);

        reef_msgs::importMatrixFromParamServer(private_nh_, zEst.xHat0, "z_x0");
        reef_msgs::importMatrixFromParamServer(private_nh_, zEst.P0, "z_P0");
        reef_msgs::importMatrixFromParamServer(private_nh_, zEst.P0forFlying, "z_P0_flying");
        reef_msgs::importMatrixFromParamServer(private_nh_, zEst.Q0, "z_Q");
        reef_msgs::importMatrixFromParamServer(private_nh_, zEst.R0, "z_R0");
        reef_msgs::importMatrixFromParamServer(private_nh_, zEst.RforFlying, "z_R_flying");
        reef_msgs::importMatrixFromParamServer(private_nh_, zEst.betaVector, "z_beta");
        zEst.Q = zEst.Q0*(zEst.dt);
        zEst.updateLinearModel();
        zEst.initialize();
        zEst.setTakeoffState(false);



        state_publisher_ = nh_.advertise<reef_msgs::XYZEstimate>("xyz_estimate", 1, true);
        pose_publisher_ =  nh_.advertise<geometry_msgs::PoseStamped>("xyz_pose", 1, true);
        delta_measurement_publisher_ =  nh_.advertise<geometry_msgs::PoseStamped>("delta_measurement", 1, true);

        if (debug_mode_) {
            debug_state_publisher_ = nh_.advertise<reef_msgs::XYZDebugEstimate>("xyz_debug_estimate", 1, true);

            //Print initial matrices in debug mode
            ROS_INFO_STREAM("xy_x0: " << std::endl << xyEst.xHat0);
            ROS_INFO_STREAM("xy_P: " << std::endl << xyEst.P);
            ROS_INFO_STREAM("xy_Q: " << std::endl << xyEst.Q);
            ROS_INFO_STREAM("xy_R: " << std::endl << xyEst.R);
            ROS_INFO_STREAM("xy_beta: " << std::endl << xyEst.betaVector << std::endl);

            ROS_INFO_STREAM("z_x0: " << std::endl << zEst.xHat0);
            ROS_INFO_STREAM("z_P: " << std::endl << zEst.P);
            ROS_INFO_STREAM("z_Q: " << std::endl << zEst.Q);
            ROS_INFO_STREAM("z_R: " << std::endl << zEst.R);
            ROS_INFO_STREAM("z_beta: " << std::endl << zEst.betaVector << std::endl);
        }
        
        is_flying_publisher_ = nh_.advertise<std_msgs::Bool>("is_flying_reef", 1, true);

        //Initialize member variables
        accSampleAverage.x = accSampleAverage.y = accSampleAverage.z = 0;
        takeoffState.data = false;
        zSigma(0) = zSigma(1) = 0;

    }

    XYZEstimator::~XYZEstimator() {}

/** This function is used to calculate the mean accelerometer bias . */
    void XYZEstimator::initializeAcc(geometry_msgs::Vector3 acc)
    {
        //Sum ACC_SAMPLE_SIZE accelerometer readings
        accSampleAverage.x += acc.x;
        accSampleAverage.y += acc.y;
        accSampleAverage.z += acc.z;
        accInitSampleCount++;

        if (accInitSampleCount == ACC_SAMPLE_SIZE) 
        {
            accSampleAverage.x /= ACC_SAMPLE_SIZE;
            accSampleAverage.y /= ACC_SAMPLE_SIZE;
            accSampleAverage.z /= ACC_SAMPLE_SIZE;

            initialAccMagnitude = getVectorMagnitude(accSampleAverage.x, accSampleAverage.y, accSampleAverage.z);
            initialAccMagnitude = 9.81;
            ROS_INFO_STREAM("initial g: " << std::endl << initialAccMagnitude << std::endl);

            accInitialized = true;
        }
    }
/** Sensor update for the IMU. */
    void XYZEstimator::sensorUpdate(sensor_msgs::Imu imu) 
    {
        //Save the stamp. This is very important for good book-keeping.
        xyzState.header.stamp = imu.header.stamp;
        xyzDebugState.header.stamp = imu.header.stamp;
        xyzPose.header.stamp = imu.header.stamp;

        if (isnan(getVectorMagnitude(imu.linear_acceleration.x, imu.linear_acceleration.y,imu.linear_acceleration.z))){
            ROS_ERROR_STREAM("IMU is giving NaNs");
            return;
        }

        //Make sure accelerometer is initialized
        if (!accInitialized) {
            initializeAcc(imu.linear_acceleration);
            last_time_stamp = imu.header.stamp.toSec();
            return;
        }
        // Compute new DT and pass it to estimators.
        xyEst.dt = imu.header.stamp.toSec() - last_time_stamp;
        zEst.dt = imu.header.stamp.toSec() - last_time_stamp;
        dt = imu.header.stamp.toSec() - last_time_stamp;
        last_time_stamp = imu.header.stamp.toSec();

        C_NED_to_body_frame = reef_msgs::quaternion_to_rotation(imu.orientation);
        //TODO incorporate rotation into range measurement.
        //zEst.H(0,0) = (-1/C_NED_to_body_frame(2,2));

        //Transform accel measurements from body to NED frame.
        accelxyz_in_body_frame << imu.linear_acceleration.x, imu.linear_acceleration.y, imu.linear_acceleration.z; //This is a column vector.
        accelxyz_in_NED_frame = C_NED_to_body_frame.transpose() * accelxyz_in_body_frame;

        //Get Gyro measurements from the IMU
        gyroxyz_in_body_frame << imu.angular_velocity.x, imu.angular_velocity.y, imu.angular_velocity.z; //This is a column vector.

        /*The specific force model is the following s = a_measured - bias_accel - noise_accel - gravity
         * The bias that is being estimated is in the inertial frame (NED)
         * */
        zEst.u(0) = accelxyz_in_NED_frame(2) + initialAccMagnitude; //We need to take avg from csv file to get a better g.
        //        zEst.u(0) = accelxyz_in_NED_frame(2) + initialAccMagnitude + zEst.xHat(2); //We need to take avg from csv file to get a better g.

        //Finally propagate.
        xyEst.nonlinearPropagation(C_NED_to_body_frame, initialAccMagnitude, accelxyz_in_body_frame, gyroxyz_in_body_frame, zEst.xHat(2));
        zEst.updateLinearModel();
        zEst.propagate();

        //Update state publisher with states, error, and three sigma bounds

        //Z estimator publisher block------------------------------------------
        if (debug_mode_)
        {
            saveMinusState();
        }

        /*We let the previous block to keep running for approximately 0.1 seconds,
            * if we haven't taken off, we reset the covariance P to its initial value and keep R big enough to
            * prevent our covariance from totally shrinking
           */

        numberOfPropagations++;
        //Reset the z estimator to its landing state every 10 propagations
        if (!takeoffState.data && numberOfPropagations >= 10)
        {
            if (enableXY)
                xyEst.resetLandingState();

            if (enableZ)
                zEst.resetLandingState();

            //Reset number of propagations
            numberOfPropagations = 0;
        }

        if (enableXY && newDeltaMeasurement) {
            if(enable_partial_update) {
                xyEst.partialUpdate();
            }
            else{
                xyEst.update();
            }


            newDeltaMeasurement = false;
        }

        if (enableZ && newSonarMeasurement) {
            //TODO adjust estimator to perform partialUpdate on z as well.
            if(enable_partial_update) {
                zEst.partialUpdate();
            }
            else {
                zEst.update();
            }
            newSonarMeasurement = false;
        }


        integrateGlobalPose();

        publishEstimates();

        if ((xyEst.XYTakeoff && (xyEst.distance > xyEst.dPoseLimit)) || (xyEst.XYTakeoff && (xyEst.xHat(xyEst.YAW) > xyEst.dYawLimit))) {

            //Save states at the time of keyframe
            C_level_keyframe_to_body_keyframe_at_k = xyEst.C_body_level_to_body_frame;
            global_pose.translation() = global_pose.translation() + C3_Yaw.transpose()*current_deltaXY;
            global_yaw = global_yaw + xyEst.xHat(xyEst.YAW);
            //Perform reset
            xyEst.relativeReset();
        }

        checkTakeoffState(accelxyz_in_body_frame.norm());
    }

    //Sonar update
    void XYZEstimator::sensorUpdate(sensor_msgs::Range range_msg) 
    {

        if (range_msg.range <= range_msg.max_range)
        {
            if (chi2Accept(range_msg.range))
            {
                zEst.z(0) = -range_msg.range; // negative size to convert to NED
                newSonarMeasurement = true;
            }
        }

    }

    //Mocap XY Pose update
    void XYZEstimator::deltaPoseUpdate(geometry_msgs::PoseStamped pose_msg){
        // TODO: Humberto please do all the coordinate frame transforms. It needs to go before chi2 rejection

        Eigen::Matrix3d delta_C;
        delta_C = reef_msgs::quaternion_to_rotation(pose_msg.pose.orientation);

        Eigen::Vector3d delta_p_camera_frame;
        delta_p_camera_frame << pose_msg.pose.position.x, pose_msg.pose.position.y, pose_msg.pose.position.z;
        Eigen::Vector3d delta_p_body_level_k_to_body_level_kplus1;
//        ROS_WARN_STREAM("BODY TO CAM P"<<xyEst.p_body_to_camera);
        delta_p_body_level_k_to_body_level_kplus1 << C_level_keyframe_to_body_keyframe_at_k.transpose() *( xyEst.p_body_to_camera + xyEst.C_body_to_camera.transpose() * ( delta_p_camera_frame + delta_C.transpose()*xyEst.C_body_to_camera*(-xyEst.p_body_to_camera)) );

        Eigen::Matrix3d delta_C_level_k_to_level_kplus1;
        delta_C_level_k_to_level_kplus1 = xyEst.C_body_level_to_body_frame.transpose() * xyEst.C_body_to_camera.transpose() * delta_C * xyEst.C_body_to_camera * C_level_keyframe_to_body_keyframe_at_k;
        double delta_yaw;
        reef_msgs::get_yaw(delta_C_level_k_to_level_kplus1, delta_yaw);

        xyEst.z << delta_p_body_level_k_to_body_level_kplus1(0), delta_p_body_level_k_to_body_level_kplus1(1), delta_yaw;
        //This message is just to perform the Chi2 rejection.
        delta_msg.header = pose_msg.header;
        delta_msg.pose.position.x =   xyEst.z(0);
        delta_msg.pose.position.y =   xyEst.z(1);
        delta_msg.pose.position.z =   xyEst.z(2);
        delta_measurement_publisher_.publish(delta_msg);

//TODO: Include the rejection scheme for the delta measurements.
        if (chi2AcceptDeltaPose(delta_msg)){
            newDeltaMeasurement = true;
        }
    }

    bool XYZEstimator::chi2Accept(float range_measurement) 
    {
        //Compute Mahalanobis distance.
        range = Eigen::VectorXd(1);
        range(0) = -range_measurement;
        expected_measurement = Eigen::VectorXd(1);
        expected_measurement = zEst.H * zEst.xHat;

        Eigen::MatrixXd S(1, 1);
        S = zEst.H * zEst.P * zEst.H.transpose() + zEst.R;

        Mahalanobis_D_hat_square(0) = (range - expected_measurement).transpose() * S.inverse() * (range - expected_measurement);
        Mahalanobis_D_hat(0) = sqrt(Mahalanobis_D_hat_square(0));

        //Value for 99% we need 6.63.
        //Value for 95% we 3.84
        if (Mahalanobis_D_hat_square(0) > mahalanobis_distance_sonar)
        {
            ROS_INFO_STREAM("Range measurement rejected");
            return false;
        } 
        else
        {
            return true;
        }
    }

    bool XYZEstimator::chi2AcceptDeltaPose(geometry_msgs::PoseStamped pose){
        //Compute Mahalanobis distance.
        Eigen::Vector3d measurement;
        measurement << pose.pose.position.x, pose.pose.position.y, pose.pose.position.z;

        Eigen::Vector3d expected_vo;
        expected_vo = xyEst.H * xyEst.xHat;
        Eigen::MatrixXd S(3, 3);

        S = xyEst.H * xyEst.P * xyEst.H.transpose() + xyEst.R;
        Mahalanobis_D_hat_square(0) = (measurement - expected_vo).transpose() * S.inverse() * (measurement - expected_vo);
        Mahalanobis_D_hat(0) = sqrt(Mahalanobis_D_hat_square(0));
        //Value for 99% we need 6.63.
        //Value for 95% we 3.84
        if (Mahalanobis_D_hat_square(0) > mahalanobis_distance_mocap_xy_)
        {
            ROS_INFO_STREAM("MOCAP XY measurement rejected");
            return false;
        } 
        else
            return true;
    }

    //Hypothesis rejector - can be used with sonar instead of chi2
    bool XYZEstimator::hypothesisAccept(float range_measurement)
    {
        float tau;
        float sigma1 = 0.03;
        float sigma2 = 0.1;
        float logLR;
        float prior_ratio;
        float meanX = zEst.xHat(0);
        float meanY = zEst.xHat(0) - 3;

        logLR = -pow(-range_measurement - meanX, 2) / pow(sigma1, 2) +
                pow(-range_measurement - meanY, 2) / pow(sigma2, 2);
        tau = 0.1 / 0.9;
        prior_ratio = 2 * log(tau / (sigma2 / sigma1));

        //ROS_INFO_STREAM(logLR);
        //ROS_INFO_STREAM(prior_ratio);

        if (logLR >= prior_ratio) 
        {
            //measurement accepted
            return true;
        } 
        else 
        {
            ROS_INFO_STREAM("Range measurement rejected");
            return false;
        }
    }

    void XYZEstimator::checkTakeoffState(double accMagnitude) 
    {
        sonarTakeoffState = zEst.z(0) <= -0.18;

        //Check variance of accelerometer vector magnitude
        if ((!accTakeoffState) || (!sonarTakeoffState && accTakeoffState))
        {
            if (numAccSamples < ACC_SAMPLE_SIZE) 
            {
                accSamples[numAccSamples++] = accMagnitude;
            } 
            else 
            {
                //shift sample in and compute the sample mean simultaneously
                accMean = 0;
                for (int i = 0; i < ACC_SAMPLE_SIZE - 1; i++) 
                {
                    accSamples[i] = accSamples[i + 1];
                    accMean += accSamples[i];
                }
                accSamples[ACC_SAMPLE_SIZE - 1] = accMagnitude;
                accMean += accMagnitude;
                accMean /= (double) ACC_SAMPLE_SIZE;

                //Finally, compute the sample variance
                accVariance = 0;
                for (int i = 0; i < ACC_SAMPLE_SIZE; i++) 
                {
                    double distFromMean = accSamples[i] - accMean;
                    accVariance += (distFromMean * distFromMean);
                }
                accVariance /= (double) (ACC_SAMPLE_SIZE - 1);
                accTakeoffState = accVariance >= ACC_TAKEOFF_VARIANCE;
            }
        } 
        else if (numAccSamples > 0) 
        {
            numAccSamples = 0;
        }

        //Publish changes in takeoff state
        if (accTakeoffState && sonarTakeoffState && !takeoffState.data) 
        {
            ROS_INFO("Takeoff!");
            zEst.setTakeoffState(true);
            xyEst.XYTakeoff = true;
            takeoffState.data = true;
            is_flying_publisher_.publish(takeoffState);
        } 
        else if (!accTakeoffState && !sonarTakeoffState && takeoffState.data) 
        {
            ROS_INFO("Landing!");
            zEst.setTakeoffState(false);
            takeoffState.data = false;
            is_flying_publisher_.publish(takeoffState);
        }
    }

    void XYZEstimator::saveMinusState()
    {
        // TODO: Humberto/Grant: Publish the integrated position as a geometry_msg/PoseStamped
        // TODO: Humberto/Grant: Publish the xyz_estimate with position and velocity here
        //XY estimator publisher block------------------------------------------

        xyzDebugState.xy_minus.x_dot = xyEst.xHat(0);
        xyzDebugState.xy_minus.y_dot = xyEst.xHat(1);
        xyzDebugState.xy_minus.pitch_bias = xyEst.xHat(xyEst.BP);
        xyzDebugState.xy_minus.roll_bias = xyEst.xHat(xyEst.BR);
        xyzDebugState.xy_minus.xa_bias = xyEst.xHat(4);
        xyzDebugState.xy_minus.ya_bias = xyEst.xHat(5);
        xyzDebugState.xy_minus.delta_x = xyEst.xHat(xyEst.PX);
        xyzDebugState.xy_minus.delta_y = xyEst.xHat(xyEst.PY);
        xyzDebugState.xy_minus.delta_yaw = xyEst.xHat(xyEst.YAW);

        xySigma = Eigen::MatrixXd(12, 1);
        xySigma(0) = 3 * sqrt(xyEst.P(0, 0));
        xySigma(1) = 3 * sqrt(xyEst.P(1, 1));
        xySigma(2) = 3 * sqrt(xyEst.P(2, 2));
        xySigma(3) = 3 * sqrt(xyEst.P(3, 3));
        xySigma(4) = 3 * sqrt(xyEst.P(4, 4));
        xySigma(5) = 3 * sqrt(xyEst.P(5, 5));
        xySigma(6) = 3 * sqrt(xyEst.P(6, 6));
        xySigma(7) = 3 * sqrt(xyEst.P(7, 7));
        xySigma(8) = 3 * sqrt(xyEst.P(8, 8));
        xySigma(9) = 3 * sqrt(xyEst.P(9, 9));
        xySigma(10) = 3 * sqrt(xyEst.P(10, 10));
        xySigma(11) = 3 * sqrt(xyEst.P(11, 11));
        xyzDebugState.xy_minus.sigma_plus[0] = xyzDebugState.xy_minus.x_dot + xySigma(0);
        xyzDebugState.xy_minus.sigma_minus[0] = xyzDebugState.xy_minus.x_dot - xySigma(0);
        xyzDebugState.xy_minus.sigma_plus[1] = xyzDebugState.xy_minus.y_dot + xySigma(1);
        xyzDebugState.xy_minus.sigma_minus[1] = xyzDebugState.xy_minus.y_dot - xySigma(1);
        xyzDebugState.xy_minus.sigma_plus[2] = xyzDebugState.xy_minus.pitch_bias + xySigma(2);
        xyzDebugState.xy_minus.sigma_minus[2] = xyzDebugState.xy_minus.pitch_bias - xySigma(2);
        xyzDebugState.xy_minus.sigma_plus[3] = xyzDebugState.xy_minus.roll_bias + xySigma(3);
        xyzDebugState.xy_minus.sigma_minus[3] = xyzDebugState.xy_minus.roll_bias - xySigma(3);
        xyzDebugState.xy_minus.sigma_plus[4] = xyzDebugState.xy_minus.xa_bias + xySigma(4);
        xyzDebugState.xy_minus.sigma_minus[4] = xyzDebugState.xy_minus.xa_bias - xySigma(4);
        xyzDebugState.xy_minus.sigma_plus[5] = xyzDebugState.xy_minus.ya_bias + xySigma(5);
        xyzDebugState.xy_minus.sigma_minus[5] = xyzDebugState.xy_minus.ya_bias - xySigma(5);

        xyzDebugState.xy_minus.sigma_plus[6] = xyzDebugState.xy_minus.ya_bias + xySigma(6);
        xyzDebugState.xy_minus.sigma_minus[6] = xyzDebugState.xy_minus.ya_bias - xySigma(6);
        xyzDebugState.xy_minus.sigma_plus[7] = xyzDebugState.xy_minus.ya_bias + xySigma(7);
        xyzDebugState.xy_minus.sigma_minus[7] = xyzDebugState.xy_minus.ya_bias - xySigma(7);
        xyzDebugState.xy_minus.sigma_plus[8] = xyzDebugState.xy_minus.ya_bias + xySigma(8);
        xyzDebugState.xy_minus.sigma_minus[8] = xyzDebugState.xy_minus.ya_bias - xySigma(8);
        xyzDebugState.xy_minus.sigma_plus[9] = xyzDebugState.xy_minus.ya_bias + xySigma(9);
        xyzDebugState.xy_minus.sigma_minus[9] = xyzDebugState.xy_minus.ya_bias - xySigma(9);
        xyzDebugState.xy_minus.sigma_plus[10] = xyzDebugState.xy_minus.ya_bias + xySigma(10);
        xyzDebugState.xy_minus.sigma_minus[10] = xyzDebugState.xy_minus.ya_bias - xySigma(10);
        xyzDebugState.xy_minus.sigma_plus[11] = xyzDebugState.xy_minus.ya_bias + xySigma(11);
        xyzDebugState.xy_minus.sigma_minus[11] = xyzDebugState.xy_minus.ya_bias - xySigma(11);

        //Z estimator publisher block------------------------------------------
        xyzDebugState.z_minus.z = zEst.xHat(0, 0);
        xyzDebugState.z_minus.z_dot = zEst.xHat(1, 0);
        xyzDebugState.z_minus.bias = zEst.xHat(2, 0);
        xyzDebugState.z_minus.u = zEst.u(0);
        reef_msgs::matrixToArray(zEst.P, xyzDebugState.z_minus.P);
        zSigma(0) = 3 * sqrt(zEst.P(0, 0));
        zSigma(1) = 3 * sqrt(zEst.P(1, 1));
        zSigma(2) = 3 * sqrt(zEst.P(2, 2));
        xyzDebugState.z_minus.sigma_plus[0] = xyzDebugState.z_minus.z + zSigma(0);
        xyzDebugState.z_minus.sigma_minus[0] = xyzDebugState.z_minus.z - zSigma(0);
        xyzDebugState.z_minus.sigma_plus[1] = xyzDebugState.z_minus.z_dot + zSigma(1);
        xyzDebugState.z_minus.sigma_minus[1] = xyzDebugState.z_minus.z_dot - zSigma(1);
        xyzDebugState.z_minus.sigma_plus[2] = xyzDebugState.z_minus.bias + zSigma(2);
        xyzDebugState.z_minus.sigma_minus[2] = xyzDebugState.z_minus.bias - zSigma(2);
    }

    void XYZEstimator::publishEstimates()
    {
        // TODO: Humberto/Grant: Publish the integrated position as a geometry_msg/PoseStamped
        // TODO: Humberto/Grant: Publish the xyz_estimate with position and velocity here
        //Publish XY and Z states

        xyzState.xy_plus.x_dot = xyEst.xHat(0);
        xyzState.xy_plus.y_dot = xyEst.xHat(1);
        xyzState.xy_plus.delta_x = xyEst.xHat(xyEst.PX);
        xyzState.xy_plus.delta_y = xyEst.xHat(xyEst.PY);
        xyzState.xy_plus.delta_yaw = xyEst.xHat(xyEst.YAW);

        xyzState.z_plus.z = zEst.xHat(0);
        xyzState.z_plus.z_dot = zEst.xHat(1);

        state_publisher_.publish(xyzState);

        xyEst.Delta.header = xyzState.header;

        if (debug_mode_)
        {

            //XY filter debug variables
            xyzDebugState.xy_plus.x_dot = xyEst.xHat(0);
            xyzDebugState.xy_plus.y_dot = xyEst.xHat(1);
            xyzDebugState.xy_plus.pitch_bias = xyEst.xHat(xyEst.BP);
            xyzDebugState.xy_plus.roll_bias = xyEst.xHat(xyEst.BR);
            xyzDebugState.xy_plus.xa_bias = xyEst.xHat(4);
            xyzDebugState.xy_plus.ya_bias = xyEst.xHat(5);
            xyzDebugState.xy_plus.delta_x = xyEst.xHat(xyEst.PX);
            xyzDebugState.xy_plus.delta_y = xyEst.xHat(xyEst.PY);
            xyzDebugState.xy_plus.delta_yaw = xyEst.xHat(xyEst.YAW);
            xyzDebugState.xy_plus.gx_bias = xyEst.xHat(xyEst.BGX);
            xyzDebugState.xy_plus.gy_bias = xyEst.xHat(xyEst.BGY);
            xyzDebugState.xy_plus.gz_bias = xyEst.xHat(xyEst.BGZ);

            xySigma = Eigen::MatrixXd(12, 1);
            xySigma(0) = 3 * sqrt(xyEst.P(0, 0));
            xySigma(1) = 3 * sqrt(xyEst.P(1, 1));
            xySigma(2) = 3 * sqrt(xyEst.P(2, 2));
            xySigma(3) = 3 * sqrt(xyEst.P(3, 3));
            xySigma(4) = 3 * sqrt(xyEst.P(4, 4));
            xySigma(5) = 3 * sqrt(xyEst.P(5, 5));
            xySigma(6) = 3 * sqrt(xyEst.P(6, 6));
            xySigma(7) = 3 * sqrt(xyEst.P(7, 7));
            xySigma(8) = 3 * sqrt(xyEst.P(8, 8));
            xySigma(9) = 3 * sqrt(xyEst.P(9, 9));
            xySigma(10) = 3 * sqrt(xyEst.P(10, 10));
            xySigma(11) = 3 * sqrt(xyEst.P(11, 11));
            xyzDebugState.xy_plus.sigma_plus[0] = xyzDebugState.xy_plus.x_dot + xySigma(0);
            xyzDebugState.xy_plus.sigma_minus[0] = xyzDebugState.xy_plus.x_dot - xySigma(0);
            xyzDebugState.xy_plus.sigma_plus[1] = xyzDebugState.xy_plus.y_dot + xySigma(1);
            xyzDebugState.xy_plus.sigma_minus[1] = xyzDebugState.xy_plus.y_dot - xySigma(1);
            xyzDebugState.xy_plus.sigma_plus[2] = xyzDebugState.xy_plus.roll_bias + xySigma(2);
            xyzDebugState.xy_plus.sigma_minus[2] = xyzDebugState.xy_plus.roll_bias - xySigma(2);
            xyzDebugState.xy_plus.sigma_plus[3] = xyzDebugState.xy_plus.pitch_bias + xySigma(3);
            xyzDebugState.xy_plus.sigma_minus[3] = xyzDebugState.xy_plus.pitch_bias - xySigma(3);
            xyzDebugState.xy_plus.sigma_plus[4] = xyzDebugState.xy_plus.xa_bias + xySigma(4);
            xyzDebugState.xy_plus.sigma_minus[4] = xyzDebugState.xy_plus.xa_bias - xySigma(4);
            xyzDebugState.xy_plus.sigma_plus[5] = xyzDebugState.xy_plus.ya_bias + xySigma(5);
            xyzDebugState.xy_plus.sigma_minus[5] = xyzDebugState.xy_plus.ya_bias - xySigma(5);

            xyzDebugState.xy_plus.sigma_plus[6] = xyzDebugState.xy_plus.delta_x + xySigma(6);
            xyzDebugState.xy_plus.sigma_minus[6] = xyzDebugState.xy_plus.delta_x - xySigma(6);
            xyzDebugState.xy_plus.sigma_plus[7] = xyzDebugState.xy_plus.delta_y + xySigma(7);
            xyzDebugState.xy_plus.sigma_minus[7] = xyzDebugState.xy_plus.delta_y - xySigma(7);
            xyzDebugState.xy_plus.sigma_plus[8] = xyzDebugState.xy_plus.delta_yaw + xySigma(8);
            xyzDebugState.xy_plus.sigma_minus[8] = xyzDebugState.xy_plus.delta_yaw - xySigma(8);
            xyzDebugState.xy_plus.sigma_plus[9] = xyzDebugState.xy_plus.gx_bias + xySigma(9);
            xyzDebugState.xy_plus.sigma_minus[9] = xyzDebugState.xy_plus.gx_bias - xySigma(9);
            xyzDebugState.xy_plus.sigma_plus[10] = xyzDebugState.xy_plus.gy_bias + xySigma(10);
            xyzDebugState.xy_plus.sigma_minus[10] = xyzDebugState.xy_plus.gy_bias - xySigma(10);
            xyzDebugState.xy_plus.sigma_plus[11] = xyzDebugState.xy_plus.gz_bias + xySigma(11);
            xyzDebugState.xy_plus.sigma_minus[11] = xyzDebugState.xy_plus.gz_bias - xySigma(11);


            //Z filter debug variables
            xyzDebugState.z_plus.z = zEst.xHat(0);
            xyzDebugState.z_plus.z_dot = zEst.xHat(1);
            xyzDebugState.z_plus.bias = zEst.xHat(2);
            xyzDebugState.z_plus.u = zEst.u(0);
            reef_msgs::matrixToArray(zEst.P, xyzDebugState.z_plus.P);
            zSigma(0) = 3 * sqrt(zEst.P(0, 0));
            zSigma(1) = 3 * sqrt(zEst.P(1, 1));
            zSigma(2) = 3 * sqrt(zEst.P(2, 2));
            xyzDebugState.z_plus.sigma_plus[0] = xyzDebugState.z_plus.z + zSigma(0);
            xyzDebugState.z_plus.sigma_minus[0] = xyzDebugState.z_plus.z - zSigma(0);
            xyzDebugState.z_plus.sigma_plus[1] = xyzDebugState.z_plus.z_dot + zSigma(1);
            xyzDebugState.z_plus.sigma_minus[1] = xyzDebugState.z_plus.z_dot - zSigma(1);
            xyzDebugState.z_plus.sigma_plus[2] = xyzDebugState.z_plus.bias + zSigma(2);
            xyzDebugState.z_plus.sigma_minus[2] = xyzDebugState.z_plus.bias - zSigma(2);

            debug_state_publisher_.publish(xyzDebugState);

        }
    }

    //Utility function
    double getVectorMagnitude(double x, double y, double z) 
    {
        return sqrt(x * x + y * y + z * z);
    }

    void XYZEstimator::integrateGlobalPose(){
        current_deltaXY << xyEst.xHat(xyEst.PX), xyEst.xHat(xyEst.PY), 0;
//        current_estimate.translation() = Eigen::Vector3d(xyEst.xHat(xyEst.PX), xyEst.xHat(xyEst.PY),0);
        C3_Yaw<< cos(global_yaw), sin(global_yaw),0,
                -sin(global_yaw), cos(global_yaw),0,
                0,                                    0,          1.0;

        //Compound current delta pose with global pose
        Eigen::Affine3d displayed_pose;
        Eigen::Vector3d displayed_global_position;
        double displayed_global_yaw;
        displayed_global_position << global_pose.translation() + C3_Yaw.transpose()*current_deltaXY;
        displayed_global_yaw = global_yaw + xyEst.xHat(xyEst.YAW);

        //Publish the accumulated pose
        xyzPose.pose.position.x = displayed_global_position.x();
        xyzPose.pose.position.y = displayed_global_position.y();
        xyzPose.pose.position.z = (displayed_global_yaw)*57.3;
        xyzPose.header = xyzState.header;
        pose_publisher_.publish(xyzPose);


    }

}

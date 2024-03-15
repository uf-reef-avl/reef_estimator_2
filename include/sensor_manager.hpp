/**
 * \class SensorManager
 *
 * \ingroup reef_estimator
 *
 * This class is meant to have all the ROS callbacks and recieve all the sensor readings
 * This class also publishes the sonar reading in the NED frame to verify it's measurements
 *
 * \author $Author: bv Humberto Ramos, William Warke, Prashant Ganesh
 *
 * \version $Revision: 1.0 $
 *
 * \date $Date: 2019/03/05 $
 *
 * Contact: prashant.ganesh@ufl.edu
 *
 */


#ifndef SENSOR_MANAGER_H
#define SENSOR_MANAGER_H

#include "rclcpp/rclcpp.hpp"
#include <std_msgs/msg/bool.hpp>
#include <sensor_msgs/msg/imu.hpp>
#include <sensor_msgs/msg/range.hpp>
#include <geometry_msgs/msg/poseWithCovarianceStamped.hpp>
#include <rosflight_msgs/msg/RCRaw.hpp>
#include <reef_msgs/msg/DeltaToVel.hpp>
#include "xyz_estimator.hpp"
#include "z_estimator.hpp"

namespace reef_estimator {
    class SensorManager {
    private:

        ros::NodeHandle private_nh_;
        ros::NodeHandle nh_;

        ros::Subscriber imu_subscriber_;
        ros::Subscriber altimeter_subscriber_;
        ros::Subscriber rc_subscriber_;
        ros::Subscriber delta_pose_subscriber_;

        ros::Publisher range_ned_publisher_;

        XYZEstimator xyzEst;
        sensor_msgs::Imu imu_msg;
        sensor_msgs::Range range_msg;
        sensor_msgs::Range range_msg_ned; //used in debug mode

        void isFlyingCallback(const std_msgs::BoolConstPtr &msg);
        void imuCallback(const sensor_msgs::ImuConstPtr &msg);
        void altimeterCallback(const sensor_msgs::RangeConstPtr &msg);
        void rcRawCallback(const rosflight_msgs::RCRawConstPtr &msg);
        void deltaPoseCallback(const geometry_msgs::PoseWithCovarianceStampedConstPtr &msg);

        bool initialized_;
        bool imuCalibrated;
        bool get_measurements;

        //Mocap override variables
        bool mocapSwitchOn;
        int mocapOverrideChannel;

    public:
        SensorManager();
        ~SensorManager();
    };
}
#endif

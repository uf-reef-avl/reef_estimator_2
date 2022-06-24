#include "sensor_manager.h"

namespace reef_estimator
{
    SensorManager::SensorManager() : private_nh_("~"), nh_(""), initialized_(false), mocapSwitchOn(false)
    {
        //Mocap override RC channel parameter
        bool enableMocapSwitch = false;
        private_nh_.param<bool>("enable_mocap_switch", enableMocapSwitch, false); ///< Enable switching between sensor readings midflight
        private_nh_.param<int>("mocap_override_channel", mocapOverrideChannel, 4); ///< RC Channel to enable switching

        if (enableMocapSwitch) {
            ROS_WARN_STREAM("Mocap override RC switch enabled on channel " << mocapOverrideChannel);
            rc_subscriber_ = nh_.subscribe("rc_raw", 1, &SensorManager::rcRawCallback, this);
        }

        std::string mocapPoseTopic;
        private_nh_.param<std::string>("mocap_pose_topic", mocapPoseTopic, "true_odom");
        delta_pose_subscriber_ = nh_.subscribe<geometry_msgs::PoseWithCovarianceStamped>("true_odom", 1, &SensorManager::deltaPoseCallback, this);

        altimeter_subscriber_ = nh_.subscribe("altimeter", 1, &SensorManager::altimeterCallback, this);
        if(xyzEst.debug_mode_)
            range_ned_publisher_ = nh_.advertise<sensor_msgs::Range>("sonar_ned", 1);

        imu_subscriber_ = nh_.subscribe("imu/data", 10, &SensorManager::imuCallback, this);

    }

    SensorManager::~SensorManager(){}

    void SensorManager::imuCallback(const sensor_msgs::ImuConstPtr &msg)
    {
        //Pass the imu message to estimator.
        xyzEst.sensorUpdate(*msg);
    }

    void SensorManager::rcRawCallback(const rosflight_msgs::RCRawConstPtr &msg) {
        //Check for toggled mocap RC switch
        if (!mocapSwitchOn && msg->values[mocapOverrideChannel] > 1500)
        {
            if (xyzEst.enableMocapXY)
            {
                xyzEst.useMocapXY = true;
                ROS_WARN("Mocap XY feedback enabled");
            }

            if (xyzEst.enableMocapZ)
            {
                xyzEst.useMocapZ = true;
                ROS_WARN("Mocap Z feedback enabled");
            }

            mocapSwitchOn = true;
        }
        else if (mocapSwitchOn && msg->values[mocapOverrideChannel] <= 1500)
        {
            if (xyzEst.enableRGBD)
            {
                xyzEst.useMocapXY = false;
                ROS_WARN("Mocap XY feedback disabled");
            }

            if (xyzEst.enableSonar)
            {
                xyzEst.useMocapZ =  false;
                ROS_WARN("Mocap Z feedback disabled");
            }

            mocapSwitchOn = false;
        }
    }

    void SensorManager::deltaPoseCallback(const geometry_msgs::PoseWithCovarianceStampedConstPtr &msg) {

        geometry_msgs::PoseStamped delta_pose;
        delta_pose.pose = msg->pose.pose;
        delta_pose.header = msg->header;
        xyzEst.deltaPoseUpdate(delta_pose);
    }

    void SensorManager::altimeterCallback(const sensor_msgs::RangeConstPtr &msg){
        xyzEst.sensorUpdate(*msg);

        if (xyzEst.debug_mode_)
        {
            //Publish the negative range measurement
            range_msg_ned = *msg;
            range_msg_ned.header.stamp = range_msg_ned.header.stamp;
            range_msg_ned.range = -range_msg_ned.range;
            range_ned_publisher_.publish(range_msg_ned);
        }
    }



}
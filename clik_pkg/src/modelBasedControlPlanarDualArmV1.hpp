#include "ros/ros.h"
#include "rosgraph_msgs/Clock.h"
#include "gazebo_msgs/SetModelConfiguration.h"
#include "sensor_msgs/JointState.h"
#include "std_msgs/Float64.h"
#include "std_srvs/Empty.h"
#include <thread>
#include <iostream>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <fstream>

class MODEL_BASED_CONTROLLER {
    private:
        ros::NodeHandle _nh;
        ros::Subscriber _gazeboClockSub, _jointStateSub;
        ros::Publisher _jLeftArmEffPub, _jRightArmEffPub;
        ros::ServiceClient _setModelConfigurationClient;
        gazebo_msgs::SetModelConfiguration _configSrv;
        ros::ServiceClient _pauseGazebo;
        double _gazeboTime;
        double _jCablesPos, _jLeftArmPos, _jRightArmPos;
        double _jCablesVel, _jLeftArmVel, _jRightArmVel;
    public:
        MODEL_BASED_CONTROLLER();
        void run();
        void setInitialState();                                 // Set system initial state
        void gazeboClockCB(rosgraph_msgs::Clock);               // Gazebo clock callback
        void jointStateCB(sensor_msgs::JointState);             // Link state callback
        void noncollocatedPFBLControl();                                     // Control Loop
        void Bfun(Eigen::MatrixXd&, Eigen::VectorXd);
        void nfun(Eigen::VectorXd&, Eigen::VectorXd, Eigen::VectorXd);
};
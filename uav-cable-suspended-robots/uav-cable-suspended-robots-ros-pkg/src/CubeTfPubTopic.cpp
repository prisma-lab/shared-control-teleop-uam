#include "ros/ros.h"
#include "boost/thread.hpp"
#include "tf/transform_broadcaster.h"
#include <tf2/LinearMath/Quaternion.h>
#include "std_msgs/Float64MultiArray.h"
#include "tf/transform_listener.h"

#define RATE 100

using namespace std;

class CUBETOPIC {
	public:
		CUBETOPIC();
		void publish_topic();
	private:
		bool _active;
		double x;
		double y;
		double z;
		double o_w;
		double o_x;
		double o_y;
		double o_z;

		ros::NodeHandle _nh;
		ros::Publisher _topic_pub;
		ros::Rate _rate;
};

CUBETOPIC::CUBETOPIC(): _rate(RATE) {
	_topic_pub = _nh.advertise< std_msgs::Float64MultiArray > ("/licasa1/cube_pos", 1);
	x = y = z = 0;
	o_w = 1;
	o_x = o_y = o_z = 0;

	boost::thread(&CUBETOPIC::publish_topic, this);
}

void CUBETOPIC::publish_topic() {
	tf::TransformListener listener;
	tf::StampedTransform cube_tf;


	while (!listener.waitForTransform("cube_goal", "cube_frame", ros::Time(0), ros::Duration(1))) {
		sleep(1);
		cout<<"WAITING FOR TRANSFORM\n";
	}

	std_msgs::Float64MultiArray cube_goal_msg;
	double roll, pitch, yaw;


	while (ros::ok()) {

		cube_goal_msg.data.clear();
		listener.lookupTransform("cube_goal","cube_frame",ros::Time(0), cube_tf);
		cube_goal_msg.data.push_back(cube_tf.getOrigin().x());
		cube_goal_msg.data.push_back(cube_tf.getOrigin().y());
		cube_goal_msg.data.push_back(cube_tf.getOrigin().z());
		tf::Matrix3x3 m( cube_tf.getRotation() ); // quaternion to RPY
		m.getRPY(roll, pitch, yaw);
		cube_goal_msg.data.push_back(roll);
		cube_goal_msg.data.push_back(pitch);
		cube_goal_msg.data.push_back(yaw);

		_topic_pub.publish(cube_goal_msg);
		
		_rate.sleep();
	}
}

int main( int argc, char** argv ) {

	//Init the ros node with ros_subscriber name
	ros::init(argc, argv, "CubeGoalTfMSg");
	ROS_INFO("CubeGoalTfMSg STARTED\n");

	//Create the sensor filter class object
	CUBETOPIC sf;
	ROS_INFO("CONSTRUCTOR WORKED\n");

	ros::spin();

	return 0;
}

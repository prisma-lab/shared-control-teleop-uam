#include "ros/ros.h"
#include "boost/thread.hpp"
#include "std_msgs/Float64MultiArray.h"
#include "std_msgs/Float64.h"
#include <ros/package.h>
#include <fstream>
#include <string>
#include <iostream>

#define RATE 100
#define Kp 1e-2
#define Ki 1e-2

using namespace std;


class PUB_TRA {
	public:
		PUB_TRA(int mode);
		void cb(std_msgs::Float64MultiArray::ConstPtr msg);
		void control_cb(std_msgs::Float64MultiArray::ConstPtr msg);
		void shared_control_cb(std_msgs::Float64MultiArray::ConstPtr msg);
		void psi_cb(std_msgs::Float64::ConstPtr msg);
		void pub_trajectory();
	private:

		int N;
		
		double* P_x;
		double* P_y;
		double* P_z;
		double* P_psi;
		double* P_d_x;
		double* P_d_y;
		double* P_d_z;
		double* P_d_psi;
		double* P_dd_x;
		double* P_dd_y;
		double* P_dd_z;
		double* P_dd_psi;

		double est[6];
		double psi;

		bool activate_control;
		double F_des; // absolute value of the desired force
		double angle; // angle of the force that the robot feels from the environment
		bool forward;

		double delta_x_send_shared;
		double delta_y_send_shared;
		double delta_z_send_shared;

		ros::NodeHandle _nh;
		ros::Publisher _topic_pub;
		ros::Publisher _delta_pub;
		ros::Subscriber _topic_sub;
		ros::Subscriber _psi_sub;
		ros::Subscriber _shared_control_sub;
		ros::Subscriber _activate_control_sub;
		ros::Rate _rate;
};

PUB_TRA::PUB_TRA(int mode): _rate(RATE) {
	_topic_sub = _nh.subscribe("/licasa1/estimate", 1, &PUB_TRA::cb, this);
	_psi_sub = _nh.subscribe("/licasa1/psi", 1, &PUB_TRA::psi_cb, this);
	_topic_pub = _nh.advertise< std_msgs::Float64MultiArray > ("/licasa1/ref_topic", 1);
	_delta_pub = _nh.advertise< std_msgs::Float64MultiArray > ("/licasa1/delta", 1);
	_activate_control_sub = _nh.subscribe("/licasa1/activate_control", 1, &PUB_TRA::control_cb, this);
	_shared_control_sub = _nh.subscribe("/licasa1/shared_control", 1, &PUB_TRA::shared_control_cb, this);

	ifstream ref;

	string pkg_loc = ros::package::getPath("hierarchical_ctrl_pkg");
	if (mode == 1) {
		ref.open(pkg_loc + "/UAV_grasp.txt", ios::in); 
	} else {
		ref.open(pkg_loc + "/UAV_push.txt", ios::in); 
	}

	if (!ref.is_open()) {
		cout<<"Error opening the file!";
		exit(1);
	}

	// trajectory length
	ref>>N;

	P_x = new double[N];
	P_y = new double[N];
	P_z = new double[N];
    P_psi = new double[N];
	P_d_x = new double[N];
	P_d_y = new double[N];
	P_d_z = new double[N];
    P_d_psi = new double[N];
	P_dd_x = new double[N];
	P_dd_y = new double[N];
	P_dd_z = new double[N];
    P_dd_psi = new double[N];

	for (int i = 0; i<N; i++) {
		ref>>P_x[i];
		ref>>P_y[i];
		ref>>P_z[i];
        ref>>P_psi[i];
		ref>>P_d_x[i];
		ref>>P_d_y[i];
		ref>>P_d_z[i];
        ref>>P_d_psi[i];
		ref>>P_dd_x[i];
		ref>>P_dd_y[i];
		ref>>P_dd_z[i];
        ref>>P_dd_psi[i];
	}

	for (int i = 0; i < 6; i++) {
		est[i] = 0;
	}

	psi = 0;

	F_des = 0;
	angle = 0;

	activate_control = false;
	forward = true;

	delta_x_send_shared = 0;
	delta_y_send_shared = 0;
	delta_z_send_shared = 0;

	ref.close();
	
	boost::thread(&PUB_TRA::pub_trajectory, this);
}

//Callback function: the input of the function is the data to read
void PUB_TRA::cb(std_msgs::Float64MultiArray::ConstPtr msg) {
	for (int i = 0; i < 6; i++) {
		est[i] = msg->data[i];
	}
}

void PUB_TRA::psi_cb(std_msgs::Float64::ConstPtr msg) {
	psi = msg->data;
}

//Callback function: the input of the function is the data to read
void PUB_TRA::control_cb(std_msgs::Float64MultiArray::ConstPtr msg) {
	if (msg->data[0] > 0) {
		activate_control = true;
	} else {
		activate_control = false;
	}

	F_des = msg->data[1];
	angle = msg->data[2];

	if (msg->data[3] > 0) {
		forward = true;
	} else {
		forward = false;
	}
}

//Callback function: the input of the function is the data to read
void PUB_TRA::shared_control_cb(std_msgs::Float64MultiArray::ConstPtr msg) {
	delta_x_send_shared = delta_x_send_shared*0.95 + 0.05*msg->data[0];
	delta_y_send_shared = delta_y_send_shared*0.95 + 0.05*msg->data[1];
	delta_z_send_shared = delta_z_send_shared*0.95 + 0.05*msg->data[2];
}

void PUB_TRA::pub_trajectory() {

	int ii = 0;
	std_msgs::Float64MultiArray ref_msg;

	double err_x = 0;
	double err_x_int = 0;
	double delta_x = 0;
	double delta_x_send = 0;
	double err_y = 0;
	double err_y_int = 0;
	double delta_y = 0;
	double delta_y_send = 0;

	while (ros::ok()) {

		if (activate_control) {
			err_x = F_des*cos(angle) - est[0];
			err_x_int = err_x_int + err_x*1/RATE;;
			delta_x = -Kp*err_x - Ki*err_x_int;
			delta_x_send = 0.9*delta_x_send + 0.1*delta_x;
			err_y = F_des*sin(angle) - est[1];
			err_y_int = err_y_int + err_y*1/RATE;;
			delta_y = -Kp*err_y - Ki*err_y_int;
			delta_y_send = 0.9*delta_y_send + 0.1*delta_y;
			cout<<"delta_x: "<<delta_x<<endl;
			cout<<"delta_x_send: "<<delta_x_send<<endl;
			cout<<"delta_y: "<<delta_y<<endl;
			cout<<"delta_y_send: "<<delta_y_send<<endl;
		} else {
			err_x = 0;
			err_x_int = 0;
			delta_x = 0;
			if (abs(delta_x_send) < 1e-4) {
				delta_x_send = 0;
			} else {
				delta_x_send = 0.995*delta_x_send;
				cout<<"delta_x_send: "<<delta_x_send<<endl;
			}
			err_y = 0;
			err_y_int = 0;
			delta_y = 0;
			if (abs(delta_y_send) < 1e-4) {
				delta_y_send = 0;
			} else {
				delta_y_send = 0.995*delta_y_send;
				cout<<"delta_y_send: "<<delta_y_send<<endl;
			}
		}

		ref_msg.data.clear();
		// ref_msg.data.push_back(P_x[ii] + delta_x_send*cos(P_psi[ii]) - delta_y_send*sin(P_psi[ii]) + delta_x_send_shared*cos(P_psi[ii]) - delta_y_send_shared*sin(P_psi[ii]));
		ref_msg.data.push_back(P_x[ii] + delta_x_send*cos(psi) - delta_y_send*sin(psi) + delta_x_send_shared);

		// cout<<"delta_x_send_shared: "<<delta_x_send_shared<<endl;
		// cout<<"delta_x_send_shared*cos(psi) - delta_y_send_shared*sin(psi): "<<delta_x_send_shared*cos(psi) - delta_y_send_shared*sin(psi)<<endl;
		// cout<<"psi: "<<psi<<endl;

		// ref_msg.data.push_back(P_y[ii] + delta_x_send*sin(P_psi[ii]) + delta_y_send*cos(P_psi[ii]) + delta_x_send_shared*sin(P_psi[ii]) + delta_y_send_shared*cos(P_psi[ii]));
		ref_msg.data.push_back(P_y[ii] + delta_x_send*sin(psi) + delta_y_send*cos(psi) + delta_y_send_shared);
		ref_msg.data.push_back(P_z[ii] + delta_z_send_shared);
		ref_msg.data.push_back(P_psi[ii]);
		ref_msg.data.push_back(P_d_x[ii]);
		ref_msg.data.push_back(P_d_y[ii]);
		ref_msg.data.push_back(P_d_z[ii]);
		ref_msg.data.push_back(P_d_psi[ii]);
		ref_msg.data.push_back(P_dd_x[ii]);
		ref_msg.data.push_back(P_dd_y[ii]);
		ref_msg.data.push_back(P_dd_z[ii]);
		ref_msg.data.push_back(P_dd_psi[ii]);

		_topic_pub.publish(ref_msg);

		std_msgs::Float64MultiArray delta_msg;
		delta_msg.data.clear();
		delta_msg.data.push_back(delta_x_send);
		delta_msg.data.push_back(delta_y_send);
		delta_msg.data.push_back(delta_x_send_shared);
		delta_msg.data.push_back(delta_y_send_shared);

		_delta_pub.publish(delta_msg);

		// cout<<"REF PUBLISHED\n";

		if (ii < N-1 && forward == true) {
			ii++;
		}

		if (ii > 0 && forward == false) {
			ii--;
		}

		_rate.sleep();
		
	}
}

int main( int argc, char** argv ) {

	int mode = 0; // 0 = push, 1 = grasp
	if (argc==2) {
        mode = atof(argv[1]);
		if (mode != 0 && mode != 1) {
			mode = 0;
			cout<<"Error in the arguments, using default: \n";
		}
    } else {
		cout<<"Error in the arguments, using default: \n";
	}
	cout<<"Selected ";
	if (mode == 1) {
		cout<<"Grasp\n";
	} else {
		cout<<"Push\n";
	}

	//Init the ros node with name
	ros::init(argc, argv, "PublishTrajectory");
	ROS_INFO("TRAJECTORY PUBLISHER STARTED\n");

	PUB_TRA pb(mode);
	
	ros::spin();

	return 0;
}

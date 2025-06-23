#include "example/modelBasedControlPlanarDualArmV1.hpp"
#include "example/utils_2.hpp"
#include "example/dynamicModelPlanarDualArmV1.hpp"
#include "example/trajectoryGeneration.hpp"


MODEL_BASED_CONTROLLER::MODEL_BASED_CONTROLLER() {
    _gazeboClockSub = _nh.subscribe("/clock",0,&MODEL_BASED_CONTROLLER::gazeboClockCB,this);                                        // Gazebo Clock Subscriber
    _jointStateSub = _nh.subscribe("/planar/joint_states",0,&MODEL_BASED_CONTROLLER::jointStateCB,this);                            // Link State subscriber
    _jLeftArmEffPub = _nh.advertise<std_msgs::Float64>("/planar/joint_left_x_effort_eff_controller/command", 0);		                        // Publisher adversisers
    _jRightArmEffPub = _nh.advertise<std_msgs::Float64>("/planar/joint_right_x_effort_eff_controller/command", 0);	                            // Publisher adversisers

    _setModelConfigurationClient = _nh.serviceClient<gazebo_msgs::SetModelConfiguration>("/gazebo/set_model_configuration");        // Set Model Configuration Client
    _pauseGazebo = _nh.serviceClient<std_srvs::Empty>("/gazebo/pause_physics");

    // Definition of the service request message
    _configSrv.request.model_name = "Planar";                               // Model name
    _configSrv.request.urdf_param_name = "robot_description";               // URDF param name (unused)
    
    // Definition of the joint names and initial position vectors
    std::vector<std::string> joints_name_vec;                               // Joint names [rad] vector
    std::vector<double> joints_pos_vec;                                     // Joint position [rad] vector
    joints_name_vec.insert(joints_name_vec.end(), {"joint_cables_x", "joint_left_x", "joint_right_x"});
    joints_pos_vec.insert(joints_pos_vec.end(), {0.3, 0, 0});                // [rad]
    // joints_pos_vec.insert(joints_pos_vec.end(), {0.0, 0, 0});             // [rad]

    // Configuration of the service
    _configSrv.request.joint_names = joints_name_vec;                       // Configuring the server with joints name list
    _configSrv.request.joint_positions = joints_pos_vec;

}



void MODEL_BASED_CONTROLLER::gazeboClockCB(rosgraph_msgs::Clock clockMsg) {
    double sec = clockMsg.clock.sec;
    double nsec = clockMsg.clock.nsec;
    double num = 1000000000;
    _gazeboTime = sec+nsec/num;
}



void MODEL_BASED_CONTROLLER::jointStateCB(sensor_msgs::JointState jState) {
    _jCablesPos = jState.position[0];
    _jLeftArmPos = jState.position[1];
    _jRightArmPos = jState.position[2];

    _jCablesVel = jState.velocity[0];
    _jLeftArmVel = jState.velocity[1];
    _jRightArmVel = jState.velocity[2];
}



void MODEL_BASED_CONTROLLER::setInitialState() {
    ros::Duration(2).sleep();                                                           // Sleep for half a second
    std::cout << "setInitialState START" << std::endl;
    _setModelConfigurationClient.call(_configSrv);
    std::cout << "setInitialState END" << std::endl << std::endl;
}


void MODEL_BASED_CONTROLLER::Bfun(Eigen::MatrixXd& B, Eigen::VectorXd q) {
  B_row1(B,q);
  B_row2(B,q);
  B_row3(B,q);
}


void MODEL_BASED_CONTROLLER::nfun(Eigen::VectorXd& n, Eigen::VectorXd q, Eigen::VectorXd qDot) {
  n_row1(n,q,qDot);
  n_row2(n,q,qDot);
  n_row3(n,q,qDot);
}



// NONCOLLOCATED PARTIAL FBL (CONTROL THE PASSIVE LINK TAKING INTO ACCOUNT DYNAMIC COUPLING WITH THE ACTIVE)
void MODEL_BASED_CONTROLLER::noncollocatedPFBLControl() {

    ros::Rate r(1000);
    
    Eigen::VectorXd qPassive(1);
    Eigen::VectorXd qDotPassive(1);
    Eigen::VectorXd qActive(2);
    Eigen::VectorXd qDotActive(2);
    Eigen::VectorXd in(1);
    Eigen::VectorXd y(2);
    Eigen::MatrixXd Kp(1,1);
    Eigen::MatrixXd Kd(1,1);
    Eigen::MatrixXd Gp(2,2);
    Eigen::MatrixXd Gd(2,2);
    Eigen::MatrixXd B(3,3);
    Eigen::VectorXd n(3);
    Eigen::MatrixXd Bb(1,1);
    Eigen::MatrixXd Bbm(1,2);
    Eigen::MatrixXd pinvBbm(2,1);
    Eigen::MatrixXd Bmb(2,1);
    Eigen::MatrixXd Bm(2,2);
    Eigen::VectorXd nb(1);
    Eigen::VectorXd nm(2);
    Eigen::MatrixXd B_tilde(2,2);
    Eigen::VectorXd n_tilde(2);
    Eigen::VectorXd B_hat(2);
    Eigen::VectorXd n_hat(2);

    Eigen::VectorXd e(1);
    Eigen::VectorXd eDot(1);
    Eigen::VectorXd tau(2);
    std_msgs::Float64 jLeftCommand;
    std_msgs::Float64 jRightCommand;

    // // Data save variables
    // std::vector<double> q1Saved, q2Saved, tauSaved;

    Kp.diagonal() << 2000;
    Kd.diagonal() << sqrt(4*Kp(0,0));

    std::vector<double> tDes, qDesPassive, qDotDesPassive, qDDotDesPassive;
    // qDesPassive.push_back(0.0);
    // qDotDesPassive.push_back(0.0);
    // qDDotDesPassive.push_back(0.0);

    // cubicVelTraj(tDes, qDesPassive, qDotDesPassive, qDDotDesPassive, 0.001, 10, 0.1, 0.0, 0.0, 0.0);
    // cubicVelTraj(tDes, qDesPassive, qDotDesPassive, qDDotDesPassive, 0.001,  1, 0.3, 0.0, 0.0, 0.0);
    // cubicVelTraj(tDes, qDesPassive, qDotDesPassive, qDDotDesPassive, 0.001,  3, 0.0, 0.0, 0.0, 0.0);
    // cubicVelTraj(tDes, qDesPassive, qDotDesPassive, qDDotDesPassive, 0.001,  2, 0.0, 0.1, 0.0, 0.0);
    // for (int i=0; i<30000; i++) {
    //     qDesPassive.push_back(qDesPassive[qDesPassive.size()-1]);
    //     qDotDesPassive.push_back(qDotDesPassive[qDotDesPassive.size()-1]);
    //     qDDotDesPassive.push_back(qDDotDesPassive[qDDotDesPassive.size()-1]);
    // }

    // dampedSinusoidalTraj(tDes, qDesPassive, qDotDesPassive, qDDotDesPassive, 0.001, 10, 0.3, 2, 0.0, 0.4);
    dampedSinusoidalTraj(tDes, qDesPassive, qDotDesPassive, qDDotDesPassive, 0.001, 20, 0.3, 0.8, 0.0, 0.6);

    for (int i = 0; i < tDes.size(); i++) {
        std::cout << tDes[i] << ", " << qDesPassive[i] << ", " << qDotDesPassive[i] << ", " << qDDotDesPassive[i] << std::endl;
    }

    std::cout << qDesPassive.size() << std::endl;

    for (int i=0; i<qDesPassive.size(); i++) {

    // while(ros::ok()) {

        if (i>30) {

            Eigen::VectorXd q(3);
            Eigen::VectorXd qDot(3);
            Eigen::VectorXd qPassive(1);
            Eigen::VectorXd qDotPassive(1);
            Eigen::VectorXd qActive(2);
            Eigen::VectorXd qDotActive(2);

            q << _jCablesPos, _jLeftArmPos, _jRightArmPos;
            qDot << _jCablesVel, _jLeftArmVel, _jRightArmVel;

            qPassive << q(0);
            qDotPassive << qDot(0);

            qActive << q(1), q(2);
            qDotActive << qDot(1), qDot(2);

            Eigen::VectorXd qDesPassiveEigen(1);
            Eigen::VectorXd qDotDesPassiveEigen(1);
            Eigen::VectorXd qDDotDesPassiveEigen(1);

            // qDesPassiveEigen << qDesPassive[0];
            // qDotDesPassiveEigen << qDotDesPassive[0];
            // qDDotDesPassiveEigen << qDDotDesPassive[0];

            qDesPassiveEigen << qDesPassive[i];
            qDotDesPassiveEigen << qDotDesPassive[i];
            qDDotDesPassiveEigen << qDDotDesPassive[i];

            // Computation auxiliar input
            e = qDesPassiveEigen - qPassive;
            std::cout << "Error: " << e.transpose() << std::endl;
            eDot = qDotDesPassiveEigen - qDotPassive;
            in = Kp*e + Kd*eDot + qDDotDesPassiveEigen;

            // Computation dynamic model
            B = Eigen::MatrixXd::Zero(3,3);
            n = Eigen::VectorXd::Zero(3);
            Bfun(B,q);
            nfun(n,q,qDot);

            // Matrix decomposition
            Bb = B.block<1,1>(0,0);
            Bbm = B.block<1,2>(0,1);
            Bmb = Bbm.transpose();
            Bm = B.block<2,2>(1,1);
            nb = n.head(1);
            nm = n.tail(2);
            pinvBbm = pseudoinverse(Bbm,0);

            // Composed matrix computation
            B_hat = Bmb-Bm*pinvBbm*Bb;
            n_hat = nm-Bm*pinvBbm*nb;

            // Control law
            tau = B_hat*in + n_hat;

            // Publishing the commands
            jLeftCommand.data = tau[0];
            jRightCommand.data = tau[1];
            _jLeftArmEffPub.publish(jLeftCommand);
            _jRightArmEffPub.publish(jRightCommand);

        //     // Saving data
        //     q1Saved.push_back(q(0));
        //     q2Saved.push_back(q(1));
        //     tauSaved.push_back(tau(0));

        }

        r.sleep();

    }

    // std::ofstream output_file("/home/giancarlo/ros_ws/src/example/src/DataExampleNoncollocated.txt", std::ofstream::out);
    // if (output_file.is_open()) {
    //     for (int i=0; i<q1Saved.size(); i++) {
    //         output_file << q1Saved[i] << " " << q2Saved[i] << " " << qDesPassive[i] << " " << tauSaved[i] << "\n";
    //     }
    //     output_file.close();
    // }
    // else std::cout << "Problem with opening file";

    std_srvs::Empty pauseSrv;
    _pauseGazebo.call(pauseSrv);
}


void MODEL_BASED_CONTROLLER::run() {
    std::thread noncollocatedPFBLControl_t(&MODEL_BASED_CONTROLLER::noncollocatedPFBLControl,this);
    ros::spin();
}



// Main loop
int main (int argc, char** argv) {
    ros::init(argc, argv, "model_based_control");
    MODEL_BASED_CONTROLLER controller;
    ros::Duration(0.5).sleep();         // Sleep for half a second
    controller.setInitialState();       // Modify the starting position of the system
    controller.run();                   // Start the identification and control
    return 0;
}
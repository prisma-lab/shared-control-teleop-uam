# Overview

This repository contains the code for the simulations on Shared-Control Teleoperation Methods for a Cable-Suspended Dual-Arm Unmanned Aerial Manipulator. The code has been tested on Ubuntu 20.04 running ROS Noetic.

# How to Launch the Simulations

## Push Wall Task

1. Open a terminal and run:
   ```bash
   roslaunch uav-cable-suspended-robots-ros-pkg push_wall.launch
   ```
2. Once Gazebo opens, unpause the simulation.
3. Wait until the joint controllers are loaded (indicated by the terminal printing `0` repeatedly).
4. Press `Enter` in the terminal to bring the arms to their initial position.
5. Once the arms have reached their position, they can be controlled with the joystick.
6. Open another terminal and run:
   ```bash
   rosrun hierarchical_ctrl_pkg sensors_hierarchical_ctrl_node
   ```
7. Wait until lift-off is completed; the drone will start hovering in place.
8. Open another terminal and run:
   ```bash
   rosrun hierarchical_ctrl_pkg publish_trajectory_node 0
   ```
9. Once the trajectory is completed, the robot will push against the wall.
10. Open another terminal and run:
    ```bash
    rqt
    ```
11. Open the **Message Publisher** plugin and create a publisher for the topic:
    ```
    /licasa1/activate_control
    ```
12. The topic requires 4 terms:
    - First term: `1` to activate force control, `0` to deactivate.
    - Second term: Absolute value of the desired force.
    - Third term: Desired angle between the reaction force and the heading direction (e.g., `3.14` to push forward).
    - Fourth term: Always set to `1` (set to `0` to make the robot follow the trajectory in reverse).

---

## Grasp Object Task

1. Open a terminal and run:
   ```bash
   roslaunch uav-cable-suspended-robots-ros-pkg grasp_obj.launch
   ```
2. Once Gazebo opens, unpause the simulation.
3. Wait until the joint controllers are loaded (indicated by the terminal printing `0` repeatedly).
4. Press `Enter` in the terminal to bring the arms to their initial position.
5. Once the arms have reached their position, they can be controlled with the joystick.
6. Open another terminal and run:
   ```bash
   rosrun hierarchical_ctrl_pkg sensors_hierarchical_ctrl_node
   ```
7. Wait until lift-off is completed; the drone will start hovering in place.
8. Open another terminal and run:
   ```bash
   rosrun hierarchical_ctrl_pkg publish_trajectory_node 1
   ```
9. Once the trajectory is completed, the robot will be a few centimeters away from the cube.
10. Press the **3rd button** on the joystick to activate world-based control for the arms and enable the shared-control algorithm for the drone.
11. Use the joystick to move the arms and drone toward the cube.
12. Once the end-effectors are on the sides of the cube, press the **1st button** on the joystick **four times** to close the arms around the cube.
13. Move the arms upwards to lift the cube.
14. Manipulate the cube as desired.

---

# How to Use the Joystick

- The sticks of the joystick move the reference frame that the end-effectors attempt to follow.
- The joystick must be set to **analog** mode.
- Reference frame orientation (from the drone's point of view):
  - **X-axis**: Forward
  - **Y-axis**: Left
  - **Z-axis**: Up

## Controls

- **Left Stick**: Moves the reference frame in the **Y-Z** plane.
- **Right Stick**: Moves the reference frame along the **X** axis.

## Buttons

- **1st Button** (X button on Logitech joystick): Decreases the distance between the arms.
- **2nd Button** (A button): Increases the distance between the arms.
- **3rd Button** (B button): Changes the arms' control mode (i.e., updates the `ref state`).

## `ref state` Behavior

- `ref state = 0`: Arm movement does **not** affect the drone.
- `ref state = 1`: Arm movement **also** controls the drone through the shared-control algorithm.

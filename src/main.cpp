#include "mecanum_controller.hpp"

int main() {
  constexpr double R = 0.20;
  constexpr double r = 0.05;
  std::array<WheelModel, 3> wheels = MecanumController::makeEquilateralTangentWheels(R);

  MecanumController robot(r, wheels);

  const double dt = 0.01;
  Pose2D pose = {0.0, 0.0, 0.0};
  std::array<double, 3> wheel_vel = {3.0, 3.0, 3.0};

  // Forward kinematics
  TwistBody twist = robot.forwardKinematics(wheel_vel);
  std::cout << "Body twist:\n";
  std::cout << "  vx    = " << twist.vx << " m/s\n";
  std::cout << "  vy    = " << twist.vy << " m/s\n";
  std::cout << "  omega = " << twist.omega << " rad/s\n";

  // Run simulation for 1 second.
  for (int k = 0; k < 100; ++k) {
    MecanumController::integrateMidpointYaw(pose, twist, dt);
  }

  std::cout << "Pose after 1.0s:\n";
  std::cout << "  x    = " << pose.x << " m\n";
  std::cout << "  y    = " << pose.y << " m\n";
  std::cout << "  theta = " << pose.theta << " rad\n";

  // Inverse kinematics
  std::array<double, 3> w = robot.computeWheelOmegas(twist);

  std::cout << "Wheel angular velocities (rad/s):\n";
  std::cout << "  w1 = " << w[0] << "\n";
  std::cout << "  w2 = " << w[1] << "\n";
  std::cout << "  w3 = " << w[2] << "\n";

  return 0;
}

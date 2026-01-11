/*
3-wheel (120 deg) omni/mecanum-like inverse kinematics (IK)

Goal:
  Given desired body twist (vx, vy, omega) in robot body frame,
  compute wheel angular velocities (rad/s) for 3 wheels.

General IK:
  s_i = t_i^T ( v + omega * (k x p_i) )
      = cos(alpha_i)*vx + sin(alpha_i)*vy
        + omega * ( -y_i*cos(alpha_i) + x_i*sin(alpha_i) )

  wheel_omega_i = s_i / r

where
  v = [vx, vy]^T  (m/s)
  omega = yaw rate (rad/s)
  p_i = [x_i, y_i]^T wheel position in body frame (m)
  t_i = [cos(alpha_i), sin(alpha_i)]^T wheel driving direction unit vector
  r = wheel radius (m)

This works for any 3-wheel geometry, not only 120 deg.

Also included:
  - Saturation (uniform scaling) to respect max wheel speed

Compile:
  g++ -O2 -std=c++17 three_wheel_ik.cpp -o ik
*/

#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>

struct WheelModel {
  // Wheel position in robot body frame (meters)
  double x;
  double y;

  // Wheel driving direction angle alpha (radians) in robot body frame.
  // t = [cos(alpha), sin(alpha)]
  double alpha;
};

struct TwistBody {
  double vx;     // m/s (forward + x)
  double vy;     // m/s (left + y)
  double omega;  // rad/s (CCW +)
};

class ThreeWheelInverseKinematics {
public:
  ThreeWheelInverseKinematics(double wheel_radius_m, const std::array<WheelModel, 3>& wheels)
      : r_(wheel_radius_m), wheels_(wheels) {
    if (r_ <= 0.0) {
      throw std::invalid_argument("wheel_radius_m must be > 0");
    }
  }

  // Compute wheel angular velocities (rad/s) from desired body twist.
  // If omega_max_rad_s > 0, wheel speeds are uniformly scaled down to satisfy |w_i| <= omega_max_rad_s.
  std::array<double, 3> computeWheelOmegas(const TwistBody& cmd, double omega_max_rad_s = 0.0) const {
    std::array<double, 3> wheel_omega{};

    for (std::size_t i = 0; i < 3; ++i) {
      const double ca = std::cos(wheels_[i].alpha);
      const double sa = std::sin(wheels_[i].alpha);

      // s_i in m/s along wheel driving direction
      const double s_i = ca * cmd.vx + sa * cmd.vy + cmd.omega * (-wheels_[i].y * ca + wheels_[i].x * sa);                                                                                                                                                 

      // wheel angular velocity rad/s
      wheel_omega[i] = s_i / r_;
    }

    if (omega_max_rad_s > 0.0) {
      wheel_omega = saturateUniform_(wheel_omega, omega_max_rad_s);
    }

    return wheel_omega;
  }

private:
  double r_;
  std::array<WheelModel, 3> wheels_;

  static std::array<double, 3> saturateUniform_(const std::array<double, 3>& w, double w_max) {
    double peak = 0.0;
    for (double wi : w) {
      const double a = std::fabs(wi);
      if (a > peak) {
        peak = a;
      }
    }
    if (peak <= w_max || peak < 1e-12) {
      return w;
    }

    const double scale = w_max / peak;
    return {w[0] * scale, w[1] * scale, w[2] * scale};
  }
};

// ---- Convenience: a common 120-deg equilateral / tangent configuration ----
// Many 3-wheel bases place wheels at beta = 0, 120, 240 degrees on a circle of radius R,
// and set the driving direction tangent: alpha = beta + 90 degrees.
//
// For that configuration, you can build WheelModel as:
//
//   p_i = R [cos(beta_i), sin(beta_i)]
//   alpha_i = beta_i + pi/2
//
static std::array<WheelModel, 3> makeEquilateralTangentWheels(double R) {
  const std::array<double, 3> beta = {0.0, 2.0 * M_PI / 3.0, 4.0 * M_PI / 3.0};

  std::array<WheelModel, 3> wheels{};
  for (std::size_t i = 0; i < 3; ++i) {
    wheels[i].x = R * std::cos(beta[i]);
    wheels[i].y = R * std::sin(beta[i]);
    wheels[i].alpha = beta[i] + M_PI / 2.0;  // tangent direction
  }
  return wheels;
}

int main() {
  // Example parameters (replace with your robot values)
  constexpr double wheel_radius = 0.05;  // m
  constexpr double R = 0.20;             // m (center -> wheel distance)
  constexpr double omega_max = 30.0;     // rad/s max wheel speed (optional)

  const auto wheels = makeEquilateralTangentWheels(R);
  ThreeWheelInverseKinematics ik(wheel_radius, wheels);

  // Desired body twist (m/s, m/s, rad/s)
  TwistBody cmd{};
  cmd.vx = -0.57735;     // forward
  cmd.vy = 0.0;     // left
  cmd.omega = 0.0;  // CCW

  const auto w = ik.computeWheelOmegas(cmd, omega_max);

  std::cout << "Wheel angular velocities (rad/s):\n";
  std::cout << "  w1 = " << w[0] << "\n";
  std::cout << "  w2 = " << w[1] << "\n";
  std::cout << "  w3 = " << w[2] << "\n";

  return 0;
}

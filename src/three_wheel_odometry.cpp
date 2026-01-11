/*
3-wheel (120 deg) omni/mecanum-like odometry

This implementation follows the general formulation:

s_i = t_i^T ( v + w * (k x p_i) )

where
- v = [vx, vy] body-frame velocity
- w = yaw rate
- p_i = [x_i, y_i] wheel position in body frame
- t_i = [cos(alpha_i), sin(alpha_i)] wheel driving (tangent) direction unit vector
- s_i = r * wheel_angular_velocity_i  (linear speed along t_i)

Stacking 3 wheels:
s = A * u,  u = [vx, vy, w]^T,  =>  u = A^{-1} s

Then integrate to world:
x += (vx*cos(theta) - vy*sin(theta)) * dt
y += (vx*sin(theta) + vy*cos(theta)) * dt
theta += w * dt

Notes:
- You MUST ensure alpha_i and wheel angular velocity signs match your hardware.
- Units: meters, seconds, radians.
- This code is self-contained and flake8 is irrelevant for C++; it is written in a clean C++ style.

Compile:
g++ -O2 -std=c++17 three_wheel_odometry.cpp -o odo
*/

#include <array>
#include <cmath>
#include <cstdint>
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
  double vx;     // m/s (forward)
  double vy;     // m/s (left)
  double omega;  // rad/s (CCW)
};

struct Pose2D {
  double x;      // m (world)
  double y;      // m (world)
  double theta;  // rad (world yaw, CCW)
};

class ThreeWheelOdometry {
public:
  ThreeWheelOdometry(double wheel_radius_m,
                     const std::array<WheelModel, 3>& wheels)
      : r_(wheel_radius_m), wheels_(wheels) {
    if (r_ <= 0.0) {
      throw std::invalid_argument("Wheel_radius_m must be > 0");
    }
    // Precompute the inverse of A (3x3) from wheel geometry.
    computeAinvOrThrow_();
  }

  // Compute body twist from wheel angular velocities (rad/s).
  TwistBody forwardKinematics(const std::array<double, 3>& wheel_omega_rad_s) const {
    // s_i = r * omega_i
    const std::array<double, 3> s = {
      r_ * wheel_omega_rad_s[0],
      r_ * wheel_omega_rad_s[1],
      r_ * wheel_omega_rad_s[2],
    };

    // u = Ainv * s
    TwistBody u{};
    u.vx = Ainv_[0][0] * s[0] + Ainv_[0][1] * s[1] + Ainv_[0][2] * s[2];
    u.vy = Ainv_[1][0] * s[0] + Ainv_[1][1] * s[1] + Ainv_[1][2] * s[2];
    u.omega = Ainv_[2][0] * s[0] + Ainv_[2][1] * s[1] + Ainv_[2][2] * s[2];

    return u;
  }

  // Integrate pose with body twist for dt seconds (mid-point yaw integration).
  static void integrateMidpointYaw(Pose2D& pose, const TwistBody& twist, double dt) {
    if (dt <= 0.0) {
      return;
    }
    const double theta_mid = pose.theta + 0.5 * twist.omega * dt;
    const double c = std::cos(theta_mid);
    const double s = std::sin(theta_mid);

    const double dx_world = (twist.vx * c - twist.vy * s) * dt;
    const double dy_world = (twist.vx * s + twist.vy * c) * dt;

    pose.x += dx_world;
    pose.y += dy_world;
    pose.theta = normalizeAngle_(pose.theta + twist.omega * dt);
  }

private:
  double r_;
  std::array<WheelModel, 3> wheels_;

  // Inverse of A (3x3)
  std::array<std::array<double, 3>, 3> Ainv_{};

  void computeAinvOrThrow_() {
    // Build A:
    // row i: [cos(ai), sin(ai), (-y_i cos(ai) + x_i sin(ai))]
    std::array<std::array<double, 3>, 3> A{};
    for (std::size_t i = 0; i < 3; ++i) {
      const double ca = std::cos(wheels_[i].alpha);
      const double sa = std::sin(wheels_[i].alpha);
      A[i][0] = ca;
      A[i][1] = sa;
      A[i][2] = (-wheels_[i].y * ca + wheels_[i].x * sa);
    }

    // Invert 3x3 via adjugate / determinant.
    const double det = det3_(A);
    const double eps = 1e-12;
    if (std::fabs(det) < eps) {
      throw std::runtime_error(
        "Wheel geometry matrix A is singular or ill-conditioned (det ~ 0). "
        "Check wheel angles/positions (alpha_i, x_i, y_i).");
    }

    const auto adj = adjugate3_(A);
    for (int r = 0; r < 3; ++r) {
      for (int c = 0; c < 3; ++c) {
        Ainv_[r][c] = adj[r][c] / det;
      }
    }
  }

  static double det3_(const std::array<std::array<double, 3>, 3>& M) {
    // |a b c|
    // |d e f| = a(ei - fh) - b(di - fg) + c(dh - eg)
    // |g h i|
    const double a = M[0][0], b = M[0][1], c = M[0][2];
    const double d = M[1][0], e = M[1][1], f = M[1][2];
    const double g = M[2][0], h = M[2][1], i = M[2][2];
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
  }

  static std::array<std::array<double, 3>, 3> adjugate3_(
      const std::array<std::array<double, 3>, 3>& M) {
    // adj(M) = cofactor(M)^T
    std::array<std::array<double, 3>, 3> C{};

    const double a = M[0][0], b = M[0][1], c = M[0][2];
    const double d = M[1][0], e = M[1][1], f = M[1][2];
    const double g = M[2][0], h = M[2][1], i = M[2][2];

    // Cofactors:
    const double C00 =  (e * i - f * h);
    const double C01 = -(d * i - f * g);
    const double C02 =  (d * h - e * g);

    const double C10 = -(b * i - c * h);
    const double C11 =  (a * i - c * g);
    const double C12 = -(a * h - b * g);

    const double C20 =  (b * f - c * e);
    const double C21 = -(a * f - c * d);
    const double C22 =  (a * e - b * d);

    // Transpose of cofactor matrix:
    C[0][0] = C00; C[0][1] = C10; C[0][2] = C20;
    C[1][0] = C01; C[1][1] = C11; C[1][2] = C21;
    C[2][0] = C02; C[2][1] = C12; C[2][2] = C22;

    return C;
  }

  static double normalizeAngle_(double a) {
    // Normalize to (-pi, pi]
    constexpr double kPi = 3.14159265358979323846;
    constexpr double kTwoPi = 2.0 * kPi;
    while (a <= -kPi) {
      a += kTwoPi;
    }
    while (a > kPi) {
      a -= kTwoPi;
    }
    return a;
  }
};

int main() {
  // Example: equilateral 3-wheel base with wheels at 0, 120, 240 degrees,
  // and wheel driving directions tangent to the circle (alpha = beta + 90deg).
  //
  // Geometry: R = distance from center to wheel (m)
  // Wheel positions p_i = R [cos(beta_i), sin(beta_i)]
  // Driving directions t_i angle alpha_i = beta_i + pi/2
  constexpr double R = 0.20;  // 20 cm
  constexpr double r = 0.05;  // 5 cm wheel radius

  const std::array<double, 3> beta = {0.0, 2.0 * M_PI / 3.0, 4.0 * M_PI / 3.0};

  std::array<WheelModel, 3> wheels{};
  for (std::size_t i = 0; i < 3; ++i) {
    wheels[i].x = R * std::cos(beta[i]);
    wheels[i].y = R * std::sin(beta[i]);
    wheels[i].alpha = beta[i] + M_PI / 2.0;  // tangent direction
  }

  ThreeWheelOdometry odo(r, wheels);

  // Suppose we read wheel angular velocities (rad/s) from encoders.
  // Replace with your measured values (with correct sign).
  std::array<double, 3> wheel_omega = {0.0, 10.0, -10.0};

  TwistBody twist = odo.forwardKinematics(wheel_omega);
  std::cout << "Body twist:\n";
  std::cout << "  vx    = " << twist.vx << " m/s\n";
  std::cout << "  vy    = " << twist.vy << " m/s\n";
  std::cout << "  omega = " << twist.omega << " rad/s\n";

  Pose2D pose{0.0, 0.0, 0.0};
  const double dt = 0.01;
  for (int k = 0; k < 100; ++k) {
    // Update twist each loop in real usage.
    ThreeWheelOdometry::integrateMidpointYaw(pose, twist, dt);
  }

  std::cout << "Pose after 1.0s:\n";
  std::cout << "  x    = " << pose.x << " m\n";
  std::cout << "  y    = " << pose.y << " m\n";
  std::cout << "  theta = " << pose.theta << " rad\n";

  return 0;
}

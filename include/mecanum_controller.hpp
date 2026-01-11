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

class MecanumController {
public:
  MecanumController( double wheel_radius_m, const std::array<WheelModel, 3>& wheels )
      : r_(wheel_radius_m), wheels_(wheels) {
    if (r_ <= 0.0) {
      throw std::invalid_argument("Wheel_radius_m must be > 0");
    }

    // Precompute the inverse of A (3x3) from wheel geometry.
    computeAinvOrThrow_();
  }

  // Forward kinematics
  /**
   * Compute body twist from wheel angular velocities (rad/s).
   */
  TwistBody forwardKinematics( const std::array<double, 3>& wheel_omega_rad_s) const;

  /*
   * Integrate pose with body twist for dt seconds (mid-point yaw integration).
   */
  static void integrateMidpointYaw( Pose2D& pose, const TwistBody& twist, double dt);

  // Inverse kinematics
  /**
   * Compute wheel angular velocities (rad/s) from desired body twist.
   * If omega_max_rad_s > 0, wheel speeds are uniformly scaled down to satisfy |w_i| <= omega_max_rad_s.
   */
  std::array<double, 3> computeWheelOmegas( const TwistBody& cmd, double omega_max_rad_s = 0.0 ) const;

  /**
   * ---- Convenience: a common 120-deg equilateral / tangent configuration ----
     Many 3-wheel bases place wheels at beta = 0, 120, 240 degrees on a circle of radius R,
     and set the driving direction tangent: alpha = beta + 90 degrees.

     For that configuration, you can build WheelModel as:

     p_i = R [cos(beta_i), sin(beta_i)]
     alpha_i = beta_i + pi/2
   */
  static std::array<WheelModel, 3> makeEquilateralTangentWheels( double R );

private:
  // Robot configurations.
  double r_;
  std::array<WheelModel, 3> wheels_;

  // Inverse of A (3x3)
  std::array<std::array<double, 3>, 3> Ainv_{};

  // Forward kinematics
  void computeAinvOrThrow_();
  static double det3_( const std::array<std::array<double, 3>, 3>& M );
  static std::array<std::array<double, 3>, 3> adjugate3_( const std::array<std::array<double, 3>, 3>& M );
  static double normalizeAngle_( double a );

  // Inverse kinematics
  static std::array<double, 3> saturateUniform_( const std::array<double, 3>& w, double w_max );
};
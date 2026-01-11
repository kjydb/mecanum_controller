#include "mecanum_controller.hpp"

TwistBody MecanumController::forwardKinematics( const std::array<double, 3>& wheel_omega_rad_s ) const {
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

void MecanumController::integrateMidpointYaw( Pose2D& pose, const TwistBody& twist, double dt ) {
  if (dt <= 0.0) {
    throw std::runtime_error("dt must be > 0");
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

std::array<double, 3> MecanumController::computeWheelOmegas( const TwistBody& cmd, double omega_max_rad_s ) const {
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

std::array<WheelModel, 3> MecanumController::makeEquilateralTangentWheels( double R ) {
  const std::array<double, 3> beta = {0.0, 2.0 * M_PI / 3.0, 4.0 * M_PI / 3.0};

  std::array<WheelModel, 3> wheels{};
  for (std::size_t i = 0; i < 3; ++i) {
    wheels[i].x = R * std::cos(beta[i]);
    wheels[i].y = R * std::sin(beta[i]);
    wheels[i].alpha = beta[i] + M_PI / 2.0;  // tangent direction
  }

  return wheels;
}

void MecanumController::computeAinvOrThrow_() {
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

double MecanumController::det3_( const std::array<std::array<double, 3>, 3>& M ) {
  // |a b c|
  // |d e f| = a(ei - fh) - b(di - fg) + c(dh - eg)
  // |g h i|
  const double a = M[0][0], b = M[0][1], c = M[0][2];
  const double d = M[1][0], e = M[1][1], f = M[1][2];
  const double g = M[2][0], h = M[2][1], i = M[2][2];

  return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
}

std::array<std::array<double, 3>, 3> MecanumController::adjugate3_( const std::array<std::array<double, 3>, 3>& M ) {
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

double MecanumController::normalizeAngle_( double a ) {
  // Normalize to (-pi, pi]
  while (a <= -M_PI) {
    a += 2 * M_PI;
  }
  while (a > M_PI) {
    a -= 2 * M_PI;
  }

  return a;
}

std::array<double, 3> MecanumController::saturateUniform_( const std::array<double, 3>& w, double w_max ) {
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

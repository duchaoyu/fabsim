// OrthotropicStVKElement.cpp
//
// Author: David Jourdan (david.jourdan@inria.fr)
// Created: 10/30/21

#include "fsim/OrthotropicStVKElement.h"

using namespace Eigen;
namespace fsim
{

OrthotropicStVKElement::OrthotropicStVKElement(const Eigen::Ref<const Mat3<double>> V,
                                                   const Eigen::Vector3i &E,
                                                   double thickness)
{
  using namespace Eigen;
  idx = E;


  // glocal to local frame
  // TODO: this is simply for test now
  Vector3d xaxis;
  xaxis << 1.0, 0.0, 0.0;   // input from the vector field
  xaxis.normalize();

  Vector3d e1 = V.row(idx(0)) - V.row(idx(2));
  Vector3d e2 = V.row(idx(1)) - V.row(idx(2));

  // construct local frame
  // origin
  Vector3d sum = V.row(idx(0)) + V.row(idx(1)) + V.row(idx(2));
  Vector3d origin = sum / 3.0;

  Vector3d zaxis = e1.cross(e2).normalized();

  double dotProduct = xaxis.dot(zaxis);
  double tolerance = 1e-6;
  if (std::abs(dotProduct) < tolerance) {
  } else {
    Vector3d v_parallel = (xaxis.dot(zaxis)) * zaxis;
    Vector3d v_projected = xaxis - v_parallel;
    xaxis << v_projected.normalized();
  }

  Vector3d yaxis = zaxis.cross(xaxis).normalized();
//    Vector3d _zaxis = xaxis.cross(yaxis).normalized();

  Matrix3d R;
  R.col(0) << xaxis;
  R.col(1) << yaxis;
  R.col(2) << zaxis;

  Matrix4d T = Matrix4d::Identity();
  T.block<3, 3>(0, 0) = R;
  T.block<3, 1>(0, 3) = origin;
  Matrix4d T_inverse = T.inverse();
  T_mul = T_inverse * Matrix4d::Identity();


  MatrixXd V_local_xy(3,2);
  V_local_xy = local_xy(V);

  // _R is the initial length
  _R << V_local_xy(1, 1) - V_local_xy(2, 1), V_local_xy(1, 0) - V_local_xy(2, 0),
        V_local_xy(2, 1) - V_local_xy(0, 1), V_local_xy(2, 0) - V_local_xy(0, 0),
        V_local_xy(0, 1) - V_local_xy(1, 1), V_local_xy(0, 0) - V_local_xy(1, 0);
//
  double d = Vector3d(V_local_xy(0, 0), V_local_xy(1, 0), V_local_xy(2, 0)).dot(_R.col(0)); // area of the triangle

  _R /= d;
  coeff = thickness / 2 * std::abs(d);



}

  Eigen::MatrixXd OrthotropicStVKElement::local_xy(const Eigen::Ref<const Mat3<double>> V) const {
    using namespace Eigen;
    MatrixXd V_local_xy(3, 3);

    for (int i = 0; i < 3; ++i) {
      // Convert each 3D point to homogeneous coordinates (Vector4d)
      Matrix<double, 1, 4> V_homogeneous_xy;
      V_homogeneous_xy << V.row(idx(i))[0], V.row(idx(i))[1], V.row(idx(i))[2], 1.0;
      Matrix<double, 1, 4>  V_transformed_xy = V_homogeneous_xy * T_mul.transpose();
      V_local_xy.row(i) << V_transformed_xy.head<3>();  // Store in the output matrix
    }

    // remove the zero z axis
    MatrixXd V_local_xy32(3,2);
    V_local_xy32 << V_local_xy.block<3, 2>(0, 0);

    // deformed geometry
  return V_local_xy32;

  }


Eigen::MatrixXd OrthotropicStVKElement::local_XY(const Eigen::Ref<const Eigen::VectorXd> X) const {
  using namespace Eigen;

  MatrixXd V_local_XY(3, 3);

  for (int i = 0; i < 3; ++i) {
    // Convert each 3D point to homogeneous coordinates (Vector4d)
    Matrix<double, 1, 4> V_homogeneous_XY;
    V_homogeneous_XY << X.segment<1>(3 * idx(i) + 0), X.segment<1>(3 * idx(i) + 1), X.segment<1>(3 * idx(i) + 2), 1.0;
    Matrix<double, 1, 4> V_transformed_XY = V_homogeneous_XY * T_mul.transpose();
    V_local_XY.row(i) << V_transformed_XY.head<3>();
  }
  return V_local_XY;
}

  Eigen::Vector3d OrthotropicStVKElement::strain(const Eigen::Ref<const Eigen::VectorXd> X) const {
    using namespace Eigen;
    Matrix3d P;
    Matrix3d V_local_XY = local_XY(X);
    P.col(0) << V_local_XY.row(0)[0], V_local_XY.row(0)[1], V_local_XY.row(0)[2];
    P.col(1) << V_local_XY.row(1)[0], V_local_XY.row(1)[1], V_local_XY.row(1)[2];
    P.col(2) << V_local_XY.row(2)[0], V_local_XY.row(2)[1], V_local_XY.row(2)[2];

    Matrix<double, 3, 2> F = P * _R;

    Vector3d res;
    res(0) = 0.5 * (F.col(0).dot(F.col(0)) - 1);
    res(1) = 0.5 * (F.col(1).dot(F.col(1)) - 1);
    res(2) = F.col(1).dot(F.col(0));

    return res;
  }

  Eigen::Vector3d OrthotropicStVKElement::stress(const Eigen::Ref<const Eigen::VectorXd> X, const Eigen::Matrix3d &_C) const {
    using namespace Eigen;
    Vector3d S_E = strain(X);
    return _C * S_E;
  }

  double
  OrthotropicStVKElement::energy(const Eigen::Ref<const Eigen::VectorXd> X, const Eigen::Matrix3d &_C, double mass) const {
    using namespace Eigen;
    Vector3d S_E = strain(X); // elastic potential energy
    return coeff *
           (0.5 * S_E.dot(_C * S_E) + 9.8 * mass * (X(3 * idx(0) + 2) + X(3 * idx(1) + 2) + X(3 * idx(2) + 2)) / 3);
  }

  Vec<double, 9>
  OrthotropicStVKElement::gradient(const Eigen::Ref<const Eigen::VectorXd> X, const Eigen::Matrix3d &_C, double mass) const {
    using namespace Eigen;
    Vector3d S = stress(X, _C);
    Matrix2d SMat = (Matrix2d(2, 2) << S(0), S(2), S(2), S(1)).finished();

    Matrix3d P;
    Matrix3d V_local_XY = local_XY(X);
    P.col(0) << V_local_XY.row(0)[0], V_local_XY.row(0)[1], V_local_XY.row(0)[2];
    P.col(1) << V_local_XY.row(1)[0], V_local_XY.row(1)[1], V_local_XY.row(1)[2];
    P.col(2) << V_local_XY.row(2)[0], V_local_XY.row(2)[1], V_local_XY.row(2)[2];

    Matrix<double, 3, 2> F = P * _R;

    Matrix3d grad = coeff * F * (SMat * _R.transpose());
    grad.row(2) += Vector3d::Constant(9.8 * coeff / 3 * mass);
    return Map<Vec < double, 9>>
    (grad.data(), 9);
  }

  Mat<double, 9, 9>
  OrthotropicStVKElement::hessian(const Eigen::Ref<const Eigen::VectorXd> X, const Eigen::Matrix3d &_C, double mass) const {
    using namespace Eigen;
    Vector3d S = stress(X, _C);

    Matrix3d P;
    Matrix3d V_local_XY = local_XY(X);
    P.col(0) << V_local_XY.row(0)[0], V_local_XY.row(0)[1], V_local_XY.row(0)[2];
    P.col(1) << V_local_XY.row(1)[0], V_local_XY.row(1)[1], V_local_XY.row(1)[2];
    P.col(2) << V_local_XY.row(2)[0], V_local_XY.row(2)[1], V_local_XY.row(2)[2];
    Matrix<double, 3, 2> F = P * _R;

    Matrix3d A = _C(0, 0) * F.col(0) * F.col(0).transpose() + S(0) * Matrix3d::Identity() +
                 _C(2, 2) * F.col(1) * F.col(1).transpose();
    Matrix3d B = _C(1, 1) * F.col(1) * F.col(1).transpose() + S(1) * Matrix3d::Identity() +
                 _C(2, 2) * F.col(0) * F.col(0).transpose();
    Matrix3d C = _C(0, 1) * F.col(0) * F.col(1).transpose() + S(2) * Matrix3d::Identity() +
                 _C(2, 2) * F.col(1) * F.col(0).transpose();

    Matrix<double, 9, 9> hess;

    for (int i = 0; i < 3; ++i)
      for (int j = i; j < 3; ++j)
        hess.block<3, 3>(3 * i, 3 * j) =
            _R(i, 0) * _R(j, 0) * A + _R(i, 1) * _R(j, 1) * B + _R(i, 0) * _R(j, 1) * C +
            _R(i, 1) * _R(j, 0) * C.transpose();

    hess.block<3, 3>(3, 0) = hess.block<3, 3>(0, 3).transpose();
    hess.block<3, 6>(6, 0) = hess.block<6, 3>(0, 6).transpose();
    return coeff * hess;

  } // namespace fsim
}
// OrthotropicStVKElement.cpp
//
// Author: David Jourdan (david.jourdan@inria.fr)
// Created: 10/30/21

#include "fsim/OrthotropicStVKElement.h"
#include <iostream>

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

  T = Matrix4d::Identity();
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
  area = 1.0 / 2 * std::abs(d);



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
//    std::cout << res << " strain"<< std::endl;
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
    if (add_pressure) {
      Vector3d e1 = X.segment<3>(3 * idx(0)) - X.segment<3>(3 * idx(2));
      Vector3d e2 = X.segment<3>(3 * idx(1)) - X.segment<3>(3 * idx(2));
      Vector3d normal = e1.cross(e2);

      Vector3d totalDisplacement = X.segment<3>(3 * idx(0)) + X.segment<3>(3 * idx(1)) + X.segment<3>(3 * idx(2));

//    Vector3d totalDisplacement = X.segment<3>(0) + X.segment<3>(3) + X.segment<3>(6);
      double avgNormalDisplacementArea = totalDisplacement.dot(normal) / 6;

      // Work done by the perpendicular force, work = F * displacement
//        double workDoneByPressure = coeff * (forcePerUnitArea * avgNormalDisplacementArea);
      double workDoneByPressure = forcePerUnitArea * avgNormalDisplacementArea;


//      Matrix3d V_local_XY = local_XY(X);
//      // Calculate area of the triangular element
//      Vector3d e1 = V_local_XY.row(0) - V_local_XY.row(2);
//      Vector3d e2 = V_local_XY.row(1) - V_local_XY.row(2);
//      Vector3d normal = e1.cross(e2);
//
//      Vector3d totalDisplacement = V_local_XY.row(0) + V_local_XY.row(1) + V_local_XY.row(2);
//      double avgNormalDisplacementArea = totalDisplacement.dot(normal) / 6;
//
//      // Work done by the perpendicular force, work = F * displacement
////        double workDoneByPressure = coeff * (forcePerUnitArea * avgNormalDisplacementArea);
//      double workDoneByPressure = forcePerUnitArea * avgNormalDisplacementArea;
      // free energy function for St. Venant-Kirchhoff Material
      double elasticEnergy = coeff * (0.5 * S_E.dot(_C * S_E));
      double workByGravity = area * 9.8 * mass * (X(3 * idx(0) + 2) + X(3 * idx(1) + 2) + X(3 * idx(2) + 2)) / 3;
//    double workByGravity = coeff * (9.8 * mass * (X(3 * idx(0) + 2) + X(3 * idx(1) + 2) + X(3 * idx(2) + 2)) / 3);
      return elasticEnergy + workByGravity - workDoneByPressure;  // Additional gravitational potential energy, +2 for the z-coordinate
    }
//    std::cout << coeff * (0.5 * S_E.dot(_C * S_E)) << " elastic energy" << std::endl;
    return coeff *
           (0.5 * S_E.dot(_C * S_E)) + area * 9.8 * mass * (X(3 * idx(0) + 2) + X(3 * idx(1) + 2) + X(3 * idx(2) + 2)) / 3;
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

    Matrix<double, 3, 2> F = P * _R; // deformation gradient

    Matrix3d grad = coeff * F * (SMat * _R.transpose()); // elastic energy gradient in the local frame
//    std::cout << grad << " gradient elastic in local frame" << std::endl;

    Matrix4d T_mul_local_global = Matrix4d::Identity().inverse() * T;
    Matrix3d grad_33_world;
    grad_33_world << grad;

    for (int i = 0; i < 3; ++i) {
      // Convert each 3D point to homogeneous coordinates (Vector4d)
      Matrix<double, 1, 4> grad14;
      grad14 << grad.col(i)[0], grad.col(i)[1], grad.col(i)[2], 1.0;
      Matrix<double, 1, 4>  grad14_xy = grad14 * T_mul_local_global.transpose();
      grad_33_world.col(i) = grad14_xy.head<3>();
    }
//    std::cout << T << " Transformation matrix" << std::endl;


    grad_33_world.row(2) += Vector3d::Constant(9.8 * area / 3 * mass); // gravity gradient
//    std::cout << grad << " total gradient in local frame" << std::endl;

    Vec<double, 9> grad_flat = Map<Vec<double, 9>>(grad_33_world.data(), 9);
//    std::cout << grad_flat << " flat total gradient in local frame" << std::endl;

//    if (add_pressure) {
//      Eigen::Vector<double, 9> pressureGrad;
//      // matlab
//      double x0_x = X(3 * idx(0) + 0);
//      double x0_y = X(3 * idx(0) + 1);
//      double x0_z = X(3 * idx(0) + 2);
//      double x1_x = X(3 * idx(1) + 0);
//      double x1_y = X(3 * idx(1) + 1);
//      double x1_z = X(3 * idx(1) + 2);
//      double x2_x = X(3 * idx(2) + 0);
//      double x2_y = X(3 * idx(2) + 1);
//      double x2_z = X(3 * idx(2) + 2);
//
//      double T_mul1_1 = T_mul(0, 0);
//      double T_mul1_2 = T_mul(0, 1);
//      double T_mul1_3 = T_mul(0, 2);
//      double T_mul1_4 = T_mul(0, 3);
//      double T_mul2_1 = T_mul(1, 0);
//      double T_mul2_2 = T_mul(1, 1);
//      double T_mul2_3 = T_mul(1, 2);
//      double T_mul2_4 = T_mul(1, 3);
//      double T_mul3_1 = T_mul(2, 0);
//      double T_mul3_2 = T_mul(2, 1);
//      double T_mul3_3 = T_mul(2, 2);
//      double T_mul3_4 = T_mul(2, 3);
//      double T_mul4_1 = T_mul(3, 0);
//      double T_mul4_2 = T_mul(3, 1);
//      double T_mul4_3 = T_mul(3, 2);
//      double T_mul4_4 = T_mul(3, 3);
//
//      double t2 = T_mul1_1*x0_x;
//      double t3 = T_mul1_1*x1_x;
//      double t4 = T_mul1_2*x0_x;
//      double t5 = T_mul1_1*x2_x;
//      double t6 = T_mul1_2*x1_x;
//      double t7 = T_mul1_3*x0_x;
//      double t8 = T_mul1_2*x2_x;
//      double t9 = T_mul1_3*x1_x;
//      double t10 = T_mul1_3*x2_x;
//      double t11 = T_mul2_1*x0_y;
//      double t12 = T_mul2_1*x1_y;
//      double t13 = T_mul2_2*x0_y;
//      double t14 = T_mul2_1*x2_y;
//      double t15 = T_mul2_2*x1_y;
//      double t16 = T_mul2_3*x0_y;
//      double t17 = T_mul2_2*x2_y;
//      double t18 = T_mul2_3*x1_y;
//      double t19 = T_mul2_3*x2_y;
//      double t20 = T_mul3_1*x0_z;
//      double t21 = T_mul3_1*x1_z;
//      double t22 = T_mul3_2*x0_z;
//      double t23 = T_mul3_1*x2_z;
//      double t24 = T_mul3_2*x1_z;
//      double t25 = T_mul3_3*x0_z;
//      double t26 = T_mul3_2*x2_z;
//      double t27 = T_mul3_3*x1_z;
//      double t28 = T_mul3_3*x2_z;
//      double t29 = T_mul4_1*3.0;
//      double t30 = T_mul4_2*3.0;
//      double t31 = T_mul4_3*3.0;
//      double t32 = -t5;
//      double t33 = -t8;
//      double t34 = -t10;
//      double t35 = -t14;
//      double t36 = -t17;
//      double t37 = -t19;
//      double t38 = -t23;
//      double t39 = -t26;
//      double t40 = -t28;
//      double t92 = t2+t3+t5+t11+t12+t14+t20+t21+t23+t29;
//      double t93 = t4+t6+t8+t13+t15+t17+t22+t24+t26+t30;
//      double t94 = t7+t9+t10+t16+t18+t19+t25+t27+t28+t31;
//      double t41 = t2+t11+t20+t32+t35+t38;
//      double t42 = t3+t12+t21+t32+t35+t38;
//      double t43 = t4+t13+t22+t33+t36+t39;
//      double t44 = t6+t15+t24+t33+t36+t39;
//      double t45 = t7+t16+t25+t34+t37+t40;
//      double t46 = t9+t18+t27+t34+t37+t40;
//      double t47 = T_mul1_2*t41;
//      double t48 = T_mul1_3*t41;
//      double t49 = T_mul1_2*t42;
//      double t50 = T_mul1_3*t42;
//      double t51 = T_mul1_1*t43;
//      double t52 = T_mul1_3*t43;
//      double t53 = T_mul1_1*t44;
//      double t54 = T_mul1_3*t44;
//      double t55 = T_mul2_2*t41;
//      double t56 = T_mul1_1*t45;
//      double t57 = T_mul2_3*t41;
//      double t58 = T_mul1_2*t45;
//      double t59 = T_mul2_2*t42;
//      double t60 = T_mul1_1*t46;
//      double t61 = T_mul2_3*t42;
//      double t62 = T_mul1_2*t46;
//      double t63 = T_mul2_1*t43;
//      double t64 = T_mul2_3*t43;
//      double t65 = T_mul2_1*t44;
//      double t66 = T_mul2_3*t44;
//      double t67 = T_mul3_2*t41;
//      double t68 = T_mul2_1*t45;
//      double t69 = T_mul3_3*t41;
//      double t70 = T_mul2_2*t45;
//      double t71 = T_mul3_2*t42;
//      double t72 = T_mul2_1*t46;
//      double t73 = T_mul3_3*t42;
//      double t74 = T_mul2_2*t46;
//      double t75 = T_mul3_1*t43;
//      double t76 = T_mul3_3*t43;
//      double t77 = T_mul3_1*t44;
//      double t78 = T_mul3_3*t44;
//      double t79 = T_mul3_1*t45;
//      double t80 = T_mul3_2*t45;
//      double t81 = T_mul3_1*t46;
//      double t82 = T_mul3_2*t46;
//      double t95 = t41*t44;
//      double t96 = t42*t43;
//      double t97 = t41*t46;
//      double t98 = t42*t45;
//      double t99 = t43*t46;
//      double t100 = t44*t45;
//      double t83 = -t51;
//      double t84 = -t56;
//      double t85 = -t58;
//      double t86 = -t63;
//      double t87 = -t68;
//      double t88 = -t70;
//      double t89 = -t75;
//      double t90 = -t79;
//      double t91 = -t80;
//      double t101 = -t96;
//      double t102 = -t98;
//      double t103 = -t100;
//      double t104 = t95+t101;
//      double t105 = t97+t102;
//      double t106 = t99+t103;
//      double t107 = (T_mul1_3*t104)/6.0;
//      double t108 = (T_mul2_3*t104)/6.0;
//      double t109 = (T_mul1_2*t105)/6.0;
//      double t110 = (T_mul3_3*t104)/6.0;
//      double t111 = (T_mul2_2*t105)/6.0;
//      double t112 = (T_mul1_1*t106)/6.0;
//      double t113 = (T_mul3_2*t105)/6.0;
//      double t114 = (T_mul2_1*t106)/6.0;
//      double t115 = (T_mul3_1*t106)/6.0;
//      double t116 = -t109;
//      double t117 = -t111;
//      double t118 = -t113;
//      pressureGrad(0) = forcePerUnitArea*(t107+t112+t116-(t94*(t49-t53))/6.0+(t93*(t50-t60))/6.0-(t92*(t54-t62))/6.0);
//      pressureGrad(1) = forcePerUnitArea*(t108+t114+t117-(t94*(t59-t65))/6.0+(t93*(t61-t72))/6.0-(t92*(t66-t74))/6.0);
//      pressureGrad(2) = forcePerUnitArea*(t110+t115+t118-(t94*(t71-t77))/6.0+(t93*(t73-t81))/6.0-(t92*(t78-t82))/6.0);
//      pressureGrad(3) = forcePerUnitArea*(t107+t112+t116+(t94*(t47+t83))/6.0-(t93*(t48+t84))/6.0+(t92*(t52+t85))/6.0);
//      pressureGrad(4) = forcePerUnitArea*(t108+t114+t117+(t94*(t55+t86))/6.0-(t93*(t57+t87))/6.0+(t92*(t64+t88))/6.0);
//      pressureGrad(5) = forcePerUnitArea*(t110+t115+t118+(t94*(t67+t89))/6.0-(t93*(t69+t90))/6.0+(t92*(t76+t91))/6.0);
//      pressureGrad(6) = forcePerUnitArea*(t107+t112+t116-(t94*(t47-t49+t53+t83))/6.0+(t93*(t48-t50+t60+t84))/6.0-(t92*(t52-t54+t62+t85))/6.0);
//      pressureGrad(7)= forcePerUnitArea*(t108+t114+t117-(t94*(t55-t59+t65+t86))/6.0+(t93*(t57-t61+t72+t87))/6.0-(t92*(t64-t66+t74+t88))/6.0);
//      pressureGrad(8) = forcePerUnitArea*(t110+t115+t118-(t94*(t67-t71+t77+t89))/6.0+(t93*(t69-t73+t81+t90))/6.0-(t92*(t76-t78+t82+t91))/6.0);
//
//      grad_flat -= pressureGrad;
//
//    }
    if (add_pressure) {
      Eigen::Vector<double, 9> pressureGrad;

      //matlab
      double x0_x = X(3 * idx(0) + 0);
      double x0_y = X(3 * idx(0) + 1);
      double x0_z = X(3 * idx(0) + 2);
      double x1_x = X(3 * idx(1) + 0);
      double x1_y = X(3 * idx(1) + 1);
      double x1_z = X(3 * idx(1) + 2);
      double x2_x = X(3 * idx(2) + 0);
      double x2_y = X(3 * idx(2) + 1);
      double x2_z = X(3 * idx(2) + 2);

      double t2 = x0_x + x1_x + x2_x;
      double t3 = x0_y + x1_y + x2_y;
      double t4 = x0_z + x1_z + x2_z;
      double t5 = -x1_x;
      double t6 = -x1_y;
      double t7 = -x2_x;
      double t8 = -x1_z;
      double t9 = -x2_y;
      double t10 = -x2_z;
      double t11 = t5 + x0_x;
      double t12 = t7 + x0_x;
      double t13 = t6 + x0_y;
      double t14 = t7 + x1_x;
      double t15 = t9 + x0_y;
      double t16 = t8 + x0_z;
      double t17 = t9 + x1_y;
      double t18 = t10 + x0_z;
      double t19 = t10 + x1_z;
      double t20 = (t12 * t17) / 6.0;
      double t21 = (t14 * t15) / 6.0;
      double t22 = (t12 * t19) / 6.0;
      double t23 = (t14 * t18) / 6.0;
      double t24 = (t15 * t19) / 6.0;
      double t25 = (t17 * t18) / 6.0;
      double t26 = -t21;
      double t27 = -t23;
      double t28 = -t25;

      pressureGrad(0) = forcePerUnitArea * (t24 + t28 + (t4 * t17) / 6.0 - (t3 * t19) / 6.0);
      pressureGrad(1) = -1 * forcePerUnitArea * (t22 + t27 + (t4 * t14) / 6.0 - (t2 * t19) / 6.0);
      pressureGrad(2) = forcePerUnitArea * (t20 + t26 + (t3 * t14) / 6.0 - (t2 * t17) / 6.0);
      pressureGrad(3) = forcePerUnitArea * (t24 + t28 - (t4 * t15) / 6.0 + (t3 * t18) / 6.0);
      pressureGrad(4) = -1 * forcePerUnitArea * (t22 + t27 - (t4 * t12) / 6.0 + (t2 * t18) / 6.0);
      pressureGrad(5) = forcePerUnitArea * (t20 + t26 - (t3 * t12) / 6.0 + (t2 * t15) / 6.0);
      pressureGrad(6) = forcePerUnitArea * (t24 + t28 + (t4 * t13) / 6.0 - (t3 * t16) / 6.0);
      pressureGrad(7) = -1 * forcePerUnitArea * (t22 + t27 + (t4 * t11) / 6.0 - (t2 * t16) / 6.0);
      pressureGrad(8) = forcePerUnitArea * (t20 + t26 + (t3 * t11) / 6.0 - (t2 * t13) / 6.0);

      grad_flat -= pressureGrad;
    }

    return grad_flat;
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


    Matrix4d T_mul_local_global = Matrix4d::Identity().inverse() * T;
    Matrix3d T_mul_local_global_R = T_mul_local_global.block<3, 3>(0, 0);

    for (int i = 0; i < 3; ++i)
      for (int j = i; j < 3; ++j){
        Matrix3d hess_block =
            _R(i, 0) * _R(j, 0) * A + _R(i, 1) * _R(j, 1) * B + _R(i, 0) * _R(j, 1) * C +
            _R(i, 1) * _R(j, 0) * C.transpose();

        Matrix3d hess_block_world;
        hess_block_world = T_mul_local_global_R.transpose() * hess_block * T_mul_local_global_R;
//        hess_block_world << hess_block;
//
//        for (int i = 0; i < 3; ++i) {
//          // Convert each 3D point to homogeneous coordinates (Vector4d)
//          Matrix<double, 1, 4> hess14;
//          hess14 << hess_block.col(i)[0], hess_block.col(i)[1], hess_block.col(i)[2], 1.0;
//          Matrix<double, 1, 4>  hess14_xy = hess14 * T_mul_local_global.transpose();
//          hess_block_world.col(i) = hess14_xy.head<3>();
//        }

        hess.block<3, 3>(3 * i, 3 * j) = hess_block_world;
      }

    hess.block<3, 3>(3, 0) = hess.block<3, 3>(0, 3).transpose();
    hess.block<3, 6>(6, 0) = hess.block<6, 3>(0, 6).transpose();
    hess *= coeff;

//    if (add_pressure) {
//      Matrix<double, 9, 9> pressureHess;
//      pressureHess.setZero();
//
//      double x0_x = X(3 * idx(0) + 0);
//      double x0_y = X(3 * idx(0) + 1);
//      double x0_z = X(3 * idx(0) + 2);
//      double x1_x = X(3 * idx(1) + 0);
//      double x1_y = X(3 * idx(1) + 1);
//      double x1_z = X(3 * idx(1) + 2);
//      double x2_x = X(3 * idx(2) + 0);
//      double x2_y = X(3 * idx(2) + 1);
//      double x2_z = X(3 * idx(2) + 2);
//
//      double T_mul1_1 = T_mul(0, 0);
//      double T_mul1_2 = T_mul(0, 1);
//      double T_mul1_3 = T_mul(0, 2);
//      double T_mul1_4 = T_mul(0, 3);
//      double T_mul2_1 = T_mul(1, 0);
//      double T_mul2_2 = T_mul(1, 1);
//      double T_mul2_3 = T_mul(1, 2);
//      double T_mul2_4 = T_mul(1, 3);
//      double T_mul3_1 = T_mul(2, 0);
//      double T_mul3_2 = T_mul(2, 1);
//      double T_mul3_3 = T_mul(2, 2);
//      double T_mul3_4 = T_mul(2, 3);
//      double T_mul4_1 = T_mul(3, 0);
//      double T_mul4_2 = T_mul(3, 1);
//      double T_mul4_3 = T_mul(3, 2);
//      double T_mul4_4 = T_mul(3, 3);
//
//      double t2 = T_mul1_1*T_mul2_2;
//      double t3 = T_mul1_2*T_mul2_1;
//      double t4 = T_mul1_1*T_mul2_3;
//      double t5 = T_mul1_3*T_mul2_1;
//      double t6 = T_mul1_2*T_mul2_3;
//      double t7 = T_mul1_3*T_mul2_2;
//      double t8 = T_mul1_1*T_mul3_2;
//      double    t9 = T_mul1_2*T_mul3_1;
//      double   t10 = T_mul1_1*T_mul3_3;
//      double   t11 = T_mul1_3*T_mul3_1;
//      double   t12 = T_mul1_2*T_mul3_3;
//      double  t13 = T_mul1_3*T_mul3_2;
//      double  t14 = T_mul2_1*T_mul3_2;
//      double  t15 = T_mul2_2*T_mul3_1;
//      double  t16 = T_mul2_1*T_mul3_3;
//      double  t17 = T_mul2_3*T_mul3_1;
//      double   t18 = T_mul2_2*T_mul3_3;
//      double  t19 = T_mul2_3*T_mul3_2;
//      double  t20 = T_mul1_1*x0_x;
//      double  t21 = T_mul1_1*x1_x;
//      double t22 = T_mul1_2*x0_x;
//      double   t23 = T_mul1_1*x2_x;
//      double  t24 = T_mul1_2*x1_x;
//      double  t25 = T_mul1_3*x0_x;
//      double  t26 = T_mul1_2*x2_x;
//      double  t27 = T_mul1_3*x1_x;
//      double  t28 = T_mul1_3*x2_x;
//      double   t29 = T_mul2_1*x0_y;
//      double  t30 = T_mul2_1*x1_y;
//      double  t31 = T_mul2_2*x0_y;
//      double  t32 = T_mul2_1*x2_y;
//      double  t33 = T_mul2_2*x1_y;
//      double  t34 = T_mul2_3*x0_y;
//      double  t35 = T_mul2_2*x2_y;
//      double  t36 = T_mul2_3*x1_y;
//      double  t37 = T_mul2_3*x2_y;
//      double  t38 = T_mul3_1*x0_z;
//      double  t39 = T_mul3_1*x1_z;
//      double  t40 = T_mul3_2*x0_z;
//      double  t41 = T_mul3_1*x2_z;
//      double  t42 = T_mul3_2*x1_z;
//      double  t43 = T_mul3_3*x0_z;
//      double  t44 = T_mul3_2*x2_z;
//      double  t45 = T_mul3_3*x1_z;
//      double  t46 = T_mul3_3*x2_z;
//      double  t47 = T_mul4_1*3.0;
//      double  t48 = T_mul4_2*3.0;
//      double  t49 = T_mul4_3*3.0;
//      double  t50 = -t3;
//      double  t51 = -t5;
//      double  t52 = -t7;
//      double  t53 = -t9;
//      double  t54 = -t11;
//      double  t55 = -t13;
//      double  t56 = -t15;
//      double  t57 = -t17;
//      double  t58 = -t19;
//      double  t59 = -t23;
//      double  t60 = -t26;
//      double  t61 = -t28;
//      double  t62 = -t32;
//      double  t63 = -t35;
//      double  t64 = -t37;
//      double  t65 = -t41;
//      double  t66 = -t44;
//      double  t67 = -t46;
//      double  t146 = t20+t21+t23+t29+t30+t32+t38+t39+t41+t47;
//      double  t147 = t22+t24+t26+t31+t33+t35+t40+t42+t44+t48;
//      double  t148 = t25+t27+t28+t34+t36+t37+t43+t45+t46+t49;
//      double  t68 = t2+t50;
//      double  t69 = t4+t51;
//      double  t70 = t6+t52;
//      double  t71 = t8+t53;
//      double  t72 = t10+t54;
//      double  t73 = t12+t55;
//      double  t74 = t14+t56;
//      double  t75 = t16+t57;
//      double  t76 = t18+t58;
//      double  t77 = t20+t29+t38+t59+t62+t65;
//      double  t78 = t21+t30+t39+t59+t62+t65;
//      double  t79 = t22+t31+t40+t60+t63+t66;
//      double  t80 = t24+t33+t42+t60+t63+t66;
//      double t81 = t25+t34+t43+t61+t64+t67;
//      double  t82 = t27+t36+t45+t61+t64+t67;
//      double  t83 = T_mul1_2*t77;
//      double  t84 = T_mul1_3*t77;
//      double  t85 = T_mul1_2*t78;
//      double  t86 = T_mul1_3*t78;
//      double  t87 = T_mul1_1*t79;
//      double  t88 = T_mul1_3*t79;
//      double  t89 = T_mul1_1*t80;
//      double  t90 = T_mul1_3*t80;
//      double  t91 = T_mul2_2*t77;
//      double  t92 = T_mul1_1*t81;
//      double t93 = T_mul2_3*t77;
//      double  t94 = T_mul1_2*t81;
//      double  t95 = T_mul2_2*t78;
//      double  t96 = T_mul1_1*t82;
//      double  t97 = T_mul2_3*t78;
//      double  t98 = T_mul1_2*t82;
//      double   t99 = T_mul2_1*t79;
//      double  t100 = T_mul2_3*t79;
//      double  t101 = T_mul2_1*t80;
//      double  t102 = T_mul2_3*t80;
//      double  t103 = T_mul3_2*t77;
//      double  t104 = T_mul2_1*t81;
//      double  t105 = T_mul3_3*t77;
//      double  t106 = T_mul2_2*t81;
//      double  t107 = T_mul3_2*t78;
//      double  t108 = T_mul2_1*t82;
//      double  t109 = T_mul3_3*t78;
//      double  t110 = T_mul2_2*t82;
//      double  t111 = T_mul3_1*t79;
//      double t112 = T_mul3_3*t79;
//      double  t113 = T_mul3_1*t80;
//      double   t114 = T_mul3_3*t80;
//      double  t115 = T_mul3_1*t81;
//      double  t116 = T_mul3_2*t81;
//      double  t117 = T_mul3_1*t82;
//      double  t118 = T_mul3_2*t82;
//      double  t149 = (t70*t146)/6.0;
//      double  t150 = (t69*t147)/6.0;
//      double  t151 = (t68*t148)/6.0;
//      double  t152 = (t73*t146)/6.0;
//      double  t153 = (t72*t147)/6.0;
//      double t154 = (t71*t148)/6.0;
//      double  t155 = (t76*t146)/6.0;
//      double  t156 = (t75*t147)/6.0;
//      double t157 = (t74*t148)/6.0;
//      double t119 = -t85;
//      double  t120 = -t86;
//      double  t121 = -t87;
//      double  t122 = -t89;
//      double t123 = -t90;
//      double t124 = -t92;
//      double  t125 = -t94;
//      double  t126 = -t95;
//      double  t127 = -t96;
//      double  t128 = -t97;
//      double  t129 = -t98;
//      double  t130 = -t99;
//      double  t131 = -t101;
//      double  t132 = -t102;
//      double  t133 = -t104;
//      double t134 = -t106;
//      double t135 = -t107;
//      double t136 = -t108;
//      double  t137 = -t109;
//      double  t138 = -t110;
//      double t139 = -t111;
//      double  t140 = -t113;
//      double  t141 = -t114;
//      double  t142 = -t115;
//      double t143 = -t116;
//      double t144 = -t117;
//      double t145 = -t118;
//      double  t158 = -t149;
//      double  t159 = -t150;
//      double  t160 = -t151;
//      double  t161 = -t152;
//      double  t162 = -t153;
//      double t163 = -t154;
//      double  t164 = -t155;
//      double  t165 = -t156;
//      double  t166 = -t157;
//      double  t167 = t83+t121;
//      double  t168 = t85+t122;
//      double t169 = t84+t124;
//      double  t170 = t86+t127;
//      double  t171 = t88+t125;
//      double  t172 = t90+t129;
//      double  t173 = t91+t130;
//      double  t174 = t95+t131;
//      double  t175 = t93+t133;
//      double  t176 = t97+t136;
//      double  t177 = t100+t134;
//      double t178 = t102+t138;
//      double t179 = t103+t139;
//      double  t180 = t107+t140;
//      double  t181 = t105+t142;
//      double  t182 = t109+t144;
//      double t183 = t112+t143;
//      double  t184 = t114+t145;
//      double t185 = (T_mul1_3*t167)/6.0;
//      double  t186 = (T_mul1_2*t169)/6.0;
//      double  t187 = (T_mul1_3*t168)/6.0;
//      double  t188 = (T_mul2_3*t167)/6.0;
//      double  t189 = (T_mul1_1*t171)/6.0;
//      double   t190 = (T_mul1_2*t170)/6.0;
//      double  t191 = (T_mul2_2*t169)/6.0;
//      double  t192 = (T_mul2_3*t168)/6.0;
//      double  t193 = (T_mul1_1*t172)/6.0;
//      double   t194 = (T_mul1_3*t173)/6.0;
//      double t195 = (T_mul3_3*t167)/6.0;
//      double  t196 = (T_mul2_1*t171)/6.0;
//      double  t197 = (T_mul2_2*t170)/6.0;
//      double  t198 = (T_mul1_2*t175)/6.0;
//      double  t199 = (T_mul1_3*t174)/6.0;
//      double  t200 = (T_mul3_2*t169)/6.0;
//      double  t201 = (T_mul3_3*t168)/6.0;
//      double t202 = (T_mul2_1*t172)/6.0;
//      double t203 = (T_mul2_3*t173)/6.0;
//      double t204 = (T_mul1_1*t177)/6.0;
//      double  t205 = (T_mul1_2*t176)/6.0;
//      double  t206 = (T_mul3_1*t171)/6.0;
//      double  t207 = (T_mul3_2*t170)/6.0;
//      double   t208 = (T_mul2_2*t175)/6.0;
//      double  t209 = (T_mul2_3*t174)/6.0;
//      double  t210 = (T_mul1_1*t178)/6.0;
//      double  t211 = (T_mul3_1*t172)/6.0;
//      double  t212 = (T_mul1_3*t179)/6.0;
//      double  t213 = (T_mul3_3*t173)/6.0;
//      double  t214 = (T_mul2_1*t177)/6.0;
//      double  t215 = (T_mul2_2*t176)/6.0;
//      double  t216 = (T_mul1_2*t181)/6.0;
//      double  t217 = (T_mul1_3*t180)/6.0;
//      double  t218 = (T_mul3_2*t175)/6.0;
//      double  t219 = (T_mul3_3*t174)/6.0;
//      double  t220 = (T_mul2_1*t178)/6.0;
//      double  t221 = (T_mul2_3*t179)/6.0;
//      double   t222 = (T_mul1_1*t183)/6.0;
//      double  t223 = (T_mul1_2*t182)/6.0;
//      double  t224 = (T_mul3_1*t177)/6.0;
//      double  t225 = (T_mul3_2*t176)/6.0;
//      double  t226 = (T_mul2_2*t181)/6.0;
//      double  t227 = (T_mul2_3*t180)/6.0;
//      double  t228 = (T_mul1_1*t184)/6.0;
//      double  t229 = (T_mul3_1*t178)/6.0;
//      double  t230 = (T_mul3_3*t179)/6.0;
//      double t231 = (T_mul2_1*t183)/6.0;
//      double  t232 = (T_mul2_2*t182)/6.0;
//      double  t233 = (T_mul3_2*t181)/6.0;
//      double  t234 = (T_mul3_3*t180)/6.0;
//      double  t235 = (T_mul2_1*t184)/6.0;
//      double  t236 = (T_mul3_1*t183)/6.0;
//      double  t237 = (T_mul3_2*t182)/6.0;
//      double  t238 = (T_mul3_1*t184)/6.0;
//      double  t281 = t89+t119+t167;
//      double  t282 = t96+t120+t169;
//      double  t283 = t98+t123+t171;
//      double  t284 = t101+t126+t173;
//      double  t285 = t108+t128+t175;
//      double  t286 = t110+t132+t177;
//      double  t287 = t113+t135+t179;
//      double t288 = t117+t137+t181;
//      double t289 = t118+t141+t183;
//      double  t239 = -t186;
//      double  t240 = -t187;
//      double t241 = -t188;
//      double  t242 = -t190;
//      double  t243 = -t191;
//      double  t244 = -t192;
//      double  t245 = -t193;
//      double   t246 = -t194;
//      double t247 = -t195;
//      double  t248 = -t196;
//      double  t249 = -t197;
//      double t250 = -t198;
//      double  t251 = -t200;
//      double  t252 = -t201;
//      double  t253 = -t202;
//      double  t254 = -t204;
//      double  t255 = -t205;
//      double  t256 = -t206;
//      double   t257 = -t207;
//      double   t258 = -t208;
//      double  t259 = -t209;
//      double  t260 = -t211;
//      double  t261 = -t212;
//      double  t262 = -t213;
//      double  t263 = -t215;
//      double  t264 = -t216;
//      double  t265 = -t218;
//      double  t266 = -t219;
//      double  t267 = -t220;
//      double  t268 = -t221;
//      double  t269 = -t222;
//      double  t270 = -t223;
//      double  t271 = -t224;
//      double  t272 = -t225;
//      double  t273 = -t226;
//      double  t274 = -t229;
//      double  t275 = -t231;
//      double  t276 = -t232;
//      double  t277 = -t233;
//      double  t278 = -t234;
//      double  t279 = -t237;
//      double  t280 = -t238;
//      double  t290 = (T_mul1_3*t281)/6.0;
//      double  t291 = (T_mul2_3*t281)/6.0;
//      double  t292 = (T_mul1_2*t282)/6.0;
//      double  t293 = (T_mul3_3*t281)/6.0;
//      double  t294 = (T_mul2_2*t282)/6.0;
//      double  t295 = (T_mul1_1*t283)/6.0;
//      double  t296 = (T_mul3_2*t282)/6.0;
//      double t297 = (T_mul2_1*t283)/6.0;
//      double t298 = (T_mul1_3*t284)/6.0;
//      double  t299 = (T_mul3_1*t283)/6.0;
//      double  t300 = (T_mul2_3*t284)/6.0;
//      double  t301 = (T_mul1_2*t285)/6.0;
//      double  t302 = (T_mul3_3*t284)/6.0;
//      double  t303 = (T_mul2_2*t285)/6.0;
//      double   t304 = (T_mul1_1*t286)/6.0;
//      double   t305 = (T_mul3_2*t285)/6.0;
//      double  t306 = (T_mul2_1*t286)/6.0;
//      double   t307 = (T_mul1_3*t287)/6.0;
//      double  t308 = (T_mul3_1*t286)/6.0;
//      double  t309 = (T_mul2_3*t287)/6.0;
//      double t310 = (T_mul1_2*t288)/6.0;
//      double t311 = (T_mul3_3*t287)/6.0;
//      double t312 = (T_mul2_2*t288)/6.0;
//      double t313 = (T_mul1_1*t289)/6.0;
//      double t314 = (T_mul3_2*t288)/6.0;
//      double t315 = (T_mul2_1*t289)/6.0;
//      double t316 = (T_mul3_1*t289)/6.0;
//      double t317 = -t290;
//      double t318 = -t292;
//      double  t319 = -t294;
//      double  t320 = -t295;
//      double  t321 = -t296;
//      double  t322 = -t298;
//      double   t323 = -t300;
//      double  t324 = -t301;
//      double  t325 = -t303;
//      double  t326 = -t304;
//      double  t327 = -t305;
//      double  t328 = -t306;
//      double  t329 = -t307;
//      double   t330 = -t309;
//      double  t331 = -t310;
//      double  t332 = -t311;
//      double  t333 = -t312;
//      double  t334 = -t313;
//      double  t335 = -t314;
//      double  t336 = -t315;
//      double  t337 = -t316;
//      double t338 = t188+t194+t196+t204+t243+t250;
//      double  t339 = t192+t199+t202+t210+t249+t255;
//      double  t340 = t195+t206+t212+t222+t251+t264;
//      double  t341 = t201+t211+t217+t228+t257+t270;
//      double  t342 = t213+t221+t224+t231+t265+t273;
//      double  t343 = t219+t227+t229+t235+t272+t276;
//      double  t344 = t185+t189+t190+t239+t240+t245;
//      double  t345 = t203+t214+t215+t258+t259+t267;
//      double  t346 = t230+t236+t237+t277+t278+t280;
//      double  t359 = t149+t151+t159+t191+t199+t210+t241+t248+t255;
//      double  t360 = t149+t151+t159+t194+t197+t204+t244+t250+t253;
//      double  t361 = t152+t154+t162+t200+t217+t228+t247+t256+t270;
//      double  t362 = t152+t154+t162+t207+t212+t222+t252+t260+t264;
//      double t363 = t155+t157+t165+t218+t227+t235+t262+t271+t276;
//      double  t364 = t155+t157+t165+t221+t225+t231+t266+t273+t274;
//      double  t347 = forcePerUnitArea*t338;
//      double t348 = forcePerUnitArea*t339;
//      double  t349 = forcePerUnitArea*t340;
//      double  t350 = forcePerUnitArea*t341;
//      double  t351 = forcePerUnitArea*t342;
//      double  t352 = forcePerUnitArea*t343;
//      double  t353 = forcePerUnitArea*t344;
//      double  t354 = forcePerUnitArea*t345;
//      double  t355 = forcePerUnitArea*t346;
//      double  t365 = forcePerUnitArea*t359;
//      double  t366 = forcePerUnitArea*t360;
//      double   t367 = forcePerUnitArea*t361;
//      double  t368 = forcePerUnitArea*t362;
//      double   t369 = forcePerUnitArea*t363;
//      double  t370 = forcePerUnitArea*t364;
//      double  t374 = t187+t193+t242+t290+t295+t318;
//      double   t375 = t209+t220+t263+t300+t306+t325;
//      double   t376 = t234+t238+t279+t311+t316+t335;
//      double   t377 = t185+t189+t239+t292+t317+t320;
//      double   t378 = t203+t214+t258+t303+t323+t328;
//      double   t379 = t230+t236+t277+t314+t332+t337;
//      double  t389 = t149+t151+t159+t192+t202+t249+t298+t304+t324;
//      double  t390 = t152+t154+t162+t201+t211+t257+t307+t313+t331;
//      double  t391 = t155+t157+t165+t219+t229+t272+t309+t315+t333;
//      double  t392 = t149+t151+t159+t198+t246+t254+t291+t297+t319;
//      double  t393 = t150+t158+t160+t199+t210+t255+t291+t297+t319;
//      double  t394 = t149+t151+t159+t188+t196+t243+t301+t322+t326;
//      double   t395 = t152+t154+t162+t216+t261+t269+t293+t299+t321;
//      double   t396 = t153+t161+t163+t217+t228+t270+t293+t299+t321;
//      double  t397 = t152+t154+t162+t195+t206+t251+t310+t329+t334;
//      double  t398 = t155+t157+t165+t226+t268+t275+t302+t308+t327;
//      double   t399 = t156+t164+t166+t227+t235+t276+t302+t308+t327;
//      double  t400 = t155+t157+t165+t213+t224+t265+t312+t330+t336;
//      double   t422 = t291+t297+t298+t304+t319+t324;
//      double   t423 = t293+t299+t307+t313+t321+t331;
//      double  t424 = t302+t308+t309+t315+t327+t333;
//      double  t356 = -t348;
//      double t357 = -t350;
//      double  t358 = -t352;
//      double  t371 = -t365;
//      double  t372 = -t367;
//      double  t373 = -t369;
//      double   t380 = forcePerUnitArea*t374;
//      double  t381 = forcePerUnitArea*t375;
//      double  t382 = forcePerUnitArea*t376;
//      double  t383 = forcePerUnitArea*t377;
//      double   t384 = forcePerUnitArea*t378;
//      double  t385 = forcePerUnitArea*t379;
//      double   t401 = forcePerUnitArea*t389;
//      double  t402 = forcePerUnitArea*t390;
//      double  t403 = forcePerUnitArea*t391;
//      double  t404 = forcePerUnitArea*t392;
//      double  t405 = forcePerUnitArea*t393;
//      double t406 = forcePerUnitArea*t394;
//      double t407 = forcePerUnitArea*t395;
//      double  t408 = forcePerUnitArea*t396;
//      double  t409 = forcePerUnitArea*t397;
//      double  t410 = forcePerUnitArea*t398;
//      double  t411 = forcePerUnitArea*t399;
//      double  t412 = forcePerUnitArea*t400;
//      double  t425 = forcePerUnitArea*t422;
//      double  t426 = forcePerUnitArea*t423;
//      double   t427 = forcePerUnitArea*t424;
//      double  t386 = -t380;
//      double   t387 = -t381;
//      double    t388 = -t382;
//      double  t413 = -t401;
//      double  t414 = -t402;
//      double  t415 = -t403;
//      double  t416 = -t404;
//      double  t417 = -t405;
//      double  t418 = -t407;
//      double  t419 = -t408;
//      double  t420 = -t410;
//      double  t421 = -t411;
//      double  t428 = -t425;
//      double  t429 = -t426;
//      double  t430 = -t427;
//
//      pressureHess.row(0)[0] = -forcePerUnitArea*((T_mul1_3*t168)/3.0-(T_mul1_2*t170)/3.0+(T_mul1_1*t172)/3.0);
//      pressureHess.row(0)[1] = t356;
//      pressureHess.row(0)[2] = t357;
//      pressureHess.row(0)[3] = t353;
//      pressureHess.row(0)[4] = t366;
//      pressureHess.row(0)[5] = t368;
//      pressureHess.row(0)[6] = t386;
//      pressureHess.row(0)[7] = t413;
//      pressureHess.row(0)[8] = t414;
//      pressureHess.row(1)[0] = t356;
//      pressureHess.row(1)[1] = -forcePerUnitArea*((T_mul2_3*t174)/3.0-(T_mul2_2*t176)/3.0+(T_mul2_1*t178)/3.0);
//      pressureHess.row(1)[2] = t358;
//      pressureHess.row(1)[3] = t371;
//      pressureHess.row(1)[4] = t354;
//      pressureHess.row(1)[5] = t370;
//      pressureHess.row(1)[6] = t417;
//      pressureHess.row(1)[7] = t387;
//      pressureHess.row(1)[8] = t415;
//      pressureHess.row(2)[0] = t357;
//      pressureHess.row(2)[1] = t358;
//      pressureHess.row(2)[2] = -forcePerUnitArea*((T_mul3_3*t180)/3.0-(T_mul3_2*t182)/3.0+(T_mul3_1*t184)/3.0);
//      pressureHess.row(2)[3] = t372;
//      pressureHess.row(2)[4] = t373;
//      pressureHess.row(2)[5] = t355;
//      pressureHess.row(2)[6] = t419;
//      pressureHess.row(2)[7] = t421;
//      pressureHess.row(2)[8] = t388;
//      pressureHess.row(3)[0] = t353;
//      pressureHess.row(3)[1] = t371;
//      pressureHess.row(3)[2] = t372;
//      pressureHess.row(3)[3] = forcePerUnitArea*((T_mul1_3*t167)/3.0-(T_mul1_2*t169)/3.0+(T_mul1_1*t171)/3.0);
//      pressureHess.row(3)[4] = t347;
//      pressureHess.row(3)[5] = t349;
//      pressureHess.row(3)[6] = t383;
//      pressureHess.row(3)[7] = t406;
//      pressureHess.row(3)[8] = t409;
//      pressureHess.row(4)[0] = t366;
//      pressureHess.row(4)[1] = t354;
//      pressureHess.row(4)[2] = t373;
//      pressureHess.row(4)[3] = t347;
//      pressureHess.row(4)[4] = forcePerUnitArea*((T_mul2_3*t173)/3.0-(T_mul2_2*t175)/3.0+(T_mul2_1*t177)/3.0);
//      pressureHess.row(4)[5] = t351;
//      pressureHess.row(4)[6] = t416;
//      pressureHess.row(4)[7] = t384;
//      pressureHess.row(4)[8] = t412;
//      pressureHess.row(5)[0] = t368;
//      pressureHess.row(5)[1] = t370;
//      pressureHess.row(5)[2] = t355;
//      pressureHess.row(5)[3] = t349;
//      pressureHess.row(5)[4] = t351;
//      pressureHess.row(5)[5] = forcePerUnitArea*((T_mul3_3*t179)/3.0-(T_mul3_2*t181)/3.0+(T_mul3_1*t183)/3.0);
//      pressureHess.row(5)[6] = t418;
//      pressureHess.row(5)[7] = t420;
//      pressureHess.row(5)[8] = t385;
//      pressureHess.row(6)[0] = t386;
//      pressureHess.row(6)[1] = t417;
//      pressureHess.row(6)[2] = t419;
//      pressureHess.row(6)[3] = t383;
//      pressureHess.row(6)[4] = t416;
//      pressureHess.row(6)[5] = t418;
//      pressureHess.row(6)[6] = -forcePerUnitArea*((T_mul1_1*t283)/3.0-(T_mul1_2*t282)/3.0+(T_mul1_3*t281)/3.0);
//      pressureHess.row(6)[7] = t428;
//      pressureHess.row(6)[8] = t429;
//      pressureHess.row(7)[0] = t413;
//      pressureHess.row(7)[1] = t387;
//      pressureHess.row(7)[2] = t421;
//      pressureHess.row(7)[3] = t406;
//      pressureHess.row(7)[4] = t384;
//      pressureHess.row(7)[5] = t420;
//      pressureHess.row(7)[6] = t428;
//      pressureHess.row(7)[7] = -forcePerUnitArea*((T_mul2_1*t286)/3.0-(T_mul2_2*t285)/3.0+(T_mul2_3*t284)/3.0);
//      pressureHess.row(7)[8] = t430;
//      pressureHess.row(8)[0] = t414;
//      pressureHess.row(8)[1] = t415;
//      pressureHess.row(8)[2] = t388;
//      pressureHess.row(8)[3] = t409;
//      pressureHess.row(8)[4] = t412;
//      pressureHess.row(8)[5] = t385;
//      pressureHess.row(8)[6] = t429;
//      pressureHess.row(8)[7] = t430;
//      pressureHess.row(8)[8] = -forcePerUnitArea*((T_mul3_1*t289)/3.0-(T_mul3_2*t288)/3.0+(T_mul3_3*t287)/3.0);
//
//      hess -= pressureHess;
//    }
    if (add_pressure) {
      Matrix<double, 9, 9> pressureHess;
      pressureHess.setZero();

      //matlab
      double x0_x = X(3 * idx(0) + 0);
      double x0_y = X(3 * idx(0) + 1);
      double x0_z = X(3 * idx(0) + 2);
      double x1_x = X(3 * idx(1) + 0);
      double x1_y = X(3 * idx(1) + 1);
      double x1_z = X(3 * idx(1) + 2);
      double x2_x = X(3 * idx(2) + 0);
      double x2_y = X(3 * idx(2) + 1);
      double x2_z = X(3 * idx(2) + 2);

      double t2 = (forcePerUnitArea * x0_x) / 2.0;
      double t3 = (forcePerUnitArea * x0_y) / 2.0;
      double t4 = (forcePerUnitArea * x1_x) / 2.0;
      double t5 = (forcePerUnitArea * x0_z) / 2.0;
      double t6 = (forcePerUnitArea * x1_y) / 2.0;
      double t7 = (forcePerUnitArea * x2_x) / 2.0;
      double t8 = (forcePerUnitArea * x1_z) / 2.0;
      double t9 = (forcePerUnitArea * x2_y) / 2.0;
      double t10 = (forcePerUnitArea * x2_z) / 2.0;
      double t11 = -t2;
      double t12 = -t3;
      double t13 = -t4;
      double t14 = -t5;
      double t15 = -t6;
      double t16 = -t7;
      double t17 = -t8;
      double t18 = -t9;
      double t19 = -t10;

      pressureHess.row(0)[4] = t10;
      pressureHess.row(0)[5] = t18;
      pressureHess.row(0)[7] = t17;
      pressureHess.row(0)[8] = t6;
      pressureHess.row(1)[3] = t19;
      pressureHess.row(1)[5] = t7;
      pressureHess.row(1)[6] = t8;
      pressureHess.row(1)[8] = t13;
      pressureHess.row(2)[3] = t9;
      pressureHess.row(2)[4] = t16;
      pressureHess.row(2)[6] = t15;
      pressureHess.row(2)[7] = t4;
      pressureHess.row(3)[1] = t19;
      pressureHess.row(3)[2] = t9;
      pressureHess.row(3)[7] = t5;
      pressureHess.row(3)[8] = t12;
      pressureHess.row(4)[0] = t10;
      pressureHess.row(4)[2] = t16;
      pressureHess.row(4)[6] = t14;
      pressureHess.row(4)[8] = t2;
      pressureHess.row(5)[0] = t18;
      pressureHess.row(5)[1] = t7;
      pressureHess.row(5)[6] = t3;
      pressureHess.row(5)[7] = t11;
      pressureHess.row(6)[1] = t8;
      pressureHess.row(6)[2] = t15;
      pressureHess.row(6)[4] = t14;
      pressureHess.row(6)[5] = t3;
      pressureHess.row(7)[0] = t17;
      pressureHess.row(7)[2] = t4;
      pressureHess.row(7)[3] = t5;
      pressureHess.row(7)[5] = t11;
      pressureHess.row(8)[0] = t6;
      pressureHess.row(8)[1] = t13;
      pressureHess.row(8)[3] = t12;
      pressureHess.row(8)[4] = t2;

      hess -= pressureHess;
//      std::cout << pressureHess << " pressure Hessian" << std::endl;
    }


    return hess;

  } // namespace fsim
}
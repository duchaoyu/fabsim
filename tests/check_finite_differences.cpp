#include "catch.hpp"
#include "fsim/HingeElement.h"
#include "fsim/util/first_fundamental_form.h"
#include "helpers.h"

#include <fsim/ElasticShell.h>
#include <fsim/ElasticRod.h>
#include <fsim/IncompressibleNeoHookeanElement.h>
#include <fsim/IncompressibleNeoHookeanMembrane.h>
#include <fsim/NeoHookeanElement.h>
#include <fsim/NeoHookeanMembrane.h>
#include <fsim/OrthotropicStVKMembrane.h>
#include <fsim/RodStencil.h>
#include <fsim/Spring.h>
#include <fsim/StVKMembrane.h>
#include <fsim/StVKElement.h>
#include <fsim/util/typedefs.h>

using namespace fsim;

// MEMBRANES

TEMPLATE_TEST_CASE("TriangleElement", "", NeoHookeanElement<>, IncompressibleNeoHookeanElement<>)
{
  using namespace Eigen;

  Mat3<double> V = GENERATE(take(2, matrix_random(3, 3)));
  double thickness = GENERATE(take(2, random(0., 1.)));
  double young_modulus = GENERATE(take(2, random(0., 1.)));
  double mass = GENERATE(take(2, random(0., 1.)));
  double poisson_ratio = GENERATE(take(2, random(0., 0.5)));
  TestType::nu = poisson_ratio;
  TestType::E = young_modulus;
  TestType::mass = mass;
  TestType e(V, Vector3i(0, 1, 2), thickness);

  SECTION("Gradient") { test_gradient(e, 1e-5); }
  SECTION("Hessian") { test_hessian(e, 1e-5); }
}

TEST_CASE("StVKElement", "[StVK]")
{
  using namespace Eigen;

  Mat2<double> V = GENERATE(take(2, matrix_random(3, 2)));
  double thickness = GENERATE(take(2, random(0., 1.)));
  double poisson_ratio = GENERATE(take(2, random(0., 0.5)));

  double E1 = GENERATE(take(2, random(0., 1.)));
  double E2 = GENERATE(take(2, random(0., 1.)));
  StVKElement<>::_C << E1, poisson_ratio * sqrt(E1 * E2), 0,
                                  poisson_ratio * sqrt(E1 * E2), E2, 0,
                                  0, 0, 0.5 * sqrt(E1 * E2) * (1 - poisson_ratio);
  StVKElement<>::_C /= (1 - std::pow(poisson_ratio, 2));
  StVKElement e(V, Vector3i(0, 1, 2), thickness);

  SECTION("Gradient") { test_gradient(e); }
  SECTION("Hessian") { test_hessian(e); }
}

// RODS

struct BundledRodStencil
{
  BundledRodStencil(const Eigen::Ref<const Mat3<double>> V,
                    const Eigen::Vector3d &t,
                    const Eigen::Vector3d &n,
                    const Eigen::Matrix<int, 5, 1> &dofs,
                    const Eigen::Vector2d &s, 
                    double m)
      : f1(t, n, t.cross(n)), f2(f1), stencil(V, f1, f2, dofs), stiffness(s), mass(m)
  {
    f2.update(V.row(dofs(1)), V.row(dofs(2)));
  }

  constexpr int nbDOFs() const { return stencil.nbDOFs(); }

  void prepare_data(const Eigen::Ref<const Eigen::VectorXd> X) const
  {
    f1.update(X.segment<3>(3 * stencil.idx(0)), X.segment<3>(3 * stencil.idx(1)));
    f2.update(X.segment<3>(3 * stencil.idx(1)), X.segment<3>(3 * stencil.idx(2)));

    stencil.updateReferenceTwist(f1, f2);
  }

  double energy(const Eigen::Ref<const Eigen::VectorXd> X) const
  {
    LocalFrame new_f1 =
        updateFrame(f1, X.segment<3>(3 * stencil.idx(0)), X.segment<3>(3 * stencil.idx(1)));
    LocalFrame new_f2 =
        updateFrame(f2, X.segment<3>(3 * stencil.idx(1)), X.segment<3>(3 * stencil.idx(2)));

    double ref_twist = stencil.getReferenceTwist();
    stencil.updateReferenceTwist(new_f1, new_f2);
    double res = stencil.energy(X, new_f1, new_f2, stiffness, mass);
    stencil.setReferenceTwist(ref_twist);
    return res;
  }
  RodStencil<>::LocalVector gradient(const Eigen::Ref<const Eigen::VectorXd> X) const
  {
    LocalFrame new_f1 =
        updateFrame(f1, X.segment<3>(3 * stencil.idx(0)), X.segment<3>(3 * stencil.idx(1)));
    LocalFrame new_f2 =
        updateFrame(f2, X.segment<3>(3 * stencil.idx(1)), X.segment<3>(3 * stencil.idx(2)));

    double ref_twist = stencil.getReferenceTwist();
    stencil.updateReferenceTwist(new_f1, new_f2);
    auto res = stencil.gradient(X, new_f1, new_f2, stiffness, mass);
    stencil.setReferenceTwist(ref_twist);
    return res;
  }
  RodStencil<>::LocalMatrix hessian(const Eigen::Ref<const Eigen::VectorXd> X) const
  {
    return stencil.hessian(X, f1, f2, stiffness, mass);
  }

  static const int NB_VERTICES = 3;
  static const int NB_DOFS = 11;

  mutable RodStencil<> stencil;
  mutable LocalFrame f1;
  mutable LocalFrame f2;
  double mass;
  Eigen::Vector2d stiffness;
};

TEST_CASE("RodStencil")
{
  using namespace Eigen;

  Mat3<double> V = GENERATE(take(2, matrix_random(3, 3)));
  Vector2d widths = GENERATE(take(2, vector_random(2, 0, 1)));
  double young_modulus = GENERATE(take(2, random(0., 1.)));
  double mass = GENERATE(take(2, random(0., 1.)));

  Vector3d n1 = GENERATE(take(2, vector_random(3))).normalized();
  Vector3d t1 = (V.row(1) - V.row(0)).normalized();
  n1 = (n1 - n1.dot(t1) * t1).normalized(); // make sure the frame is orthogonal

  VectorXi dofs(5);
  dofs << 0, 1, 2, 9, 10;

  Vector2d stiffnesses(pow(widths(0), 3) * widths(1), pow(widths(1), 3) * widths(0));
  stiffnesses *= young_modulus / 12;

  BundledRodStencil rod(V, t1, n1, dofs, stiffnesses, mass);

  SECTION("Gradient") { test_gradient(rod, 1e-5); }
  SECTION("Hessian") 
  { 
    VectorXd var(rod.nbDOFs());
    MatrixXd hessian_computed = MatrixXd::Zero(rod.nbDOFs(), rod.nbDOFs());
    MatrixXd hessian_numerical = MatrixXd::Zero(rod.nbDOFs(), rod.nbDOFs());

    for(int i = 0; i < 10; ++i)
    {
      var = VectorXd::NullaryExpr(rod.nbDOFs(), RandomRange(-1.0, 1.0));
      rod.prepare_data(var);

      hessian_computed += MatrixXd(rod.hessian(var));
      hessian_numerical += sym(MatrixXd(finite_differences([&rod](const Eigen::VectorXd &X) { return rod.gradient(X); }, var)));
    }

    hessian_computed /= 10;
    hessian_numerical /= 10;
    MatrixXd diff = hessian_computed - hessian_numerical;

    INFO("Numerical hessian\n" << hessian_numerical);
    INFO("Computed hessian\n" << hessian_computed);
    INFO("Difference\n" << diff);
    REQUIRE(diff.norm() / hessian_numerical.norm() == Approx(0.0).margin(1e-5));
  }
}

// SHELLS

bool check_normals(const Eigen::VectorXd &X)
{
  using namespace Eigen;
  Vector3d n0 = (X.segment<3>(0) - X.segment<3>(6)).cross(X.segment<3>(3) - X.segment<3>(6));
  Vector3d n1 = (X.segment<3>(3) - X.segment<3>(9)).cross(X.segment<3>(0) - X.segment<3>(9));

  return (n0.norm() > 1e-5 && n1.norm() > 1e-5 && (n0.normalized() + n1.normalized()).norm() > 1e-2);
};

TEMPLATE_TEST_CASE("HingeElement", "", SquaredAngleFormulation, TanAngleFormulation)
{
  using namespace Eigen;

  double coeff = GENERATE(take(2, random(0., 1.)));
  VectorXd X = GENERATE(take(2, filter(check_normals, vector_random(12))));
  Map<Mat3<double>> V(X.data(), 4, 3);
  HingeElement<TestType, true> e1(V, Vector4i(0, 1, 2, 3), coeff);

  SECTION("Gradient") { test_gradient(e1, 5e-5, check_normals); }
  SECTION("Hessian") { test_hessian(e1, 5e-5, check_normals); }
}

// SPRINGS

TEST_CASE("Spring")
{
  using namespace Eigen;

  double rest_length = GENERATE(take(2, random(0., 1.)));
  Spring s(0, 1, rest_length);

  SECTION("Gradient") 
  { 
    test_gradient(
      [&s](const auto &X) { return s.energy(X); },
      [&s](const auto &X) { 
        Vec<double, 6> res;
        res.segment<3>(0) = -s.force(X);
        res.segment<3>(3) = s.force(X);
        return res;
       }, 6); 
  }
  SECTION("Hessian") 
  { 
    test_hessian(
      [&s](const auto &X) { 
        Vec<double, 6> res;
        res.segment<3>(0) = -s.force(X);
        res.segment<3>(3) = s.force(X);
        return res;
       },
       [&s](const auto &X) { 
          MatrixXd res(6,6);
          Matrix3d h = s.hessian(X);
          res.block<3, 3>(0, 0) = -h;
          res.block<3, 3>(3, 3) = -h;
          res.block<3, 3>(0, 3) = h;
          res.block<3, 3>(3, 0) = h;
          return res;
       }, 6); 
  }
}

TEST_CASE("First Fundamental form")
{
  using namespace Eigen;

  Vector3i face(0, 1, 2);
  auto i = GENERATE(0, 1, 2, 3);

  SECTION("Gradient") 
  { 
    test_gradient(
      [&](const auto &X) { 
        return first_fundamental_form(Map<Mat3<double>>(const_cast<double*>(X.data()), 3, 3), face)(i); 
      }, [&](const auto &X) { 
        Matrix<double, 4, 9> deriv;
        first_fundamental_form(Map<Mat3<double>>(const_cast<double*>(X.data()), 3, 3), face, &deriv);
        return deriv.row(i);
      }, 9, 1e-6
    ); 
  }
  SECTION("Hessian") { 
    test_hessian(
      [&](const auto &X) { 
        Matrix<double, 4, 9> deriv;
        first_fundamental_form(Map<Mat3<double>>(const_cast<double*>(X.data()), 3, 3), face, &deriv);
        return deriv.row(i);
      }, [&](const auto &X) { 
        Matrix<double, 36, 9> dderiv;
        first_fundamental_form(Map<Mat3<double>>(const_cast<double*>(X.data()), 3, 3), face, nullptr, &dderiv);
        return dderiv.block<9,9>(9 * i, 0);
      }, 9, 1e-6
    ); 
  }
}

TEST_CASE("bendAngleGradient")
{
  using namespace Eigen;

  Vec<int, 4> indices;
  indices << 0, 1, 2, 3;

  SECTION("Gradient")
  { 
    test_gradient(
      [&](const Eigen::VectorXd &X) { 
        Vector3d n0 = (X.segment<3>(6) - X.segment<3>(0)).cross(X.segment<3>(6) - X.segment<3>(3));
        Vector3d n1 = (X.segment<3>(9) - X.segment<3>(3)).cross(X.segment<3>(9) - X.segment<3>(0));
        return signed_angle(n0, n1, X.segment<3>(3) - X.segment<3>(0));
      }, [&](const auto &X) { 
        return HingeElement<>::bendAngleGradient(Map<Mat3<double>>(const_cast<double*>(X.data()), 4, 3), indices);
      }, 12, 1e-6
    ); 
  }
  SECTION("Hessian") { 
    test_hessian(
      [&](const auto &X) { 
        return HingeElement<>::bendAngleGradient(Map<Mat3<double>>(const_cast<double*>(X.data()), 4, 3), indices);
      }, [&](const auto &X) { 
        Mat<double, 12, 12> dderiv;
        HingeElement<>::bendAngleGradient(Map<Mat3<double>>(const_cast<double*>(X.data()), 4, 3), indices, &dderiv);
        return dderiv;
      }, 12, 1e-6
    ); 
  }
}
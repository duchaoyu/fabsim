// ElasticRod.cpp
//
// Implementation of the Discrete Elastic Rods model of Bergou et al. as presented in
// "Discrete Viscous Threads" (https://doi.org/10.1145/1778765.1778853),
// see also  "A discrete, geometrically exact method for simulating nonlinear, elastic or
// non-elastic beams"  (https://hal.archives-ouvertes.fr/hal-02352879v1)
//
// Please note that, at the moment this implementation assumes a rectangular cross-section whose
// dimensions are given by the normal and binormal widths variables
//
// Author: David Jourdan (david.jourdan@inria.fr)
// Created: 01/09/20

#include "fsim/ElasticRod.h"

#include "fsim/RodStencil.h"
#include "fsim/util/vector_utils.h"

namespace fsim
{

template <CrossSection c>
ElasticRod<c>::ElasticRod(const Eigen::Ref<const Mat3<double>> V,
                       const Eigen::Ref<const Eigen::VectorXi> indices,
                       const Eigen::Vector3d &n,
                       double thickness,
                       double width,
                       double young_modulus,
                       double mass)
    : _stretch_modulus(1e3 * young_modulus * thickness * width)
{
  using namespace Eigen;

  RodStencil::mass = mass;

  nV = V.rows();
  int nR = indices.size();
  MatrixX2i extremal_edges(nR, 2); // rod extremal edge indices

  Mat3<double> D1, D2;
  Map<VectorXi> E(const_cast<int *>(indices.data()), indices.size());
  ElasticRod<>::bishopFrame(V, E, n, D1, D2);
  for(int j = 0; j < E.size() - 1; ++j)
  {
    _frames.emplace_back((V.row(E(j + 1)) - V.row(E(j))).normalized(), D1.row(j), D2.row(j));
    _springs.emplace_back(E(j), E(j + 1), (V.row(E(j)) - V.row(E(j + 1))).norm());
  }

  for(int j = 1; j < E.size() - 1; ++j)
  {
    Matrix<int, 5, 1> dofs;
    dofs << E(j - 1), E(j), E(j + 1), 3 * nV + j - 1, 3 * nV + j;
    _stencils.emplace_back(V, _frames[j - 1], _frames[j], dofs, Vector2d(thickness, width), young_modulus);
  }

  assert(_springs.size() == _frames.size());
}

template <CrossSection c>
ElasticRod<c>::ElasticRod(const Eigen::Ref<const Mat3<double>> V,
                       const Eigen::Vector3d &n,
                       double thickness,
                       double width,
                       double young_modulus,
                       double mass)
    : ElasticRod(V,
                 Eigen::VectorXi::LinSpaced(V.rows(), 0, V.rows() - 1),
                 n,
                 thickness,
                 width,
                 young_modulus,
                 mass)
{}

template <CrossSection c>
double ElasticRod<c>::energy(const Eigen::Ref<const Eigen::VectorXd> X) const
{
  using namespace Eigen;

  double result = 0;
  for(RodStencil &e: _stencils)
  {
    double old_twist = e.getReferenceTwist();
    // update properties
    LocalFrame f1 = getFrame(X, e.idx(0), e.idx(1), e.idx(3));
    f1.update(X.segment<3>(3 * e.idx(0)), X.segment<3>(3 * e.idx(1)));
    LocalFrame f2 = getFrame(X, e.idx(1), e.idx(2), e.idx(4));
    f2.update(X.segment<3>(3 * e.idx(1)), X.segment<3>(3 * e.idx(2)));

    e.updateReferenceTwist(f1, f2);

    result += e.energy(X, f1, f2);

    // restore previous state
    e.setReferenceTwist(old_twist);
  }

  for(const auto &s: _springs)
    result += _stretch_modulus * s.energy(X);

  return result;
}

template <CrossSection c>
void ElasticRod<c>::gradient(const Eigen::Ref<const Eigen::VectorXd> X, Eigen::Ref<Eigen::VectorXd> Y) const
{
  using namespace Eigen;

  for(RodStencil &e: _stencils)
  {
    LocalFrame f1 = getFrame(X, e.idx(0), e.idx(1), e.idx(3));
    LocalFrame f2 = getFrame(X, e.idx(1), e.idx(2), e.idx(4));
    auto grad = e.gradient(X, f1, f2);

    int n = e.nbVertices();
    for(int j = 0; j < n; ++j)
      Y.segment<3>(3 * e.idx(j)) += grad.segment<3>(3 * j);

    for(int j = 0; j < e.idx.size() - n; ++j)
      Y(e.idx(n + j)) += grad(3 * n + j);
  }

  for(const auto &s: _springs)
  {
    Vector3d force = _stretch_modulus * s.force(X);
    Y.segment<3>(3 * s.i) -= force;
    Y.segment<3>(3 * s.j) += force;
  }
}

template <CrossSection c>
Eigen::VectorXd ElasticRod<c>::gradient(const Eigen::Ref<const Eigen::VectorXd> X) const
{
  using namespace Eigen;

  VectorXd Y = VectorXd::Zero(X.size());
  gradient(X, Y);
  return Y;
}

template <CrossSection c>
Eigen::SparseMatrix<double> ElasticRod<c>::hessian(const Eigen::Ref<const Eigen::VectorXd> X) const
{
  using namespace Eigen;

  std::vector<Triplet<double>> triplets = hessianTriplets(X);

  SparseMatrix<double> hess(X.size(), X.size());
  hess.setFromTriplets(triplets.begin(), triplets.end());
  return hess;
}

template <CrossSection c>
std::vector<Eigen::Triplet<double>> ElasticRod<c>::hessianTriplets(const Eigen::Ref<const Eigen::VectorXd> X) const
{
  using namespace Eigen;

  std::vector<Triplet<double>> triplets;

  triplets.reserve(11 * 11 * _stencils.size() + 6 * 6 * _springs.size());
  for(auto &e: _stencils)
  {
    LocalFrame f1 = getFrame(X, e.idx(0), e.idx(1), e.idx(3));
    LocalFrame f2 = getFrame(X, e.idx(1), e.idx(2), e.idx(4));
    auto hess = e.hessian(X, f1, f2);

    int n = e.nbVertices();
    for(int j = 0; j < n; ++j)
      for(int k = 0; k < n; ++k)
        for(int l = 0; l < 3; ++l)
          for(int m = 0; m < 3; ++m)
            triplets.emplace_back(3 * e.idx(j) + l, 3 * e.idx(k) + m, hess(3 * j + l, 3 * k + m));

    for(int j = 0; j < e.idx.size() - n; ++j)
      for(int k = 0; k < n; ++k)
        for(int l = 0; l < 3; ++l)
        {
          triplets.emplace_back(e.idx(n + j), 3 * e.idx(k) + l, hess(3 * n + j, 3 * k + l));
          triplets.emplace_back(3 * e.idx(k) + l, e.idx(n + j), hess(3 * k + l, 3 * n + j));
        }

    for(int j = 0; j < e.idx.size() - n; ++j)
      for(int k = 0; k < e.idx.size() - n; ++k)
        triplets.emplace_back(e.idx(n + j), e.idx(n + k), hess(3 * n + j, 3 * n + k));
  }

  for(const auto &s: _springs)
  {
    Matrix3d h = _stretch_modulus * s.hessian(X);

    for(int k = 0; k < 3; ++k)
      for(int l = 0; l < 3; ++l)
      {
        triplets.emplace_back(3 * s.i + k, 3 * s.j + l, h(k, l));
        triplets.emplace_back(3 * s.j + k, 3 * s.i + l, h(k, l));
        triplets.emplace_back(3 * s.i + k, 3 * s.i + l, -h(k, l));
        triplets.emplace_back(3 * s.j + k, 3 * s.j + l, -h(k, l));
      }
  }

  return triplets;
}

template <CrossSection c>
void ElasticRod<c>::updateProperties(const Eigen::Ref<const Eigen::VectorXd> X)
{
  int k = 0;
  for(auto &spring: _springs)
  {
    _frames[k++].update(X.segment<3>(3 * spring.i), X.segment<3>(3 * spring.j));
  }

  for(auto &e: _stencils)
  {
    LocalFrame f1 = getFrame(X, e.idx(0), e.idx(1), e.idx(3));
    LocalFrame f2 = getFrame(X, e.idx(1), e.idx(2), e.idx(4));
    e.updateReferenceTwist(f1, f2);
  }
}

template <CrossSection c>
void ElasticRod<c>::getReferenceDirectors(Mat3<double> &D1, Mat3<double> &D2) const
{
  using namespace Eigen;
  D1.resize(_frames.size(), 3);
  D2.resize(_frames.size(), 3);
  for(int i = 0; i < _frames.size(); ++i)
  {
    D1.row(i) = _frames[i].d1;
    D2.row(i) = _frames[i].d2;
  }
}

template <CrossSection c>
void ElasticRod<c>::getRotatedDirectors(const Eigen::Ref<const Eigen::VectorXd> theta,
                                     Mat3<double> &P1,
                                     Mat3<double> &P2) const
{
  using namespace Eigen;
  P1.resize(_frames.size(), 3);
  P2.resize(_frames.size(), 3);
  for(int i = 0; i < _frames.size(); ++i)
  {
    P1.row(i) = cos(theta(i)) * _frames[i].d1 + sin(theta(i)) * _frames[i].d2;
    P2.row(i) = -sin(theta(i)) * _frames[i].d1 + cos(theta(i)) * _frames[i].d2;
  }
}

template <CrossSection c>
LocalFrame ElasticRod<c>::getFrame(const Eigen::Ref<const Eigen::VectorXd> X, int x0, int x1, int id) const
{
  LocalFrame f = _frames[id - 3 * nV];
  if(f.t.dot(X.segment<3>(3 * x1) - X.segment<3>(3 * x0)) < 0)
  {
    f.t *= -1;
    f.d2 *= -1;
  }
  return f;
}

template <CrossSection c>
void ElasticRod<c>::bishopFrame(const Eigen::Ref<const Mat3<double>> V,
                             const Eigen::Ref<const Eigen::VectorXi> E,
                             const Eigen::Vector3d &n,
                             Mat3<double> &P1,
                             Mat3<double> &P2)
{
  using namespace Eigen;
  int nE = E.rows() - 1; // number of edges

  P1.resize(nE, 3);
  P2.resize(nE, 3);

  // define first frame
  Vector3d t = (V.row(E(1)) - V.row(E(0))).normalized();
  Vector3d u;
  // if n isn't orthogonal to the first edge, the direction chosen to be the normal will be as close as possible to n
  if(abs(t.dot(n)) < 1e-6)
    u = n.normalized();
  else
    u = (n - n.dot(t) * t).normalized();
  Vector3d v = t.cross(u);
  P1.row(0) = u;
  P2.row(0) = v;

  for(int j = 1; j < nE; ++j)
  {
    u = parallel_transport(u, t, (V.row(E(j + 1)) - V.row(E(j))).normalized());
    t = (V.row(E(j + 1)) - V.row(E(j))).normalized();

    u = (u - u.dot(t) * t).normalized();
    P1.row(j) = u;
    P2.row(j) = t.cross(u);
  }
}

template <CrossSection c>
Mat3<double> ElasticRod<c>::curvatureBinormals(const Eigen::Ref<const Mat3<double>> P,
                                            const Eigen::Ref<const Eigen::VectorXi> E)
{
  using namespace Eigen;
  Mat3<double> KB(E.size() - 2, 3);

  for(int i = 1; i < E.size() - 1; ++i)
  {
    Vector3d t0 = (P.row(E(i)) - P.row(E(i - 1))).normalized();
    Vector3d t1 = (P.row(E(i + 1)) - P.row(E(i))).normalized();
    KB.row(i - 1) = 2 * t0.cross(t1) / (1 + t0.dot(t1));
  }
  return KB;
}

// instantiation
template class ElasticRod<CrossSection::Square>;
template class ElasticRod<CrossSection::Circle>;

} // namespace fsim

// ElasticMembraneModel.h
//
// Author: David Jourdan (david.jourdan@inria.fr)
// Created: 01/19/20

#pragma once

#include "ModelBase.h"
#include "TriangleElement.h"
#include "util/typedefs.h"

#include <iostream>

/**
 * template class for any FEM-based membrane material model (e.g. StVK, neohookean...)
 * @tparam Element  triangle stencil class such as StVKElement, NeoHookeanElement, etc.
 */
template <class Element>
class ElasticMembraneModel : public ModelBase<Element>
{
public:
  /**
   * constructor for ElasticMembraneModel
   * @param V  nV by 3 list of vertex positions
   * @param F  nF by 3 list of face indices
   * @param young_modulus  membrane's Young's modulus, controls the resistance to bending
   * @param poisson_ratio  membrane's Poisson's ratio
   * @param thickness  membrane's thickness
   * @param stress  membrane's amount of stretch in length (e.g. 1.3 if it has been stretched from 10cm long to 13cm
   * long)
   */
  ElasticMembraneModel(const Eigen::Ref<const fsim::Mat3<double>> V,
                       const Eigen::Ref<const fsim::Mat3<int>> F,
                       double young_modulus,
                       double poisson_ratio,
                       double thickness,
                       double stress = 1);

  ElasticMembraneModel(const Eigen::Ref<const fsim::Mat3<double>> V,
                       const Eigen::Ref<const fsim::Mat3<int>> F,
                       const std::vector<double> &young_moduli,
                       const std::vector<double> &thicknesses,
                       double poisson_ratio,
                       double stress = 1);

  // ElasticMembraneModel(const Eigen::Ref<const fsim::Mat3<double>> V,
  //                      const Eigen::Ref<const fsim::Mat3<int>> F,
  //                      const Eigen::Ref<const Eigen::VectorXd> stressX,
  //                      const Eigen::Ref<const Eigen::VectorXd> stressY,
  //                      const Eigen::Ref<const fsim::Mat2<double>> stressDirX,
  //                      const Eigen::Ref<const fsim::Mat2<double>> stressDirY,
  //                      const std::vector<double> &young_moduli,
  //                      const std::vector<double> &thicknesses,
  //                      double poisson_ratio);

  int nb_dofs() const { return 3 * nV; }
  void set_poisson_ratio(double poisson_ratio);
  void set_young_modulus(double E);
  void set_thickness(double t);
  void set_stress(double stress);

  double get_lame_alpha() { return Element::alpha; }
  double get_lame_beta() { return Element::beta; }
  double get_thickness() { return _thickness; }
  double get_stress() { return _stress; }

private:
  int nV, nF;
  double _stress;
  double _thickness;
  double _young_modulus;
};

// the ids are there to disambiguate between different instances so that they don't have the same Lamé coefficients
// and thicknesses (which are stored as static variables in each TriangleElement)
// If you only want to declare one Membrane instance (or if you're using several with the same Lamé coefficents),
// you can safely use the provided typedefs (e.g. StVKMembrane). However if you declared several instances with
// different Lamé coefficents, please declare them as e.g. StVKMembrane<0>, StVKMembrane<1>, etc.

using StVKMembrane = ElasticMembraneModel<StVKElement<0>>;
using NeoHookeanMembrane = ElasticMembraneModel<NeoHookeanElement<0>>;
using NHIncompressibleMembrane = ElasticMembraneModel<NHIncompressibleElement<0>>;

template <int id = 0>
using StVKMembraneModel = ElasticMembraneModel<StVKElement<id>>;
template <int id = 0>
using NeoHookeanMembraneModel = ElasticMembraneModel<NeoHookeanElement<id>>;
template <int id = 0>
using NHIncompressibleMembraneModel = ElasticMembraneModel<NHIncompressibleElement<id>>;

template <class Element>
ElasticMembraneModel<Element>::ElasticMembraneModel(const Eigen::Ref<const fsim::Mat3<double>> V,
                                                    const Eigen::Ref<const fsim::Mat3<int>> F,
                                                    double young_modulus,
                                                    double poisson_ratio,
                                                    double thickness,
                                                    double stress)
    : ElasticMembraneModel(V,
                           F,
                           std::vector<double>(F.rows(), young_modulus),
                           std::vector<double>(F.rows(), thickness),
                           poisson_ratio,
                           stress)
{}

template <class Element>
ElasticMembraneModel<Element>::ElasticMembraneModel(const Eigen::Ref<const fsim::Mat3<double>> V,
                                                    const Eigen::Ref<const fsim::Mat3<int>> F,
                                                    const std::vector<double> &young_moduli,
                                                    const std::vector<double> &thicknesses,
                                                    double poisson_ratio,
                                                    double stress)
    : _stress{stress}, _thickness{thicknesses[0]}, _young_modulus{young_moduli[0]}
{
  using namespace Eigen;

  nV = V.rows();
  nF = F.rows();

  if(Element::alpha != 0 || Element::beta != 0)
    std::cerr << "Warning: overwriting Lamé coefficents. Please declare your different instances as e.g. "
                 "StVKMembraneModel<0>, StVKMembraneModel<1>, etc.\n";
  Element::alpha = poisson_ratio / (1.0 - pow(poisson_ratio, 2));
  Element::beta = 0.5 / (1.0 + poisson_ratio);

  MatrixX3d smallerV = V / stress;

  this->_elements.reserve(nF);
  for(int i = 0; i < nF; ++i)
  {
    this->_elements.emplace_back(smallerV, F.row(i), thicknesses[i], young_moduli[i]);
  }
}

// template <class Element>
// ElasticMembraneModel<Element>::ElasticMembraneModel(const Eigen::Ref<const fsim::Mat3<double>> V,
//                                                     const Eigen::Ref<const fsim::Mat3<int>> F,
//                                                     const Eigen::Ref<const Eigen::VectorXd> stressX,
//                                                     const Eigen::Ref<const Eigen::VectorXd> stressY,
//                                                     const Eigen::Ref<const fsim::Mat2<double>> stressDirX,
//                                                     const Eigen::Ref<const fsim::Mat2<double>> stressDirY,
//                                                     const std::vector<double> &young_moduli,
//                                                     const std::vector<double> &thicknesses,
//                                                     double poisson_ratio)
//     : _stress{1}, _thickness{thicknesses[0]}, _young_modulus{young_moduli[0]}
// {
//   using namespace Eigen;

//   nV = V.rows();
//   nF = F.rows();

//   Element::alpha = poisson_ratio / (1.0 - pow(poisson_ratio, 2));
//   Element::beta = 0.5 / (1.0 + poisson_ratio);

//   this->_elements.reserve(nF);
//   for(int i = 0; i < nF; ++i)
//   {
//     this->_elements.emplace_back(V, F.row(i), thicknesses[i], young_moduli[i]);

//     Matrix2d E;
//     E.col(0) = V.row(F(i, 1)) - V.row(F(i, 0));
//     E.col(1) = V.row(F(i, 2)) - V.row(F(i, 0));
//     DiagonalMatrix<double, 2, 2> S(pow(stressX(i), 2), pow(stressY(i), 2));
//     Matrix2d R;
//     R.col(0) = stressDirX.row(i);
//     R.col(1) = stressDirY.row(i);

//     Matrix2d abar = E.transpose() * R * S * R.transpose() * E;
//     this->_elements[i].coeff = young_moduli[i] * thicknesses[i] / 2 * sqrt(abar.determinant());
//     this->_elements[i].abar_inv = abar.inverse();
//   }
// }

template <class Element>
void ElasticMembraneModel<Element>::set_poisson_ratio(double poisson_ratio)
{
  Element::alpha = poisson_ratio / (1 - pow(poisson_ratio, 2));
  Element::beta = 0.5 / (1 + poisson_ratio);
}

template <class Element>
void ElasticMembraneModel<Element>::set_young_modulus(double young_modulus)
{
  for(auto &e: this->_elements)
  {
    e.coeff *= young_modulus / _young_modulus;
  }
  _young_modulus = young_modulus;
}

template <class Element>
void ElasticMembraneModel<Element>::set_stress(double stress)
{
  assert(stress != 0);
  for(auto &e: this->_elements)
  {
    e.set_stretch_factor(stress, _stress);
  }
  _stress = stress;
}

template <class Element>
void ElasticMembraneModel<Element>::set_thickness(double t)
{
  assert(_thickness != 0);
  for(auto &elem: this->_elements)
  {
    elem.coeff *= t / _thickness;
  }
  _thickness = t;
}
// OrthotropicStVKElement.h
//
// StVK version of "Stable Orthotropic Materials" by Li and Barbič (https://doi.org/10.1109/tvcg.2015.2448105)
// Parameterizes the elasticity tensor with 2 Young's moduli and 1 Poisson's ratio
//
// Author: David Jourdan (david.jourdan@inria.fr)
// Created: 30/10/21

#include "ModelBase.h"
#include "OrthotropicStVKElement.h"

#include <array>

namespace fsim
{

template <int id = 0>
class OrthotropicStVKModel : public ModelBase<OrthotropicStVKElement<id>>
{
public:
  /**
   * constructor for OrthotropicStVKModel
   * @param V  nV by 3 list of vertex positions
   * @param F  nF by 3 list of face indices
   * @param thicknesses  membrane's thickness (per-triangle value)
   * @param E1  0 degree Young's modulus
   * @param E2  90 degrees Young's modulus
   * @param poisson_ratio  membrane's Poisson's ratio
   */
  OrthotropicStVKModel(const Eigen::Ref<const Mat2<double>> V,
                       const Eigen::Ref<const Mat3<int>> F,
                       double thickness,
                       double E1,
                       double E2,
                       double poisson_ratio);

  OrthotropicStVKModel(const Eigen::Ref<const Mat2<double>> V,
                       const Eigen::Ref<const Mat3<int>> F,
                       const std::vector<double> &thicknesses,
                       double E1,
                       double E2,
                       double poisson_ratio);

  int nb_dofs() const { return 3 * nV; }
  void set_poisson_ratio(double poisson_ratio);
  void set_young_moduli(double E1, double E2);
  void set_thickness(double t);

  double get_poisson_ratio() { return _poisson_ratio; }
  std::array<double, 2> get_young_moduli() { return {_E1, _E2}; }
  double get_thickness() { return _thickness; }

private:
  int nV;
  double _thickness = -1;
  double _poisson_ratio;
  double _E1;
  double _E2;
};

using OrthotropicStVKMembrane = OrthotropicStVKModel<0>;

template <int id>
OrthotropicStVKModel<id>::OrthotropicStVKModel(const Eigen::Ref<const Mat2<double>> V,
                                               const Eigen::Ref<const Mat3<int>> F,
                                               double thickness,
                                               double E1,
                                               double E2,
                                               double poisson_ratio)
    : OrthotropicStVKModel(V, F, std::vector<double>(F.rows(), thickness), E1, E2, poisson_ratio)
{
  _thickness = thickness;
}

template <int id>
OrthotropicStVKModel<id>::OrthotropicStVKModel(const Eigen::Ref<const Mat2<double>> V,
                                               const Eigen::Ref<const Mat3<int>> F,
                                               const std::vector<double> &thicknesses,
                                               double E1,
                                               double E2,
                                               double poisson_ratio)
  : _poisson_ratio{poisson_ratio}, _E1{E1}, _E2{E2}
{
  using namespace Eigen;

  nV = V.rows();

  if(OrthotropicStVKElement<id>::_C.norm() == 0)
    std::cerr << "Warning: overwriting elasticity tensor. Please declare your different instances as "
                 "OrthotropicStVKModel<0>, OrthotropicStVKModel<1>, etc.\n";

  OrthotropicStVKElement<id>::_C << 
    E1, poisson_ratio * sqrt(E1 * E2), 0, 
    poisson_ratio * sqrt(E1 * E2), E2, 0, 0, 
    0, 0.5 * sqrt(E1 * E2) * (1 - poisson_ratio);

  OrthotropicStVKElement<id>::_C /= (1 - std::pow(poisson_ratio, 2));

  int nF = F.rows();
  this->_elements.reserve(nF);
  for(int i = 0; i < nF; ++i)
    this->_elements.emplace_back(V, F.row(i), thicknesses[i]);
}

template <int id>
void OrthotropicStVKModel<id>::set_poisson_ratio(double poisson_ratio)
{
  _poisson_ratio = poisson_ratio;
  
  OrthotropicStVKElement<id>::_C << 
    _E1, _poisson_ratio * sqrt(_E1 * _E2), 0, 
    _poisson_ratio * sqrt(_E1 * _E2), _E2, 0, 0, 
    0, 0.5 * sqrt(_E1 * _E2) * (1 - _poisson_ratio);

  OrthotropicStVKElement<id>::_C /= (1 - std::pow(_poisson_ratio, 2));
}

template <int id>
void OrthotropicStVKModel<id>::set_young_moduli(double E1, double E2)
{
  _E1 = E1;
  _E2 = E2;

  OrthotropicStVKElement<id>::_C << 
    _E1, _poisson_ratio * sqrt(_E1 * _E2), 0, 
    _poisson_ratio * sqrt(_E1 * _E2), _E2, 0, 0, 
    0, 0.5 * sqrt(_E1 * _E2) * (1 - _poisson_ratio);

  OrthotropicStVKElement<id>::_C /= (1 - std::pow(_poisson_ratio, 2));
}

template <int id>
void OrthotropicStVKModel<id>::set_thickness(double t)
{
  if(_thickness <= 0)
    throw std::runtime_error(
        "Warning: membrane may have a locally varying thickness\nCan't set it to a constant value\n");
  for(auto &elem: this->_elements)
  {
    elem.coeff *= t / _thickness;
  }
  _thickness = t;
}

} // namespace fsim

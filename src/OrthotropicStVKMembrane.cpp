// OrthotropicStVKMembrane.cpp
//
// StVK version of "Stable Orthotropic Materials" by Li and Barbič (https://doi.org/10.1109/tvcg.2015.2448105)
// Parameterizes the elasticity tensor with 2 Young's moduli and 1 Poisson's ratio
//
// Author: David Jourdan (david.jourdan@inria.fr)
// Created: 10/30/21

#include "fsim/OrthotropicStVKMembrane.h"

#include <array>
#include <iostream>

namespace fsim
{

OrthotropicStVKMembrane::OrthotropicStVKMembrane(const Eigen::Ref<const Mat3<double>> V,
                                                     const Eigen::Ref<const Mat3<int>> F,
                                                     double thickness,
                                                     double E1,
                                                     double E2,
                                                     double poisson_ratio,
                                                     double mass,
                                                     double pressure)
    : OrthotropicStVKMembrane(V, F, std::vector<double>(F.rows(), thickness), E1, E2, poisson_ratio, mass, pressure)
{
  _thickness = thickness;
  // TODO: local frame, e1, e2
}

OrthotropicStVKMembrane::OrthotropicStVKMembrane(const Eigen::Ref<const Mat3<double>> V,
                                                     const Eigen::Ref<const Mat3<int>> F,
                                                     const std::vector<double> &thicknesses,
                                                     double E1,
                                                     double E2,
                                                     double poisson_ratio,
                                                     double mass,
                                                     double pressure)
    : _poisson_ratio{poisson_ratio}, _E1{E1}, _E2{E2}, _mass{mass}, _nu12{poisson_ratio}, _nu21{poisson_ratio*E2/E1}, _pressure{pressure}
{
  using namespace Eigen;

  nV = V.rows();

  _nu21 = _nu12 * _E2 / _E1;
  _G12 = _E1 * _E2 / (_E1 + _E2 + 2 * _E2 * _nu12); // https://help.solidworks.com/2011/english/solidworks/cworks/legacyhelp/simulation/materials/material_models/linear_elastic_orthotropic_model.html


  int nF = F.rows();
  this->_elements.reserve(nF);
  for(int i = 0; i < nF; ++i) {
    this->_elements.emplace_back(V, F.row(i), thicknesses[i]);
    if(pressure > 0){this->_elements[i].addPressure(pressure);}
  }
}

double OrthotropicStVKMembrane::energy(const Eigen::Ref<const Eigen::VectorXd> X) const
{
  Eigen::Matrix3d C;  // compliance matric / elastic material matrix

  C << _E1, _poisson_ratio * sqrt(_E1 * _E2), 0, 
       _poisson_ratio * sqrt(_E1 * _E2), _E2, 0, 
       0, 0, 0.5 * sqrt(_E1 * _E2) * (1 - _poisson_ratio);
  C /= (1 - std::pow(_poisson_ratio, 2));
  return ModelBase<OrthotropicStVKElement>::energy(X, C, _mass);
}

void OrthotropicStVKMembrane::gradient(const Eigen::Ref<const Eigen::VectorXd> X, Eigen::Ref<Eigen::VectorXd> Y) const
{
  Eigen::Matrix3d C;
  C << _E1, _poisson_ratio * sqrt(_E1 * _E2), 0, 
       _poisson_ratio * sqrt(_E1 * _E2), _E2, 0, 
       0, 0, 0.5 * sqrt(_E1 * _E2) * (1 - _poisson_ratio);
  C /= (1 - std::pow(_poisson_ratio, 2));
  ModelBase<OrthotropicStVKElement>::gradient(X, Y, C, _mass);
}

Eigen::VectorXd OrthotropicStVKMembrane::gradient(const Eigen::Ref<const Eigen::VectorXd> X) const
{
  using namespace Eigen;

  VectorXd Y = VectorXd::Zero(X.size());
  gradient(X, Y);
  return Y;
}

Eigen::SparseMatrix<double> OrthotropicStVKMembrane::hessian(const Eigen::Ref<const Eigen::VectorXd> X) const
{
  Eigen::Matrix3d C;
  C << _E1, _poisson_ratio * sqrt(_E1 * _E2), 0, 
       _poisson_ratio * sqrt(_E1 * _E2), _E2, 0, 
       0, 0, 0.5 * sqrt(_E1 * _E2) * (1 - _poisson_ratio);
  C /= (1 - std::pow(_poisson_ratio, 2));
  return ModelBase<OrthotropicStVKElement>::hessian(X, C, _mass);
}

std::vector<Eigen::Triplet<double>>
OrthotropicStVKMembrane::hessianTriplets(const Eigen::Ref<const Eigen::VectorXd> X) const
{
  Eigen::Matrix3d C;
  C << _E1, _poisson_ratio * sqrt(_E1 * _E2), 0, 
       _poisson_ratio * sqrt(_E1 * _E2), _E2, 0, 
       0, 0, 0.5 * sqrt(_E1 * _E2) * (1 - _poisson_ratio);
  C /= (1 - std::pow(_poisson_ratio, 2));
  return ModelBase<OrthotropicStVKElement>::hessianTriplets(X, C, _mass);
}


void OrthotropicStVKMembrane::setPoissonRatio(double poisson_ratio)
{
  _poisson_ratio = poisson_ratio;
  _nu12 = poisson_ratio;
  _nu21 = poisson_ratio *_E2/_E1;
}

void OrthotropicStVKMembrane::setPressure(double pressure)
    {
      _pressure = pressure;
      for(auto &elem: this->_elements){
        elem.addPressure(pressure);
      }
    }

void OrthotropicStVKMembrane:: setE1(double E1){
  _E1 = E1;
  _nu21 = _nu12 * _E2 / _E1;
}

void OrthotropicStVKMembrane:: setE2(double E2){
  _E2 = E2;
  _nu21 = _nu12 * _E2 / _E1;
}

void OrthotropicStVKMembrane::setThickness(double t)
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

void OrthotropicStVKMembrane::setMass(double mass)
{
  _mass = mass;
}

} // namespace fsim

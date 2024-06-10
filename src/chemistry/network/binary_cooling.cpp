//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file binary_cooling.cpp
//! \brief implementation of functions in class ChemNetwork
//! cooling function for binary disk simulation


// this class header
#include "binary_cooling.hpp"

// C headers

// C++ header
#include <iostream>   // endl
#include <limits>     // inf
#include <sstream>    // stringstream

// Athena++ header
#include "../../defs.hpp"
#include "../../eos/eos.hpp"
#include "../../hydro/hydro.hpp"
#include "../../mesh/mesh.hpp"
#include "../../parameter_input.hpp"
#include "../../scalars/scalars.hpp"
#include "../../units/units.hpp"
#include "network.hpp"

// species names
// This species is included to keep the code running but the reaction rate is zero.
const std::array<std::string, NSPECIES> ChemNetwork::species_names = {"spec0"};

//----------------------------------------------------------------------------------------
//! \brief ChemNetwork constructor

ChemNetwork::ChemNetwork(MeshBlock *pmb, ParameterInput *pin) {
  pmy_spec_ = pmb->pscalars;
  pmy_mb_ = pmb;
}

//----------------------------------------------------------------------------------------
//! \brief ChemNetwork destructor

ChemNetwork::~ChemNetwork() {
}

//----------------------------------------------------------------------------------------
//! \fn void ChemNetwork::InitializeNextStep(const int k, const int j, const int i)
//! \brief Set the rates of chemical reactions, eg. through density and radiation field.
//!
//! k, j, i are the corresponding index of the grid

void ChemNetwork::InitializeNextStep(const int k, const int j, const int i) {
  //TODO(Munan Gong): we need to change this likely
  Real rho, rho_floor;
  // density
  rho = pmy_mb_->phydro->w(IDN, k, j, i);
  // apply density floor
  rho_floor = pmy_mb_->peos->GetDensityFloor();
  rho = (rho > rho_floor) ?  rho : rho_floor;
  // 2D density -> surface density
  sigma_ =  rho;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ChemNetwork::RHS(const Real t, const Real *y, const Real ED,
//!                       Real *ydot)
//! \brief RHS: right-hand-side of ODE.
//!
//! dy/dt = ydot(t, y). Here y are the abundance
//! of species. details see CVODE package documentation.
//! all input/output variables are in code units

void ChemNetwork::RHS(const Real t, const Real *y, const Real ED, Real *ydot) {
  for (int i=0; i<NSPECIES; i++) {
    ydot[i] = 0.;
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn Real ChemNetwork::Edot(const Real t, const Real *y, const Real ED)
//! \brief energy equation dED/dt
//!
//! all input/output variables are in code units (ED is the energy density)

Real ChemNetwork::Edot(const Real t, const Real *y, const Real ED) {
  // isothermal
  if (!NON_BAROTROPIC_EOS) {
    return 0;
  }
  const Real T_floor = 1.; // temperature floor for cooling
  //TODO(Munan Gong): need to implement cooling
  // ernergy per hydrogen atom
  //const Real E_ergs = ED * pmy_mb_->pmy_mesh->punit->code_energydensity_cgs / nH_;
  //Real T = E_ergs / Thermo::CvCold(x_H2, x_He, x_e);
  //if (T < T_floor) {
  //  return 0;
  //}
  //Real dEdt = - Thermo::alpha_GD_ * nH_ * std::sqrt(T) * T;
  // return in code units
  //Real dEDdt = (dEdt * nH_ / pmy_mb_->pmy_mesh->punit->code_energydensity_cgs)
  //              * pmy_mb_->pmy_mesh->punit->code_time_cgs;
  //return dEDdt;
  return 0;
}

//----------------------------------------------------------------------------------------
//! \fn Real GetKappa(const Real temp)
//! \brief opacity from dust per gas surface density
//!
//! input:
//!   temp: dust temperature in Kelvin
//! output:
//!   kappa: opacity from dust per gas surface density, in cm2/g

Real ChemNetwork::GetKappa(const Real temp) {
  Real kappa = 0.;
  const Real temp_sub = 1500. //dust sublimation temperature in K
  if (temp > temp_sub) {
    return 0.;
  }
  if (temp < 150.) {
    kappa = 2.0e-4 * temp*temp;
  } else {
    kappa = 4.5;
  }
  return kappa;
}

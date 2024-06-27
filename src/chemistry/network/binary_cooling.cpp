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
#include "../utils/chemistry_utils.hpp"
#include "network.hpp"

// species names
// This species is included to keep the code running but the reaction rate is zero.
const std::array<std::string, NSPECIES> ChemNetwork::species_names = {"spec0"};

//----------------------------------------------------------------------------------------
//! \brief ChemNetwork constructor

ChemNetwork::ChemNetwork(MeshBlock *pmb, ParameterInput *pin) {
  pmy_spec_ = pmb->pscalars;
  pmy_mb_ = pmb;
  gm1_ = pin->GetReal("hydro", "gamma") - 1.;
  muH_ = pin->GetReal("problem", "muH");
  rstar_rsun_ =  pin->GetReal("problem", "rstar_rsun");
  rstar_ = rstar_rsun_ * Constants::solar_radius_cgs
            / pmy_mb_->pmy_mesh->punit->code_length_cgs;
  alpha_vis_ =  pin->GetReal("problem", "alpha_vis");
  f_lacc_ =  pin->GetReal("problem", "f_lacc");

  //calculate viscosity
  const Real omega_cgs = 1. / pmy_mb_->pmy_mesh->punit->code_time_cgs;

  //softening radius
  rsoft_cgs_ = 0.1 * pmy_mb_->pmy_mesh->punit->code_length_cgs;
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
  Real rho, rho_floor;
  // density
  rho = pmy_mb_->phydro->w(IDN, k, j, i);
  // apply density floor
  rho_floor = pmy_mb_->peos->GetDensityFloor();
  rho = (rho > rho_floor) ?  rho : rho_floor;
  // 2D density -> surface density
  sigma_ =  rho;
  sigma_cgs_ = sigma_ * pmy_mb_->pmy_mesh->punit->code_density_cgs
                * pmy_mb_->pmy_mesh->punit->code_length_cgs;

  // distance to binary stars
  const Real x1p = pmy_mb_->ruser_meshblock_data[1](0);
  const Real x2p = pmy_mb_->ruser_meshblock_data[1](1);
  const Real x3p = pmy_mb_->ruser_meshblock_data[1](2);
  const Real x1s = pmy_mb_->ruser_meshblock_data[2](0);
  const Real x2s = pmy_mb_->ruser_meshblock_data[2](1);
  const Real x3s = pmy_mb_->ruser_meshblock_data[2](2);
  const Real x1 = pmy_mb_->pcoord->x1v(i);
  const Real x2 = pmy_mb_->pcoord->x2v(j);
  const Real x3 = pmy_mb_->pcoord->x3v(k);
  rdiskp_cgs_ = sqrt(SQR(x1-x1p)+SQR(x2-x2p)+SQR(x3-x3p))
                 * pmy_mb_->pmy_mesh->punit->code_length_cgs;
  rdisks_cgs_ = sqrt(SQR(x1-x1s)+SQR(x2-x2s)+SQR(x3-x3s))
                 * pmy_mb_->pmy_mesh->punit->code_length_cgs;
  // mass accretion rates
  const Real mdotp = pmy_mb_->ruser_meshblock_data[0](2);
  const Real mdots = pmy_mb_->ruser_meshblock_data[0](3);
  mdotp_cgs_ = mdotp * pmy_mb_->pmy_mesh->punit->code_mass_cgs
                 / pmy_mb_->pmy_mesh->punit->code_time_cgs;
  mdots_cgs_ = mdots * pmy_mb_->pmy_mesh->punit->code_mass_cgs
                 / pmy_mb_->pmy_mesh->punit->code_time_cgs;
  // TODO(Munan Gong): for test in steady-state. Remove for time-dependent run
  // mdotp_cgs_ = 2.e20;
  // mdots_cgs_ = 1.e20;
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
  // sound speed cs^2 in cgs
  const Real cs_sq_cgs = ED/sigma_ * gm1_
                           * SQR(pmy_mb_->pmy_mesh->punit->code_velocity_cgs);
  const Real T = cs_sq_cgs *  muH_ * Constants::hydrogen_mass_cgs
                  / Constants::k_boltzmann_cgs;
  if (T < T_floor) {
    return 0;
  }
  const Real tau = sigma_cgs_ * GetKappa(T);
  // calculate accretion luminosity and irradiation flux
  const Real lum_acc_p = f_lacc_ * (1./rstar_)
                          * SQR(pmy_mb_->pmy_mesh->punit->code_velocity_cgs)
                          * mdotp_cgs_;
  const Real lum_acc_s = f_lacc_ * (1./rstar_)
                          * SQR(pmy_mb_->pmy_mesh->punit->code_velocity_cgs)
                          * mdots_cgs_;
  const Real f_firr = 0.127456; // factor 0.1/0.5**0.35
  const Real flux_irr_p = f_firr * lum_acc_p/( 4.*PI*SQR(rdiskp_cgs_+rsoft_cgs_) );
  const Real flux_irr_s = f_firr * lum_acc_s/( 4.*PI*SQR(rdisks_cgs_+rsoft_cgs_) );
  const Real flux_irr = flux_irr_p + flux_irr_s;
  // calculate dust cooling
  const Real flux_cool = ( Constants::stefan_boltzmann_cgs * SQR(T)*SQR(T)
                           - flux_irr ) / (tau + 1./tau);
  // return in code units
  const Real dEDdt = - flux_cool / ( pmy_mb_->pmy_mesh->punit->code_energydensity_cgs
                      * pmy_mb_->pmy_mesh->punit->code_length_cgs
                      / pmy_mb_->pmy_mesh->punit->code_time_cgs);
  return dEDdt;
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
  const Real temp_sub = 1500.; //dust sublimation temperature in K
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

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
  cfl_cool = pin->GetReal("chemistry", "cfl_cool");
  gm1_ = pin->GetReal("hydro", "gamma") - 1.;
  muH_ = pin->GetReal("problem", "muH");
  rstar_rsun_ =  pin->GetReal("problem", "rstar_rsun");
  rstar_ = rstar_rsun_ * Constants::solar_radius_cgs
            / pmy_mb_->pmy_mesh->punit->code_length_cgs;
  alpha_vis_ =  pin->GetReal("problem", "alpha_vis");
  f_lacc_ =  pin->GetReal("problem", "f_lacc");

  //calculate viscosity
  const Real omega_cgs = 1. / pmy_mb_->pmy_mesh->punit->code_time_cgs;
  nu_ = ChemistryUtility::GetViscosityBinaryDisk(alpha_vis_, omega_cgs);

  //TODO(Munan Gong): these need to be modified for binary disk run
  rdisk_ = 1.; //need to be distance to either binary star for calculating L_star(s)


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
  // sound speed cs^2 in cgs
  const Real cs_sq_cgs = ED/sigma_ * gm1_
                           * SQR(pmy_mb_->pmy_mesh->punit->code_velocity_cgs);
  const Real T = cs_sq_cgs *  muH_ * Constants::hydrogen_mass_cgs
                  / Constants::k_boltzmann_cgs;
  if (T < T_floor) {
    return 0;
  }
  const Real sigma_cgs = sigma_ * pmy_mb_->pmy_mesh->punit->code_density_cgs
                          * pmy_mb_->pmy_mesh->punit->code_length_cgs;
  const Real tau = sigma_cgs * GetKappa(T);
  //TODO(Munan Gong): need to read mdot from binary accretion
  mdot_cgs_ = 3. * PI * nu_ * sigma_cgs;
  // calculate accretion luminosity and irradiation flux
  // TODO(Munan Gong): need to account for luminosity from both stars
  const Real lum_acc = f_lacc_ * (1./rstar_)
                        * SQR(pmy_mb_->pmy_mesh->punit->code_velocity_cgs)
                        * mdot_cgs_;
  const Real f_firr = 0.127456; // factor 0.1/0.5**0.35
  const Real rdisk_cgs = rdisk_ * pmy_mb_->pmy_mesh->punit->code_length_cgs;
  const Real flux_irr = f_firr * lum_acc/( 4.*PI*SQR(rdisk_cgs) );
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

//----------------------------------------------------------------------------------------
//! \fn static Real CoolingTimeStep(MeshBlock *pmb)
//! \brief calculate cooling timestep for constraining timestep
//!   return cooling timestep in code units
Real ChemNetwork::CoolingTimeStep(MeshBlock *pmb) {
  Real real_max = std::numeric_limits<Real>::max();
  Real min_dt = real_max;
  return min_dt;
}

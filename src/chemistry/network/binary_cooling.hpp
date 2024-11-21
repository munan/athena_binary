#ifndef CHEMISTRY_NETWORK_BINARY_COOLING_HPP_
#define CHEMISTRY_NETWORK_BINARY_COOLING_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file binary_cooling.hpp
//! \brief definitions cooling function in binary disk simulations

// C headers

// C++ headers
#include <array>
#include <string>

// Athena++ classes headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "network.hpp"

class Units;

//! \class ChemNetwork
//! \brief binary cooling function
//!
//!  Note: species reaction rate set to zero.
class ChemNetwork : public NetworkWrapper {
  friend class MeshBlock;
 public:
  ChemNetwork(MeshBlock *pmb, ParameterInput *pin);
  ~ChemNetwork();

  // a list of species name, used in output
  static const std::array<std::string, NSPECIES> species_names;

  void InitializeNextStep(const int k, const int j, const int i);
  void RHS(const Real t, const Real *y, const Real ED, Real *ydot);
  Real Edot(const Real t, const Real *y, const Real ED);

 private:
  PassiveScalars *pmy_spec_;
  MeshBlock *pmy_mb_;

  Real GetKappa(const Real temp); //opacity from dust per gas surface density

  //input variables
  Real gm1_;  //adiabatic index (gamma - 1)
  Real muH_;  //mean molecular weight per hydrogen nuclei
  Real rstar_rsun_;     //stellar radius in solar radius
  Real rstar_;      //stellar radius in code units
  Real f_lacc_; // efficiency for accretion luminosity
  Real T_floor_; //temperature floor below which dEdt=0
  Real tau_floor_; //dust optical depth floor
  Real Tirr_ceiling_; //dust optical depth floor

  //variables updated at InitializeNextStep from hydro variable
  Real sigma_;     // surface density in code units
  Real sigma_cgs_; // surface density in cgs

  //variables for cooling
  Real rsoft_cgs_;      // softening radius for calcuating radiation
  Real rdiskp_cgs_;     // distance to primary star in code units
  Real rdisks_cgs_;     // distance to secondary star in code units
  Real mdotp_cgs_;  // mass accretion rate for primary in cgs
  Real mdots_cgs_;  // mass accretion rate for secondary in cgs
  Real mp_;         // mass of primary in code units
  Real ms_;         // mass of secondary in code units
};
#endif // CHEMISTRY_NETWORK_BINARY_COOLING_HPP_

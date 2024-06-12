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

  //calculate cooling timestep
  static Real CoolingTimeStep(MeshBlock *pmb);

  //cfl number for cooling time
  Real cfl_cool;

 private:
  PassiveScalars *pmy_spec_;
  MeshBlock *pmy_mb_;

  Real GetKappa(const Real temp); //opacity from dust per gas surface density

  //input variables
  Real gm1_;  //adiabatic index (gamma - 1)
  Real muH_;  //mean molecular weight per hydrogen nuclei
  Real rstar_rsun_;     //stellar radius in solar radius
  Real rstar_;      //stellar radius in code units
  Real alpha_vis_;      //alpha for calculating viscosity
  Real nu_;   //viscosity in cgs, calculated from alpha_vis
  Real f_lacc_; // efficiency for accretion luminosity

  //variables updated at InitializeNextStep from hydro variable
  Real sigma_;     // surface density in code units

  //variables for 0D cooling
  //TODO(Munan Gong): need to be modified for 2D disk simulation
  Real rdisk_;     // circumbinary disk radius in code units
  Real mdot_cgs_;  // mass accretion rate in cgs
};
#endif // CHEMISTRY_NETWORK_BINARY_COOLING_HPP_

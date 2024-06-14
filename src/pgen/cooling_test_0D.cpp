//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file cooling_test_0D.cpp
//! \brief problem generator, uniform 2D mesh with cooling, for test binary disk cooling
//======================================================================================

// C headers

// C++ headers
#include <algorithm>  // std::find()
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // std::runtime_error()
#include <string>     // c_str()
#include <vector>     // vector container

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../chem_rad/chem_rad.hpp"
#include "../chem_rad/integrators/rad_integrators.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"
#include "../units/units.hpp"

namespace {
Real HistoryT(MeshBlock *pmb, int iout); // average temperature output
Real muH; //mean molecular weight per hydrogen nuclei
} // namespace

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//========================================================================================
void Mesh::InitUserMeshData(ParameterInput *pin) {
  //mean molecular weight per hydrogen nuclei
  muH = pin->GetReal("problem", "muH");

  AllocateUserHistoryOutput(1);
  EnrollUserHistoryOutput(0, HistoryT, "T");
  return;
}

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief initialize problem with uniform 2D mesh and cooling parameters
//======================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // read initial gas surface density and sound speed
  const Real sigma0 = pin->GetReal("problem", "sigma0");
  const Real iso_cs = pin->GetReal("hydro", "iso_sound_speed");
  const Real pres = sigma0*SQR(iso_cs);
  const Real gm1  = peos->GetGamma() - 1.0;
  const Real r_init = 0;
  const Real vx = 0;

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        // density
        phydro->u(IDN, k, j, i) = sigma0;
        // velocity, x direction
        phydro->u(IM1, k, j, i) = sigma0*vx;
        // energy
        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN, k, j, i) = pres/gm1 + 0.5*sigma0*SQR(vx);
        }
      }
    }
  }

  // intialize chemical species
  if (NSPECIES > 0) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          for (int ispec=0; ispec < NSPECIES; ++ispec) {
            pscalars->s(ispec, k, j, i) = r_init*sigma0;
            if (CHEMISTRY_ENABLED) {
              Real s_ispec = pin->GetOrAddReal("problem",
                  "r_init_"+pscalars->chemnet.species_names[ispec], -1);
              if (s_ispec >= 0.) {
                pscalars->s(ispec, k, j, i) = s_ispec*sigma0;
              }
            }
          }
        }
      }
    }
  }
  return;
}

namespace {

Real HistoryT(MeshBlock *pmb, int iout) {
  int is = pmb->is, ie = pmb->ie, js = pmb->js, je = pmb->je, ks = pmb->ks, ke = pmb->ke;
  Real T = 0;
  AthenaArray<Real> volume; // 1D array of volumes
  // allocate 1D array for cell volume used in usr def history
  volume.NewAthenaArray(pmb->ncells1);
  // total volume
  Real vol_tot= (pmb->pmy_mesh->mesh_size.x1max - pmb->pmy_mesh->mesh_size.x1min)
    *(pmb->pmy_mesh->mesh_size.x2max - pmb->pmy_mesh->mesh_size.x2min)
    *(pmb->pmy_mesh->mesh_size.x3max - pmb->pmy_mesh->mesh_size.x3min);

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      pmb->pcoord->CellVolume(k,j,pmb->is,pmb->ie,volume);
      for (int i=is; i<=ie; i++) {
        T += volume(i)
             * (pmb->phydro->w(IPR,k,j,i)*pmb->pmy_mesh->punit->code_pressure_cgs)
             / (pmb->phydro->w(IDN,k,j,i)*pmb->pmy_mesh->punit->code_density_cgs)
             * muH * Constants::hydrogen_mass_cgs / Constants::k_boltzmann_cgs;
      }
    }
  }
  T /= vol_tot;
  return T;
}

} // namespace

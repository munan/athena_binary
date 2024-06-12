//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================
//! \file cooling_test_1D.cpp
//! \brief problem generator, 1D shock with cooling for binary disk test
//======================================================================================

// C headers

// C++ headers
#include <iostream>   // endl

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief initialize problem with 1D shock tube
//======================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // read initial gas surface density and velocity
  const Real sigma0 = pin->GetReal("problem", "sigma0");
  const Real vxp_kms = pin->GetReal("problem", "vx_kms_p"); // for x > 0
  const Real vxm_kms = pin->GetReal("problem", "vx_kms_m"); // for x < 0
  const Real iso_cs = pin->GetReal("hydro", "iso_sound_speed");
  const Real pres = sigma0*SQR(iso_cs);
  const Real gm1  = peos->GetGamma() - 1.0;
  const Real r_init = 0; // passive scalar

  // calculate initial velocity in code units
  const Real vxp = vxp_kms * 1e5 / pmy_mesh->punit->code_velocity_cgs;
  const Real vxm = vxm_kms * 1e5 / pmy_mesh->punit->code_velocity_cgs;

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        if (pcoord->x1v(i) < 0.) {
          // density
          phydro->u(IDN, k, j, i) = sigma0;
          // velocity, x direction
          phydro->u(IM1, k, j, i) = sigma0*vxm;
          // energy
          if (NON_BAROTROPIC_EOS) {
            phydro->u(IEN, k, j, i) = pres/gm1 + 0.5*sigma0*SQR(vxm);
          }
        } else {
          // density
          phydro->u(IDN, k, j, i) = sigma0;
          // velocity, x direction
          phydro->u(IM1, k, j, i) = sigma0*vxp;
          // energy
          if (NON_BAROTROPIC_EOS) {
            phydro->u(IEN, k, j, i) = pres/gm1 + 0.5*sigma0*SQR(vxp);
          }
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

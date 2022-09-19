/*!
 * \file trans_sources.hpp
 * \brief Numerics classes for integration of source terms in transition problems.
 * \author R. Roos
 * \version 7.4.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once
#include "../../../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../../scalar/scalar_sources.hpp"


/*!
 * \class CSourcePieceWise_TranEN
 * \brief Class for integrating the source terms of the e^N transition model equations.
 * \ingroup SourceDiscr
 * \author R. Roos
 */
template <class FlowIndices>
class CSourcePieceWise_TransEN final : public CNumerics {
 private:
  const FlowIndices idx; /*!< \brief Object to manage the access to the flow primitives. */

  su2double g_eff_i,
  g_eff_j,
  g_sep_i,
  g_sep_j;

  /*--- e^N Closure constants ---*/
  const su2double sigma_n = 1.0;

  su2double Vorticity;
  su2double Residual[2];
  su2double* Jacobian_i[2];
  su2double Jacobian_Buffer[4];// Static storage for the Jacobian (which needs to be pointer for return type).

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourcePieceWise_TransEN(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config)
      : CNumerics(val_nDim, 2, config),
        idx(val_nDim, config->GetnSpecies()) {
    
    /*--- "Allocate" the Jacobian using the static buffer. ---*/
    Jacobian_i[0] = Jacobian_Buffer;
    Jacobian_i[1] = Jacobian_Buffer + 2;
  }

  /*!
   * \brief Residual for source term integration.
   * \param[in] config - Definition of the particular problem.
   * \return A lightweight const-view (read-only) of the residual/flux and Jacobians.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override {
  /*--- ScalarVar[0] = k, ScalarVar[0] = w, TransVar[0] = gamma, and TransVar[0] = ReThetaT ---*/
  /*--- dU/dx = PrimVar_Grad[1][0] ---*/
  AD::StartPreacc();
  AD::SetPreaccIn(StrainMag_i);
  AD::SetPreaccIn(ScalarVar_i, nVar);
  AD::SetPreaccIn(ScalarVar_Grad_i, nVar, nDim);
  AD::SetPreaccIn(TransVar_i, nVar);
  AD::SetPreaccIn(TransVar_Grad_i, nVar, nDim);
  AD::SetPreaccIn(Volume); AD::SetPreaccIn(dist_i);
  AD::SetPreaccIn(&V_i[idx.Velocity()], nDim);
  AD::SetPreaccIn(PrimVar_Grad_i, nDim+idx.Velocity(), nDim);
  AD::SetPreaccIn(Vorticity_i, 3);

  su2double VorticityMag = sqrt(Vorticity_i[0]*Vorticity_i[0] +
                                Vorticity_i[1]*Vorticity_i[1] +
                                Vorticity_i[2]*Vorticity_i[2]);
  
  const su2double vel_u = V_i[idx.Velocity()];
  const su2double vel_v = V_i[1+idx.Velocity()];
  const su2double vel_w = (nDim ==3) ? V_i[2+idx.Velocity()] : 0.0;

  const su2double Velocity_Mag = sqrt(vel_u*vel_u + vel_v*vel_v + vel_w*vel_w);

  AD::SetPreaccIn(V_i[idx.Density()], V_i[idx.LaminarViscosity()], V_i[idx.EddyViscosity()]);

  Density_i = V_i[idx.Density()];
  Laminar_Viscosity_i = V_i[idx.LaminarViscosity()];
  Eddy_Viscosity_i = V_i[idx.EddyViscosity()];

  Residual[0] = 0.0;       Residual[1] = 0.0;
  Jacobian_i[0][0] = 0.0;  Jacobian_i[0][1] = 0.0;
  Jacobian_i[1][0] = 0.0;  Jacobian_i[1][1] = 0.0;
  
  if (dist_i > 1e-10) {

    /*-- Gradient of velocity magnitude ---*/

    su2double dU_dx = 0.5/Velocity_Mag*( 2.*vel_u*PrimVar_Grad_i[1][0]
                              +2.*vel_v*PrimVar_Grad_i[2][0]);
    if (nDim==3)
      dU_dx += 0.5/Velocity_Mag*( 2.*vel_w*PrimVar_Grad_i[3][0]);

    su2double dU_dy = 0.5/Velocity_Mag*( 2.*vel_u*PrimVar_Grad_i[1][1]
                              +2.*vel_v*PrimVar_Grad_i[2][1]);
    if (nDim==3)
      dU_dy += 0.5/Velocity_Mag*( 2.*vel_w*PrimVar_Grad_i[3][1]);

    su2double dU_dz = 0.0;
    if (nDim==3)
      dU_dz = 0.5/Velocity_Mag*( 2.*vel_u*PrimVar_Grad_i[1][2]
                                +2.*vel_v*PrimVar_Grad_i[2][2]
                                +2.*vel_w*PrimVar_Grad_i[3][2]);

    su2double du_ds = vel_u/Velocity_Mag*dU_dx + vel_v/Velocity_Mag*dU_dy;
    if (nDim==3)
      du_ds += vel_w/Velocity_Mag * dU_dz;

    /*--- Source ---*/
    Residual[0] += 1.0*Volume;
    Residual[1] += 1.0*Volume;

    /*--- Implicit part ---*/   
    Jacobian_i[0][0] = (1.0)*Volume;
    Jacobian_i[0][1] = 0.0;
    Jacobian_i[1][0] = 0.0;
    Jacobian_i[1][1] = 1.0*Volume;
  }  
  
  AD::SetPreaccOut(Residual, nVar);
  AD::EndPreacc();

  return ResidualType<>(Residual, Jacobian_i, nullptr);
  }
};


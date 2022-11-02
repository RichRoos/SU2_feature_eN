/*!
 * \file trans_diffusion.hpp
 * \brief Declarations of numerics classes for discretization of
 *        viscous fluxes in transition problems.
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


#include "../../scalar/scalar_diffusion.hpp"

/*!
 * \class CAvgGrad_TransEN
 * \brief Class for computing viscous term using average of gradient with correction (e^N transition model).
 * \ingroup ViscDiscr
 * \author R. Roos
 */
template <class FlowIndices>
class CAvgGrad_TransEN final : public CAvgGrad_Scalar<FlowIndices> {
private:
  using Base = CAvgGrad_Scalar<FlowIndices>;
  using Base::Laminar_Viscosity_i;
  using Base::Laminar_Viscosity_j;
  using Base::Eddy_Viscosity_i;
  using Base::Eddy_Viscosity_j;
  using Base::Density_i;
  using Base::Density_j;
  using Base::ScalarVar_i;
  using Base::ScalarVar_j;
  using Base::Proj_Mean_GradScalarVar;
  using Base::proj_vector_ij;
  using Base::implicit;
  using Base::Flux;
  using Base::Jacobian_i;
  using Base::Jacobian_j;

  const su2double sigma_n = 1.0;

  /*!
   * \brief Adds any extra variables to AD
   */
  void ExtraADPreaccIn() override {}

  /*!
   * \brief SA specific steps in the ComputeResidual method
   * \param[in] config - Definition of the particular problem.
   */
  void FinishResidualCalc(const CConfig* config) override {

	/*--- Compute mean effective dynamic viscosity ---*/
	const su2double diff_i_amplification = (Laminar_Viscosity_i + Eddy_Viscosity_i)/sigma_n;
	const su2double diff_j_amplification = (Laminar_Viscosity_j + Eddy_Viscosity_j)/sigma_n;

	const su2double diff_amplification = 0.5*(diff_i_amplification + diff_j_amplification);

    Flux[0] = diff_amplification*Proj_Mean_GradScalarVar[0];

    /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/

    if (implicit) {
      Jacobian_i[0][0] = (0.5*Proj_Mean_GradScalarVar[0]-diff_amplification*proj_vector_ij);
      Jacobian_j[0][0] = (0.5*Proj_Mean_GradScalarVar[0]+diff_amplification*proj_vector_ij);
    }

    //cout<<"EN Diffusion jacobian = "<<Jacobian_i[0][0]<<", "<<Jacobian_j[0][1]<<endl;
  }

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] correct_grad - Whether to correct gradient for skewness.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_TransEN(unsigned short val_nDim, unsigned short val_nVar,
                  bool correct_grad, const CConfig* config)
    : CAvgGrad_Scalar<FlowIndices>(val_nDim, val_nVar, correct_grad, config) {}
};

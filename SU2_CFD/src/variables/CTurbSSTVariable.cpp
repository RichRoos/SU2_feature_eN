/*!
 * \file CTurbSSTVariable.cpp
 * \brief Definition of the solution fields.
 * \author F. Palacios, A. Bueno
 * \version 7.0.3 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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


#include "../../include/variables/CTurbSSTVariable.hpp"


CTurbSSTVariable::CTurbSSTVariable(su2double kine, su2double omega, su2double mut, unsigned long npoint, unsigned long ndim, unsigned long nvar, const su2double* constants, CConfig *config)
  : CTurbVariable(npoint, ndim, nvar, config) {

  nPrimVar = nVar;

  Primitive.resize(nPoint, nPrimVar) = su2double(0.0);

  for(unsigned long iPoint=0; iPoint<nPoint; ++iPoint)
  {
    Primitive(iPoint,0) = kine;
    Primitive(iPoint,1) = omega;
    Solution(iPoint,0) = mut*omega;
    Solution(iPoint,1) = mut*omega*omega/kine;
  }

  Solution_Old = Solution;

  sigma_om2 = constants[3];
  beta_star = constants[6];

  F1.resize(nPoint) = su2double(1.0);
  F2.resize(nPoint) = su2double(0.0);
  CDkw.resize(nPoint) = su2double(0.0);

  muT.resize(nPoint) = mut;

  Max_Lambda_Inv.resize(nPoint) = su2double(0.0);
  Max_Lambda_Visc.resize(nPoint) = su2double(0.0);
  Delta_Time.resize(nPoint) = su2double(0.0);
}

void CTurbSSTVariable::SetCrossDiff(unsigned long iPoint, const su2double val_density) {
  AD::StartPreacc();
  AD::SetPreaccIn(val_density);
  AD::SetPreaccIn(Primitive(iPoint,1));
  AD::SetPreaccIn(Gradient[iPoint], nVar, nDim);

  /*--- Cross diffusion ---*/

  CDkw(iPoint) = 0.0;
  for (unsigned long iDim = 0; iDim < nDim; iDim++)
    CDkw(iPoint) += Gradient(iPoint,0,iDim)*Gradient(iPoint,1,iDim);
  CDkw(iPoint) *= 2.0*val_density*sigma_om2/Primitive(iPoint,1);
  // CDkw(iPoint) = max(CDkw(iPoint), CDKW_MIN);

  AD::SetPreaccOut(CDkw(iPoint));
  AD::EndPreacc();
}

void CTurbSSTVariable::SetBlendingFunc(unsigned long iPoint, const su2double val_viscosity,
                                       const su2double val_dist, const su2double val_density, const su2double cdkw_max) {
  su2double arg2, arg2A, arg2B, arg1;

  AD::StartPreacc();
  AD::SetPreaccIn(val_viscosity);  
  AD::SetPreaccIn(val_dist);
  AD::SetPreaccIn(val_density);
  AD::SetPreaccIn(Primitive[iPoint], nVar);
  AD::SetPreaccIn(cdkw_max);
  AD::SetPreaccIn(CDkw(iPoint));

  /*--- Note: Hellsten uses max CDkw in domain times 10^-8 as CDkw_min in arg1 for
              numerical stability. ---*/
  const su2double eps = numeric_limits<passivedouble>::epsilon();
  // const su2double cdkw_min = cdkw_max*1.0e-8;

  /*--- F1 ---*/

  arg2A = sqrt(Primitive(iPoint,0))/(beta_star*Primitive(iPoint,1)*val_dist+eps);
  arg2B = 500.0*val_viscosity / (val_density*val_dist*val_dist*Primitive(iPoint,1)+eps);
  arg2 = max(arg2A, arg2B);
  arg1 = min(arg2, 4.0*val_density*sigma_om2*Primitive(iPoint,0) / (max(CDkw(iPoint),CDKW_MIN)*val_dist*val_dist+eps));
  F1(iPoint) = tanh(pow(arg1, 4.0));

  /*--- F2 ---*/

  arg2 = max(2.0*arg2A, arg2B);
  F2(iPoint) = tanh(pow(arg2, 2.0));

  AD::SetPreaccOut(F1(iPoint)); AD::SetPreaccOut(F2(iPoint));
  AD::EndPreacc();

}

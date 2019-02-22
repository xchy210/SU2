/*!
 * \file sgs_model.cpp
 * \brief Subroutines of the <i>sgs_model.hpp</i> file.
 * \author E. van der Weide, T. Economon, P. Urbanczyk
 * \version 6.2.0 "Falcon"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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

#include "../include/sgs_model.hpp"

void CSGSModel::ComputeEddyViscosity_2D(const int       nEntities,
                                        const su2double *rho,
                                        const su2double *dudx,
                                        const su2double *dudy,
                                        const su2double *dvdx,
                                        const su2double *dvdy,
                                        const su2double lenScale,
                                        const su2double *distToWall,
                                        su2double       *muTurb) {
  for(int i=0; i<nEntities; ++i)
    muTurb[i] = 0.0;
}

void CSGSModel::ComputeEddyViscosity_3D(const int       nEntities,
                                        const su2double *rho,
                                        const su2double *dudx,
                                        const su2double *dudy,
                                        const su2double *dudz,
                                        const su2double *dvdx,
                                        const su2double *dvdy,
                                        const su2double *dvdz,
                                        const su2double *dwdx,
                                        const su2double *dwdy,
                                        const su2double *dwdz,
                                        const su2double lenScale,
                                        const su2double *distToWall,
                                        su2double       *muTurb) {
  for(int i=0; i<nEntities; ++i)
    muTurb[i] = 0.0;
}

void CSGSModel::ComputeGradEddyViscosity_2D(const int       nEntities,
                                            const su2double *rho,
                                            const su2double *drhodx,
                                            const su2double *drhody,
                                            const su2double *dudx,
                                            const su2double *dudy,
                                            const su2double *dvdx,
                                            const su2double *dvdy,
                                            const su2double *d2udx2,
                                            const su2double *d2udy2,
                                            const su2double *d2udxdy,
                                            const su2double *d2vdx2,
                                            const su2double *d2vdy2,
                                            const su2double *d2vdxdy,
                                            const su2double lenScale,
                                            const su2double *distToWall,
                                                  su2double *dMuTdx,
                                                  su2double *dMuTdy) {
  for(int i=0; i<nEntities; ++i)
    dMuTdx[i] = dMuTdy[i] = 0.0;
}

void CSGSModel::ComputeGradEddyViscosity_3D(const int       nEntities,
                                            const su2double *rho,
                                            const su2double *drhodx,
                                            const su2double *drhody,
                                            const su2double *drhodz,
                                            const su2double *dudx,
                                            const su2double *dudy,
                                            const su2double *dudz,
                                            const su2double *dvdx,
                                            const su2double *dvdy,
                                            const su2double *dvdz,
                                            const su2double *dwdx,
                                            const su2double *dwdy,
                                            const su2double *dwdz,
                                            const su2double *d2udx2,
                                            const su2double *d2udy2,
                                            const su2double *d2udz2,
                                            const su2double *d2udxdy,
                                            const su2double *d2udxdz,
                                            const su2double *d2udydz,
                                            const su2double *d2vdx2,
                                            const su2double *d2vdy2,
                                            const su2double *d2vdz2,
                                            const su2double *d2vdxdy,
                                            const su2double *d2vdxdz,
                                            const su2double *d2vdydz,
                                            const su2double *d2wdx2,
                                            const su2double *d2wdy2,
                                            const su2double *d2wdz2,
                                            const su2double *d2wdxdy,
                                            const su2double *d2wdxdz,
                                            const su2double *d2wdydz,
                                            const su2double lenScale,
                                            const su2double *distToWall,
                                                  su2double *dMuTdx,
                                                  su2double *dMuTdy,
                                                  su2double *dMuTdz) {
  for(int i=0; i<nEntities; ++i)
    dMuTdx[i] = dMuTdy[i] = dMuTdz[i] = 0.0;
}

void CSmagorinskyModel::ComputeEddyViscosity_2D(const int       nEntities,
                                                const su2double *rho,
                                                const su2double *dudx,
                                                const su2double *dudy,
                                                const su2double *dvdx,
                                                const su2double *dvdy,
                                                const su2double lenScale,
                                                const su2double *distToWall,
                                                su2double       *muTurb) {

  /* Constant coefficient Smagorinsky SGS is calculated:
   * ( C_s * L_c )^2 * |S(x,t)|
   * C_s = Smagorinsky constant
   * L_c = Filter width
   * S(x,t) = Rate of Strain Tensor ( 1/2 [ du_i/dx_j + du_j/dx_i] )
   */
  const su2double C_s_filter_width  = const_smag*filter_mult*lenScale;
  const su2double C_s_filter_width2 = C_s_filter_width*C_s_filter_width;

  /* Loop over the number of entities for which the eddy viscosity
     must be computed. */
  for(int i=0; i<nEntities; ++i) {

    const su2double S12          = 0.5*(dudy[i] + dvdx[i]);
    const su2double strain_rate2 = 2.0*(dudx[i]*dudx[i] + dvdy[i]*dvdy[i]
                                 + 2.0*S12*S12);

    muTurb[i] = rho[i]*C_s_filter_width2*sqrt(strain_rate2);
  }
}

void CSmagorinskyModel::ComputeEddyViscosity_3D(const int       nEntities,
                                                const su2double *rho,
                                                const su2double *dudx,
                                                const su2double *dudy,
                                                const su2double *dudz,
                                                const su2double *dvdx,
                                                const su2double *dvdy,
                                                const su2double *dvdz,
                                                const su2double *dwdx,
                                                const su2double *dwdy,
                                                const su2double *dwdz,
                                                const su2double lenScale,
                                                const su2double *distToWall,
                                                su2double       *muTurb) {

  /* Constant coefficient Smagorinsky SGS is calculated:
   * ( C_s * L_c )^2 * |S(x,t)|
   * C_s = Smagorinsky constant
   * L_c = Filter width
   * S(x,t) = Rate of Strain Tensor ( 1/2 [ du_i/dx_j + du_j/dx_i] )
   */
  const su2double C_s_filter_width = const_smag*filter_mult*lenScale;
  const su2double C_s_filter_width2 = C_s_filter_width*C_s_filter_width;

  /* Loop over the number of entities for which the eddy viscosity
     must be computed. */
  for(int i=0; i<nEntities; ++i) {

    const su2double S12 = 0.5*(dudy[i] + dvdx[i]);
    const su2double S13 = 0.5*(dudz[i] + dwdx[i]);
    const su2double S23 = 0.5*(dvdz[i] + dwdy[i]);

    const su2double strain_rate2 = 2.0*(dudx[i]*dudx[i] + dvdy[i]*dvdy[i]
                                 +      dwdz[i]*dwdz[i]
                                 +      2.0*(S12*S12 + S13*S13 + S23*S23));

    /* Return the SGS dynamic viscosity. */
    muTurb[i] = rho[i]*C_s_filter_width2*sqrt(strain_rate2);
  }
}

void CSmagorinskyModel::ComputeGradEddyViscosity_2D(const int       nEntities,
                                                    const su2double *rho,
                                                    const su2double *drhodx,
                                                    const su2double *drhody,
                                                    const su2double *dudx,
                                                    const su2double *dudy,
                                                    const su2double *dvdx,
                                                    const su2double *dvdy,
                                                    const su2double *d2udx2,
                                                    const su2double *d2udy2,
                                                    const su2double *d2udxdy,
                                                    const su2double *d2vdx2,
                                                    const su2double *d2vdy2,
                                                    const su2double *d2vdxdy,
                                                    const su2double lenScale,
                                                    const su2double *distToWall,
                                                          su2double *dMuTdx,
                                                          su2double *dMuTdy) {

  cout << "CSmagorinskyModel::ComputeGradEddyViscosity_2D: Not implemented yet" << endl;
  exit(1);
}

void CSmagorinskyModel::ComputeGradEddyViscosity_3D(const int       nEntities,
                                                    const su2double *rho,
                                                    const su2double *drhodx,
                                                    const su2double *drhody,
                                                    const su2double *drhodz,
                                                    const su2double *dudx,
                                                    const su2double *dudy,
                                                    const su2double *dudz,
                                                    const su2double *dvdx,
                                                    const su2double *dvdy,
                                                    const su2double *dvdz,
                                                    const su2double *dwdx,
                                                    const su2double *dwdy,
                                                    const su2double *dwdz,
                                                    const su2double *d2udx2,
                                                    const su2double *d2udy2,
                                                    const su2double *d2udz2,
                                                    const su2double *d2udxdy,
                                                    const su2double *d2udxdz,
                                                    const su2double *d2udydz,
                                                    const su2double *d2vdx2,
                                                    const su2double *d2vdy2,
                                                    const su2double *d2vdz2,
                                                    const su2double *d2vdxdy,
                                                    const su2double *d2vdxdz,
                                                    const su2double *d2vdydz,
                                                    const su2double *d2wdx2,
                                                    const su2double *d2wdy2,
                                                    const su2double *d2wdz2,
                                                    const su2double *d2wdxdy,
                                                    const su2double *d2wdxdz,
                                                    const su2double *d2wdydz,
                                                    const su2double lenScale,
                                                    const su2double *distToWall,
                                                          su2double *dMuTdx,
                                                          su2double *dMuTdy,
                                                          su2double *dMuTdz) {

  cout << "CSmagorinskyModel::ComputeGradEddyViscosity_3D: Not implemented yet" << endl;
  exit(1);
}

void CWALEModel::ComputeEddyViscosity_2D(const int       nEntities,
                                         const su2double *rho,
                                         const su2double *dudx,
                                         const su2double *dudy,
                                         const su2double *dvdx,
                                         const su2double *dvdy,
                                         const su2double lenScale,
                                         const su2double *distToWall,
                                         su2double       *muTurb) {

  /* Compute the length scale in WALE. */
  const su2double lenScaleWale  = const_WALE*lenScale;
  const su2double lenScaleWale2 = lenScaleWale*lenScaleWale;

  /* Loop over the number of entities for which the eddy viscosity
     must be computed. */
  for(int i=0; i<nEntities; ++i) {

    /* Compute the strain rate tensor, which is symmetric. */
    const su2double S11 = dudx[i], S22 = dvdy[i];
    const su2double S12 = 0.5*(dudy[i] + dvdx[i]);

    /* Compute the values of the Sd tensor. First without the trace
       correction of the diagonal terms. */
    su2double Sd11 = dudx[i]*dudx[i] + dudy[i]*dvdx[i];
    su2double Sd22 = dvdx[i]*dudy[i] + dvdy[i]*dvdy[i];

    const su2double Sd12 = 0.5*(dudx[i]*dudy[i] + dudy[i]*dvdy[i]
                         +      dvdx[i]*dudx[i] + dvdy[i]*dvdx[i]);

    /* Correct the diagonal elements, such that the trace of the Sd tensor is zero
       Note that this comes from the 3D formulation. */
    const su2double thirdTrace = (Sd11 + Sd22)/3.0;

    Sd11 -= thirdTrace;
    Sd22 -= thirdTrace;

    /* Compute the summation of both tensors. */
    const su2double sumS  = S11 *S11  + S22 *S22  + 2.0*S12 *S12;
    const su2double sumSd = Sd11*Sd11 + Sd22*Sd22 + 2.0*Sd12*Sd12;

    /* Compute the kinematic eddy viscosity. */
    const su2double sumSdPow3_2 = sumSd*sqrt(sumSd);
    const su2double sumSdPow5_4 = sqrt(sumSdPow3_2*sumSd);
    const su2double sumSPow5_2  = sumS*sumS*sqrt(sumS);
    const su2double denom       = sumSPow5_2 + sumSdPow5_4;

    const su2double nuEddy = lenScaleWale2*sumSdPow3_2/max(denom, 1.e-20);

    /* Compute the SGS dynamic viscosity. */
    muTurb[i] = rho[i]*nuEddy;
  }
}

void CWALEModel::ComputeEddyViscosity_3D(const int       nEntities,
                                         const su2double *rho,
                                         const su2double *dudx,
                                         const su2double *dudy,
                                         const su2double *dudz,
                                         const su2double *dvdx,
                                         const su2double *dvdy,
                                         const su2double *dvdz,
                                         const su2double *dwdx,
                                         const su2double *dwdy,
                                         const su2double *dwdz,
                                         const su2double lenScale,
                                         const su2double *distToWall,
                                         su2double       *muTurb) {

  /* Compute the length scale in WALE. */
  const su2double lenScaleWale  = const_WALE*lenScale;
  const su2double lenScaleWale2 = lenScaleWale*lenScaleWale;

  /* Loop over the number of entities for which the eddy viscosity
     must be computed. */
  for(int i=0; i<nEntities; ++i) {

    /* Compute the strain rate tensor, which is symmetric. */
    const su2double S11 = dudx[i], S22 = dvdy[i], S33 = dwdz[i];
    const su2double S12 = 0.5*(dudy[i] + dvdx[i]);
    const su2double S13 = 0.5*(dudz[i] + dwdx[i]);
    const su2double S23 = 0.5*(dvdz[i] + dwdy[i]);

    /* Compute the values of the Sd tensor. First without the trace
       correction of the diagonal terms. */
    su2double Sd11 = dudx[i]*dudx[i] + dudy[i]*dvdx[i] + dudz[i]*dwdx[i];
    su2double Sd22 = dvdx[i]*dudy[i] + dvdy[i]*dvdy[i] + dvdz[i]*dwdy[i];
    su2double Sd33 = dwdx[i]*dudz[i] + dwdy[i]*dvdz[i] + dwdz[i]*dwdz[i];

    const su2double Sd12 = 0.5*(dudx[i]*dudy[i] + dudy[i]*dvdy[i] + dudz[i]*dwdy[i]
                         +      dvdx[i]*dudx[i] + dvdy[i]*dvdx[i] + dvdz[i]*dwdx[i]);
    const su2double Sd13 = 0.5*(dudx[i]*dudz[i] + dudy[i]*dvdz[i] + dudz[i]*dwdz[i]
                         +      dwdx[i]*dudx[i] + dwdy[i]*dvdx[i] + dwdz[i]*dwdx[i]);
    const su2double Sd23 = 0.5*(dvdx[i]*dudz[i] + dvdy[i]*dvdz[i] + dvdz[i]*dwdz[i]
                         +      dwdx[i]*dudy[i] + dwdy[i]*dvdy[i] + dwdz[i]*dwdy[i]);

    /* Correct the diagonal elements, such that the trace of the Sd tensor is zero. */
    const su2double thirdTrace = (Sd11 + Sd22 + Sd33)/3.0;

    Sd11 -= thirdTrace;
    Sd22 -= thirdTrace;
    Sd33 -= thirdTrace;

    /* Compute the summation of both tensors. */
    const su2double sumS  = S11*S11 + S22*S22 + S33*S33
                          + 2.0*(S12*S12 + S13*S13 + S23*S23);
    const su2double sumSd = Sd11*Sd11 + Sd22*Sd22 + Sd33*Sd33
                          + 2.0*(Sd12*Sd12 + Sd13*Sd13 + Sd23*Sd23);

    /* Compute the kinematic eddy viscosity. */
    const su2double sumSdPow3_2 = sumSd*sqrt(sumSd);
    const su2double sumSdPow5_4 = sqrt(sumSdPow3_2*sumSd);
    const su2double sumSPow5_2  = sumS*sumS*sqrt(sumS);
    const su2double denom       = sumSPow5_2 + sumSdPow5_4;

    const su2double nuEddy = lenScaleWale2*sumSdPow3_2/max(denom, 1.e-20);

    /* Compute the SGS dynamic viscosity. */
    muTurb[i] = rho[i]*nuEddy;
  }
}

void CWALEModel::ComputeGradEddyViscosity_2D(const int       nEntities,
                                             const su2double *rho,
                                             const su2double *drhodx,
                                             const su2double *drhody,
                                             const su2double *dudx,
                                             const su2double *dudy,
                                             const su2double *dvdx,
                                             const su2double *dvdy,
                                             const su2double *d2udx2,
                                             const su2double *d2udy2,
                                             const su2double *d2udxdy,
                                             const su2double *d2vdx2,
                                             const su2double *d2vdy2,
                                             const su2double *d2vdxdy,
                                             const su2double lenScale,
                                             const su2double *distToWall,
                                                   su2double *dMuTdx,
                                                   su2double *dMuTdy) {


  cout << "CWALEModel::ComputeGradEddyViscosity_2D: Not implemented yet" << endl;
  exit(1);
}

void CWALEModel::ComputeGradEddyViscosity_3D(const int       nEntities,
                                             const su2double *rho,
                                             const su2double *drhodx,
                                             const su2double *drhody,
                                             const su2double *drhodz,
                                             const su2double *dudx,
                                             const su2double *dudy,
                                             const su2double *dudz,
                                             const su2double *dvdx,
                                             const su2double *dvdy,
                                             const su2double *dvdz,
                                             const su2double *dwdx,
                                             const su2double *dwdy,
                                             const su2double *dwdz,
                                             const su2double *d2udx2,
                                             const su2double *d2udy2,
                                             const su2double *d2udz2,
                                             const su2double *d2udxdy,
                                             const su2double *d2udxdz,
                                             const su2double *d2udydz,
                                             const su2double *d2vdx2,
                                             const su2double *d2vdy2,
                                             const su2double *d2vdz2,
                                             const su2double *d2vdxdy,
                                             const su2double *d2vdxdz,
                                             const su2double *d2vdydz,
                                             const su2double *d2wdx2,
                                             const su2double *d2wdy2,
                                             const su2double *d2wdz2,
                                             const su2double *d2wdxdy,
                                             const su2double *d2wdxdz,
                                             const su2double *d2wdydz,
                                             const su2double lenScale,
                                             const su2double *distToWall,
                                                   su2double *dMuTdx,
                                                   su2double *dMuTdy,
                                                   su2double *dMuTdz) {

  cout << "CWALEModel::ComputeGradEddyViscosity_3D: Not implemented yet" << endl;
  exit(1);
}

void CVremanModel::ComputeEddyViscosity_2D(const int       nEntities,
                                           const su2double *rho,
                                           const su2double *dudx,
                                           const su2double *dudy,
                                           const su2double *dvdx,
                                           const su2double *dvdy,
                                           const su2double lenScale,
                                           const su2double *distToWall,
                                           su2double       *muTurb) {

  /* Compute the square of the length scale. */
  const su2double lenScale2 = lenScale*lenScale;

  /* Loop over the number of entities for which the eddy viscosity
     must be computed. */
  for(int i=0; i<nEntities; ++i) {

    su2double alpha11 = dudx[i];
    su2double alpha22 = dvdy[i];

    //Check if it is necessary to remove the trace.
    const su2double tmp = (alpha11 + alpha22)/3.0;
    alpha11 -= tmp;
    alpha22 -= tmp;

    const su2double alpha12 = dudy[i];
    const su2double alpha21 = dvdx[i];

    const su2double beta11  = lenScale2*(alpha11*alpha11 + alpha12*alpha12);
    const su2double beta12  = lenScale2*(alpha11*alpha21 + alpha12*alpha22);
    const su2double beta22  = lenScale2*(alpha21*alpha21 + alpha22*alpha22);

    su2double B = beta11*beta22 - beta12*beta12;
              B = (B + fabs(B))*0.5;
    const su2double denon = alpha11*alpha11 + alpha22*alpha22
                          + alpha12*alpha12 + alpha21*alpha21;

    const su2double nuEddy_Vreman = sqrt(B/(denon+1.0E-20));

    /* Compute the SGS dynamic viscosity. */
    muTurb[i] = rho[i]*const_Vreman*nuEddy_Vreman;
  }
}

void CVremanModel::ComputeEddyViscosity_3D(const int       nEntities,
                                           const su2double *rho,
                                           const su2double *dudx,
                                           const su2double *dudy,
                                           const su2double *dudz,
                                           const su2double *dvdx,
                                           const su2double *dvdy,
                                           const su2double *dvdz,
                                           const su2double *dwdx,
                                           const su2double *dwdy,
                                           const su2double *dwdz,
                                           const su2double lenScale,
                                           const su2double *distToWall,
                                           su2double       *muTurb) {
  
  /* Compute the square of the length scale. */
  const su2double lenScale2 = lenScale*lenScale;

  /* Loop over the number of entities for which the eddy viscosity
     must be computed. */
  for(int i=0; i<nEntities; ++i) {

    su2double alpha11 = dudx[i];
    su2double alpha22 = dvdy[i];
    su2double alpha33 = dwdz[i];

    //Check if it is necessary to remove the trace.
    const su2double tmp = (alpha11 + alpha22 + alpha33)/3.0;
    alpha11 -= tmp;
    alpha22 -= tmp;
    alpha33 -= tmp;

    const su2double alpha12 = dudy[i];
    const su2double alpha13 = dudz[i];
    const su2double alpha23 = dvdz[i];

    const su2double alpha21 = dvdx[i];
    const su2double alpha31 = dwdx[i];
    const su2double alpha32 = dwdy[i];

    const su2double beta11 = lenScale2*(alpha11*alpha11 + alpha12*alpha12 + alpha13*alpha13);
    const su2double beta12 = lenScale2*(alpha11*alpha21 + alpha12*alpha22 + alpha13*alpha23);
    const su2double beta13 = lenScale2*(alpha11*alpha31 + alpha12*alpha32 + alpha13*alpha33);
    const su2double beta22 = lenScale2*(alpha21*alpha21 + alpha22*alpha22 + alpha23*alpha23);
    const su2double beta23 = lenScale2*(alpha21*alpha31 + alpha22*alpha32 + alpha23*alpha33);
    const su2double beta33 = lenScale2*(alpha31*alpha31 + alpha32*alpha32 + alpha33*alpha33);

    su2double B = beta11*beta22-beta12*beta12+beta11*beta33-beta13*beta13+beta22*beta33-beta23*beta23;
              B = (B + fabs(B))*0.5;
    const su2double denon = alpha11*alpha11+alpha22*alpha22+alpha33*alpha33 +
                            alpha12*alpha12+alpha13*alpha13+alpha23*alpha23 +
                            alpha21*alpha21+alpha31*alpha31+alpha32*alpha32;

    const su2double nuEddy_Vreman = sqrt(B/(denon+1.0E-20));

    /* Compute the SGS dynamic viscosity. */
    muTurb[i] = rho[i]*const_Vreman*nuEddy_Vreman;
  }
}

void CVremanModel::ComputeGradEddyViscosity_2D(const int       nEntities,
                                               const su2double *rho,
                                               const su2double *drhodx,
                                               const su2double *drhody,
                                               const su2double *dudx,
                                               const su2double *dudy,
                                               const su2double *dvdx,
                                               const su2double *dvdy,
                                               const su2double *d2udx2,
                                               const su2double *d2udy2,
                                               const su2double *d2udxdy,
                                               const su2double *d2vdx2,
                                               const su2double *d2vdy2,
                                               const su2double *d2vdxdy,
                                               const su2double lenScale,
                                               const su2double *distToWall,
                                                     su2double *dMuTdx,
                                                     su2double *dMuTdy) {

  cout << "CWALEModel::ComputeGradEddyViscosity_2D: Not implemented yet" << endl;
  exit(1);
}

void CVremanModel::ComputeGradEddyViscosity_3D(const int       nEntities,
                                               const su2double *rho,
                                               const su2double *drhodx,
                                               const su2double *drhody,
                                               const su2double *drhodz,
                                               const su2double *dudx,
                                               const su2double *dudy,
                                               const su2double *dudz,
                                               const su2double *dvdx,
                                               const su2double *dvdy,
                                               const su2double *dvdz,
                                               const su2double *dwdx,
                                               const su2double *dwdy,
                                               const su2double *dwdz,
                                               const su2double *d2udx2,
                                               const su2double *d2udy2,
                                               const su2double *d2udz2,
                                               const su2double *d2udxdy,
                                               const su2double *d2udxdz,
                                               const su2double *d2udydz,
                                               const su2double *d2vdx2,
                                               const su2double *d2vdy2,
                                               const su2double *d2vdz2,
                                               const su2double *d2vdxdy,
                                               const su2double *d2vdxdz,
                                               const su2double *d2vdydz,
                                               const su2double *d2wdx2,
                                               const su2double *d2wdy2,
                                               const su2double *d2wdz2,
                                               const su2double *d2wdxdy,
                                               const su2double *d2wdxdz,
                                               const su2double *d2wdydz,
                                               const su2double lenScale,
                                               const su2double *distToWall,
                                                     su2double *dMuTdx,
                                                     su2double *dMuTdy,
                                                     su2double *dMuTdz) {

  cout << "CWALEModel::ComputeGradEddyViscosity_3D: Not implemented yet" << endl;
  exit(1);
}

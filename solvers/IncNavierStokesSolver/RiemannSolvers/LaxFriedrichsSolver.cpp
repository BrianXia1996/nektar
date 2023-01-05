///////////////////////////////////////////////////////////////////////////////
//
// File: LaxFriedrichsSolver.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: LaxFriedrichs Riemann solver for incNS
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/RiemannSolvers/LaxFriedrichsSolver.h>

namespace Nektar
{
std::string LaxFriedrichsSolver::solverName =
    SolverUtils::GetRiemannSolverFactory().RegisterCreatorFunction(
        "LaxFriedrichs", LaxFriedrichsSolver::create,
        "Lax-Friedrichs Riemann solver");

LaxFriedrichsSolver::LaxFriedrichsSolver(
    const LibUtilities::SessionReaderSharedPtr &pSession)
    : IncompressibleSolver(pSession)
{
}

/**
 * @brief Lax-Friedrichs Riemann solver
 *
 * @param uL     x-velocity component left state.
 * @param uR     x-velocity component right state.
 * @param vL     y-velocity component left state.
 * @param vR     y-velocity component right state.
 * @param wL     z-velocity component left state.
 * @param wR     z-velocity component right state.
 * @param EL        Energy left state.
 * @param ER        Energy right state.
 * @param uf     Computed Riemann flux for x-velocity component
 * @param vf     Computed Riemann flux for y-velocity component
 * @param wf     Computed Riemann flux for z-velocity component
 * @param Ef     Computed Riemann flux for energy.
 */
void LaxFriedrichsSolver::v_PointSolve( NekDouble uL, NekDouble vL, NekDouble wL, NekDouble EL,
                                        NekDouble uR, NekDouble vR, NekDouble wR, NekDouble ER,
                                        NekDouble nx, NekDouble ny, NekDouble nz,
                                        NekDouble &uf, NekDouble &vf, NekDouble &wf, NekDouble &Ef)
{
    boost::ignore_unused(EL, ER, Ef);   // currently energy eqn not supportred     
    // maximum eigenvalues
     auto LambdaL = 2.0*fabs(nx*uL+ny*vL+nz*wL);
     auto LambdaR = 2.0*fabs(nx*uR+ny*vR+nz*wR);
     if (LambdaL<LambdaR){
        LambdaL = LambdaR;
     }
     // lax flux
     uf = 0.5*(uL*uL*nx + uL*vL*ny + uL*wL*nz 
             + uR*uR*nx + uR*vR*ny + uR*wR*nz 
             + LambdaL*(uL-uR));
     vf = 0.5*(vL*uL*nx + vL*vL*ny + vL*wL*nz 
             + vR*uR*nx + vR*vR*ny + vR*wR*nz 
             + LambdaL*(vL-vR));
     wf = 0.5*(wL*uL*nx + wL*vL*ny + wL*wL*nz 
             + wR*uR*nx + wR*vR*ny + wR*wR*nz 
             + LambdaL*(wL-wR));

}
} // namespace Nektar

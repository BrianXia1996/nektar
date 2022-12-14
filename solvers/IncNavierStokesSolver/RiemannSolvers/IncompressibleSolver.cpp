///////////////////////////////////////////////////////////////////////////////
//
// File: CompressibleSolver.cpp
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
// Description: Incompressible Riemann solver.
//
///////////////////////////////////////////////////////////////////////////////

#include "IncompressibleSolver.h"
#include <boost/algorithm/string/predicate.hpp>
#include <boost/core/ignore_unused.hpp>

namespace Nektar
{
IncompressibleSolver::IncompressibleSolver(
    const LibUtilities::SessionReaderSharedPtr &pSession)
    : RiemannSolver(pSession), m_pointSolve(true)
{
    m_requiresRotation = false;

    // Create equation of state object
    // std::string eosType;
    // pSession->LoadSolverInfo("EquationOfState", eosType, "IdealGas");
    // m_eos = GetEquationOfStateFactory().CreateInstance(eosType, pSession);
    // Check if using ideal gas
    // m_idealGas = boost::iequals(eosType, "IdealGas");
}

IncompressibleSolver::IncompressibleSolver()
    : RiemannSolver(), m_pointSolve(true)
{
    m_requiresRotation = false;
}

void IncompressibleSolver::v_Solve(
    const int nDim, const Array<OneD, const Array<OneD, NekDouble>> &Fwd,
    const Array<OneD, const Array<OneD, NekDouble>> &Bwd,
    Array<OneD, Array<OneD, NekDouble>> &flux)
{
    if (m_pointSolve)
    {
        size_t expDim = nDim;
        NekDouble vf{}, wf{}, Ef{}; // Ef(energy) has not been implemented

        if (expDim == 1)
        {
            for (size_t i = 0; i < Fwd[0].size(); ++i)
            {
                v_PointSolve(Fwd[0][i], 0.0, 0.0, 0.0,
                             Bwd[0][i], 0.0, 0.0, 0.0,
                             m_vectors["N"]()[0][i], 0.0, 0.0,
                             flux[0][i], vf, wf, Ef);
            }
        }
        else if (expDim == 2)
        {
            for (size_t i = 0; i < Fwd[0].size(); ++i)
            {
                v_PointSolve(Fwd[0][i], Fwd[1][i], 0.0, 0.0,
                             Bwd[0][i], Bwd[1][i], 0.0, 0.0,
                             m_vectors["N"]()[0][i], m_vectors["N"]()[1][i], 0.0,
                             flux[0][i], flux[1][i], wf, Ef);
            }
        }
        else if (expDim == 3)
        {
            for (size_t i = 0; i < Fwd[0].size(); ++i)
            {
                v_PointSolve(Fwd[0][i], Fwd[1][i], Fwd[2][i], 0.0,
                             Bwd[0][i], Bwd[1][i], Bwd[2][i], 0.0,
                             m_vectors["N"]()[0][i], m_vectors["N"]()[1][i], m_vectors["N"]()[2][i], 
                             flux[0][i], flux[1][i],flux[2][i], Ef);
            }
        }
    }
    else
    {
        v_ArraySolve(Fwd, Bwd, flux);
    }
}



} // namespace Nektar

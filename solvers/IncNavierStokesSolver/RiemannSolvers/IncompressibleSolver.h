///////////////////////////////////////////////////////////////////////////////
//
// File: CompressibleSolver.h
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

#ifndef NEKTAR_SOLVERS_INCNAVIERSTOKESSOLVER_RIEMANNSOLVER_INCOMPRESSIBLESOLVER
#define NEKTAR_SOLVERS_INCNAVIERSTOKESSOLVER_RIEMANNSOLVER_INCOMPRESSIBLESOLVER

#include <boost/core/ignore_unused.hpp>

// #include <CompressibleFlowSolver/Misc/EquationOfState.h>
#include <SolverUtils/RiemannSolvers/RiemannSolver.h>

using namespace Nektar::SolverUtils;

namespace Nektar
{
class IncompressibleSolver : public RiemannSolver
{
protected:
    bool m_pointSolve;
    // EquationOfStateSharedPtr m_eos;
    // bool m_idealGas;

    /// Session ctor
    IncompressibleSolver(const LibUtilities::SessionReaderSharedPtr &pSession);

    /// Programmatic ctor
    IncompressibleSolver();

    using ND = NekDouble;

    virtual void v_Solve(const int nDim,
                         const Array<OneD, const Array<OneD, ND>> &Fwd,
                         const Array<OneD, const Array<OneD, ND>> &Bwd,
                         Array<OneD, Array<OneD, ND>> &flux) override;

    virtual void v_ArraySolve(const Array<OneD, const Array<OneD, ND>> &Fwd,
                              const Array<OneD, const Array<OneD, ND>> &Bwd,
                              Array<OneD, Array<OneD, ND>> &flux)
    {
        boost::ignore_unused(Fwd, Bwd, flux);
        NEKERROR(ErrorUtil::efatal,
                 "This function should be defined by subclasses.");
    }
    // 有什么用？
    virtual void v_PointSolve(ND uL, ND vL, ND wL, ND EL,
                              ND uR, ND vR, ND wR, ND ER,
                              ND nx, ND ny, ND nz,
                              ND &uf, ND &vf, ND &wf, ND &Ef)
    {
        boost::ignore_unused(uL, vL, wL, EL, uR, vR,
                             wR, ER, uf, vf, wf, Ef);
        NEKERROR(ErrorUtil::efatal,
                 "This function should be defined by subclasses.");
    }

    // ND GetRoeSoundSpeed(ND rhoL, ND pL, ND eL, ND HL, ND srL, ND rhoR, ND pR,
    //                     ND eR, ND HR, ND srR, ND HRoe, ND URoe2, ND srLR);
};
} // namespace Nektar

#endif

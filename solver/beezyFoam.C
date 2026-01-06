/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    multiphaseInterFoam

Description
    Solver for n incompressible fluids which captures the interfaces and
    includes surface-tension and contact-angle effects for each phase, with
    optional mesh motion and mesh topology changes.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// ============================================================
// EXTRA BLOCK
// ============================================================   
// #include <limits>

#include "multiphaseMixture.H"
#include "incompressibleMomentumTransportModels.H"
#include "pimpleControl.H"
#include "pressureReference.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "CorrectPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"

    turbulence->validate();

    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    // const surfaceScalarField& rhoPhi(mixture.rhoPhi());                   // removed 02.01.26 (1.)

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (pimple.run(runTime))
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "alphaCourantNo.H"                                         
        #include "setDeltaT.H"

        fvModels.preUpdateMesh();

        mesh.update();

        runTime++;

        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        while (pimple.loop())
        {
            if (pimple.firstPimpleIter() || moveMeshOuterCorrectors)
            {
                scalar timeBeforeMeshUpdate = runTime.elapsedCpuTime();

                mesh.move();

                if (mesh.changing())
                {
                    Info<< "Execution time for mesh.update() = "
                        << runTime.elapsedCpuTime() - timeBeforeMeshUpdate
                        << " s" << endl;

                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;

                    MRF.update();

                    if (correctPhi)
                    {
                        #include "correctPhi.H"
                    }

                    mixture.correct();

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            fvModels.correct();
 
            // 1.) multiphase: alpha solving (MULES)
            mixture.solve();
            mixture.correct();                                                       // added 02.01.26 (1.)
            // rho = mixture.rho();                                                  // removed 02.01.26 (1.)

            const int nEnergyCorr =                                                  // added 02.01.26 (1.)
            pimple.dict().lookupOrDefault<int>("nEnergyCorrectors", 2);              // added 02.01.26 (1.)
             
            // --- 2) Inner coupling loop: energy <-> phase change <-> properties    // added 02.01.26 (1.)
           for (int iEnergy=0; iEnergy<nEnergyCorr; ++iEnergy)                       // added 02.01.26 (1.)
           {                                                                         // added 02.01.26 (1.)
              // Update density from current alpha 
              rho = mixture.rho();                                                   // added 02.01.26 (1.)

              // Update mass flux
              mixture.correct();                                                     // added 02.01.26 (1.)
              const surfaceScalarField& rhoPhi = mixture.rhoPhi();                   // added 02.01.26 (1.)

              // Solve enthalpy + phase change
              #include "hEqn.H"                                                      // added 02.01.26 (1.)
            }
            // 2. ) energy equtaion h and T update (contains updateFs and updateKEff)
            // 2.1) calculate temperature by enthalpy
            // 2.2) calculate solid fraction --- #include "updateFs.H"
            // 2.3) calculate thermal conductivity --- #include "updateKEff.H"
            // 2.4) include latent heat source term
            // 2.5) enthalpy equation
            // 2.6) calculate temperature by enthalpy
           //  #include "hEqn.H"
            
            // --- DEBUG: log T and h on plattform patch ---
            // const label plattformPatchID =
            // mesh.boundaryMesh().findPatchID("plattform");

            // if (plattformPatchID != -1)
            // {
            // const scalar TminP =
            // gMin(T.boundaryField()[plattformPatchID]);
            // const scalar TmaxP =
            // gMax(T.boundaryField()[plattformPatchID]);

            // const scalar hminP =
            // gMin(h.boundaryField()[plattformPatchID]);
            // const scalar hmaxP =
            // gMax(h.boundaryField()[plattformPatchID]);

            // Info<< "plattform patch | "
            //     << "T[min,max]=(" << TminP << ", " << TmaxP << ") "
            //     << "h[min,max]=(" << hminP << ", " << hmaxP << ")"
            //     << nl;
            // }

            // 3. ) Momentum equation U
            // 3.1) Mushy zone declaration --- darcy equation
            // 3.2) momentum equation contaning darcy damping 

            rho = mixture.rho();                                                     // added 02.01.26 (1.)
            const surfaceScalarField& rhoPhi = mixture.rhoPhi();                     // added 02.01.26 (1.)
            #include "UEqn.H"

            // 4.) pressure equation/ correction p_rgh
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

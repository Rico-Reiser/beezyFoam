MRF.correctBoundaryVelocity(U);

// =====================
// Sicherheit fs
// =====================
fs.max(0.0);
fs.min(1.0);

// --- Mushy-zone Darcy damping ---
dimensionedScalar Amushy
(
    "Amushy",
    dimensionSet(0, 0, -1, 0, 0, 0, 0),  // 1/s
    1e4
);
// ============================================================
// EXTRA BLOCK
// ============================================================
// fsDarcy: nur obere Begrenzung (keine Untergrenze!)
// volScalarField fsDarcy
//(
//    "fsDarcy",
//    min(max(fs, scalar(0.0)), scalar(0.95))
//);

volScalarField DarcyCoeff
(
    "DarcyCoeff",
    Amushy * sqr(fs) / (pow3(1.0 - fs) + scalar(1e-3))
    // "DarcyCoeff",
    // Amushy * sqr(fsDarcy) / (pow3(1.0 - fsDarcy) + scalar(1e-3))
);

   Info<< "fs[min,max]=(" << gMin(fs) << "," << gMax(fs)
    << ") Darcy[max]=" << gMax(DarcyCoeff) << endl;

// ============================================================
// EXTRA BLOCK  --- Kopplung des Darcy Terms nur an alpha.water      
// ============================================================
// Nur Wasser ---- 
if (mesh.foundObject<volScalarField>("alpha.water"))
{
    const volScalarField& alphaW =
        mesh.lookupObject<volScalarField>("alpha.water");

    DarcyCoeff *= pos(alphaW - SMALL);
}

// Darcy erst aktivieren, wenn wirklich Mushy (fs > 0.05)
   DarcyCoeff *= pos(fs - 0.05);

// Darc Begrenzungen (CLAMP)
// DarcyCoeff.max(1e7);
   DarcyCoeff.min(1e7);

// ============================================================
// EXTRA BLOCK
// ============================================================
// Debug
// Info<< "fs[min,max]=(" << gMin(fs) << "," << gMax(fs) << nl << endl;
//
//     << ") fsDarcy[min,max]=(" << gMin(fsDarcy) << "," << gMax(fsDarcy)
//     << ") Darcy[max]=" << gMax(DarcyCoeff) << nl << endl;

// =====================
// Momentum equation
// =====================
fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(rhoPhi, U)
  + MRF.DDt(rho, U)
  + turbulence->divDevTau(rho, U)
  // + fvm::Sp(rho * DarcyCoeff, U) // impliziter Darcy-Term
  // + rho * DarcyCoeff * U
  - fvm::Sp(rho*DarcyCoeff, U)                                       // added 02.01.25 (1.)

 ==
    fvModels.source(rho, U)
);

UEqn.relax();
fvConstraints.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve
    (
        UEqn
     ==
        fvc::reconstruct
        (
            (
                mixture.surfaceTensionForce()
              - ghf*fvc::snGrad(rho)
              - fvc::snGrad(p_rgh)
            ) * mesh.magSf()
        )
    );

    fvConstraints.constrain(U);
}

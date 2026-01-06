// --- 1) Temperatur konsistent aus h
T = Tref + h/cp;
T.correctBoundaryConditions();

// ============================================================
// EXTRA BLOCK
// ============================================================
// ADDITIONS: platform patch nur einmal am Ende (02.01.26)

// const label pidPlat = mesh.boundaryMesh().findPatchID("plattform");
// ---2) Patch-Temperatur explizit aus h setzen ===
// const label pid = mesh.boundaryMesh().findPatchID("plattform");
// if (pid != -1)
// {
//     T.boundaryFieldRef()[pid] =
//         Tref.value() + h.boundaryField()[pid]/cp.value();
// }
// T.correctBoundaryConditions();
// ============================================================
// EXTRA BLOCK
// ============================================================
// #include "updateKEff.H"

// --- 3) fs_old explizit sichern 
volScalarField fsOld
(
    IOobject
    (
        "fsOld",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fs
);
// ============================================================
// EXTRA BLOCK
// ============================================================
// 1) h ohne Latentwärme
//{
//    fvScalarMatrix hEqn
//    (
//        fvm::ddt(rho, h)
//      + fvm::div(rhoPhi, h)
//      - fvm::laplacian(kEff/cp, h)
//    );
//
//    hEqn.relax();
//    hEqn.solve();
//}

// --- 4) Solid fraction aus T aktualisieren
#include "updateFs.H"

fs.max(0.0);                               // added 02.01.26 (1.)
fs.min(1.0);                               // added 02.01.26 (1.)

// fs.correctBoundaryConditions();

// Sicherheit
if (mesh.foundObject<volScalarField>("alpha.water"))
{
    const volScalarField& alphaW =
        mesh.lookupObject<volScalarField>("alpha.water");

    forAll(fs, celli)
    {
        fs[celli] = min(fs[celli], alphaW[celli]);
    }
}

fs.correctBoundaryConditions();

// --- 5) kEff aus fs
#include "updateKEff.H"

// --- 6) d(fs)/dt sauber berechnen
// const dimensionedScalar dt = runTime.deltaT();              // removed 02.01.26 (1.)
volScalarField dFsdt
(
    IOobject
    (
        "dFsdt",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    (fs - fsOld)/runTime.deltaT()
);
const volScalarField& alphaW =
    mesh.lookupObject<volScalarField>("alpha.water");

dFsdt *= pos(alphaW - SMALL);                                  // added 02.01.26 (1.)

// === 7: dFsdt nur zulassen, wenn Wasser vorhanden ===        // removed 02.01.26 (1.)
// if (mesh.foundObject<volScalarField>("alpha.water"))        // removed 02.01.26 (1.)
// {
//     const volScalarField& alphaW =                          // removed 02.01.26 (1.)
//        mesh.lookupObject<volScalarField>("alpha.water");    // removed 02.01.26 (1.)
// 
//     dFsdt *= pos(alphaW - SMALL);                           // removed 02.01.26 (1.)
// }
// ============================================================
// EXTRA BLOCK
// ============================================================
//T = Tref + h/cp;
//T.correctBoundaryConditions();



// ---- WICHTIG: DIMENSIONIERTE LIMITER!
dimensionedScalar dFsdtMax("dFsdtMax", dimless/dimTime, 1e2);  // eher 1e2 statt 1e3
dFsdt.min(dFsdtMax);
dFsdt.max(-dFsdtMax);

//dimensionedScalar zeroDf("zeroDf", dimless/dimTime, 0.0);
//dimensionedScalar dFsdtMax("dFsdtMax", dimless/dimTime, 1.0);

//dFsdt.max(zeroDf);                                           // kein Schmelzen
//dFsdt.min(dFsdtMax);                                         // Stabilitaet

// --- 8) Latent heat: verhindert unphysiaklische Aufheizen

volScalarField mushyMask
(
    IOobject
    (
        "mushyMask",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    // pos(T - Tsolid) * pos(Tliquid - T)                      // removed 02.01.26 (1.)
    pos(T - Tsolid) * pos(Tliquid - T) * pos(alphaW - SMALL)   // added 02.01.26 (1.)
);

volScalarField Sh
(
    IOobject
    (
        "Sh",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho * L * dFsdt * mushyMask
);
Info<< "Tsolid=" << Tsolid.value()
    << " Tliquid=" << Tliquid.value() << nl;

Info<< "dFsdt[min,max]=(" << gMin(dFsdt) << "," << gMax(dFsdt) << ")\n";
Info<< "Sh[min,max]=(" << gMin(Sh) << "," << gMax(Sh) << ")\n";



// ============================================================
// EXTRA BLOCK
// ============================================================
// latent heat source ohne Einbindung mushy-zone
//volScalarField Sh
//(
//    IOobject
//    (
//        "Sh",
//        runTime.timeName(),
//        mesh,
//        IOobject::NO_READ,
//        IOobject::NO_WRITE
//    ),
//    rho * L * dFsdt
//);

// --- 6) Gleichung fuer sensible h = cp*T
{
    fvScalarMatrix hEqn
    (
        fvm::ddt(rho, h)
      + fvm::div(rhoPhi, h)
      - fvm::laplacian(kEff/cp, h)
     ==
        Sh
    );

    hEqn.relax();
    hEqn.solve();
}

// --- 8) Temperatur aktualisieren
T = Tref + h/cp;
T.correctBoundaryConditions();

const label pid2 = mesh.boundaryMesh().findPatchID("bottom");
if (pid2 != -1)
{
    T.boundaryFieldRef()[pid2] =
        Tref.value() + h.boundaryField()[pid2]/cp.value();
}

// ---X) Debugs
Info<< "h[min,max]=(" << gMin(h) << "," << gMax(h) << ")\n";
Info<< "T[min,max]=(" << gMin(T) << "," << gMax(T) << ")\n";
Info<< "fs[min,max]=(" << gMin(fs) << "," << gMax(fs) << ")\n";
// ============================================================
// EXTRA BLOCK
// ============================================================
// Debugs
//Info<< "FS DEBUG\n";
//Info<< "  T(min,max) = (" << gMin(T) << ", " << gMax(T) << ")\n";
//Info<< "  fs(min,max) = (" << gMin(fs) << ", " << gMax(fs) << ")\n";
//Info<< "  dFsdt(max) = " << gMax(dFsdt)
//    << "  Sh(max) = " << gMax(Sh) << nl;

//if (pidPlat != -1)
//{
//    const fvPatchScalarField& Tp =
//        T.boundaryField()[pidPlat];
//    Info<< "plattform TcellNear(min,max)=("
//        << gMin(Tp.patchInternalField()) << ", "
//        << gMax(Tp.patchInternalField()) << ")\n";
//}

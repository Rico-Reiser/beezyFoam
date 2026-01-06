Info<< "Reading field p_rgh\n" << endl;
// 1.) definition pressure field
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
// 2.) definition velocity field
Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
// 3.) include flux calculation
#include "createPhi.H"

// 4.) include fluid interactions
multiphaseMixture mixture(U, phi);

// 5.) density field
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mixture.rho()
);
rho.oldTime();

// 6.) Construct incompressible turbulence model
autoPtr<incompressible::momentumTransportModel> turbulence
(
    incompressible::momentumTransportModel::New(U, phi, mixture)
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

// 7.) calculate pressure by p_rgh
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

pressureReference pressureReference(p, p_rgh, pimple.dict());

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pressureReference.refValue()
      - getRefCellValue(p, pressureReference.refCell())
    );
    p_rgh = p - rho*gh;
}

mesh.schemes().setFluxRequired(p_rgh.name());

#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"

// --- addition --- thermophysical ---
// ============================================================
// EXTRA BLOCK
// ============================================================ 
// 8.) Temperature --- added 08.12.25
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("T0", dimTemperature, 290.0)
);


// ============================================================
// EXTRA BLOCK
// ============================================================
// Info<< "Creating Twater diagnostic field\n" << endl;
// volScalarField Twater
// (
//     IOobject
//     (
//         "Twater",
//          runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("zero", T.dimensions(), 0.0)
// );


// 9.) Solid fraction fs
volScalarField fs
(
    IOobject
    (
        "fs",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);


// 10.) solidification Properties --- added 08.12.25
Info<< "Reading solidificationProperties\n" << endl;

IOdictionary solidProps
(
    IOobject
    (
        "solidificationProperties",
        "constant",
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// 10 alternative
//dimensionedScalar Tsolid  ("Tsolid", dimTemperature, solidProps.lookup("Tsolid"));
//dimensionedScalar Tliquid ("Tliquid", dimTemperature, solidProps.lookup("Tliquid"));
//dimensionedScalar kSolid  ("kSolid", dimPower/dimLength/dimTemperature, solidProps.lookup("kSolid"));
//dimensionedScalar kLiquid ("kLiquid", dimPower/dimLength/dimTemperature, solidProps.lookup("kLiquid"));
//dimensionedScalar cp      ("cp", dimEnergy/dimMass/dimTemperature, solidProps.lookup("cp"));
//dimensionedScalar L       ("L", dimEnergy/dimMass, solidProps.lookup("latentHeat"));

//10.1) temperature --- solid phase
dimensionedScalar Tsolid
(
    "Tsolid",
    dimTemperature,
    solidProps.lookup("Tsolid")
);
// 10.2) temeprature --- liquid phase
dimensionedScalar Tliquid
(
    "Tliquid",
    dimTemperature,
    solidProps.lookup("Tliquid")
);
// 10.3) heat conduction --- solid phase
dimensionedScalar kSolid // [1 1 -3 -1 0 0]
(
    "kSolid",
    dimensionSet(1,1,-3,-1,0,0,0),
    solidProps.lookup("kSolid")
);
// 10.4) heat conduction --- liquid phase
dimensionedScalar kLiquid // [1 1 -3 -1 0 0]
(
    "kLiquid",
    dimensionSet(1,1,-3,-1,0,0,0),
    solidProps.lookup("kLiquid")
);
// 10.5) heat capacity --- dim. mentioned in solProps
dimensionedScalar cp      // [0 2 -2 -1 0 0]
(
    "cp",
    dimensionSet(0,2,-2,-1,0,0,0),
    solidProps.lookup("cp")
);

// 10.6) latent heat 
dimensionedScalar L     // [0 2 -2 0 0 0 0]
(
    "L",
    dimensionSet(0,2,-2,0,0,0,0),
    solidProps.lookup("latentHeat")
);

// 11.) use script updateFs.H --- solid fraction 
//if (mesh.foundObject<volScalarField>("alpha.water"))
//{
//    #include "updateFs.H"
//
//    const volScalarField& alphaW =
//        mesh.lookupObject<volScalarField>("alpha.water");
//
//    fs.max(0.0);
//
//    forAll(fs, celli)
//    {
//        fs[celli] = min(fs[celli], alphaW[celli]);
//    }
//
//    fs.correctBoundaryConditions();
//}
//else
//{
//    fs = scalar(0);
//    fs.correctBoundaryConditions();
//}
if (mesh.foundObject<volScalarField>("alpha.water"))
{
    const volScalarField& alphaW =
        mesh.lookupObject<volScalarField>("alpha.water");

    // initial fs from initial T (no storeOldTime, no relaxation)
    forAll(fs, celli)
    {
        const scalar Ti = T[celli];
        scalar fi;

        if (Ti <= Tsolid.value())      fi = 1.0;
        else if (Ti >= Tliquid.value()) fi = 0.0;
        else
        {
            const scalar xi = (Ti - Tsolid.value())/((Tliquid.value() - Tsolid.value()) + VSMALL);
            fi = 0.5*(1.0 + Foam::cos(constant::mathematical::pi*min(max(xi, 0.0), 1.0)));
        }

        fs[celli] = alphaW[celli]*fi;
    }

    fs.max(0.0);
    forAll(fs, celli)
    {
        fs[celli] = min(fs[celli], alphaW[celli]);
    }

    fs.correctBoundaryConditions();
}
else
{
    fs = scalar(0);
    fs.correctBoundaryConditions();
}


// compute fs from initial temperature --- added 08.12.25
//forAll(fs, i)
//{
//    const scalar Ti = T[i];
//    scalar f = (Tsolid.value() - Ti)/(Tsolid.value() - Tliquid.value());
//    fs[i] = min(max(f, 0.0), 1.0);
//}
//fs.correctBoundaryConditions();

// 12.) Enthalpy --- added 08.12.25
// 12.1) defintion volScalarField
dimensionedScalar Tref
(
    "Tref",
    dimTemperature,
    Tliquid.value()   // oder 273.15
);

volScalarField h
(
    IOobject
    (
        "h",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimEnergy/dimMass, 0)   // [0 2 -2 0 0 0 0]
);

//if (mesh.time().timeIndex() == 0)
//{
//    // sensible enthalpy relative zu Tref
//    h = cp*(T - Tref);
//    h.correctBoundaryConditions();
//}

//if (mesh.time().timeIndex() == 0 && !IOobject("h", runTime.timeName(), mesh).headerOk())
//{
//    h = cp*(T - Tref);
//    h.correctBoundaryConditions();
//}
// Initialisation only if h does not exist
if (runTime.timeIndex() == 0 &&
    !IOobject("h", runTime.timeName(), mesh).headerOk())
{
    h = cp*(T - Tref);
    h.correctBoundaryConditions();
}
T = Tref + h/cp;
T.correctBoundaryConditions();

// 13.)  effective conductivity --- alternatives : test both
// 13.1) volScalarField
volScalarField kEff
(
    IOobject
    (
        "kEff",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("kEff", kSolid.dimensions(), 0.0)
);

if (mesh.foundObject<volScalarField>("alpha.water"))                         // added 02.01.26 (1.)
{                                                                            // added 02.01.26 (1.)
    const volScalarField& alphaW =                                           // added 02.01.26 (1.)
        mesh.lookupObject<volScalarField>("alpha.water");                    // added 02.01.26 (1.)
 
    kEff =                                                                   // added 02.01.26 (1.)
        alphaW*(fs*kSolid + (scalar(1) - fs)*kLiquid)                        // added 02.01.26 (1.)
      + (scalar(1) - alphaW)*kLiquid;   // Luft: kLiquid oder kAir           // added 02.01.26 (1.)
}                                                                            // added 02.01.26 (1.)
else                                                                         // added 02.01.26 (1.)
{                                                                            // added 02.01.26 (1.)
    kEff = kLiquid;                                                          // added 02.01.26 (1.)
}                                                                            // added 02.01.26 (1.)

kEff.correctBoundaryConditions();                                            // added 02.01.26 (1.)


// if (mesh.time().timeIndex() == 0)                                          // removed 02.01.26 (1.)
// {
//     kEff = fs*kSolid + (scalar(1) - fs)*kLiquid;                           // removed 02.01.26 (1.)
//     kEff.correctBoundaryConditions();                                      // removed 02.01.26 (1.)
// }                                                                          // removed 02.01.26 (1.)

// 13.3) alternative through includes 
// #include updateKEff.H
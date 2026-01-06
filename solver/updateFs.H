{
    const volScalarField& alphaW =
        mesh.lookupObject<volScalarField>("alpha.water");

    const scalar Tsol = Tsolid.value();
    const scalar Tliq = Tliquid.value();
    const scalar dT   = (Tliq - Tsol) + VSMALL;

    // ============================================================
    // 1) fs aus Temperatur berechnen (nur Wasser)
    // ============================================================
    fs.storeOldTime();

    forAll(fs, celli)
    {
        const scalar Ti = T[celli];
        scalar fi;

        if (Ti <= Tsol)    
           {
           fi = 1.0;             // 100 % solid
           }
        else if (Ti >= Tliq)
           {       
           fi = 0.0;             // 100% liquid
           }
        else
        {
            scalar xi = (Ti - Tsol)/dT;                                 // mushy Zone
            xi = min(max(xi, 0.0), 1.0);
            fi = 0.5*(1.0 + Foam::cos(constant::mathematical::pi*xi));
        }

        fs[celli] = alphaW[celli] * fi;
    }

    // bounds
    fs.max(0.0);
    fs.min(1.0);                                                       // added 02.01.26 (1.)


    forAll(fs, celli)
    {
        fs[celli] = min(fs[celli], alphaW[celli]);
    }

    // Relaxation (fsOld kommt aus hEqn.H!)
    const scalar fsRelax = 0.5;                       // alternativ größer  // 02.01.26 (1.): 0.05 - langsamer Phasenwechsel // 02.01.26 (2.): 0.25 // 05.01.26 (3.): 0.5
    fs = (1.0 - fsRelax)*fsOld + fsRelax*fs;

    fs.correctBoundaryConditions();

    // ============================================================
    // EXTRA BLOCK
    // ============================================================
    // Info<< "FS DEBUG\n"
    // << "  T(min,max) = (" << gMin(T) << ", " << gMax(T) << ")\n"
    // << "  fs(min,max) = (" << gMin(fs) << ", " << gMax(fs) << ")\n"
    // << "  alphaW(min,max) = ("
    // << gMin(alphaW) << ", " << gMax(alphaW) << ")\n"
    // << endl;

}



{
    
    // kEff = fs*kSolid + (scalar(1) - fs)*kLiquid;
    
    
    // ============================================================
    // EXTRA BLOCK
    // ============================================================    
    // Effektive Wärmeleitfähigkeit abhängig vom Aggregatzustand

    // kEff =
    //     fs*kSolid
    //   + (scalar(1) - fs)*kLiquid;

    // minimale Leitfähigkeit gegen numerische Degeneration
    // kEff.max
    // (
    //     dimensionedScalar
    //     (
    //         "kMin",
    //         kEff.dimensions(),
    //         1e-6
    //     )
    // );

    // kEff.correctBoundaryConditions();
        if (mesh.foundObject<volScalarField>("alpha.water"))
    {
        const volScalarField& alphaW =
            mesh.lookupObject<volScalarField>("alpha.water");

        // kEff:
        // - Eis/Wasser nur im Wasseranteil
        // - Luft bekommt keine feste Leitfähigkeit
        kEff =
            alphaW*(fs*kSolid + (scalar(1) - fs)*kLiquid)
          + (scalar(1) - alphaW)*kLiquid;   // optional: kAir
    }
    else
    {
        // Fallback: kein VOF
        kEff = fs*kSolid + (scalar(1) - fs)*kLiquid;
    }

    kEff.correctBoundaryConditions();
}

    // ============================================================
    // EXTRA BLOCK
    // ============================================================   

    Info<< "kEff[min,max]=("
    << gMin(kEff) << ", "
    << gMax(kEff) << ")\n";
    // Debug
    // Info<< "kEff[min,max]=(" << gMin(kEff) << "," << gMax(kEff) << ")\n";

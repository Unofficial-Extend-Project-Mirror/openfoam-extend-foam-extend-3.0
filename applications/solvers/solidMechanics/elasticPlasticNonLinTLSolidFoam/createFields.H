    Info<< "Reading field DU\n" << endl;
    volVectorField DU
    (
        IOobject
        (
            "DU",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volTensorField gradDU = fvc::grad(DU);

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
	mesh,
	dimensionedVector("zero", dimLength, vector::zero)
    );

    volTensorField gradU
    (
        IOobject
        (
	 "grad(U)",
	 runTime.timeName(),
	 mesh,
	 IOobject::READ_IF_PRESENT,
	 IOobject::AUTO_WRITE
	 ),
        mesh,
        dimensionedTensor("zero", dimless, tensor::zero)
     );

    //- Increment of Green finite strain tensor
    volSymmTensorField DEpsilon
    (
        IOobject
        (
            "DEpsilon",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
	mesh,
	dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    );

    volSymmTensorField epsilon
    (
        IOobject
        (
            "epsilon",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
	mesh,
	dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    );

    //- plastic strain
    volSymmTensorField epsilonP
    (
        IOobject
        (
            "epsilonP",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    );


    //- Increment of 2nd Piola-Kirchhoff stress tensor
    volSymmTensorField DSigma
    (
        IOobject
        (
            "DSigma",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
    );

    //- 2nd Piola-Kirchhoff stress tensor
    volSymmTensorField sigma
    (
        IOobject
        (
            "sigma",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
	mesh,
	dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
    );

      volVectorField divDSigmaExp
      (
        IOobject
        (
	 "divDSigmaExp",
	 runTime.timeName(),
	 mesh,
	 IOobject::NO_READ,
	 IOobject::NO_WRITE
	 ),
        mesh,
        dimensionedVector("zero", dimensionSet(1,-2,-2,0,0,0,0), vector::zero)
       );

      volVectorField divDSigmaNonLinExp
      (
        IOobject
        (
	 "divDSigmaNonLinExp",
	 runTime.timeName(),
	 mesh,
	 IOobject::NO_READ,
	 IOobject::NO_WRITE
	 ),
        mesh,
        dimensionedVector("zero", dimensionSet(1,-2,-2,0,0,0,0), vector::zero)
       );

    constitutiveModel rheology(sigma, DU);
    
    volScalarField rho = rheology.rho();

    volScalarField mu = rheology.mu();
    volScalarField lambda = rheology.lambda();
    surfaceScalarField muf = fvc::interpolate(mu, "mu");
    surfaceScalarField lambdaf = fvc::interpolate(lambda, "lambda");

    surfaceVectorField n = mesh.Sf()/mesh.magSf();

   // plastic strain increment
   const volSymmTensorField& DEpsilonP = rheology.DEpsilonP();

   // for aitken relaxation
   volVectorField aitkenDelta
    (
        IOobject
        (
            "aitkenDelta",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
	dimensionedVector("zero", dimLength, vector::zero)
    );
// aitken relaxation factor
scalar aitkenInitialRes = 1.0;
scalar aitkenTheta = 0.1;

//    volVectorField resid
//    (
//         IOobject
//         (
//          "resid",
//          runTime.timeName(),
//          mesh,
//          IOobject::NO_READ,
//          IOobject::AUTO_WRITE
//          ),
//         mesh,
//         dimensionedVector("zero", dimless, vector::zero)
//     );

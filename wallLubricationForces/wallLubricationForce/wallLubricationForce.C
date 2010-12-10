/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "wallLubricationForce.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

namespace multiphase
{

namespace surfaceForces
{

    defineTypeNameAndDebug(wallLubricationForce, 0);
    defineRunTimeSelectionTable(wallLubricationForce, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

wallLubricationForce::wallLubricationForce
(
    const dictionary& multiphaseSurfaceForcesDict,
    const multiphase::transport& mtm,
    const PtrList<volScalarField>& alpha,
    const volScalarField& beta
)
:
	multiphaseSurfaceForcesDict_(multiphaseSurfaceForcesDict),
    alpha_(alpha),
    beta_(beta),
    mtm_(mtm),
    g_
    (
    	IOobject
        (
            "g",
            beta.time().constant(),
            beta.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    y_(beta.mesh()),
    yr_(beta.mesh()),
    Eotvos_
    (
        IOobject
        (
            "Eotvos",
            beta.time().timeName(),
            beta.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        beta.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)
    ),
	wallForce_
	(
		IOobject
		(
			"wallForce",
			beta.time().timeName(),
			beta.mesh(),
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
        beta.mesh(),
        dimensionedVector("zero", dimensionSet(0, 0, 0, 0, 0), vector(0.0, 0.0, 0.0))
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

wallLubricationForce::~wallLubricationForce()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceForces

} // End namespace multiphase

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //

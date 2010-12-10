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

#include "Antal.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

namespace multiphase
{

namespace surfaceForces
{

    defineTypeNameAndDebug(Antal, 0);

    addToRunTimeSelectionTable
    (
    	wallLubricationForce,
    	Antal,
        dictionary
    );


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Antal::Antal
(
    const dictionary& multiphaseSurfaceForcesDict,
    const multiphase::transport& mtm,
    const PtrList<volScalarField>& alpha,
    const volScalarField& beta
)
:
	wallLubricationForce(multiphaseSurfaceForcesDict, mtm, alpha, beta),
	parametersDict_(multiphaseSurfaceForcesDict.subDict("wallLubricationParameters")),
	modelParametersDict_(parametersDict_.subDict(typeName+"Parameters")),
    Cw1_(modelParametersDict_.lookup("Cw1")),
    Cw2_(modelParametersDict_.lookup("Cw2")),
    Cwl_
    (
        IOobject
        (
            "Cwl",
            beta.time().timeName(),
            beta.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        beta.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)
    ),
	Urw_
	(
		IOobject
		(
			"Urw",
			beta.time().timeName(),
			beta.mesh(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		beta.mesh(),
		dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)
	)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Antal::~Antal()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/*!
 * Calculate and return the wall lubrication force.
 *
 * \param[in] Ur Relative Velocity
 * \param[in] sauterdiameter/diameter Diameter of bubbles/droplets.
 */
Foam::tmp<Foam::volVectorField> Antal::wallForce
(
    const volVectorField& Ur,
    const volScalarField& diameter
)
{

	update(diameter);

	// Relative velocity in wall plane
	Urw_ = Foam::mag(Ur^yr_.n());

	wallForce_ = alpha_[0]*mtm_.rhoc()*Foam::sqr(Urw_)*Cwl_*yr_.n();
	wallForce_.correctBoundaryConditions();

	return wallForce_;

}

/*!
 * Update the variables used to calculate the wall lubrication force.
 */
void Antal::update(const volScalarField& diameter_)
{

	Cwl_ = (Cw1_/diameter_)+(Cw2_/y_.y());
	Cwl_.max(0.0);

    Cwl_.correctBoundaryConditions();

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceForces

} // End namespace multiphase

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //

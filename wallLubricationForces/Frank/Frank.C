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

#include "Frank.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

namespace multiphase
{

namespace surfaceForces
{

    defineTypeNameAndDebug(Frank, 0);

    addToRunTimeSelectionTable
    (
    	wallLubricationForce,
    	Frank,
        dictionary
    );


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Frank::Frank
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
    p_(modelParametersDict_.lookup("p")),
    Cw3_
    (
        IOobject
        (
            "Cw3",
            beta.time().timeName(),
            beta.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        beta.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)
    ),
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

Frank::~Frank()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/*!
 * Calculate and return the wall lubrication force.
 *
 * \param[in] Ur Relative Velocity
 * \param[in] sauterdiameter/diameter Diameter of bubbles/droplets.
 */
Foam::tmp<Foam::volVectorField> Frank::wallForce
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
void Frank::update(const volScalarField& diameter_)
{

	Eotvos_ = (mtm_.rhoc()-mtm_.rhod(0))*Foam::mag(g_)*Foam::sqr(diameter_)/mtm_.sigmad(0);

	Cw3_ = 0.00599*Eotvos_-0.0187;

	forAll(Eotvos_,celli)
	{
	    if (Eotvos_[celli] <= 5)
	    {
	        Cw3_[celli] = Foam::exp(-0.933*Eotvos_[celli]+0.179);
	    }
	    else if (Eotvos_[celli] > 33)
	    {
	        Cw3_[celli] = 0.179;
	    }
	}

	forAll(Eotvos_.boundaryField(),patchi)
	{
		forAll(Eotvos_.boundaryField()[patchi],facei)
		{
			if (Eotvos_.boundaryField()[patchi][facei] <= 5)
			{
				Cw3_.boundaryField()[patchi][facei] = Foam::exp(-0.933*Eotvos_.boundaryField()[patchi][facei]+0.179);
			}
			else if (Eotvos_.boundaryField()[patchi][facei] > 33)
			{
				Cw3_.boundaryField()[patchi][facei] = 0.179;
			}
		}
	}

    Cwl_ = Cw3_*(1/Cw2_)*( (1-(y_.y()/ (Cw1_*diameter_)) )  / (y_.y()*Foam::pow(y_.y()/(Cw1_*diameter_),   (p_-1))) );
    Cwl_.max(0.0);

    Cwl_.correctBoundaryConditions();

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceForces

} // End namespace multiphase

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //

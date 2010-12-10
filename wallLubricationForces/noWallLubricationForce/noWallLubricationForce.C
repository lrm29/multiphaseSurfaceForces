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

#include "noWallLubricationForce.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

namespace multiphase
{

namespace surfaceForces
{

    defineTypeNameAndDebug(noWallLubricationForce, 0);

    addToRunTimeSelectionTable
    (
    	wallLubricationForce,
    	noWallLubricationForce,
        dictionary
    );


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

noWallLubricationForce::noWallLubricationForce
(
    const dictionary& multiphaseSurfaceForcesDict,
    const multiphase::transport& mtm,
    const PtrList<volScalarField>& alpha,
    const volScalarField& beta
)
:
	wallLubricationForce(multiphaseSurfaceForcesDict, mtm, alpha, beta)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

noWallLubricationForce::~noWallLubricationForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/*!
 * Return zero for the wall lubrication force.
 *
 * \param[in] Ur Relative Velocity
 * \param[in] sauterdiameter/diameter Diameter of bubbles/droplets.
 */
Foam::tmp<Foam::volVectorField> noWallLubricationForce::wallForce
(
    const volVectorField& Ur,
    const volScalarField& diameter
)
{

	return 0.0*wallForce_;

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceForces

} // End namespace multiphase

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //

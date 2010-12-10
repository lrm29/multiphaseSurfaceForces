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

Class
    Foam::multiphase::surfaceForces::noWallLubricationForce

Description

SourceFiles
    noWallLubricationForce.C

\*---------------------------------------------------------------------------*/

#ifndef noWallLubricationForce_H
#define noWallLubricationForce_H

#include "wallLubricationForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace multiphase
{

namespace surfaceForces
{

/*---------------------------------------------------------------------------*\
                           Class noWallLubricationForce Declaration
\*---------------------------------------------------------------------------*/

/*!
 * \ingroup wallLubricationForces
 * \brief noWallLubricationForce wall lubrication force model.
 */
class noWallLubricationForce
:
    public wallLubricationForce
{

public:

    //- Runtime type information
    TypeName("noWallLubricationForce");


    // Constructors

        //- Construct from components
		noWallLubricationForce
        (
            const dictionary& multiphaseSurfaceForcesDict,
            const multiphase::transport& mtm,
            const PtrList<volScalarField>& alpha,
            const volScalarField& beta
        );


    //- Destructor
    virtual ~noWallLubricationForce();


    // Member Functions

        //! Wall Lubrication Force.
        tmp<volVectorField> wallForce
        (
            const volVectorField& Ur,
            const volScalarField& sauterDiameter
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceForces

} // End namespace multiphase

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

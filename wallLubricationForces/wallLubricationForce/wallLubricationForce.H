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
    Foam::multiphase::surfaceForces::wallLubricationForce

Description

SourceFiles
    wallLubricationForce.C
    newWallLubricationForce.C

\*---------------------------------------------------------------------------*/

#ifndef wallLubricationForce_H
#define wallLubricationForce_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dictionary.H"
#include "multiphaseTransport.H"
#include "runTimeSelectionTables.H"
#include "wallDistReflection.H"
#include "wallDist.H"
#include "uniformDimensionedFields.H"

namespace Foam
{

namespace multiphase
{

namespace surfaceForces
{

/*---------------------------------------------------------------------------*\
                           Class wallLubricationForce Declaration
\*---------------------------------------------------------------------------*/

//! \ingroup multiphaseSurfaceForces
//! \defgroup wallLubricationForces Multiphase Wall Lubrication Force Model.
//! \brief This module contains documentation of the multiphase wallLubricationForce models.

//! \brief wallLubricationForce class.
class wallLubricationForce
{

protected:

    // Protected data

        const dictionary& multiphaseSurfaceForcesDict_;
        const PtrList<volScalarField>& alpha_;
        const volScalarField& beta_;
        const multiphase::transport& mtm_;

        const uniformDimensionedVectorField g_;

        //  Distance to nearest wall
        wallDist y_;
        //  Wall normal unit vector
        wallDistReflection yr_;

        //! Eotvos dimensionless number.
        volScalarField Eotvos_;

        volVectorField wallForce_;

public:

    // Runtime type information
    TypeName("wallLubricationForce");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            wallLubricationForce,
            dictionary,
            (
                const dictionary& multiphaseSurfaceForcesDict,
                const multiphase::transport& mtm,
                const PtrList<volScalarField>& alpha,
                const volScalarField& beta
            ),
            (multiphaseSurfaceForcesDict, mtm, alpha, beta)
        );


    // Constructors

        wallLubricationForce
        (
            const dictionary& multiphaseSurfaceForcesDict,
            const multiphase::transport& mtm,
            const PtrList<volScalarField>& alpha,
            const volScalarField& beta
        );


    //! Destructor
    virtual ~wallLubricationForce();


    // Selectors

        static autoPtr<wallLubricationForce> New
        (
                const dictionary& multiphaseSurfaceForcesDict,
                const multiphase::transport& mtm,
                const PtrList<volScalarField>& alpha,
                const volScalarField& beta
        );


    // Member Functions

        //! The wall force
        virtual tmp<volVectorField> wallForce
        (
            const volVectorField& Ur,
            const volScalarField& sauterDiameter
        ) = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceForces

} // End namespace multiphase

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

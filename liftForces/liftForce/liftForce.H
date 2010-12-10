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
    Foam::multiphase::surfaceForces::liftForce

Description

SourceFiles
    liftForce.C
    newLiftForce.C

\*---------------------------------------------------------------------------*/

#ifndef liftForce_H
#define liftForce_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dictionary.H"
#include "multiphaseTransport.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

namespace multiphase
{

namespace surfaceForces
{

/*---------------------------------------------------------------------------*\
                           Class liftForce Declaration
\*---------------------------------------------------------------------------*/

//! \ingroup multiphaseSurfaceForces
//! \defgroup liftForces Multiphase Lift Force Model.
//! \brief This module contains documentation of the multiphase liftForce models.

//! \brief liftForce class.
class liftForce
{

protected:

    // Protected data

        const dictionary& multiphaseSurfaceForcesDict_;
        const PtrList<volScalarField>& alpha_;
        const volScalarField& beta_;
        const multiphase::transport& mtm_;

        volVectorField liftForce_;

public:

    // Runtime type information
    TypeName("liftForce");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            liftForce,
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

        liftForce
        (
            const dictionary& multiphaseSurfaceForcesDict,
            const multiphase::transport& mtm,
            const PtrList<volScalarField>& alpha,
            const volScalarField& beta
        );


    //! Destructor
    virtual ~liftForce();


    // Selectors

        static autoPtr<liftForce> New
        (
                const dictionary& multiphaseSurfaceForcesDict,
                const multiphase::transport& mtm,
                const PtrList<volScalarField>& alpha,
                const volScalarField& beta
        );


    // Member Functions

        //! The lift force used in the momentum eq.
        virtual tmp<volVectorField> lift
        (
            const volVectorField& Uc,
            const volVectorField& Ud
        ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceForces

} // End namespace multiphase

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

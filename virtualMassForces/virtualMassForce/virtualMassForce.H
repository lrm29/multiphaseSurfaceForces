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
    Foam::multiphase::surfaceForces::virtualMassForce

Description

SourceFiles
    virtualMassForce.C
    newVirtualMassForce.C

\*---------------------------------------------------------------------------*/

#ifndef virtualMassForce_H
#define virtualMassForce_H

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
                           Class virtualMassForce Declaration
\*---------------------------------------------------------------------------*/

//! \ingroup multiphaseSurfaceForces
//! \defgroup virtualMassForces Multiphase Virtual Mass Force Model.
//! \brief This module contains documentation of the multiphase virtualMassForce models.

//! \brief virtualMassForce class.
class virtualMassForce
{

protected:

    // Protected data

        const dictionary& multiphaseSurfaceForcesDict_;
        const PtrList<volScalarField>& alpha_;
        const volScalarField& beta_;
        const multiphase::transport& mtm_;

        volVectorField virtualMassForce_;

public:

    // Runtime type information
    TypeName("virtualMassForce");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            virtualMassForce,
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

        virtualMassForce
        (
            const dictionary& multiphaseSurfaceForcesDict,
            const multiphase::transport& mtm,
            const PtrList<volScalarField>& alpha,
            const volScalarField& beta
        );


    //! Destructor
    virtual ~virtualMassForce();


    // Selectors

        static autoPtr<virtualMassForce> New
        (
                const dictionary& multiphaseSurfaceForcesDict,
                const multiphase::transport& mtm,
                const PtrList<volScalarField>& alpha,
                const volScalarField& beta
        );


    // Member Functions

        //! The virtual mass force used in the momentum equation.
        virtual tmp<volVectorField> virtualMass
        (
            const volVectorField& DDtUc,
            const volVectorField& DDtUd
        ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceForces

} // End namespace multiphase

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

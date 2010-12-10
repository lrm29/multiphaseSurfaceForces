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
    Foam::multiphase::surfaceForces::dragForce

Description

SourceFiles
    dragForce.C
    newDragForce.C

\*---------------------------------------------------------------------------*/

#ifndef dragForce_H
#define dragForce_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dictionary.H"
#include "multiphaseTransport.H"
#include "runTimeSelectionTables.H"

//! \defgroup multiphaseSurfaceForces Multiphase Surface Forces Models.
//! multiphase::surfaceForces calculates interfacial forces.

namespace Foam
{

namespace multiphase
{

//! \brief Namespace containing surface force models to be used with multiphase solvers.
namespace surfaceForces
{

/*---------------------------------------------------------------------------*\
                           Class dragForce Declaration
\*---------------------------------------------------------------------------*/

//! \ingroup multiphaseSurfaceForces
//! \defgroup dragForces Multiphase Drag Force Model.
//! \brief This module contains documentation of the multiphase dragForce models.

//! \brief dragForce class.
class dragForce
{

protected:

    // Protected data

        const dictionary& multiphaseSurfaceForcesDict_;
        const PtrList<volScalarField>& alpha_;
        const volScalarField& beta_;
        const multiphase::transport& mtm_;

        volScalarField dragForce_;

public:

    // Runtime type information
    TypeName("dragForce");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            dragForce,
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

        dragForce
        (
            const dictionary& multiphaseSurfaceForcesDict,
            const multiphase::transport& mtm,
            const PtrList<volScalarField>& alpha,
            const volScalarField& beta
        );


    //! Destructor
    virtual ~dragForce();


    // Selectors

        static autoPtr<dragForce> New
        (
                const dictionary& multiphaseSurfaceForcesDict,
                const multiphase::transport& mtm,
                const PtrList<volScalarField>& alpha,
                const volScalarField& beta
        );


    // Member Functions

        //! The drag function K used in the momentum eq.
        virtual tmp<volScalarField> K
        (
            const volScalarField& Ur,
            const volScalarField& sauterDiameter,
            const int& dispersedPhase
        ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceForces

} // End namespace multiphase

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

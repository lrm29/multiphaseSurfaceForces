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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::multiphase::surfaceForces::wallLubricationForce>
Foam::multiphase::surfaceForces::wallLubricationForce::New
(
	const dictionary& multiphaseSurfaceForcesDict,
	const multiphase::transport& mtm,
	const PtrList<volScalarField>& alpha,
	const volScalarField& beta
)
{
    word wallLubricationForceType
    (
    	multiphaseSurfaceForcesDict.lookup("wallLubricationForce")
    );

    Info << "Selecting wallLubricationForce "
        << ": "
        << wallLubricationForceType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(wallLubricationForceType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "wallLubricationForce::New : " << endl
                << "    unknown wallLubricationForceType type "
                << wallLubricationForceType
                << ", constructor not in hash table" << endl << endl
                << "    Valid wallLubricationForce types are : " << endl;
        Info << dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return cstrIter()(multiphaseSurfaceForcesDict, mtm, alpha, beta);
}


// ************************************************************************* //

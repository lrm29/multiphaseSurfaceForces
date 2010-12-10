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

#include "SchillerNaumann.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

namespace multiphase
{

namespace surfaceForces
{

    defineTypeNameAndDebug(SchillerNaumann, 0);

    addToRunTimeSelectionTable
    (
    	dragForce,
        SchillerNaumann,
        dictionary
    );


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

SchillerNaumann::SchillerNaumann
(
    const dictionary& multiphaseSurfaceForcesDict,
    const multiphase::transport& mtm,
    const PtrList<volScalarField>& alpha,
    const volScalarField& beta
)
:
    dragForce(multiphaseSurfaceForcesDict, mtm, alpha, beta)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

SchillerNaumann::~SchillerNaumann()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> SchillerNaumann::K
(
    const volScalarField& Ur,
    const volScalarField& sauterDiameter,
    const int& dispersedPhase
)
{

    volScalarField Re = max( (1/mtm_.nuc()) * sauterDiameter * ( mag(Ur) ) , scalar(1.0e-3) );
    volScalarField Cds = 24.0*(scalar(1) + 0.15*pow(Re, 0.687))/Re;

    forAll(Re, celli)
    {
        if(Re[celli] > 1000.0)
        {
            Cds[celli] = 0.44;
        }
    }

    dragForce_ = 0.75*alpha_[dispersedPhase]*beta_*mtm_.rhoc()*(1/sauterDiameter)*Cds*Foam::mag(Ur);
    dragForce_.correctBoundaryConditions();

    return dragForce_;

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceForces

} // End namespace multiphase

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //

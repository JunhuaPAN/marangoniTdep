/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

#include "marangoniTdepFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"

#include "volFields.H"
#include "transformFvPatchFields.H"
#include "symmTransformField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

marangoniTdepFvPatchVectorField::marangoniTdepFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    transformFvPatchField<vector>(p, iF),
    //transformFvPatchVectorField(p, iF),
    gamma_(0),
    mu0_(0),
    A_(0)
{
    //this->checkVolField();
}


marangoniTdepFvPatchVectorField::marangoniTdepFvPatchVectorField
(
    const marangoniTdepFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    transformFvPatchField<vector>(ptf, p, iF, mapper),
    //transformFvPatchVectorField(ptf, p, iF, mapper),
    gamma_(ptf.gamma_),
    mu0_(ptf.mu0_),
    A_(ptf.A_)
{
    //this->checkVolField();
}


marangoniTdepFvPatchVectorField::marangoniTdepFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    transformFvPatchField<vector>(p, iF)
    //transformFvPatchVectorField(p, iF)
{
    //this->checkVolField();
    gamma_=readScalar(dict.lookup("gamma"));
    mu0_=readScalar(dict.lookup("mu0"));
    A_=readScalar(dict.lookup("A"));
    evaluate();
}


marangoniTdepFvPatchVectorField::marangoniTdepFvPatchVectorField
(
    const marangoniTdepFvPatchVectorField& ptf
)
:
    transformFvPatchField<vector>(ptf),
    //transformFvPatchVectorField(ptf),
    gamma_(ptf.gamma_),
    mu0_(ptf.mu0_),
    A_(ptf.A_)
{
    //this->checkVolField();
}


marangoniTdepFvPatchVectorField::marangoniTdepFvPatchVectorField
(
    const marangoniTdepFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    transformFvPatchField<vector>(ptf, iF),
    //transformFvPatchVectorField(ptf, iF),
    gamma_(ptf. gamma_),
    mu0_(ptf. mu0_),
    A_(ptf. A_)
{
    //this->checkVolField();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void marangoniTdepFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    vectorField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void marangoniTdepFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    transformFvPatchField<vector>::rmap(ptf, addr);
    //transformFvPatchVectorField::rmap(ptf, addr);

    const marangoniTdepFvPatchVectorField& dmptf =
        refCast<const marangoniTdepFvPatchVectorField >(ptf);

    gamma_=dmptf.gamma_;
    mu0_=dmptf.mu0_;
    A_=dmptf.A_;
}


// Return gradient at boundary
tmp<vectorField > marangoniTdepFvPatchVectorField::snGrad() const
{
    //Info << "entering  marangoniTdepFvPatchVectorField::snGrad()" << endl;
    vectorField nHat = this->patch().nf();
    vectorField pif = this->patchInternalField();
    vectorField result;
    if(!db().foundObject<vectorField>("gradT")) {
 	Info << " marangoniTdepFvPatchVectorField::snGrad(): object gradT not found!" << endl;
 	Foam::error theError(" marangoniTdepFvPatchVectorField::evaluate(): object gradT not found!");
	return transform(I - sqr(nHat),pif);
    }
    else if(!db().foundObject<scalarField>("T")) {
 	Info << " marangoniTdepFvPatchVectorField::snGrad(): object T not found!" << endl;
 	Foam::error theError(" marangoniTdepFvPatchVectorField::evaluate(): object T not found!");
	return transform(I - sqr(nHat),pif);
    }
    else if(!db().foundObject<scalarField>("mu")) {
 	Info << " marangoniTdepFvPatchVectorField::snGrad(): object mu not found!" << endl;
 	Foam::error theError(" marangoniTdepFvPatchVectorField::evaluate(): object mu not found!");
	return transform(I - sqr(nHat),pif);
    }
    else {
	
 	fvPatchField<vector> gradT =this->patch().lookupPatchField<volVectorField, vector>("gradT");
 	fvPatchField<scalar> T =this->patch().lookupPatchField<volScalarField, scalar>("T");
 	vectorField  gradT_internal = gradT.patchInternalField();
 	vectorField gradTplane= transform(I - sqr(nHat),gradT_internal);
 	vectorField pifplane= transform(I - sqr(nHat),pif);
        scalarField T_internal = T.patchInternalField();

        // Get mu from the solver patch field
 	fvPatchField<scalar> mu =this->patch().lookupPatchField<volScalarField, scalar>("mu");
        scalarField mu_internal = mu.patchInternalField();
        // Compute marangoni coeffient using value from patch field
        scalarField marangonicoeff = gamma_/mu_internal;

        /*
        // Get mu by computing it here directly and imposing the same constraint
        // as in the solver (messier)
        scalarField Tvisc = Foam::max(1689.15, T_internal);
        scalarField marangonicoeff = gamma_/(mu0_*exp(A_/Tvisc));
        */

        // Use marangonicoeff as boundary condition
	result=pifplane+marangonicoeff*gradTplane/this->patch().deltaCoeffs();
	return (result-pif)*this->patch().deltaCoeffs(); // Original result from downloaded code

        // BSM: I believe the returned result is a fancy way of setting normal
        // velocity to zero on boundary. If normal velocity is zero and cell
        // normal velocity is un, then gradient in normal direction is (0-un)/d,
        // which is basically what the (result-pif) part does
        //
        // The rest of the gradient is just the marangoniceoff*gradTplane
        //
        // I wonder if this way of setting 0 velocity on boundary, which is
        // different than fixing cell velocity to be zero, could be causing
        // problems?
    }
    //Info << "leaving  marangoniTdepFvPatchVectorField::snGrad()" << endl;    
    return result;
}


// Evaluate the field on the patch
void marangoniTdepFvPatchVectorField::evaluate()
{
    //Info << "entering  marangoniTdepFvPatchVectorField::evaluate()" << endl;
    if (!this->updated())
    {
	//Info << "marangoniTdepFvPatchVectorField::evaluate(): calling updatecoeffs" << endl;
        this->updateCoeffs();
    }

     vectorField nHat = this->patch().nf();
     vectorField pif = this->patchInternalField();
     scalarField deltas=this->patch().deltaCoeffs();

     if(!db().foundObject<vectorField>("gradT")) {
 	Info << " marangoniTdepFvPatchVectorField::snGrad(): object gradT not found!" << endl;
 	Foam::error theError(" marangoniTdepFvPatchVectorField::evaluate(): object gradT not found!");
 	vectorField::operator=
 	    (
 	      	transform(I - sqr(nHat),pif)
 	     );
 	// theError.exit();
     }
     else if(!db().foundObject<scalarField>("T")) {
 	Info << " marangoniTdepFvPatchVectorField::snGrad(): object T not found!" << endl;
 	Foam::error theError(" marangoniTdepFvPatchVectorField::evaluate(): object T not found!");
 	vectorField::operator=
 	    (
 	      	transform(I - sqr(nHat),pif)
 	     );
 	// theError.exit();
     }
     else if(!db().foundObject<scalarField>("mu")) {
 	Info << " marangoniTdepFvPatchVectorField::snGrad(): object mu not found!" << endl;
 	Foam::error theError(" marangoniTdepFvPatchVectorField::evaluate(): object mu not found!");
 	vectorField::operator=
 	    (
 	      	transform(I - sqr(nHat),pif)
 	     );
 	// theError.exit();
     }
     else {

 	fvPatchField<vector> gradT =this->patch().lookupPatchField<volVectorField, vector>("gradT");
 	fvPatchField<scalar> T =this->patch().lookupPatchField<volScalarField, scalar>("T");
 	vectorField  gradT_internal = gradT.patchInternalField();
 	vectorField gradTplane= transform(I - sqr(nHat),gradT_internal);
 	vectorField pifplane= transform(I - sqr(nHat),pif);
        scalarField T_internal = T.patchInternalField();

        // Get mu from the solver patch field
 	fvPatchField<scalar> mu =this->patch().lookupPatchField<volScalarField, scalar>("mu");
        scalarField mu_internal = mu.patchInternalField();
        // Compute marangoni coeffient using value from patch field
        scalarField marangonicoeff = gamma_/mu_internal;

        /*
        // Get mu by computing it here directly and imposing the same constraint
        // as in the solver (messier)
        scalarField Tvisc = Foam::max(1689.15, T_internal);
        scalarField marangonicoeff = gamma_/(mu0_*exp(A_/Tvisc));
        */

        // Use coefficient in final result
	vectorField result=pifplane+marangonicoeff*gradTplane/deltas;

 	vectorField::operator=
 	    (
	     result
 	     );
     }

    transformFvPatchField<vector>::evaluate();
    //transformFvPatchVectorField::evaluate();
}


// Return defining fields
tmp<vectorField > marangoniTdepFvPatchVectorField::snGradTransformDiag() const
{
    vectorField nHat = this->patch().nf();
    vectorField diag(nHat.size());

    diag.replace(vector::X, mag(nHat.component(vector::X)));
    diag.replace(vector::Y, mag(nHat.component(vector::Y)));
    diag.replace(vector::Z, mag(nHat.component(vector::Z)));

    return transformFieldMask<vector>(pow<vector, pTraits<vector>::rank>(diag));
}


// Write
void marangoniTdepFvPatchVectorField::write(Ostream& os) const
{
    transformFvPatchField<vector>::write(os);
    //transformFvPatchVectorField::write(os);
    os.writeKeyword("gamma") << gamma_ << token::END_STATEMENT << nl;
    os.writeKeyword("mu0") << mu0_ << token::END_STATEMENT << nl;
    os.writeKeyword("A") << A_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, marangoniTdepFvPatchVectorField);


} // End namespace Foam

// ************************************************************************* //


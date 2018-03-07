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

Class
    marangoniTdepFvPatchVectorField

Description

SourceFiles
    marangoniTdepFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef marangoniTdepFvPatchVectorField_H
#define marangoniTdepFvPatchVectorField_H
#include "fvPatchFields.H"
#include "transformFvPatchFields.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class marangoniTdepFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class marangoniTdepFvPatchVectorField
:
    public transformFvPatchField<vector>
    //public transformFvPatchVectorField
{
    // Private data

    // proportionality factor btween thermal and momntum gradient
    scalar gamma_;
    // Constants for Arrhenius temperature-dependent viscosity
    scalar mu0_;
    scalar A_;

public:

    //- Runtime type information
    TypeName("marangoniTdep");

    // Constructors

        //- Construct from patch and internal field
        marangoniTdepFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        marangoniTdepFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given marangoniTdepFvPatchVectorField
        //  onto a new patch
        marangoniTdepFvPatchVectorField
        (
            const marangoniTdepFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        marangoniTdepFvPatchVectorField
        (
            const marangoniTdepFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField > clone() const
        {
            return tmp<fvPatchVectorField >
            (
                new marangoniTdepFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        marangoniTdepFvPatchVectorField
        (
            const marangoniTdepFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField > clone(const DimensionedField<vector, volMesh>& iF) const
        {
            return tmp<fvPatchVectorField >
            (
                new marangoniTdepFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            );


    // Return defining fields

    // virtual vectorField& fixedValue()
    //       {
    //          return fixedValue_;
    //      }

    //      virtual const vectorField& fixedValue() const
    //      {
    //          return fixedValue_;
    //      }


        // Evaluation functions

            //- Return gradient at boundary
            virtual tmp<vectorField > snGrad() const;

            //- Evaluate the patch field
            virtual void evaluate();

            //- Return face-gradient transform diagonal
            virtual tmp<vectorField > snGradTransformDiag() const;


        //- Write
        virtual void write(Ostream&) const;


    // Member operators

    virtual void operator=(const UList<vector>&) {}

    virtual void operator=(const fvPatchVectorField&) {}
      virtual void operator+=(const fvPatchVectorField&) {}
      virtual void operator-=(const fvPatchVectorField&) {}
      virtual void operator*=(const fvPatchField<scalar>&) {}
      virtual void operator/=(const fvPatchField<scalar>&) {}
      virtual void operator+=(const vectorField&) {}
      virtual void operator-=(const vectorField&) {}

      virtual void operator*=(const Field<scalar>&) {}
      virtual void operator/=(const Field<scalar>&) {}

            virtual void operator=(const vector&) {}
      virtual void operator+=(const vector&) {}
      virtual void operator-=(const vector&) {}
      virtual void operator*=(const scalar) {}
      virtual void operator/=(const scalar) {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//#   include "marangoniTdepFvPatchVectorField.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

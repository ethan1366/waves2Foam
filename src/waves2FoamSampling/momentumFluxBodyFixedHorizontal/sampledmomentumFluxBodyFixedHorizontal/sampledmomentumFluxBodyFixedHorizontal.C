/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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

#include "sampledmomentumFluxBodyFixedHorizontal.H"
#include "dictionary.H"
//#if EXTBRANCH==1 && OFVERSION>310
//    #include "foamTime.H"
//#else
    #include "Time.H"
//#endif
#include "volFields.H"
#include "ListListOps.H"
#include "SortableList.H"
#include "volPointInterpolation.H"
#include "mixedFvPatchField.H"
#include "convexPolyhedral.H"
#include "waveTheory.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledmomentumFluxBodyFixedHorizontal, 0);
}

bool Foam::sampledmomentumFluxBodyFixedHorizontal::verbose_ = false;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


bool Foam::sampledmomentumFluxBodyFixedHorizontal::checkFieldTypes()
{
    wordList fieldTypes(fieldNames_.size());

    // check files for a particular time
    if (loadFromFiles_)
    {
        forAll (fieldNames_, fieldi)
        {
            IOobject io
            (
                fieldNames_[fieldi],
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            );

#if OFPLUSBRANCH == 1
    #if OFVERSION<1606
            if (io.headerOk())
            {
                fieldTypes[fieldi] = io.headerClassName();
            }
            else
            {
                fieldTypes[fieldi] = "(notFound)";
            }
    #else
            // Circumventing the check and only check for scalarField, since 
            // we know that we do not need to check for anything else.
            if (io.typeHeaderOk<volScalarField>(true))
            {
                fieldTypes[fieldi] = io.headerClassName();
            }
            // below add line 
            else
            {
               if (io.typeHeaderOk<volVectorField>(true)) 
               {
                  fieldTypes[fieldi] = io.headerClassName();
               }
            
               else
              {
                fieldTypes[fieldi] = "(notFound)";
              }

            }
           // above add line

           /* else
            {
                fieldTypes[fieldi] = "(notFound)";
            }*/
            
    #endif
#else
            if (io.headerOk())
            {
                fieldTypes[fieldi] = io.headerClassName();
            }
            else
            {
                fieldTypes[fieldi] = "(notFound)";
            }
#endif
        }
    }
    else
    {
        // check objectRegistry
        forAll (fieldNames_, fieldi)
        {
            objectRegistry::const_iterator iter =
                mesh_.find(fieldNames_[fieldi]);

            if (iter != mesh_.objectRegistry::end())
            {
                fieldTypes[fieldi] = iter()->type();
            }
            else
            {
                fieldTypes[fieldi] = "(notFound)";
            }
        }
    }


    label nFields = 0;

    // classify fieldTypes
    nFields += grep(scalarFields_, fieldTypes); // grep means search file; when two fields exist, i.e. alpha.water and U, this line search alpha.water,give this value to scalarFields_;
    nFields += grep(vectorFields_, fieldTypes);
//    nFields += grep(sphericalTensorFields_, fieldTypes);
//    nFields += grep(symmTensorFields_, fieldTypes);
//    nFields += grep(tensorFields_, fieldTypes);

    if (Pstream::master())
    {
        if (debug)
        {
            Pout<< "timeName = " << mesh_.time().timeName() << nl
                << "scalarFields    " << scalarFields_ << nl
                << "vectorFields    " << vectorFields_ << nl;
//                << "sphTensorFields " << sphericalTensorFields_ << nl
//                << "symTensorFields " << symmTensorFields_ <<nl
//                << "tensorFields    " << tensorFields_ <<nl;
        }

//        if (nFields > 0)
//        {
//            if (debug)
//            {
//                Pout<< "Creating directory "
//                    << outputPath_/mesh_.time().timeName()
//                    << nl << endl;
//            }
//
//            mkDir(outputPath_/mesh_.time().timeName());
//        }
    }

    return nFields > 0;
}


void Foam::sampledmomentumFluxBodyFixedHorizontal::combineSampledSets
(
    PtrList<coordSet>& masterSampledSets,
    labelListList& indexSets
)
{
    // Combine sampleSets from processors. Sort by curveDist. Return
    // ordering in indexSets.
    // Note: only master results are valid

    masterSampledSets_.clear();
    masterSampledSets_.setSize(size());
    indexSets_.setSize(size());

    const PtrList<sampledSet>& sampledSets = *this;

    forAll (sampledSets, seti)
    {
        const sampledSet& samplePts = sampledSets[seti];

        // Collect data from all processors
        List<List<point> > gatheredPts(Pstream::nProcs());
        gatheredPts[Pstream::myProcNo()] = samplePts;
        Pstream::gatherList(gatheredPts);

        List<labelList> gatheredSegments(Pstream::nProcs());
        gatheredSegments[Pstream::myProcNo()] = samplePts.segments();
        Pstream::gatherList(gatheredSegments);

        List<scalarList> gatheredDist(Pstream::nProcs());
        gatheredDist[Pstream::myProcNo()] = samplePts.curveDist();
        Pstream::gatherList(gatheredDist);


        // Combine processor lists into one big list.
        List<point> allPts
        (
            ListListOps::combine<List<point> >
            (
                gatheredPts, accessOp<List<point> >()
            )
        );
        labelList allSegments
        (
            ListListOps::combine<labelList>
            (
                gatheredSegments, accessOp<labelList>()
            )
        );
        scalarList allCurveDist
        (
            ListListOps::combine<scalarList>
            (
                gatheredDist, accessOp<scalarList>()
            )
        );

        // Sort curveDist and use to fill masterSamplePts
        SortableList<scalar> sortedDist(allCurveDist);
        indexSets[seti] = sortedDist.indices();


        // The constructor for coordSet has changed as of version 2.0.
        // This is taken care of using these pre-processor statements.
#if OFVERSION < 200 || EXTBRANCH==1
        // Get reference point (note: only master has all points)
        point refPt;

        if (allPts.size())
        {
            refPt = samplePts.getRefPoint(allPts);
        }
        else
        {
            refPt = vector::zero;
        }

        masterSampledSets.set
        (
            seti,
            new coordSet
            (
                samplePts.name(),
                samplePts.axis(),
                UIndirectList<point>(allPts, indexSets[seti]),
                refPt
            )
        );
#else
        masterSampledSets.set
        (
            seti,
            new coordSet
            (
                samplePts.name(),
                samplePts.axis(),
                List<point>(UIndirectList<point>(allPts, indexSets[seti])),
                allCurveDist
            )
        );
#endif
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledmomentumFluxBodyFixedHorizontal::sampledmomentumFluxBodyFixedHorizontal
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    PtrList<sampledSet>(),
    name_(name),
    mesh_(refCast<const fvMesh>(obr)),
    loadFromFiles_(loadFromFiles),
    outputPath_(fileName::null),
#if OFVERSION < 210
    searchEngine_(mesh_, true),
#else
    searchEngine_(mesh_),
#endif
    fieldNames_(),
    interpolationScheme_(word::null),
    writeFormat_(word::null),
    momentumFluxBodyFixedHorizontalFilePtr_( NULL )
{
    startTime_ = dict.lookupOrDefault<scalar>("samplingStartTime", 0.0);
    nextSampleTime_ = startTime_;
    momentumFluxBodyFixedHorizontalSampleDeltaT_ =
        dict.lookupOrDefault<scalar>("momentumFluxBodyFixedSampleDeltaT", -1);
   // pitchAngle_ = dict.lookupOrDefault<scalar>("pitchAngle", 0.0);

    if (Pstream::parRun())
    {
        outputPath_ = mesh_.time().path()/".."/name_;
    }
    else
    {
        outputPath_ = mesh_.time().path()/name_;
    }
    if (mesh_.name() != fvMesh::defaultRegion)
    {
        outputPath_ = outputPath_/mesh_.name();
    }

//    mkDir(outputPath_ + "/" + mesh_.time().timeName() );

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::sampledmomentumFluxBodyFixedHorizontal::~sampledmomentumFluxBodyFixedHorizontal()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::sampledmomentumFluxBodyFixedHorizontal::verbose(const bool verbosity)
{
    verbose_ = verbosity;
}


void Foam::sampledmomentumFluxBodyFixedHorizontal::execute()
{
    // Do nothing - only valid on write
}


void Foam::sampledmomentumFluxBodyFixedHorizontal::end()
{
    // Do nothing - only valid on write
}


void Foam::sampledmomentumFluxBodyFixedHorizontal::write()
{
    if (size() && checkFieldTypes())
    {
        sampleIntegrateAndWrite(scalarFields_,vectorFields_);
    }
}

bool Foam::sampledmomentumFluxBodyFixedHorizontal::performAction()
{
    if (momentumFluxBodyFixedHorizontalSampleDeltaT_ <= 10*SMALL)
    {
        // This line is needed to update the locations of the interpolation
        // lines.
    	// Note, that performAction() is only called in case the upstream
    	// time controls passes, i.e. a given timeIndex or at outputTime().
        if (mesh_.moving())
        {
    	    this->correct();
        }

        return mesh_.time().value() >= startTime_;
    }
    else
    {
        if (mesh_.time().value() < nextSampleTime_)
        {
            return false;
        }
        else
        {
            while (mesh_.time().value() > nextSampleTime_)
            {
                nextSampleTime_ += momentumFluxBodyFixedHorizontalSampleDeltaT_;
            }

            // This line is needed to update the locations of the interpolation
            // lines.
            if (mesh_.moving())
            {
        	    this->correct();
            }

            return true;
        }
    }
}


void Foam::sampledmomentumFluxBodyFixedHorizontal::sampleIntegrateAndWrite
(
    fieldGroup<scalar>& fieldsScalar,
    fieldGroup<vector>& fieldsVector
)
{
    if (fieldsScalar.size() && fieldsVector.size() && performAction())
    {
        scalarField result(0);
        sampleAndIntegrate(scalarFields_, vectorFields_, result);

        if (Pstream::master())
        {
            // create file if not already there, notice: this shall be
            // done on master node only
            if (momentumFluxBodyFixedHorizontalFilePtr_.empty())
            {
                mkDir( outputPath_ + "/" + mesh_.time().timeName() );
                momentumFluxBodyFixedHorizontalFilePtr_.reset
                (
                    new OFstream
                    (
                        outputPath_ + "/" + mesh_.time().timeName()
                      + "/momentumFluxBodyFixedHorizontal.dat"
                    )
                );

                // write header
                if (momentumFluxBodyFixedHorizontalFilePtr_.valid())
                {
                    momentumFluxBodyFixedHorizontalFilePtr_() << "Time";

                    forAll (masterSampledSets_, seti)
                    {
                        momentumFluxBodyFixedHorizontalFilePtr_() << tab
                            << masterSampledSets_[seti].name();
                    }
                    momentumFluxBodyFixedHorizontalFilePtr_() << endl;

                    for (int coordi = 0; coordi < 3; coordi++)
                    {
                        momentumFluxBodyFixedHorizontalFilePtr_() << -1 - coordi;

                        forAll (masterSampledSets_, seti)
                        {
                            momentumFluxBodyFixedHorizontalFilePtr_() << tab
                                 << masterSampledSets_[seti][0].component(coordi);
                        }
                        momentumFluxBodyFixedHorizontalFilePtr_() << endl;
                    }
                }
                else
                {
                    FatalErrorIn
                    (
                       "void Foam::sampledmomentumFlux::sampleIntegrateAndWrite( ... )"
                    )
                    << "Output file could not be opened in " << outputPath_
                    << "/" << mesh_.time().timeName() << endl << endl
                    << exit(FatalError);
                }
            }

            if (momentumFluxBodyFixedHorizontalFilePtr_.valid())
            {
                momentumFluxBodyFixedHorizontalFilePtr_() << mesh_.time().value();

                forAll (result, seti)
                {
                    momentumFluxBodyFixedHorizontalFilePtr_() << tab << result[seti];
                }

                momentumFluxBodyFixedHorizontalFilePtr_() << endl;
            }
        }
    }
}


void Foam::sampledmomentumFluxBodyFixedHorizontal::sampleAndIntegrate
(
    fieldGroup<scalar>& fieldsScalar,
    fieldGroup<vector>& fieldsVector,
    Field<scalar>& result
)
{
    result.setSize(0);

    if (fieldsScalar.size() && fieldsVector.size())
    {
        bool interpolate = interpolationScheme_ != "cell";

        // Create or use existing writer
        if (fieldsScalar.formatter.empty())
        {
            fieldsScalar.formatter = writer<scalar>::New(writeFormat_);
        }
        if (fieldsVector.formatter.empty())
        {
            fieldsVector.formatter = writer<vector>::New(writeFormat_);
        }        

        // Storage for interpolated values
        PtrList<volFieldSampler<scalar> > sampledFieldsScalar(fieldsScalar.size());
        PtrList<volFieldSampler<vector> > sampledFieldsVector(fieldsVector.size());
          // for scalar field
        forAll (fieldsScalar, fieldi)
        {
            if (Pstream::master() && verbose_)
            {
                Pout<< "sampledSets::sampleAndWrite: "
                    << fieldsScalar[fieldi] << endl;
            }

            if (loadFromFiles_)
            {
                GeometricField<scalar, fvPatchField, volMesh> vf
                (
                    IOobject
                    (
                        fieldsScalar[fieldi],
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    mesh_
                );

                if (interpolate)
                {
                    sampledFieldsScalar.set
                    (
                        fieldi,
                        new volFieldSampler<scalar>
                        (
                            interpolationScheme_,
                            vf,
                            *this
                        )
                    );
                }
                else
                {
                    sampledFieldsScalar.set
                    (
                        fieldi,
                        new volFieldSampler<scalar>(vf, *this)
                    );
                }
            }
            else
            {
                if (interpolate)
                {
                    sampledFieldsScalar.set
                    (
                        fieldi,
                        new volFieldSampler<scalar>
                        (
                            interpolationScheme_,
                            mesh_.lookupObject
                            <GeometricField<scalar, fvPatchField, volMesh> >
                            (fieldsScalar[fieldi]),
                            *this
                        )
                    );
                }
                else
                {
                    sampledFieldsScalar.set
                    (
                        fieldi,
                        new volFieldSampler<scalar>
                        (
                            mesh_.lookupObject
                            <GeometricField<scalar, fvPatchField, volMesh> >
                            (fieldsScalar[fieldi]),
                            *this
                        )
                    );
                }
            }
        }
            
              // for vector field
        forAll (fieldsVector, fieldi)
        {
            if (Pstream::master() && verbose_)
            {
                Pout<< "sampledSets::sampleAndWrite: "
                    << fieldsVector[fieldi] << endl;
            }

            if (loadFromFiles_)
            {
                GeometricField<vector, fvPatchField, volMesh> vf
                (
                    IOobject
                    (
                        fieldsVector[fieldi],
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    mesh_
                );

                if (interpolate)
                {
                    sampledFieldsVector.set
                    (
                        fieldi,
                        new volFieldSampler<vector>
                        (
                            interpolationScheme_,
                            vf,
                            *this
                        )
                    );
                }
                else
                {
                    sampledFieldsVector.set
                    (
                        fieldi,
                        new volFieldSampler<vector>(vf, *this)
                    );
                }
            }
            else
            {
                if (interpolate)
                {
                    sampledFieldsVector.set
                    (
                        fieldi,
                        new volFieldSampler<vector>
                        (
                            interpolationScheme_,
                            mesh_.lookupObject
                            <GeometricField<vector, fvPatchField, volMesh> >
                            (fieldsVector[fieldi]),
                            *this
                        )
                    );
                }
                else
                {
                    sampledFieldsVector.set
                    (
                        fieldi,
                        new volFieldSampler<vector>
                        (
                            mesh_.lookupObject
                            <GeometricField<vector, fvPatchField, volMesh> >
                            (fieldsVector[fieldi]),
                            *this
                        )
                    );
                }
            }
        }





        // Combine sampled fields from processors.
        // Note: only master results are valid
        PtrList<volFieldSampler<scalar> > masterFieldsScalar(sampledFieldsScalar.size());
        combineSampledValues(sampledFieldsScalar, indexSets_, masterFieldsScalar);

        PtrList<volFieldSampler<vector> > masterFieldsVector(sampledFieldsVector.size());
        combineSampledValues(sampledFieldsVector, indexSets_, masterFieldsVector);

        result.setSize(masterSampledSets_.size(), 0.0);

        if (Pstream::master())
        {
            forAll (masterSampledSets_, seti)
            {
                const coordSet & cs( masterSampledSets_[seti] );
                // scalar fields
                List< const Field<scalar>*> valueSetsScalar(masterFieldsScalar.size());
                valueSetsScalar[0] = &masterFieldsScalar[0][seti];  // & take the address of a variable

                List<const Field<scalar>*> columnsScalar(valueSetsScalar.size());
                columnsScalar[0] = valueSetsScalar[0];

                const Field<scalar>& alpha = *columnsScalar[0];  // * take the value corresponding to an address
                // vector fields

                List< const Field<vector>*> valueSetsVector(masterFieldsVector.size());
                valueSetsVector[0] = &masterFieldsVector[0][seti];

                List<const Field<vector>*> columnsVector(valueSetsVector.size());
                columnsVector[0] = valueSetsVector[0];

                const Field<vector>& U = *columnsVector[0];



                scalar tolerance(0.0001);

                if (alpha.size() < 2)
                {
                    result[seti] = -GREAT;
                }
                else if
                (
                    (alpha[0] < tolerance && alpha[alpha.size()-1] < tolerance)
                    ||
                    (
                        alpha[0] > 1.0 - tolerance &&
                        alpha[alpha.size()-1] > 1.0 - tolerance
                    )
                )
                {
                    result[seti] = -GREAT;
                }
                else
                {
                    scalar value(0);
 //                   scalar signi(1);
 //                   scalar signiplus1(1);
 //                   scalar minScalarCoord(cs.scalarCoord(0));
                    
                   vector coordVector0=cs.vectorCoord (0);    // a set of probe, the first point, component of a class using [] while component of class member using ()
                   vector coordVector1=cs.vectorCoord (1);    // a set of probe, the second point
                   scalar pitchAngle=Foam::atan2 (coordVector1.y()-coordVector0.y(), coordVector1.x()-coordVector0.x())-constant::mathematical::pi/2; // calculate the pitch angle , here the coordinate system is x right, y upward and z pointed outside the plane; atan2(x1,x2), x1 faces to the angle while x2 is adjacent to the angle.

                    for (int pointi=0; pointi < alpha.size() - 1; pointi++)
                    {
 /*                       if (U[pointi + 1].component(0) >= 0 )
                        {
                              signiplus1=1;
                        }
                        else 
                        {
                              signiplus1=-1;
                        }
                        if (U[pointi].component(0) >= 0 )
                        {
                              signi=1;
                        }
                        else 
                        {
                              signi=-1;
                        }    */                    
                        value +=
                            (
                                cs.scalarCoord(pointi + 1)            // note in probeDefinitions, axis is set as y, so here cs.scalarCoord refer to y coordinate.
                              - cs.scalarCoord(pointi)
                             ) 
                             *( alpha[pointi + 1]*(U[pointi + 1].component(0)*cos(pitchAngle)+U[pointi + 1].component(1)*sin(pitchAngle))
                                 * (U[pointi + 1].component(0)*cos(pitchAngle)+U[pointi + 1].component(1)*sin(pitchAngle))
                                 + alpha[pointi]*(U[pointi].component(0)*cos(pitchAngle)+U[pointi].component(1)*sin(pitchAngle))
                                 * (U[pointi].component(0)*cos(pitchAngle)+U[pointi].component(1)*sin(pitchAngle))
                              );

  /*                      minScalarCoord =
                            Foam::min
                            (
                                minScalarCoord,
                                cs.scalarCoord(pointi + 1)
                            );*/
                    }
                    value= value*0.5/cos(pitchAngle);     // because previous calculation is deltay, here divided by cos(pitchAngle), it is distance between two points.

                    result[seti] = value;
                }
            }
        }
    }
}


void Foam::sampledmomentumFluxBodyFixedHorizontal::read(const dictionary& dict)
{
    dict_ = dict;

    fieldNames_ = wordList(dict_.lookup("fields"));

    interpolationScheme_ = "cell";
    dict_.readIfPresent("interpolationScheme", interpolationScheme_);

    dict_.lookup("setFormat") >> writeFormat_;

    //pitchAngle_ = dict.lookupOrDefault<scalar>("pitchAngle", 0.0);

    scalarFields_.clear();

    PtrList<sampledSet> newList
    (
        dict_.lookup("sets"),
        sampledSet::iNew(mesh_, searchEngine_)
    );
    transfer(newList);
    combineSampledSets(masterSampledSets_, indexSets_);

    if (Pstream::master() && debug)
    {
        Pout<< "sample fields:" << fieldNames_ << nl
            << "sample sets:" << nl << "(" << nl;

        forAll (*this, si)
        {
            Pout << "  " << operator[](si) << endl;
        }
        Pout << ")" << endl;
    }
}


void Foam::sampledmomentumFluxBodyFixedHorizontal::correct()
{
    // Reset interpolation
	// These two lines make the moving mesh algorithms crash
	// (tested: velocityLaplacian)
	// NGJ: 16/03/2015
//    pointMesh::Delete(mesh_);
//    volPointInterpolation::Delete(mesh_);

    searchEngine_.correct();

    // A quick test has shown that this takes a lot of time on moving meshes
    // Potentially return to improve - if possible.
    // NGJ: 16/03/2015.
    PtrList<sampledSet> newList
    (
        dict_.lookup("sets"),
        sampledSet::iNew(mesh_, searchEngine_)
    );

    transfer(newList);

    combineSampledSets(masterSampledSets_, indexSets_);
}


void Foam::sampledmomentumFluxBodyFixedHorizontal::updateMesh(const mapPolyMesh&)
{
    correct();
}


#if OFVERSION<220 || EXTBRANCH==1
void Foam::sampledmomentumFluxBodyFixedHorizontal::movePoints(const pointField&)
{
    correct();
}
#else
void Foam::sampledmomentumFluxBodyFixedHorizontal::movePoints(const polyMesh&)
{
    correct();
}

#if OFVERSION > 220 && EXTBRANCH==0
    bool Foam::sampledmomentumFluxBodyFixedHorizontal::timeSet()
    {
        // Do nothing
        return true;
    }
#elif XVERSION && EXTBRANCH==0
    bool Foam::sampledmomentumFluxBodyFixedHorizontal::timeSet()
    {
        // Do nothing
        return true;
    }
#endif

#endif


void Foam::sampledmomentumFluxBodyFixedHorizontal::readUpdate
(
    const polyMesh::readUpdateState state
)
{
    if (state != polyMesh::UNCHANGED)
    {
        correct();
    }
}


// ************************************************************************* //

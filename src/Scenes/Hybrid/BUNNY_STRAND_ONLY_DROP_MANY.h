#ifndef BUNNY_STRAND_ONLY_DROP_MANY_H
#define BUNNY_STRAND_ONLY_DROP_MANY_H

#include "SIMULATION_SCENE.h"
#include "Timestepper/Volume/BDF_1.h"
#include "Timestepper/Volume/BDF_2.h"
#include "Timestepper/Hybrid/TIMESTEPPER.h"
#include "Timestepper/Hybrid/VOLUME_BDF_1.h"
#include "Geometry/STRAND_MESH_FASTER.h"
#include "Geometry/STRAND_MESH.h"
#include "Geometry/BOWL.h"
#include "Geometry/TRIANGLE_MESH.h"
#include "Geometry/TRIANGLE_MESH_FASTER.h"
#include "Hyperelastic/Strand/QUADRATIC_STRETCHING.h"
#include "Hyperelastic/Strand/BERGOU_2010.h"
#include "Hyperelastic/Strand/ISOTROPIC_BENDING.h"
#include "Hyperelastic/Strand/QUADRATIC_BENDING.h"
#include "Hyperelastic/Strand/QUADRATIC_UNIT_BENDING.h"
#include "Hyperelastic/Strand/ISOTROPIC_THETA.h"
#include "Hyperelastic/Shell/STVK.h"
#include "Hyperelastic/Shell/ARAP.h"
#include "Hyperelastic/Shell/BW_STRETCH.h"
#include "Hyperelastic/Shell/BW_SHEAR.h"
#include "Hyperelastic/Shell/BARAFF_WITKIN.h"
#include "Hyperelastic/Shell/BENDING_SPRING.h"
#include "Hyperelastic/Shell/DIHEDRAL.h"

namespace HOBAK {

class BUNNY_STRAND_ONLY_DROP_MANY : public SIMULATION_SCENE {
public:

  virtual void printSceneDescription() override
  {
    cout << "=====================================================================" << endl;
    cout << " Dropping a bunny down an obstacle course to test out both kinematic" << endl;
    cout << " and self-collisions. Both VF and EE collisions are enabled." << endl;
    cout << "=====================================================================" << endl;
  }

  virtual bool buildScene() override
  {
    // this will determine the MOV and JSON filenames
    _sceneName = "bunny_strand_only_drop_many";
    _writeToFile = true;

    // read in the tet mesh file
    //_tetMeshFilename = string("../data/scorpion_0_125.tobj");
    _tetMeshFilename = string("../data/bunny/bunny_5.tobj");
    _strandMeshFilename = "../data/singleCurl.sobj";
    // _triangleMeshFilename = string("../data/objs/ribbon_4x44.obj");
    _triangleMeshFilename = string("../data/objs/curtain_10x10.obj");
    //_tetMeshFilename = string("../data/Brachiosaurus_15.tobj");
    //_tetMeshFilename = string("../data/dragon.tobj");

    if(_writeToFile){
      if(!fs::exists("../data/renders"))
        fs::create_directory("../data/renders");
      _outDir = string("../data/renders/") + _sceneName + string("/");
      if(!fs::exists(_outDir))
        fs::create_directory(_outDir);
      if(!fs::exists(_outDir + string("volume")))
        fs::create_directory(_outDir + string("volume"));
      if(!fs::exists(_outDir + string("strand")))
        fs::create_directory(_outDir + string("strand"));
      if(!fs::exists(_outDir + string("shell")))
        fs::create_directory(_outDir + string("shell"));
      // if(!fs::exists("../data/renders/shell"))
      //   fs::create_directory("../data/renders/shell");

      _volumeFilePath = _outDir + string("volume/");
      _strandFilePath = _outDir + string("strand/");
      _shellFilePath = _outDir + string("shell/");
    }

    // read tet
    vector<VECTOR3> vertices;
    vector<VECTOR4I> tets;
    bool success = TET_MESH::readTobjFile(_tetMeshFilename, vertices, tets);
    if (!success)
    {
      cout << " Failed to open file " << _tetMeshFilename.c_str() << endl;
      return false;
    }
    vertices = TET_MESH::normalizeVertices(vertices);
    _normalizedVertices = true;


    // strand settings
    using namespace Eigen;
    // _E = 1e8;
    // _G = 1e8;
    // _spacing = 1.0;
    // _density = 1.5;
    // _baseRadius = 0.005;
    // _tipRadius = 0.005;
    // _drawVertices = true;
    // _totalStrands = 1;

    _E = 3.9e6;
    _G = 3.9e6;
    _spacing = 1.0;
    _density = 1.3;
    _baseRadius = 0.1;
    _tipRadius = 0.1;
    _drawVertices = true;
    _totalStrands = 1;

    // build strand
    // read strand from file
    // vector<vector<int> > strandIndices;
    // std::vector<VECTOR3> restVerticesStrand;
    // if (!readSOBJFile(_strandMeshFilename, restVerticesStrand, strandIndices)) { 
    //   cout << " Failed to open file " << _strandMeshFilename << endl;
    //   std::exit(1); }
    // // _initialAStrand = 0.15 * MATRIX3::Identity() * AngleAxisd(0.5 * M_PI, VECTOR3::UnitZ());
    // _initialAStrand = MATRIX3::Identity() * AngleAxisd(0.5 * M_PI, VECTOR3::UnitZ());
    // // _initialAStrand = MATRIX3::Identity();
    // _initialTranslationStrand = VECTOR3(2.0, -4.0, 0.0);
    // // _initialTranslationStrand = VECTOR3(2.0, -1.0, 0.0);
    // for (unsigned int x = 0; x < restVerticesStrand.size(); x++)
    //   restVerticesStrand[x] = _initialAStrand * restVerticesStrand[x] + _initialTranslationStrand;
    // _strandMesh = new TET_WISP_MESH(restVerticesStrand, strandIndices,
    //                       _E, _G, _density, _baseRadius, _tipRadius);

    // construct hair tie
    // ---- all transformation --------
    // ----- strand ------
    _initialAStrand = MATRIX3::Identity();
    // _initialAStrand = MATRIX3::Identity();
    _initialTranslationStrand = VECTOR3(-4.5, 2.0, 0.0);
    // _initialTranslationStrand = VECTOR3(2.0, -1.0, 0.0);

    // ----- shell ------
    _initialAShell = 13 * MATRIX3::Identity()* AngleAxisd(0, VECTOR3::UnitX())
                                            // * AngleAxisd(0.5 * M_PI,  VECTOR3::UnitY())
                                            * AngleAxisd(0,  VECTOR3::UnitY())
                                            * AngleAxisd(0, VECTOR3::UnitZ());;
    // _initialTranslationShell = VECTOR3(6.0, -10.0, 0.0);
    _initialTranslationShell = VECTOR3(0.0, 8.2, 0.0);

    // ----- volume ------
    MATRIX3 M;
    // M =   AngleAxisd(-0.5 * M_PI, VECTOR3::UnitX())
    // //M =   AngleAxisd(0 , VECTOR3::UnitX())
    //     * AngleAxisd(0,  VECTOR3::UnitY())
    //     * AngleAxisd(0, VECTOR3::UnitZ());
    REAL bunnyScale = 5.2;
    M =  bunnyScale * MATRIX3::Identity() * AngleAxisd(-0.5 * M_PI, VECTOR3::UnitX())
    //M =   AngleAxisd(0 , VECTOR3::UnitX())
        * AngleAxisd(0,  VECTOR3::UnitY())
        * AngleAxisd(0, VECTOR3::UnitZ());
    VECTOR3 half(0.5, 0.5, 0.5);
   
    _initialA           = M;
    _initialTranslation = half - M * half + VECTOR3(-1.5,-5.3,1.5);


    buildHairTieMany();

    _strandMesh -> setCollisionEps(0.01);
    _strandMesh->bendingForceFilterEnabled() = true;

    const REAL k = 1.0;
    _stretchingEnergyStrand = new STRAND::QUADRATIC_STRETCHING(k);

    // read and build triangle mesh
    vector<VECTOR3> verticesShell;
    vector<VECTOR3I> trianglesShell;
    success = TRIANGLE_MESH::readObjFile(_triangleMeshFilename, verticesShell, trianglesShell);
    if (!success)
    {
      cout << " Failed to open file " << _triangleMeshFilename.c_str() << endl;
      return false;
    }
    // transformation
    // _initialAShell = 12 * MATRIX3::Identity()* AngleAxisd(0, VECTOR3::UnitX())
    //                                         // * AngleAxisd(0.5 * M_PI,  VECTOR3::UnitY())
    //                                         * AngleAxisd(0,  VECTOR3::UnitY())
    //                                         * AngleAxisd(0, VECTOR3::UnitZ());;
    // // _initialTranslationShell = VECTOR3(6.0, -10.0, 0.0);
    // _initialTranslationShell = VECTOR3(0.0, 0.2, 0.0);

    for (unsigned int x = 0; x < verticesShell.size(); x++)
      verticesShell[x] = _initialAShell * verticesShell[x] + _initialTranslationShell;
    // _triangleMesh = new TRIANGLE_MESH(vertices, triangles);
    _triangleMesh = new TRIANGLE_MESH_FASTER(verticesShell, trianglesShell);
    _triangleMesh -> setCollisionEps(0.1);
    //_strechingEnergy = new SHELL::STVK(10.0, 1.0);
    _strechingEnergyShell = new SHELL::ARAP(80.0, 0.0);
    //_strechingEnergy = new SHELL::BW_STRETCH(1.0, 0.0);
    //_strechingEnergy = new SHELL::BW_SHEAR(1.0, 0.0);
    //_strechingEnergy = new SHELL::BARAFF_WITKIN(1.0, 0.0);
    //_bendingEnergy = new SHELL::BENDING_SPRING(20);
    //_bendingEnergy = new SHELL::BENDING_SPRING(10);
    // _bendingEnergyShell = new SHELL::DIHEDRAL(0.01);
    _bendingEnergyShell = new SHELL::DIHEDRAL(5.0);
    // _bendingEnergy = new SHELL::DIHEDRAL(1);
    // _bendingEnergy = new SHELL::DIHEDRAL(10);
    //_bendingEnergy = new SHELL::DIHEDRAL(100);


    // build tet mesh
    // _initialTranslation = half - M * half;


    for (unsigned int x = 0; x < vertices.size(); x++)
      vertices[x] = _initialA * vertices[x] + _initialTranslation;

    // get more bunnies
    const int numOtherBunnies = 4;
    vector<VECTOR3> otherBunnyTranslations(numOtherBunnies);
    const int bunnyVertices = (int)vertices.size();
    const int bunnyTets = (int)tets.size();
    otherBunnyTranslations[0] = VECTOR3(-bunnyScale, -1.0, 0.0);
    otherBunnyTranslations[1] = VECTOR3(bunnyScale, -1.0, 0.0);
    otherBunnyTranslations[2] = VECTOR3(0.0, -1.0, -bunnyScale);
    otherBunnyTranslations[3] = VECTOR3(0.0, -1.0, bunnyScale);
    for(unsigned int y = 0; y < numOtherBunnies; y++) {
      for(unsigned int x = 0; x < bunnyVertices; x++) {
        vertices.push_back(vertices[x] + otherBunnyTranslations[y]);
      }
      for(unsigned int z = 0; z < bunnyTets; z++) {
        tets.push_back(tets[z] + VECTOR4I(bunnyVertices * (y+1), bunnyVertices * (y+1),
                                          bunnyVertices * (y+1), bunnyVertices * (y+1)));
      }
    }

    _gravity = VECTOR3(0, -1.0, 0);

    // REAL E = 3.0;
    // //REAL nu = 0.3; // lambda \approx 10
    // REAL nu = 0.45; // lambda \approx 10
    REAL E = 200.0;
    REAL nu = 0.45; 

    REAL mu     = VOLUME::HYPERELASTIC::computeMu(E, nu);
    REAL lambda = VOLUME::HYPERELASTIC::computeLambda(E, nu);
    cout << " mu:     " << mu << endl;
    cout << " lambda: " << lambda << endl;


    // build the tet mesh object
    _tetMesh = new TET_MESH_FASTER(vertices, tets);
    _hyperelastic = new VOLUME::SNH(mu, lambda);
    //_hyperelastic = new VOLUME::ARAP(mu * 4, lambda);

    const vector<REAL>& areas = _tetMesh->surfaceTriangleAreas();
    REAL smallest = areas[0];
    REAL largest = areas[0];
    for (unsigned int x = 1; x < areas.size(); x++)
    {
      if (areas[x] > largest) largest  = areas[x];
      if (areas[x] < largest) smallest = areas[x];
    }
    cout << "Largest triangle area:  "  << largest << endl;
    cout << "Smallest triangle area: "  << smallest << endl;

    // build the time integrator
    // _hybridSolver = new VOLUME::BACKWARD_EULER_POSITION(*_tetMesh, *_hyperelastic);
    cout<<"strand vertices: "<<_strandMesh->totalVertices()<<endl;
    cout<<"volume vertices: "<<_tetMesh->totalVertices()<<endl;
    cout<<"strand vertices: "<<_triangleMesh->totalVertices()<<endl;

    _hybridSolver = new HYBRID::TIMESTEPPER(*_strandMesh, *_stretchingEnergyStrand, 
                                            *_tetMesh, *_hyperelastic,
                                            *_triangleMesh, *_strechingEnergyShell, *_bendingEnergyShell);

    const int maxNewtonIterations = 3;  // works with warm start
    _hybridSolver->maxNewtonIterations() = maxNewtonIterations;
    //_hybridSolver = new VOLUME::BDF_2(*_tetMesh, *_hyperelastic);
    //
    // FLT_EPSILON in VOLUME::TIMESTEPPER::findNewSurfaceConstraints does not play nice here;
    // it seems to like it when it's set to zero. Needs further investigation.
    // Plus the conditioning of the matrix seems much worse than position-based
    //_hybridSolver = new VOLUME::BACKWARD_EULER_VELOCITY(*_tetMesh, *_hyperelastic);
    // _hybridSolver->setDt(1.0 / 30.0);
    // _subStep = 1;
    _hybridSolver->setDt(1.0 / 90.0);
    _subStep = 3;

    _kinematicShapes.reserve(10);
    vector<VECTOR3> centers;
    centers.reserve(10);

    // VECTOR3 center(0.0, -10, 0.0);
    // centers.push_back(center);
    // _kinematicShapes.push_back(new CUBE(centers.back(), 10.0));
    // _hybridSolver->addKinematicCollisionObject(_kinematicShapes.back());
    
    // ground
    // VECTOR3 center(0.0, -35, 0.0);
    // centers.push_back(center);
    // _kinematicShapes.push_back(new CUBE(centers.back(), 40.0));
    // _hybridSolver->addKinematicCollisionObject(_kinematicShapes.back());

    // bowl
    VECTOR3 center(0.0, -3.0, 0.0);
    centers.push_back(center);
    BOWL* bowlPtr = new BOWL(centers.back(), 0.9,  12.0);
    _kinematicShapes.push_back(bowlPtr);
    _hybridSolver->addKinematicCollisionObject(_kinematicShapes.back());
    int q = 50; // subdivisions of longitudinal arc
    int p = 100; // subdivisions of latitudinal circles 
    string bowlPath = _outDir + string("bowl.obj");
    cout<<"bowl path: "<<bowlPath<<endl;
    bowlPtr->writeToObj(bowlPath, p, q);

    // container
    center = VECTOR3(0.0, 8.0, 0.0);
    centers.push_back(center);
    _container = new CUBE(centers.back(), 14.0);
    _kinematicShapes.push_back(_container);
    _hybridSolver->attachKinematicSurfaceConstraints(_kinematicShapes.back(), true);

    // center = VECTOR3(-2.5, 1.0, 0.0);
    // centers.push_back(center);
    // _kinematicShapes.push_back(new CUBE(centers.back(), 6.0));
    // _hybridSolver->attachKinematicSurfaceConstraints(_kinematicShapes.back(), true);
    // _hybridSolver->addKinematicCollisionObject(_kinematicShapes.back());

    // top container                                
    // center = VECTOR3(0.0, 25, 0.0);
    // centers.push_back(center);
    // _kinematicShapes.push_back(new CUBE(centers.back(), 40.0));
    // _hybridSolver->attachKinematicSurfaceConstraints(_kinematicShapes.back(), true);
    // _hybridSolver->addKinematicCollisionObject(_kinematicShapes.back());
    
    // center = VECTOR3(0.25, 0.0, 0.25);
    // centers.push_back(center);
    // _kinematicShapes.push_back(new CUBE(centers.back(), 1.0));
    // _kinematicShapes.back()->rotation() = Eigen::AngleAxisd(M_PI * 0.25, VECTOR3::UnitZ());
    // _hybridSolver->addKinematicCollisionObject(_kinematicShapes.back());

    // center = VECTOR3(2.0, -0.75, 0.25);
    // centers.push_back(center);
    // _kinematicShapes.push_back(new CUBE(centers.back(), 1.0));
    // _kinematicShapes.back()->rotation() = Eigen::AngleAxisd(M_PI * 0.25, VECTOR3::UnitZ());
    // _hybridSolver->addKinematicCollisionObject(_kinematicShapes.back());

    // center = VECTOR3(0.25, -1.5, 0.25);
    // centers.push_back(center);
    // _kinematicShapes.push_back(new CUBE(centers.back(), 1.0));
    // _kinematicShapes.back()->rotation() = Eigen::AngleAxisd(M_PI * 0.25, VECTOR3::UnitZ());
    // _hybridSolver->addKinematicCollisionObject(_kinematicShapes.back());

    // center = VECTOR3(2.0, -2.25, 0.25);
    // centers.push_back(center);
    // _kinematicShapes.push_back(new CUBE(centers.back(), 1.0));
    // _kinematicShapes.back()->rotation() = Eigen::AngleAxisd(M_PI * 0.25, VECTOR3::UnitZ());
    // _hybridSolver->addKinematicCollisionObject(_kinematicShapes.back());

    // center = VECTOR3(0.25, -3.0, 0.25);
    // centers.push_back(center);
    // _kinematicShapes.push_back(new CUBE(centers.back(), 1.0));
    // _kinematicShapes.back()->rotation() = Eigen::AngleAxisd(M_PI * 0.25, VECTOR3::UnitZ());
    // _hybridSolver->addKinematicCollisionObject(_kinematicShapes.back());

    // center = VECTOR3(2.0, -3.75, 0.25);
    // centers.push_back(center);
    // _kinematicShapes.push_back(new CUBE(centers.back(), 1.0));
    // _kinematicShapes.back()->rotation() = Eigen::AngleAxisd(M_PI * 0.25, VECTOR3::UnitZ());
    // _hybridSolver->addKinematicCollisionObject(_kinematicShapes.back());

    // center = VECTOR3(0.25, -4.5, 0.25);
    // centers.push_back(center);
    // _kinematicShapes.push_back(new CUBE(centers.back(), 1.0));
    // _kinematicShapes.back()->rotation() = Eigen::AngleAxisd(M_PI * 0.25, VECTOR3::UnitZ());
    // _hybridSolver->addKinematicCollisionObject(_kinematicShapes.back());

    // center = VECTOR3(2.0, -5.25, 0.25);
    // centers.push_back(center);
    // _kinematicShapes.push_back(new CUBE(centers.back(), 1.0));
    // _kinematicShapes.back()->rotation() = Eigen::AngleAxisd(M_PI * 0.25, VECTOR3::UnitZ());
    // _hybridSolver->addKinematicCollisionObject(_kinematicShapes.back());
    // collision constants
    // const REAL collisionMu = 1000.0;
    _hybridSolver->collisionsEnabled() = true;
    _hybridSolver->pcgEnabled() = false;
    _hybridSolver->hessianClampingEnabled() = true;
    _hybridSolver->collisionStiffnessVolume() = 1000;
    _hybridSolver->collisionDampingBetaVolume() = 0.01;
    _hybridSolver->collisionStiffnessStrand() = 1e5;
    _hybridSolver->collisionDampingBetaStrand() = 0.0;
    _hybridSolver->collisionStiffnessShell() = 2e3;
    _hybridSolver->collisionDampingBetaShell() = 0.01;
    _hybridSolver->collisionDampingBetaAll() = 0.0001;
    _hybridSolver->collisionStiffnessAll() = 1e4;
    _hybridSolver->collisionEpsEdgeEdge() = 0.02;
    _hybridSolver->collisionEpsVertexFace() = 0.02;
    _hybridSolver->vertexFaceSelfCollisionsOn() = true;
    _hybridSolver->edgeEdgeSelfCollisionsOn() = true;
    _hybridSolver->setRayelighStrand(0.1, 0.1);
    _hybridSolver->setRayelighShell(0.1, 0.1);

    _strandMesh->bendingForceFilterEnabled() = false;

    // _eye    = VECTOR3(1.7, -2.25, 8.5);
    // _lookAt = VECTOR3(1.6, -2.25, 7.5);
    // _up     = VECTOR3(0.0, 1.0, 0.0);

    _eye    = VECTOR3(0.7, 2.5, 40.5);
    _lookAt = VECTOR3(0.6, -3.5, 7.5);
    _up     = VECTOR3(0.0, 1.0, 0.0);

    _worldCenter = VECTOR3( 0.497, 0.785889, 0.452556);
    // _pauseFrame = 500;
    _pauseFrame = 400;
    return true;
  }

  // void buildHairTieMany() 
  // {
  //   cout<<"building hair tie."<<endl;
  //   vector<VECTOR3> restVertices;
  //   vector<vector<int> > strands;

  //   // const REAL amplitude = 0.5;
  //   _spacing = 0.1;
  //   const REAL radiusCurl = 0.50;
  //   const REAL radiusRing = 2.0;
  //   _totalPoints = 100;
  //   int frequency = 5;
  //   REAL ringStep = 2 * M_PI / _totalPoints;
  //   REAL curlStep = 2 * M_PI / frequency;
  //   REAL thetaRing = 0.0;
  //   REAL thetaCurl = 0.0;
  //   int idxCurl = 0;
  //   vector<int> strand0;
  //   for(unsigned int i = 0; i < _totalPoints; i++) 
  //   {
  //     VECTOR3 ringPosition(radiusRing * cos(thetaRing), 0.0, - radiusRing * sin(thetaRing));
  //     VECTOR3 curlBasis0 = ringPosition.normalized() * radiusCurl;
  //     VECTOR3 curlBasis1(0.0, radiusCurl, 0.0);
  //     VECTOR3 curlPosition = cos(thetaCurl) * curlBasis0 + sin(thetaCurl) * curlBasis1;
  //     restVertices.push_back(ringPosition + curlPosition);
  //     strand0.push_back(i);
  //     thetaRing += ringStep;
  //     idxCurl = (idxCurl + 1) % frequency;
  //     thetaCurl = -idxCurl * curlStep;
  //   }
  //   strand0.push_back(0);
  //   strands.push_back(strand0);

  //   for (unsigned int x = 0; x < restVertices.size(); x++)
  //     restVertices[x] = _initialAStrand * restVertices[x] + _initialTranslationStrand;
    
  //   _strandMesh = new TET_WISP_MESH(restVertices, strands,
  //                                   _E, _G, _density, _baseRadius, _tipRadius);
  //   cout<<"done."<<endl;
  // }

  void buildHairTieMany()
  {
    cout<<"building hair tie."<<endl;
    vector<VECTOR3> restVertices;
    vector<vector<int> > strands;

    // const REAL amplitude = 0.5;
    _spacing = 5.8;
    const REAL radiusCurl = 0.50;
    const REAL radiusRing = 2.0;
    _totalPoints = 100;
    int frequency = 5;
    REAL ringStep = 2 * M_PI / _totalPoints;
    REAL curlStep = 2 * M_PI / frequency;
    REAL thetaRing = 0.0;
    REAL thetaCurl = 0.0;
    int idxCurl = 0;
    vector<int> strand0;
    const int numStrands = 4;
    vector<VECTOR3> centerPositions(numStrands);
    centerPositions[0] = VECTOR3(0.0, 0.0, 0.0);
    centerPositions[1] = VECTOR3(_spacing, 0.0, 0.0);
    centerPositions[2] = VECTOR3(_spacing * 0.5, _spacing * 0.5, -_spacing * 0.5);
    centerPositions[3] = VECTOR3(_spacing * 0.5, _spacing * 0.5, _spacing * 0.5);
    
    // vector<VECTOR3> centerPositions(numStrands);
    // centerPositions[0] = VECTOR3(_spacing * 0.5, 0.0, 0.0);
    // centerPositions[1] = VECTOR3(_spacing, 0.0, 0.0);
    // centerPositions[2] = VECTOR3(_spacing * 0.5, _spacing * 0.5, -_spacing * 0.5);
    // centerPositions[3] = VECTOR3(_spacing * 0.5, _spacing * 0.5, _spacing * 0.5);


    for(unsigned c = 0; c < numStrands; c++){
      strand0.clear();
      for(unsigned int i = 0; i < _totalPoints; i++) 
      {
        VECTOR3& centerPosition = centerPositions[c];
        VECTOR3 ringPosition(radiusRing * cos(thetaRing), 0.0, - radiusRing * sin(thetaRing));
        VECTOR3 curlBasis0 = ringPosition.normalized() * radiusCurl;
        VECTOR3 curlBasis1(0.0, radiusCurl, 0.0);
        VECTOR3 curlPosition = cos(thetaCurl) * curlBasis0 + sin(thetaCurl) * curlBasis1;
        restVertices.push_back(ringPosition + curlPosition + centerPosition);
        strand0.push_back(i + c * _totalPoints);
        thetaRing += ringStep;
        idxCurl = (idxCurl + 1) % frequency;
        thetaCurl = -idxCurl * curlStep;
      }
      strand0.push_back(c * _totalPoints);
      strands.push_back(strand0);
    }

    for (unsigned int x = 0; x < restVertices.size(); x++)
      restVertices[x] = _initialAStrand * restVertices[x] + _initialTranslationStrand;
    
    _strandMesh = new TET_WISP_MESH(restVertices, strands,
                                    _E, _G, _density, _baseRadius, _tipRadius);
    cout<<"done."<<endl;
  }

  virtual void stepSimulation(const bool verbose = true) override
  {
    // write out the initial frame
    if (_writeToFile && _frameNumber == 0)
    {
      char buffer[256];
      sprintf(buffer, "%04i", _frameNumber);
      string filenameVolume = _volumeFilePath + string("initial.obj");
      string filenameStrand = _strandFilePath + string("initial.sobj");
      string filenameShell = _shellFilePath + string("initial.obj");
      // strand
      cout << " Writing file " << filenameStrand.c_str() << " ... " << flush;
      writeSOBJFile(filenameStrand.c_str(), _strandMesh->vertices(), _strandMesh->strandIndices());
      cout << "done." << endl;

      // volume
      cout << " Writing file " << filenameVolume.c_str() << " ... " << flush;
      vector<VECTOR3> surfaceV;
      const vector<int>& surfaceVIndices = _tetMesh->surfaceVertices();
      const vector<VECTOR3>& tetVertices = _tetMesh->vertices(); 
      for(unsigned int i = 0; i < surfaceVIndices.size(); i++)
      {
        surfaceV.push_back(tetVertices[surfaceVIndices[i]]);
      }
      writeOBJFile(filenameVolume.c_str(), surfaceV, _tetMesh->surfaceTrianglesIntoSurfaceVertices());
      cout << "done." << endl;

      // shell
      cout << " Writing file " << filenameShell.c_str() << " ... " << flush;
      writeOBJFile(filenameShell.c_str(), _triangleMesh->vertices(), _triangleMesh->triangles());
      cout << "done." << endl;

      TIMER gzipper("Gzipping");
      string gzip = string("gzip -f -9 ");
      cout << " Gzipping strand mesh ..." << flush;
      system((gzip + filenameStrand).c_str());
      cout << " Gzipping volume mesh ..." << flush;
      system((gzip + filenameVolume).c_str());
      cout << " Gzipping shell mesh ..." << flush;
      system((gzip + filenameShell).c_str());
      cout << "done." << endl;
      gzipper.stop();
    }
    for(unsigned int i = 0; i < _subStep; i++){
      _hybridSolver->externalForces().setZero();
      _hybridSolver->addGravity(_gravity);
      _hybridSolver->solveDynamics(verbose);
    }
    // _hybridSolver->solveNewton(verbose);

    if (_writeToFile) {
      char buffer[256];
      sprintf(buffer, "%04i", _frameNumber);
      string filenameVolume = _volumeFilePath + _sceneName + string("_frame_") + string(buffer) + string(".obj");
      string filenameStrand = _strandFilePath + _sceneName + string("_frame_") + string(buffer) + string(".sobj");
      string filenameShell = _shellFilePath + _sceneName + string("_frame_") + string(buffer) + string(".obj");

      // strand
      cout << " Writing file " << filenameStrand.c_str() << " ... " << flush;
      writeSOBJFile(filenameStrand.c_str(), _strandMesh->vertices(), _strandMesh->strandIndices());
      cout << "done." << endl;

      // volume
      cout << " Writing file " << filenameVolume.c_str() << " ... " << flush;
      vector<VECTOR3> surfaceV;
      const vector<int>& surfaceVIndices = _tetMesh->surfaceVertices();
      const vector<VECTOR3>& tetVertices = _tetMesh->vertices(); 
      for(unsigned int i = 0; i < surfaceVIndices.size(); i++)
      {
        surfaceV.push_back(tetVertices[surfaceVIndices[i]]);
      }
      writeOBJFile(filenameVolume.c_str(), surfaceV, _tetMesh->visibleSurfaceTrianglesIntoSurfaceVertices());
      cout << "done." << endl;

      // shell
      cout << " Writing file " << filenameShell.c_str() << " ... " << flush;
      writeOBJFile(filenameShell.c_str(), _triangleMesh->vertices(), _triangleMesh->triangles());
      cout << "done." << endl;

      TIMER gzipper("Gzipping");
      string gzip = string("gzip -f -9 ");
      cout << " Gzipping strand mesh ..." << flush;
      system((gzip + filenameStrand).c_str());
      cout << " Gzipping volume mesh ..." << flush;
      system((gzip + filenameVolume).c_str());
      cout << " Gzipping shell mesh ..." << flush;
      system((gzip + filenameShell).c_str());
      cout << "done." << endl;
      gzipper.stop();
    }


    if(_frameNumber == 80){
      _hybridSolver->removeKinematicSurfaceConstraints();
      _hybridSolver->setStrandKinematicConstraint(false);
      _hybridSolver->attachKinematicSurfaceConstraints(_container, true);
    }

    // if(_frameNumber == 130)
    //   _hybridSolver->popBackKinematicSurfaceConstraints();

    _frameNumber++;
  }
#ifndef GL_DISABLED
  virtual void drawScene()
  {
    //drawAxes();
    drawStrandMesh(*_strandMesh);
    drawSurfaceTriangles(*_tetMesh, false, VECTOR3(0.0, 0.5, 0.5), VECTOR3(0,0,0));
    // drawTriangleMesh(*_triangleMesh, false);
    //drawStrandMeshOld(*_strandMesh);

    glEnable(GL_DEPTH_TEST);
    for (unsigned int x = 0; x < _kinematicShapes.size(); x++)
      drawKinematicShape(*_kinematicShapes[x]);

    //drawCollisionsOld(*_strandMesh);
    // drawCollisions(*_strandMesh);
    // drawKinematicConstraints(_triangleMesh, _hybridSolver);
    // drawVertex(*_triangleMesh, _arrowCounter);
    /*
    drawKinematicConstraints(_strandMesh, _strandSolver);
    //drawPlaneConstraints(_triangleMesh, _strandSolver);
    drawStrandTwistFreeFrames(*_strandMesh);
    drawStrandTwistFrames(*_strandMesh);
    drawStrandTwistForces(*_strandMesh);
    */
  };
#endif
protected:
  // strand generation parameters
  string _volumeFilePath;
  string _strandFilePath;
  string _shellFilePath;
  string _outDir;
  bool _writeToFile;
  char* _strandMeshFilename;
  string _triangleMeshFilename;
  STRAND_MESH* _strandMesh;
  STRAND::STRETCHING* _stretchingEnergyStrand;
  TRIANGLE_MESH* _triangleMesh;
  SHELL::STRETCHING* _strechingEnergyShell;
  SHELL::BENDING* _bendingEnergyShell;
  // strand building parameters
  unsigned int _totalPoints;
  unsigned int _totalStrands;
  int _subStep;
  REAL _spacing;
  REAL _E;
  REAL _nu;
  REAL _density;
  REAL _baseRadius;
  REAL _tipRadius;

  HYBRID::TIMESTEPPER* _hybridSolver;
  KINEMATIC_SHAPE* _container;

  // GL drawing params
  bool _drawVertices;

  // initial rotation-scale and translation of tet mesh
  MATRIX3 _initialAStrand;
  VECTOR3 _initialTranslationStrand;
  MATRIX3 _initialAShell;
  VECTOR3 _initialTranslationShell;
};

}

#endif

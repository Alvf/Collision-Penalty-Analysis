#ifndef BUNNY_STRAND_DROP_DEBUG_H
#define BUNNY_STRAND_DROP_DEBUG_H

#include "SIMULATION_SCENE.h"
#include "Timestepper/Volume/BDF_1.h"
#include "Timestepper/Volume/BDF_2.h"
#include "Timestepper/Strand_Volume/TIMESTEPPER.h"
#include "Timestepper/Strand/VOLUME_TIMESTEPPER.h"
#include "Timestepper/Strand/VOLUME_BDF_1.h"
#include "Geometry/STRAND_MESH_FASTER.h"
#include "Geometry/STRAND_MESH.h"
#include "Hyperelastic/Strand/QUADRATIC_STRETCHING.h"
#include "Hyperelastic/Strand/BERGOU_2010.h"
#include "Hyperelastic/Strand/ISOTROPIC_BENDING.h"
#include "Hyperelastic/Strand/QUADRATIC_BENDING.h"
#include "Hyperelastic/Strand/QUADRATIC_UNIT_BENDING.h"
#include "Hyperelastic/Strand/ISOTROPIC_THETA.h"

namespace HOBAK {

class BUNNY_STRAND_DROP_DEBUG : public SIMULATION_SCENE {
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
    _sceneName = "bunny_strand_drop_debug";

    // read in the tet mesh file
    //_tetMeshFilename = string("../data/scorpion_0_125.tobj");
    // _tetMeshFilename = string("../data/bunny/bunny_5.tobj");
    _strandMeshFilename = "../data/singleCurl.sobj";
    //_tetMeshFilename = string("../data/Brachiosaurus_15.tobj");
    //_tetMeshFilename = string("../data/dragon.tobj");

    // read tet
    // vector<VECTOR3> vertices;
    // vector<VECTOR4I> tets;
    // bool success = TET_MESH::readTobjFile(_tetMeshFilename, vertices, tets);
    // if (!success)
    // {
    //   cout << " Failed to open file " << _tetMeshFilename.c_str() << endl;
    //   return false;
    // }
    // vertices = TET_MESH::normalizeVertices(vertices);
    // _normalizedVertices = true;

    // read strand
    vector<vector<int> > strandIndices;
    std::vector<VECTOR3> restVerticesStrand;
    if (!readSOBJFile(_strandMeshFilename, restVerticesStrand, strandIndices)) { 
      cout << " Failed to open file " << _strandMeshFilename << endl;
      std::exit(1); }

    using namespace Eigen;
    _E = 3.9e9;
    _G = 3.9e9;
    // _nu = 500;
    _spacing = 1.0;
    _density = 1.3;
    _baseRadius = 0.1;
    _tipRadius = 0.1;
    _drawVertices = true;
    _totalStrands = 1;

    // _initialAStrand = 0.15 * MATRIX3::Identity() * AngleAxisd(0.5 * M_PI, VECTOR3::UnitZ());
    _initialAStrand = MATRIX3::Identity() * AngleAxisd(0.5 * M_PI, VECTOR3::UnitZ());
    // _initialAStrand = MATRIX3::Identity();
    _initialTranslationStrand = VECTOR3(2.0, -4.0, 0.0);
    // _initialTranslationStrand = VECTOR3(2.0, -1.0, 0.0);
    for (unsigned int x = 0; x < restVerticesStrand.size(); x++)
      restVerticesStrand[x] = _initialAStrand * restVerticesStrand[x] + _initialTranslationStrand;

    _strandMesh = new TET_WISP_MESH(restVerticesStrand, strandIndices,
                          _E, _G, _density, _baseRadius, _tipRadius);
    _strandMesh -> setCollisionEps(0.01);
    _strandMesh->bendingForceFilterEnabled() = true;

    const REAL k = 1.0;
    _stretchingEnergy = new STRAND::QUADRATIC_STRETCHING(k);

    // MATRIX3 M;
    // M =   AngleAxisd(-0.5 * M_PI, VECTOR3::UnitX())
    // //M =   AngleAxisd(0 , VECTOR3::UnitX())
    //     * AngleAxisd(0,  VECTOR3::UnitY())
    //     * AngleAxisd(0, VECTOR3::UnitZ());
    // VECTOR3 half(0.5, 0.5, 0.5);
   
    // _initialA           = M;
    // _initialTranslation = half - M * half + VECTOR3(-1.0,-5.0,0.0);
    // // _initialTranslation = half - M * half;

    // for (unsigned int x = 0; x < vertices.size(); x++)
    //   vertices[x] = _initialA * vertices[x] + _initialTranslation;

    _gravity = VECTOR3(0, -1.0, 0);

    // REAL E = 3.0;
    //REAL nu = 0.3; // lambda \approx 10
    // REAL nu = 0.45; // lambda \approx 10

    // REAL mu     = VOLUME::HYPERELASTIC::computeMu(E, nu);
    // REAL lambda = VOLUME::HYPERELASTIC::computeLambda(E, nu);
    // cout << " mu:     " << mu << endl;
    // cout << " lambda: " << lambda << endl;

    // build the tet mesh object
    // _tetMesh = new TET_MESH_FASTER(vertices, tets);
    // _hyperelastic = new VOLUME::SNH(mu, lambda);
    //_hyperelastic = new VOLUME::ARAP(mu * 4, lambda);

    // const vector<REAL>& areas = _tetMesh->surfaceTriangleAreas();
    // REAL smallest = areas[0];
    // REAL largest = areas[0];
    // for (unsigned int x = 1; x < areas.size(); x++)
    // {
    //   if (areas[x] > largest) largest  = areas[x];
    //   if (areas[x] < largest) smallest = areas[x];
    // }
    // cout << "Largest triangle area:  "  << largest << endl;
    // cout << "Smallest triangle area: "  << smallest << endl;

    // build the time integrator
    // _strandSolver = new VOLUME::BACKWARD_EULER_POSITION(*_tetMesh, *_hyperelastic);
    // _strandSolver = new STRAND::VOLUME_BDF_1(*_strandMesh, *_stretchingEnergy);
    // cout<<"------------using volume bdf 1.---------"<<endl;
    _strandSolver = new STRAND::VOLUME_TIMESTEPPER(*_strandMesh, *_stretchingEnergy);
    cout<<"------------using volume timestepper.---------"<<endl;

    const int maxNewtonIterations = 3;  // works with warm start
    _strandSolver->maxNewtonIterations() = maxNewtonIterations;
    //_strandSolver = new VOLUME::BDF_2(*_tetMesh, *_hyperelastic);
    //
    // FLT_EPSILON in VOLUME::TIMESTEPPER::findNewSurfaceConstraints does not play nice here;
    // it seems to like it when it's set to zero. Needs further investigation.
    // Plus the conditioning of the matrix seems much worse than position-based
    //_strandSolver = new VOLUME::BACKWARD_EULER_VELOCITY(*_tetMesh, *_hyperelastic);
    _strandSolver->setDt(1.0 / 30.0);

    _kinematicShapes.reserve(10);
    vector<VECTOR3> centers;
    centers.reserve(10);

    VECTOR3 center(0.0, -25, 0.0);
    centers.push_back(center);
    _kinematicShapes.push_back(new CUBE(centers.back(), 40.0));
    _strandSolver->addKinematicCollisionObject(_kinematicShapes.back());
    
    // VECTOR3 center(0.0, -25, 0.0);
    // centers.push_back(center);
    // _kinematicShapes.push_back(new CUBE(centers.back(), 40.0));
    // _strandSolver->addKinematicCollisionObject(_kinematicShapes.back());
    // center = VECTOR3(0.25, 0.0, 0.25);
    // centers.push_back(center);
    // _kinematicShapes.push_back(new CUBE(centers.back(), 1.0));
    // _kinematicShapes.back()->rotation() = Eigen::AngleAxisd(M_PI * 0.25, VECTOR3::UnitZ());
    // _strandSolver->addKinematicCollisionObject(_kinematicShapes.back());

    // center = VECTOR3(2.0, -0.75, 0.25);
    // centers.push_back(center);
    // _kinematicShapes.push_back(new CUBE(centers.back(), 1.0));
    // _kinematicShapes.back()->rotation() = Eigen::AngleAxisd(M_PI * 0.25, VECTOR3::UnitZ());
    // _strandSolver->addKinematicCollisionObject(_kinematicShapes.back());

    // center = VECTOR3(0.25, -1.5, 0.25);
    // centers.push_back(center);
    // _kinematicShapes.push_back(new CUBE(centers.back(), 1.0));
    // _kinematicShapes.back()->rotation() = Eigen::AngleAxisd(M_PI * 0.25, VECTOR3::UnitZ());
    // _strandSolver->addKinematicCollisionObject(_kinematicShapes.back());

    // center = VECTOR3(2.0, -2.25, 0.25);
    // centers.push_back(center);
    // _kinematicShapes.push_back(new CUBE(centers.back(), 1.0));
    // _kinematicShapes.back()->rotation() = Eigen::AngleAxisd(M_PI * 0.25, VECTOR3::UnitZ());
    // _strandSolver->addKinematicCollisionObject(_kinematicShapes.back());

    // center = VECTOR3(0.25, -3.0, 0.25);
    // centers.push_back(center);
    // _kinematicShapes.push_back(new CUBE(centers.back(), 1.0));
    // _kinematicShapes.back()->rotation() = Eigen::AngleAxisd(M_PI * 0.25, VECTOR3::UnitZ());
    // _strandSolver->addKinematicCollisionObject(_kinematicShapes.back());

    // center = VECTOR3(2.0, -3.75, 0.25);
    // centers.push_back(center);
    // _kinematicShapes.push_back(new CUBE(centers.back(), 1.0));
    // _kinematicShapes.back()->rotation() = Eigen::AngleAxisd(M_PI * 0.25, VECTOR3::UnitZ());
    // _strandSolver->addKinematicCollisionObject(_kinematicShapes.back());

    // center = VECTOR3(0.25, -4.5, 0.25);
    // centers.push_back(center);
    // _kinematicShapes.push_back(new CUBE(centers.back(), 1.0));
    // _kinematicShapes.back()->rotation() = Eigen::AngleAxisd(M_PI * 0.25, VECTOR3::UnitZ());
    // _strandSolver->addKinematicCollisionObject(_kinematicShapes.back());

    // center = VECTOR3(2.0, -5.25, 0.25);
    // centers.push_back(center);
    // _kinematicShapes.push_back(new CUBE(centers.back(), 1.0));
    // _kinematicShapes.back()->rotation() = Eigen::AngleAxisd(M_PI * 0.25, VECTOR3::UnitZ());
    // _strandSolver->addKinematicCollisionObject(_kinematicShapes.back());
    // collision constants
    const REAL collisionMu = 1000.0;
    _strandSolver->collisionsEnabled() = false;
    _strandSolver->pcgEnabled() = true;
    _strandSolver->hessianClampingEnabled() = true;
    // _strandSolver->collisionStiffnessVolume() = collisionMu;
    // _strandSolver->collisionDampingBetaVolume() = 0.01;
    _strandSolver->collisionStiffness() = 1000;
    _strandSolver->collisionDampingBeta() = 0.0;
    // _strandSolver->vertexFaceSelfCollisionsOn() = true;
    // _strandSolver->edgeEdgeSelfCollisionsOn() = true;

    _strandMesh->bendingForceFilterEnabled() = false;
    cout<<"------------bending force filtering enabled: "<< _strandMesh->bendingForceFilterEnabled()<<"---------"<<endl;

    // _eye    = VECTOR3(1.7, -2.25, 8.5);
    // _lookAt = VECTOR3(1.6, -2.25, 7.5);
    // _up     = VECTOR3(0.0, 1.0, 0.0);

    _eye    = VECTOR3(1.7, -2.25, 30.5);
    _lookAt = VECTOR3(1.6, -2.25, 7.5);
    _up     = VECTOR3(0.0, 1.0, 0.0);

    _worldCenter = VECTOR3( 0.497, 0.785889, 0.452556);
    //_pauseFrame = 800;
    _pauseFrame = 400;
    return true;
  }

  virtual void stepSimulation(const bool verbose = true) override
  {
    _strandSolver->externalForces().setZero();
    _strandSolver->addGravity(_gravity);
    _strandSolver->solveDynamics(verbose);
    // _strandSolver->solveNewton(verbose);

    if (_writeToFile) {
      if (!writeFrameToFile()) {
        // Do something here that isn't too noisy.
      }
    }

    _frameNumber++;
  }
#ifndef GL_DISABLED
  virtual void drawScene()
  {
    //drawAxes();
    drawStrandMesh(*_strandMesh);
    // drawSurfaceTriangles(*_tetMesh, true);
    //drawStrandMeshOld(*_strandMesh);

    glEnable(GL_DEPTH_TEST);
    for (unsigned int x = 0; x < _kinematicShapes.size(); x++)
      drawKinematicShape(*_kinematicShapes[x]);

    //drawCollisionsOld(*_strandMesh);
    // drawCollisions(*_strandMesh);

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
  char* _strandMeshFilename;
  TET_WISP_MESH* _strandMesh;
  STRAND::STRETCHING* _stretchingEnergy;
  unsigned int _totalPoints;
  unsigned int _totalStrands;
  REAL _spacing;
  REAL _E;
  REAL _nu;
  REAL _density;
  REAL _baseRadius;
  REAL _tipRadius;
  STRAND::TIMESTEPPER* _strandSolver;

  // GL drawing params
  bool _drawVertices;

  // initial rotation-scale and translation of tet mesh
  MATRIX3 _initialAStrand;
  VECTOR3 _initialTranslationStrand;
};

}

#endif

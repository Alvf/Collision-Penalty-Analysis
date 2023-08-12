#ifndef BUNNY_DROP_H
#define BUNNY_DROP_H

#include "Scenes/SIMULATION_SCENE.h"
#include "Timestepper/Volume/BDF_1.h"
#include "Timestepper/Volume/BDF_2.h"
#include "Collision/C_PLANES_EE.h"

namespace HOBAK {

class BUNNY_DROP : public SIMULATION_SCENE {
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
    _sceneName = "bunny_drop";

    // read in the tet mesh file
    //_tetMeshFilename = string("../data/scorpion_0_125.tobj");
    _tetMeshFilename = string("../data/bunny/bunny_5.tobj");
    //_tetMeshFilename = string("../data/Brachiosaurus_15.tobj");
    //_tetMeshFilename = string("../data/dragon.tobj");
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

    using namespace Eigen;
    MATRIX3 M;
    M =   AngleAxisd(-0.5 * M_PI, VECTOR3::UnitX())
    //M =   AngleAxisd(0 , VECTOR3::UnitX())
        * AngleAxisd(0,  VECTOR3::UnitY())
        * AngleAxisd(0, VECTOR3::UnitZ());
    VECTOR3 half(0.5, 0.5, 0.5);
   
    _initialA           = M;
    //_initialTranslation = half - M * half + VECTOR3(0,0,-0.5);
    _initialTranslation = half - M * half;

    for (unsigned int x = 0; x < vertices.size(); x++)
      vertices[x] = _initialA * vertices[x] + _initialTranslation;

    _gravity = VECTOR3(0, -1.0, 0);

    REAL E = 3.0;
    //REAL nu = 0.3; // lambda \approx 10
    REAL nu = 0.45; // lambda \approx 10

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
    
    // setCollisions(1000.0, 0.01);
    setUnsignedCollisions(1000.0, 0.02);
    // build the time integrator
    _volumeSolver = new VOLUME::BACKWARD_EULER_POSITION(*_tetMesh, *_hyperelastic, *_vfGeneral, *_eeGeneral);
    //_volumeSolver = new VOLUME::BDF_2(*_tetMesh, *_hyperelastic);
    //
    // FLT_EPSILON in VOLUME::TIMESTEPPER::findNewSurfaceConstraints does not play nice here;
    // it seems to like it when it's set to zero. Needs further investigation.
    // Plus the conditioning of the matrix seems much worse than position-based
    //_volumeSolver = new VOLUME::BACKWARD_EULER_VELOCITY(*_tetMesh, *_hyperelastic);
    _volumeSolver->setDt(1.0 / 60.0);

    _kinematicShapes.reserve(10);
    vector<VECTOR3> centers;
    centers.reserve(10);

    VECTOR3 center(0.0, -10, 0.0);
    centers.push_back(center);
    _kinematicShapes.push_back(new CUBE(centers.back(), 10.0));
    _volumeSolver->addKinematicCollisionObject(_kinematicShapes.back());
    
    center = VECTOR3(0.25, 0.0, 0.25);
    centers.push_back(center);
    _kinematicShapes.push_back(new CUBE(centers.back(), 1.0));
    _kinematicShapes.back()->rotation() = Eigen::AngleAxisd(M_PI * 0.25, VECTOR3::UnitZ());
    _volumeSolver->addKinematicCollisionObject(_kinematicShapes.back());

    center = VECTOR3(2.0, -0.75, 0.25);
    centers.push_back(center);
    _kinematicShapes.push_back(new CUBE(centers.back(), 1.0));
    _kinematicShapes.back()->rotation() = Eigen::AngleAxisd(M_PI * 0.25, VECTOR3::UnitZ());
    _volumeSolver->addKinematicCollisionObject(_kinematicShapes.back());

    center = VECTOR3(0.25, -1.5, 0.25);
    centers.push_back(center);
    _kinematicShapes.push_back(new CUBE(centers.back(), 1.0));
    _kinematicShapes.back()->rotation() = Eigen::AngleAxisd(M_PI * 0.25, VECTOR3::UnitZ());
    _volumeSolver->addKinematicCollisionObject(_kinematicShapes.back());

    center = VECTOR3(2.0, -2.25, 0.25);
    centers.push_back(center);
    _kinematicShapes.push_back(new CUBE(centers.back(), 1.0));
    _kinematicShapes.back()->rotation() = Eigen::AngleAxisd(M_PI * 0.25, VECTOR3::UnitZ());
    _volumeSolver->addKinematicCollisionObject(_kinematicShapes.back());

    center = VECTOR3(0.25, -3.0, 0.25);
    centers.push_back(center);
    _kinematicShapes.push_back(new CUBE(centers.back(), 1.0));
    _kinematicShapes.back()->rotation() = Eigen::AngleAxisd(M_PI * 0.25, VECTOR3::UnitZ());
    _volumeSolver->addKinematicCollisionObject(_kinematicShapes.back());

    center = VECTOR3(2.0, -3.75, 0.25);
    centers.push_back(center);
    _kinematicShapes.push_back(new CUBE(centers.back(), 1.0));
    _kinematicShapes.back()->rotation() = Eigen::AngleAxisd(M_PI * 0.25, VECTOR3::UnitZ());
    _volumeSolver->addKinematicCollisionObject(_kinematicShapes.back());

    center = VECTOR3(0.25, -4.5, 0.25);
    centers.push_back(center);
    _kinematicShapes.push_back(new CUBE(centers.back(), 1.0));
    _kinematicShapes.back()->rotation() = Eigen::AngleAxisd(M_PI * 0.25, VECTOR3::UnitZ());
    _volumeSolver->addKinematicCollisionObject(_kinematicShapes.back());

    center = VECTOR3(2.0, -5.25, 0.25);
    centers.push_back(center);
    _kinematicShapes.push_back(new CUBE(centers.back(), 1.0));
    _kinematicShapes.back()->rotation() = Eigen::AngleAxisd(M_PI * 0.25, VECTOR3::UnitZ());
    _volumeSolver->addKinematicCollisionObject(_kinematicShapes.back());

    // collision constants
    //const REAL collisionMu = 1000.0;
    //_volumeSolver->collisionStiffness() = collisionMu;
    _volumeSolver->collisionDampingBeta() = 0.01;
    
    _volumeSolver->vertexFaceSelfCollisionsOn() = true;
    _volumeSolver->edgeEdgeSelfCollisionsOn() = true;

    _eye    = VECTOR3(1.7, -2.25, 8.5);
    _lookAt = VECTOR3(1.6, -2.25, 7.5);
    _up     = VECTOR3(0.0, 1.0, 0.0);

    _worldCenter = VECTOR3( 0.497, 0.785889, 0.452556);
    //_pauseFrame = 800;
    _pauseFrame = 400;

    return true;
  }
};

}

#endif

#ifndef NET_SCENE_DROP_H
#define NET_SCENE_DROP_H

#include "Scenes/SIMULATION_SCENE.h"
#include "Geometry/STRAND_MESH.h"
#include "Geometry/STRAND_MESH_FASTER.h"
#include "Geometry/STRAND_NET_MESH.h"
#include "Timestepper/Strand/TIMESTEPPER.h"
#include "Hyperelastic/Strand/QUADRATIC_STRETCHING.h"
#include "Hyperelastic/Strand/ISOTROPIC_BENDING.h"
#include "Hyperelastic/Strand/QUADRATIC_BENDING.h"
#include "Hyperelastic/Strand/QUADRATIC_UNIT_BENDING.h"
#include "Hyperelastic/Strand/ISOTROPIC_THETA.h"

namespace HOBAK {

class NET_SCENE_DROP : public SIMULATION_SCENE {
public:
  virtual void printSceneDescription() override
  {
    cout << "=====================================================================" << endl;
    cout << " Strand grid collision scene" << endl;
    cout << "=====================================================================" << endl;
  }


  // for testing stretch
  virtual bool buildScene() override
  {
    // this will determine the MOV and JSON filenames
    _sceneName = "net_drop";
    _n0 = 21; _n1 = 21; _subSegment = 2;
    _netCorner = VECTOR3(0.75, 1.75, -0.25);
    _netScale = 0.8;

 
    //const REAL bottomStart = 2.0 - 0.5 * _bottomStrandSpacing * _bottomStrands;

    //const int totalTopPoints = 5;
    //const int totalBottomPoints = 8;
    //const int totalTopPoints = 100;
    //const int totalTopPoints = 50;
    //const int totalTopPoints = 30;

    buildNet();

    _gravity = VECTOR3(0,-1.0,0);

    const REAL E = 1e4; // Taz settings, seems to stack well
    //const REAL E = 1e6;
    //const REAL E = 1e5; // things seem to misbehave here
    const REAL nu = 0.36;
    const REAL density = 1.32;
    const REAL radiusA = 0.1;
    const REAL radiusB = 0.1;
    const REAL collisionEps = 0.001;

    _strandMesh = new STRAND_NET_MESH(_restVertices, _strandIndices, E, nu, density, radiusA, radiusB);

    // _strandMesh = new STRAND_MESH_FASTER(restVertices, strandIndices, E, nu, density, radiusA, radiusB);

    _strandMesh->setCollisionEps(collisionEps);
    _strandMesh->bendingForceFilterEnabled() = false;

    // create the integrator
    const REAL k = E * M_PI * radiusA * radiusB;
    cout <<" Stretching stiffness: " << k << endl;
    STRAND::QUADRATIC_STRETCHING* stretchingEnergy = new STRAND::QUADRATIC_STRETCHING(k);
    _strandSolver = new STRAND::TIMESTEPPER(*_strandMesh, *stretchingEnergy);
    _strandSolver->setDt(1.0 / 30.0); // lots of interpenetrations
    //_strandSolver->setDt(1.0 / 100.0);// passes through
    //_strandSolver->setDt(1.0 / 200.0);  // some pass-throughs, looks like a good CCD debugging scenario
    // _strandSolver->setDt(1.0 / 300.0);  // works fine

    _strandSolver->hessianClampingEnabled() = true;
    _strandSolver->pcgEnabled() = false;
    //_strandSolver->collisionStiffness() = 1e5;  // jittery
    //_strandSolver->collisionStiffness() = 100;  // smooth
    _strandSolver->collisionStiffness() = 10;
    _strandSolver->collisionDampingBeta() = 0;

    //const REAL cubeScale = 0.5;
    // REAL sphereScale = 0.3;
    // VECTOR3 center0(0.5, -1.0, 0.5);
    // _kinematicShapes.push_back(new SPHERE(center0, sphereScale));
    // _strandSolver->addKinematicCollisionObject(_kinematicShapes[0]);

    _kinematicShapes.reserve(10);
    vector<VECTOR3> centers;
    centers.reserve(10);

    VECTOR3 center(0.0, -10, 0.0);
    centers.push_back(center);
    _kinematicShapes.push_back(new CUBE(centers.back(), 10.0));
    _strandSolver->addKinematicCollisionObject(_kinematicShapes.back());
    
    center = VECTOR3(0.25, 0.0, 0.25);
    centers.push_back(center);
    _kinematicShapes.push_back(new CUBE(centers.back(), 1.2));
    _kinematicShapes.back()->rotation() = Eigen::AngleAxisd(M_PI * 0.25, VECTOR3::UnitZ());
    _strandSolver->addKinematicCollisionObject(_kinematicShapes.back());

    center = VECTOR3(2.0, -0.75, 0.25);
    centers.push_back(center);
    _kinematicShapes.push_back(new CUBE(centers.back(), 1.2));
    _kinematicShapes.back()->rotation() = Eigen::AngleAxisd(M_PI * 0.25, VECTOR3::UnitZ());
    _strandSolver->addKinematicCollisionObject(_kinematicShapes.back());

    center = VECTOR3(0.25, -1.5, 0.25);
    centers.push_back(center);
    _kinematicShapes.push_back(new CUBE(centers.back(), 1.2));
    _kinematicShapes.back()->rotation() = Eigen::AngleAxisd(M_PI * 0.25, VECTOR3::UnitZ());
    _strandSolver->addKinematicCollisionObject(_kinematicShapes.back());

    center = VECTOR3(2.0, -2.25, 0.25);
    centers.push_back(center);
    _kinematicShapes.push_back(new CUBE(centers.back(), 1.2));
    _kinematicShapes.back()->rotation() = Eigen::AngleAxisd(M_PI * 0.25, VECTOR3::UnitZ());
    _strandSolver->addKinematicCollisionObject(_kinematicShapes.back());

    center = VECTOR3(0.25, -3.0, 0.25);
    centers.push_back(center);
    _kinematicShapes.push_back(new CUBE(centers.back(), 1.2));
    _kinematicShapes.back()->rotation() = Eigen::AngleAxisd(M_PI * 0.25, VECTOR3::UnitZ());
    _strandSolver->addKinematicCollisionObject(_kinematicShapes.back());

    center = VECTOR3(2.0, -3.75, 0.25);
    centers.push_back(center);
    _kinematicShapes.push_back(new CUBE(centers.back(), 1.2));
    _kinematicShapes.back()->rotation() = Eigen::AngleAxisd(M_PI * 0.25, VECTOR3::UnitZ());
    _strandSolver->addKinematicCollisionObject(_kinematicShapes.back());

    center = VECTOR3(0.25, -4.5, 0.25);
    centers.push_back(center);
    _kinematicShapes.push_back(new CUBE(centers.back(), 1.2));
    _kinematicShapes.back()->rotation() = Eigen::AngleAxisd(M_PI * 0.25, VECTOR3::UnitZ());
    _strandSolver->addKinematicCollisionObject(_kinematicShapes.back());

    center = VECTOR3(2.0, -5.25, 0.25);
    centers.push_back(center);
    _kinematicShapes.push_back(new CUBE(centers.back(), 1.2));
    _kinematicShapes.back()->rotation() = Eigen::AngleAxisd(M_PI * 0.25, VECTOR3::UnitZ());
    _strandSolver->addKinematicCollisionObject(_kinematicShapes.back());

    
    _eye    = VECTOR3(1.7, -2.25, 8.5);
    _lookAt = VECTOR3(1.6, -2.25, 7.5);
    _up     = VECTOR3(0.0, 1.0, 0.0);
    // _eye    = VECTOR3(2.0, 3.0, 4.0);
    // _lookAt = VECTOR3(0.3, 0.0, 0.0);
    // _up     = VECTOR3(0.0, 1.0, 0.0);

    // _eye    = VECTOR3(0.327816, 0.174481, 1.04042);
    // _lookAt = VECTOR3(-0.309947, 0.414465, 1.77231);
    // _up     = VECTOR3(-0.108162, -0.96871, 0.223384);

    // _eye    = VECTOR3(2.98866, -1.08057, -2.52685);
    // _lookAt = VECTOR3(2.41212, -0.820819, -1.75217);
    // _up     = VECTOR3(-0.119018, -0.964707, 0.234895);

    // _pauseFrame = 10;
    _pauseFrame = 300;
    //_pauseFrame = 315;
    //_pauseFrame = 91;

    return true;
  }

  // simulation loop
  virtual void stepSimulation(const bool verbose = true) override
  {
    /*
    //const VECTOR3 translation(0, -0.005, 0.0); // 30 explodes
    const VECTOR3 translation(0, -0.01, 0.0); // 30 explodes
    //const VECTOR3 translation(0, -0.02, 0.0); // 20 explodes

    for (unsigned int x = 0; x < 2 * _topStrands; x++)
      _kinematicShapes[x]->translation() -= translation;

    cout << " Translation: " << _kinematicShapes[0]->translation().transpose() << endl;
    */

    _strandSolver->externalForces().setZero();
    _strandSolver->addGravity(_gravity);
    _strandSolver->solveDynamics(verbose);
    //_strandSolver->solveNewton(verbose);

    _frameNumber++;
  };

  void buildNet()
  {
    // procedurally generate a net

    // vertices
    const REAL spacing = _netScale/(REAL(_n0) - 1.0);
    const REAL subSpacing = spacing / REAL(_subSegment);
    // _strandIndices.resize(_n0 + _n1);
    // for(int x = 0; x < (_n0 + _n1); x++)
    //   _strandIndices[x].clear();
    // // cout<<"spacing"<< spacing<<endl;
    // const int num0 = 1 + (_n1 - 1) * _subSegment, num1 = 1 + (_n0 - 1) * _subSegment;
    // int vId = 0;
    // for(int z = 0; z < _n1 - 1; z++){
    //   for(int x = 0; x < _n0 - 1; x++){
    //     _restVertices.push_back(VECTOR3(x * _subSegment * spacing, 0.0, z * _subSegment * spacing));
    //     _strandIndices[z].push_back(vId);
    //     _strandIndices[_n1]
    //     vId ++;
    //     for(int i = 0; i < _subSegment - 1; i++){
    //       _restVertices.push_back(VECTOR3(x * _subSegment * spacing + i * spacing, 0.0, z * _subSegment * spacing));
    //       _strandIndices[z].push_back(vId);
    //       vId ++;
    //     }
    //   }

    // }
    
    // vertex
    for(int z = 0; z < _n1; z++){
      for(int x = 0; x < _n0; x++){
        _restVertices.push_back(VECTOR3(x * spacing, 0.0 , z * spacing) + _netCorner);
      }
    }


    // strands
    // along x
    for(int z = 0; z < _n1; z++){
      vector<int> strand;
      for(int x = 0; x < _n0; x++){
        strand.push_back(_n0 * z + x);
      }
      _strandIndices.push_back(strand);
    }
    // along z
    for(int x = 0; x < _n0; x++){
      vector<int> strand;
      for(int z = 0; z < _n1; z++){
        strand.push_back(_n0 * z + x);
      }
      _strandIndices.push_back(strand);
    }
  }

// drawScene MUST have this guard, so it can be deactivated for regression testing
#ifndef GL_DISABLED
  // what to draw?
  virtual void drawScene()
  {
    const bool drawOld = false;

    drawAxes();
    if (drawOld)
      drawStrandMeshOld(*_strandMesh);
    else
      drawStrandMesh(*_strandMesh, true, false);

    glEnable(GL_DEPTH_TEST);
    for (unsigned int x = 0; x < _kinematicShapes.size(); x++)
      drawKinematicShape(*_kinematicShapes[x]);

#if 0
    if (drawOld)
      drawCollisionsOld(*_strandMesh);
    else
      drawCollisions(*_strandMesh);
#endif

    /*
    drawKinematicConstraints(_strandMesh, _strandSolver);
    //drawPlaneConstraints(_triangleMesh, _strandSolver);
    drawStrandTwistFreeFrames(*_strandMesh);
    drawStrandTwistFrames(*_strandMesh);
    drawStrandTwistForces(*_strandMesh);
    */

    char buffer[512];
    sprintf(buffer, "%i", _frameNumber);

    string drawString("Frame ");
    drawString = drawString + string(buffer);
    printGlString(drawString);
  };
#endif

protected:
  //STRAND_MESH* _strandMesh;
  STRAND_NET_MESH* _strandMesh;
  STRAND::TIMESTEPPER* _strandSolver;
  STRAND::STRETCHING* _strechingEnergy;

  //shape parameters
  int _n0, _n1, _subSegment; // _n0 along x; _n1 along z; 
  VECTOR3 _netCorner;
  REAL _netScale;
  std::vector<VECTOR3> _restVertices;
  vector<vector<int> > _strandIndices;
};

}

#endif

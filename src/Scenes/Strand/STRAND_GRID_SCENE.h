#ifndef STRAND_GRID_SCENE_H
#define STRAND_GRID_SCENE_H

#include "Scenes/SIMULATION_SCENE.h"
#include "Geometry/STRAND_MESH.h"
#include "Geometry/STRAND_MESH_FASTER.h"
#include "Timestepper/Strand/TIMESTEPPER.h"
#include "Hyperelastic/Strand/QUADRATIC_STRETCHING.h"
#include "Hyperelastic/Strand/ISOTROPIC_BENDING.h"
#include "Hyperelastic/Strand/QUADRATIC_BENDING.h"
#include "Hyperelastic/Strand/QUADRATIC_UNIT_BENDING.h"
#include "Hyperelastic/Strand/ISOTROPIC_THETA.h"

namespace HOBAK {

class STRAND_GRID_SCENE : public SIMULATION_SCENE {
public:
  virtual void printSceneDescription() override
  {
    cout << "=====================================================================" << endl;
    cout << " Strand grid collision scene" << endl;
    cout << "=====================================================================" << endl;
  }

  static void printMatlab(const VECTOR& p, const string name)
  {
    // print out some Matlab-ready output
    cout << " " << name.c_str() << " = [";
    for (unsigned int x = 0; x < p.size(); x++)
      cout << p[x] << " ";
    cout << "];" << endl;
  }

  static void printMatlab(const vector<VECTOR3>& p, const string name)
  {
    // print out some Matlab-ready output
    cout << " " << name.c_str() << " = [";
    for (unsigned int y = 0; y < 3; y++)
    {
      for (unsigned int x = 0; x < p.size(); x++)
        cout << p[x][y] << " ";
      if (y != 2)
        cout << ";" << endl << "\t";
      else
        cout << "];" << endl;
    }
  }

  // for testing stretch
  virtual bool buildScene() override
  {
    // this will determine the MOV and JSON filenames
    _sceneName = "strand_grid_collision";

    _topStrands = 10;
    _bottomStrands = 10;
    //_topStrands = 1;
    //_bottomStrands = 8;

    //_topStrands = 1;
    //_topStrands = 3;
    _topStrandSpacing = 0.2;
    //_bottomStrands = 3;
    _bottomStrandSpacing = 0.2;

    // how far down the bottom strand is
    const REAL bottomHeight = 0.2;  // close up test
    //const REAL bottomHeight = 1.5;  // more difficult
    const REAL collisionEps = 0.1;

    const REAL topStart = -0.5 * _topStrandSpacing * _topStrands;
    const REAL bottomStart = 3.0 - 0.5 * _bottomStrandSpacing * _bottomStrands;
    //const REAL bottomStart = 2.0 - 0.5 * _bottomStrandSpacing * _bottomStrands;

    //const int totalTopPoints = 5;
    //const int totalBottomPoints = 8;
    //const int totalTopPoints = 100;
    //const int totalTopPoints = 50;
    //const int totalTopPoints = 30;
    const unsigned int totalTopPoints = 20;
    //const int totalTopPoints = 25;
    const unsigned int totalBottomPoints = 25;
    const REAL deltaTop = 3.0 / (REAL)(totalTopPoints - 1);
    const REAL deltaBottom = 3.0 / (REAL)(totalBottomPoints - 1);
    assert(totalTopPoints >= 4);

    // procedurally generate a scene
    std::vector<VECTOR3> restVertices;

    // vector of all the strands
    vector<vector<int> > strandIndices;

    for (unsigned int y = 0; y < _topStrands; y++)
    {
      for (unsigned int x = 0; x < totalTopPoints; x++)
        restVertices.push_back(VECTOR3(topStart  + y * _topStrandSpacing,0, deltaTop * x));

      // store the indices of individual strands
      vector<int> strand0;
      const unsigned int nextStart = y * totalTopPoints;
      //for (unsigned int x = 0; x < totalTopPoints; x++)
      for (unsigned int x = nextStart; x < nextStart + totalTopPoints; x++)
        strand0.push_back(x);
      strandIndices.push_back(strand0);
    } 

    for (unsigned int y = 0; y < _bottomStrands; y++)
    {
      const unsigned int nextStart = restVertices.size();
      for (unsigned int x = 0; x < totalBottomPoints; x++)
        restVertices.push_back(VECTOR3(-1.5 + x * deltaBottom,bottomHeight, bottomStart - y * _bottomStrandSpacing));
      
      // store the indices of individual strands
      vector<int> strand1;
      for (unsigned int x = nextStart ; x < nextStart + totalBottomPoints; x++)
        strand1.push_back(x);
      strandIndices.push_back(strand1);
    }

    _gravity = VECTOR3(0,981,0);

    const REAL E = 1e7; // Taz settings, seems to stack well
    //const REAL E = 1e6;
    //const REAL E = 1e5; // things seem to misbehave here
    const REAL nu = 0.36;
    const REAL density = 1.32;
    const REAL radiusA = 0.005;
    const REAL radiusB = 0.005;

    //_strandMesh = new STRAND_MESH(restVertices, strandIndices, E, nu, density, radiusA, radiusB);
    _strandMesh = new STRAND_MESH_FASTER(restVertices, strandIndices, E, nu, density, radiusA, radiusB);

    setUnsignedCollisions(10.0, collisionEps);
    //_strandMesh->setCollisionEps(collisionEps);
    _strandMesh->bendingForceFilterEnabled() = false;

    //setUnsignedCollisions(1000.0, 0.1);

    // create the integrator
    const REAL k = E * M_PI * radiusA * radiusB;
    cout <<" Stretching stiffness: " << k << endl;
    STRAND::QUADRATIC_STRETCHING* stretchingEnergy = new STRAND::QUADRATIC_STRETCHING(k);
    _strandSolver = new STRAND::TIMESTEPPER(*_strandMesh, *stretchingEnergy, *_eeGeneral);
    //_strandSolver->setDt(1.0 / 30.0); // lots of interpenetrations
    //_strandSolver->setDt(1.0 / 100.0);// passes through
    //_strandSolver->setDt(1.0 / 200.0);  // some pass-throughs, looks like a good CCD debugging scenario
    _strandSolver->setDt(1.0 / 300.0);  // works fine

    _strandSolver->hessianClampingEnabled() = true;
    //_strandSolver->collisionStiffness() = 1e5;  // jittery
    //_strandSolver->collisionStiffness() = 100;  // smooth
    //_strandSolver->collisionStiffness() = 10;
    _strandSolver->collisionDampingBeta() = 0;

    //const REAL cubeScale = 0.5;
    const REAL cubeScale = 2.1;

    // constrain the first vertex
    for (unsigned int y = 0; y < _topStrands; y++)
    {
      VECTOR3 center1(topStart + y * _topStrandSpacing,0,3);
      _kinematicShapes.push_back(new CUBE(center1, cubeScale * deltaTop));
    }
    for (unsigned int y = 0; y < _bottomStrands; y++)
    {
      VECTOR3 center2(-1.5,bottomHeight,bottomStart - y * _bottomStrandSpacing);
      _kinematicShapes.push_back(new CUBE(center2, cubeScale * deltaBottom));
      VECTOR3 center3( 1.5,bottomHeight,bottomStart - y * _bottomStrandSpacing);
      _kinematicShapes.push_back(new CUBE(center3, cubeScale * deltaBottom));
    }

    for (unsigned int x = 0; x < _kinematicShapes.size(); x++)
      _strandSolver->attachKinematicSurfaceConstraints(_kinematicShapes[x]);

    _eye    = VECTOR3(2.30806, -0.305108, -2.46584);
    _lookAt = VECTOR3(1.84931, -0.102653, -1.60064);
    _up     = VECTOR3(-0.0820219, -0.97919, 0.185639);

    _eye    = VECTOR3(0.327816, 0.174481, 1.04042);
    _lookAt = VECTOR3(-0.309947, 0.414465, 1.77231);
    _up     = VECTOR3(-0.108162, -0.96871, 0.223384);

    _eye    = VECTOR3(2.98866, -1.08057, -2.52685);
    _lookAt = VECTOR3(2.41212, -0.820819, -1.75217);
    _up     = VECTOR3(-0.119018, -0.964707, 0.234895);

    _pauseFrame = 10;
    //_pauseFrame = 99;
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
      drawStrandMesh(*_strandMesh);

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
  STRAND_MESH* _strandMesh;
  STRAND::TIMESTEPPER* _strandSolver;
  STRAND::STRETCHING* _strechingEnergy;

  unsigned int _topStrands;
  REAL _topStrandSpacing;
  
  unsigned int _bottomStrands;
  REAL _bottomStrandSpacing;
};

}

#endif

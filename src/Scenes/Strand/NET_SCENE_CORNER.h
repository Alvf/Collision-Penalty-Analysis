#ifndef NET_SCENE_CORNER_H
#define NET_SCENE_CORNER_H

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
#include "Collision/C_PLANES_EE.h"
#include <filesystem>

namespace HOBAK {

#define USING_FILESYSTEM 0

#if USING_FILESYSTEM
#if __APPLE__
namespace fs = std::__fs::filesystem;
#else
namespace fs = std::filesystem;
#endif
#endif

class NET_SCENE_CORNER : public SIMULATION_SCENE {
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
    // _sceneName = "net_moving_corners";
    _n0 = 41; _n1 = 41; _subSegment = 2;
    // string method = "ours";
    string method = "tan";
    // string method = "brute";
    // string method = "GN";
    // string method = "unFIl";

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

    // file output
    // _writeToFile = true;
    char buffer[100];
    sprintf(buffer, "_E%.1f", E);
    _sceneName = "net_moving_corners" + method + buffer;
    _writeToFile = true;

#if USING_FILESYSTEM
    if(_writeToFile){
      if(!fs::exists("../data/renders"))
        fs::create_directory("../data/renders");
      string outDir = string("../data/renders/") + _sceneName + string("/");
      if(!fs::exists(outDir))
        fs::create_directory(outDir);
      if(!fs::exists(outDir + string("strand")))
        fs::create_directory(outDir + string("strand"));

      _strandFilePath = outDir + string("strand/");
    }
#endif

    _strandMesh = new STRAND_NET_MESH(_restVertices, _strandIndices, E, nu, density, radiusA, radiusB);

    // _strandMesh = new STRAND_MESH_FASTER(restVertices, strandIndices, E, nu, density, radiusA, radiusB);

    _strandMesh->setCollisionEps(collisionEps);
    _strandMesh->bendingForceFilterEnabled() = false;

    //Build the collision penalties
    vector<REAL> springArgs;
    springArgs.push_back(10.0);
    springArgs.push_back(0.001);

    SIGNED_LEN_PLANES* lengthFuncEE = new SIGNED_LEN_PLANES();
    ENERGY_1D* normalSpringEE = new ENERGY_1D(springArgs);
    C_PLANES_EE* cFuncEE = new C_PLANES_EE();

    ENERGY_12D* eeGeneral = new ENERGY_12D(normalSpringEE, cFuncEE, lengthFuncEE);

    // create the integrator
    const REAL k = E * M_PI * radiusA * radiusB;
    cout <<" Stretching stiffness: " << k << endl;
    STRAND::QUADRATIC_STRETCHING* stretchingEnergy = new STRAND::QUADRATIC_STRETCHING(k);
    _strandSolver = new STRAND::TIMESTEPPER(*_strandMesh, *stretchingEnergy, *eeGeneral);
    _strandSolver->setDt(1.0 / 30.0); // lots of interpenetrations
    _strandSolver->collisionsEnabled() = false;
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

    REAL sphereScale = 0.02;
    VECTOR3 center0(0.0, 0.0, 0.0);
    VECTOR3 center1(1.0, 0.0, 0.0);
    VECTOR3 center2(0.0, 0.0, 1.0);
    VECTOR3 center3(1.0, 0.0, 1.0);
    _kinematicShapes.push_back(new SPHERE(center0, sphereScale));
    _kinematicShapes.push_back(new SPHERE(center1, sphereScale));
    _kinematicShapes.push_back(new SPHERE(center2, sphereScale));
    _kinematicShapes.push_back(new SPHERE(center3, sphereScale));
    _strandSolver->attachKinematicSurfaceConstraints(_kinematicShapes[0]);
    _strandSolver->attachKinematicSurfaceConstraints(_kinematicShapes[1]);
    _strandSolver->attachKinematicSurfaceConstraints(_kinematicShapes[2]);
    _strandSolver->attachKinematicSurfaceConstraints(_kinematicShapes[3]);

    

    _eye    = VECTOR3(2.0, 3.0, 4.0);
    _lookAt = VECTOR3(0.3, 0.0, 0.0);
    _up     = VECTOR3(0.0, 1.0, 0.0);

    // _eye    = VECTOR3(0.327816, 0.174481, 1.04042);
    // _lookAt = VECTOR3(-0.309947, 0.414465, 1.77231);
    // _up     = VECTOR3(-0.108162, -0.96871, 0.223384);

    // _eye    = VECTOR3(2.98866, -1.08057, -2.52685);
    // _lookAt = VECTOR3(2.41212, -0.820819, -1.75217);
    // _up     = VECTOR3(-0.119018, -0.964707, 0.234895);

    // _pauseFrame = 10;
    _pauseFrame = 400;
    //_pauseFrame = 315;
    //_pauseFrame = 91;

    return true;
  }

  // simulation loop
  virtual void stepSimulation(const bool verbose = true) override
  {
    if (_writeToFile && _frameNumber == 0)
    {
      char buffer[256];
      sprintf(buffer, "%04i", _frameNumber);
      string filenameStrand = _strandFilePath + _sceneName + string("_frame_") + string(buffer) + string(".sobj");
      // strand
      cout << " Writing file " << filenameStrand.c_str() << " ... " << flush;
      writeSOBJFile(filenameStrand.c_str(), _strandMesh->vertices(), _strandMesh->strandIndices());
      cout << "done." << endl;
    }

    _strandSolver->externalForces().setZero();
    _strandSolver->addGravity(_gravity);
    _strandSolver->solveDynamics(verbose);
    //_strandSolver->solveNewton(verbose);

    if (_writeToFile)
    {
      char buffer[256];
      sprintf(buffer, "%04i", _frameNumber + 1);
      string filenameStrand = _strandFilePath + _sceneName + string("_frame_") + string(buffer) + string(".sobj");
      // strand
      cout << " Writing file " << filenameStrand.c_str() << " ... " << flush;
      writeSOBJFile(filenameStrand.c_str(), _strandMesh->vertices(), _strandMesh->strandIndices());
      cout << "done." << endl;
    }

    // move corners
    _kinematicShapes[0]->translation() = _kinematicShapes[0]->translation() + VECTOR3(0.001, 0.0, 0.001);
    _kinematicShapes[1]->translation() = _kinematicShapes[1]->translation() + VECTOR3(-0.001, 0.0, 0.001);
    _kinematicShapes[2]->translation() = _kinematicShapes[2]->translation() + VECTOR3(0.001, 0.0, -0.001);
    _kinematicShapes[3]->translation() = _kinematicShapes[3]->translation() + VECTOR3(-0.001, 0.0, -0.001);

    _frameNumber++;
  };

  void buildNet()
  {
    // procedurally generate a net

    // vertices
    const REAL spacing = 1.0/(REAL(_n0) - 1.0);
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
        _restVertices.push_back(VECTOR3(x * spacing, 0.0 , z * spacing));
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
  std::vector<VECTOR3> _restVertices;
  vector<vector<int> > _strandIndices;

  bool _writeToFile;
  string _strandFilePath;
};

}

#endif

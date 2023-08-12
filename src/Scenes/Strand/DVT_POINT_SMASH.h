#ifndef DVT_POINT_SMASH_H
#define DVT_POINT_SMASH_H

#include "Scenes/SIMULATION_SCENE.h"
#include "Geometry/WISP_MESH.h"
#include "Geometry/STRAND_MESH_FASTER.h"
#include "Timestepper/Strand/TIMESTEPPER.h"
#include "Hyperelastic/Strand/QUADRATIC_STRETCHING.h"

namespace HOBAK {

class DVT_POINT_SMASH : public SIMULATION_SCENE {
public:
  virtual void printSceneDescription() override
  {
    cout << "=====================================================================" << endl;
    cout << " Hair wisp scene" << endl;
    cout << "=====================================================================" << endl;
  }

  DVT_POINT_SMASH()
  {
    _E = 3.9e9;
    //_G = _E;
    REAL nu = 0.48;
    _G = _E / (2.0 * (1 + nu));
    _density = 1.3; // 5% of the usual density
    //_baseRadius = 0.0037;
    //_tipRadius = 0.0037;
    _baseRadius = 0.1;
    _tipRadius = 0.1;

    //_smash = 8e-1; // recovers
    //_smash = 7.5e-1; // recovers
    //_smash = 7e-1;  // recovers
    //_smash = 6e-1;  // diverges, Gauss-Newton recovers
    //_smash = 1e-1;  // diverges, Gauss-Newton diverges
    //_smash = 5e-1;  // Gauss-Newton converges
    //_smash = 4e-1;  // diverges
    //_smash = 0.25;  // ???
    //_smash = 1e-5;  // diverges
   
    //_smash = 7e-1;  // recovers (DVT)
    //_smash = 6e-1;  // diverges, Gauss-Newton recovers (DVT)
    _smash = 4e-1;  // diverges, Gauss-Newton diverges (DVT)

    //_dt = 1.0 / 300000.0;   // recovers
    //_dt = 1.0 / 30000.0;
    //_dt = 1.0 / 3000.0;
    //_dt = 1.0 / 300.0;
    _dt = 1.0 / 30.0;
    
    _writeToFile = true;
    //_writeToFile = false;
  }

  // for testing stretch
  virtual bool buildScene() override
  {
    _autoplay = false;

    // this will determine the MOV and JSON filenames
    char buffer[512];
    //snprintf(buffer,256, "smash=%.3f", _smash);
    snprintf(buffer,256, "smash=%4.2e", _smash);
    string smashString(buffer);
    _sceneName = "dvt_" + smashString;

    //_totalPoints = 50;
    _spacing = 1.0;

    //_gravity = VECTOR3(0,-981,0);
    _gravity = VECTOR3(0,0,0);

    _totalPoints = 100; // decent for type 4C

    _drawVertices = true;

    // which type of hair?
    //buildType1();
    //buildType4B();
    
    buildType4C();
    //buildType4CAlternate();
    //buildSampleMatch();

    cout << " Strand length: " << _strandMesh->strandLength(0) << endl;

    const REAL k = 1.0;
    STRAND::QUADRATIC_STRETCHING* stretchingEnergy = new STRAND::QUADRATIC_STRETCHING(k);
    //_strandSolver = new STRAND::TIMESTEPPER(*_strandMesh, *stretchingEnergy);
    _strandSolver = new STRAND::BDF_1(*_strandMesh, *stretchingEnergy);

    _substeps = 1;
    //const REAL dt = 1.0 / 3000.0;
    //const REAL dt = 1.0 / 300.0;
    //const REAL dt = 1.0 / 30.0;
    _strandSolver->setDt(_dt);

    const REAL radius = 15.25;
    VECTOR3 sphereCenter0(0,-radius,0); 
    _kinematicShapes.push_back(new SPHERE(sphereCenter0, radius));

    // attach hairs to the sphere
    for (unsigned int x = 0; x < _kinematicShapes.size(); x++)
      _strandSolver->attachKinematicSurfaceConstraints(_kinematicShapes[x]);

    // make the sphere a collision object
    //_strandSolver->addKinematicCollisionObject(_kinematicShapes.back());

    // apply the smash
    smashVertices();
    //scrambleVertices();
    //_strandSolver->position() = _strandMesh->getDisplacement();

    //_strandSolver->collisionsEnabled() = true;
    _strandSolver->collisionsEnabled() = false;
    //_strandSolver->pcgEnabled() = false;
    _strandSolver->pcgEnabled() = true;

    // full filtering
#if 0
    _strandSolver->hessianClampingEnabled() = true;
    _strandSolver->gaussNewtonEnabled() = false;
    _strandSolver->gaussNewtonFilteredEnabled() = false;
    _sceneName = _sceneName + string("_filtered_");
#endif
    // full Gauss-Newton
#if 0
    _strandSolver->hessianClampingEnabled() = false;
    _strandSolver->gaussNewtonEnabled() = true;
    _strandSolver->gaussNewtonFilteredEnabled() = false;
    _sceneName = _sceneName + string("_GN_");
#endif
    // filtered Gauss-Newton
#if 0
    _strandSolver->hessianClampingEnabled() = false;
    _strandSolver->gaussNewtonEnabled() = false;
    _strandSolver->gaussNewtonFilteredEnabled() = true;
    _sceneName = _sceneName + string("_GNF_");
#endif
    // unfiltered
#if 1
    _strandSolver->hessianClampingEnabled() = false;
    _strandSolver->gaussNewtonEnabled() = false;
    _strandSolver->gaussNewtonFilteredEnabled() = false;
    _sceneName = _sceneName + string("_unfiltered_");
#endif
    _strandSolver->collisionStiffness() = 10;
    _strandSolver->collisionDampingBeta() = 0;
    _strandSolver->disablePreconditioner() = true;
    _strandSolver->maxNewtonIterations() = 3;

    _strandMesh->bendingForceFilterEnabled() = false;
    //_strandMesh->bendingForceFilterEnabled() = true;

    //_eye    = VECTOR3(16.6521, 13.8222, 42.1509);
    //_lookAt = VECTOR3(16.3798, 13.7082, 41.1954);
    //_up     = VECTOR3(0.0296061, 0.991492, -0.126738);

    _eye    = VECTOR3(3.21116, 6.56164, 17.5514);
    _lookAt = VECTOR3(3.06594, 6.42116, 16.572);
    _up     = VECTOR3(0.017671, 0.989341, -0.144537);

    _worldCenter = VECTOR3(0,0,0);

    _pauseFrame = 100;
    //_pauseFrame = 200;
    //_pauseFrame = -1;

    //_writeToFile = true;
    //_writeToFile = false;

    snprintf(buffer,256, "E=%4.2e", _E);
    string Estring(buffer);
    snprintf(buffer,256, "r0=%.3f_r1=%.3f", _baseRadius, _tipRadius);
    string Rstring(buffer);

    int totalStrands = _strandMesh->totalStrands();
    cout << " total strands: " << totalStrands << endl;

    snprintf(buffer,256, "s=%d", _totalPoints);
    string Sstring(buffer);
    
    snprintf(buffer,256, "dt=%4.2e", _dt);
    string Tstring(buffer);

    _sceneName = _sceneName + string("_") + Estring + string("_") + Rstring + string("_") + Tstring;
    cout << " Scene name: " << _sceneName.c_str() << endl;
    
    _filePath = string("../data/render/") + _sceneName + string("/");

    if (_writeToFile)
    {
      string mkdir("mkdir ");
      mkdir = mkdir + _filePath;
      system(mkdir.c_str());
    }

    return true;
  }

  // mix up the vertices in space
  void scrambleVertices()
  {
    std::mt19937 gen(123);
    std::uniform_real_distribution<REAL> rng(-15.0, 15.0); 

    vector<VECTOR3>& vertices = _strandMesh->vertices();
    for (unsigned int x = 0; x < vertices.size(); x++)
    {
      // skip constrained vertices
      if (_strandSolver->isKinematicallyConstrained(x)) continue;

      vertices[x][0] = rng(gen);
      vertices[x][1] = rng(gen) + 15;
      vertices[x][2] = rng(gen);
    }
    _strandMesh->updateProperties();
  }

  // mix up the vertices in space
  void smashVertices()
  {
    vector<VECTOR3>& vertices = _strandMesh->vertices();

    for (unsigned int x = 0; x < vertices.size(); x++)
    {
      // skip constrained vertices
      if (_strandSolver->isKinematicallyConstrained(x)) continue;

      vertices[x][0] *= _smash;
      //vertices[x][1] *= smash;
      //vertices[x][2] *= _smash;
    }
    _strandMesh->updateProperties();
  }

  // Type O hair on LOIS scale
  void buildType4C()
  {
    _drawVertices = false;

    vector<VECTOR3> restVertices;
    vector<vector<int> > strands;

    const REAL frequency = 1;
    const REAL amplitude = 0.5;

    //_spacing = 0.714621;
    //_spacing = 1.0;
    _spacing = 0.1;

    const REAL radius = -0.25;
    const int stemPoints = 1;

    VECTOR3 scalpVertex;

    restVertices.push_back(VECTOR3(-0.1, radius - stemPoints, 0.0));
    for (unsigned int x = 0; x < _totalPoints; x++)
    {
      if (x <= stemPoints)
      {
        restVertices.push_back(VECTOR3(0, x + radius - stemPoints, 0.0));
        if (x == stemPoints)
          scalpVertex = restVertices.back();
        continue;
      }

      const REAL yReal = _spacing * (x - stemPoints + 1) + radius;
      const REAL theta = (_spacing * (x - stemPoints + 1) * 4.0 * M_PI);
      const REAL xReal = -cos(frequency * theta) * amplitude;
      const REAL zReal = sin(frequency * theta) * amplitude;

      restVertices.push_back(VECTOR3(xReal, yReal, zReal));
    }

    const REAL curlLength = (scalpVertex - restVertices.back()).norm();
    cout << " Curl length: " << curlLength << endl;

    vector<int> strand0;
    for (unsigned int x = 0; x < _totalPoints; x++)
      strand0.push_back(x);
    strands.push_back(strand0);
      
    //_strandMesh = new STRAND_MESH(restVertices, strands,
    //                              _E, _G, _density, _radiusA, _radiusB);
    _strandMesh = new WISP_MESH(restVertices, strands,
                                _E, _G, _density, _baseRadius, _tipRadius);
  }

  // simulation loop
  virtual void stepSimulation(const bool verbose = true) override
  {
    // write out the initial frame
    if (_writeToFile && _frameNumber == 0)
    {
      char buffer[256];
      snprintf(buffer,256, "%04i", _frameNumber);
      string filename = _filePath + string("initial.sobj");

      cout << " Writing file " << filename.c_str() << " ... " << flush;
      writeSOBJFile(filename.c_str(), _strandMesh->vertices(), _strandMesh->strandIndices());
      cout << "done." << endl;

      TIMER gzipper("Gzipping");
      cout << " Gzipping ..." << flush;
      string gzip = string("gzip -f -9 ") + filename;
      system(gzip.c_str());
      cout << "done." << endl;
      gzipper.stop();
    }

    for (unsigned int x = 0; x < _substeps; x++)
    {
      _strandSolver->externalForces().setZero();
      _strandSolver->addGravity(_gravity);
      _strandSolver->solveDynamics(verbose);
      //_strandSolver->solveDynamics(false);
    }

    //if (_frameNumber > 0 && (_frameNumber % 10 == 0))
    if (_frameNumber > 0 && (_frameNumber % 5 == 0))
    {
      //TIMER::printTimings();
      TIMER::printTimingsPerFrame(_frameNumber);
    }

    if (_writeToFile)
    {
      char buffer[256];
      snprintf(buffer,256, "%04i", _frameNumber);
      string frameName = string("_frame_") + string(buffer) + string(".sobj");
      string filename = _filePath + _sceneName + frameName;

      cout << " Writing file " << filename.c_str() << " ... " << flush;
      writeSOBJFile(filename.c_str(), _strandMesh->vertices(), _strandMesh->strandIndices());
      cout << "done." << endl;

      TIMER gzipper("Gzipping");
      cout << " Gzipping ..." << flush;
      string gzip = string("gzip -f -9 ") + filename;
      system(gzip.c_str());
      cout << "done." << endl;
      gzipper.stop();
    }

    _frameNumber++;
  };

// drawScene MUST have this guard, so it can be deactivated for regression testing
#ifndef GL_DISABLED
  // what to draw?
  virtual void drawScene()
  {
    drawAxes();
    //drawStrand(*_strandMesh, 250);
    //drawStrandMeshOld(*_strandMesh, 250);
    //drawStrandMesh(*_strandMesh, -1, true);

    WISP_MESH* wisp = (WISP_MESH*)_strandMesh;
    drawWispMesh(*wisp, -1, _drawVertices);
    //drawStrandMeshOld(*_strandMesh);

    glEnable(GL_DEPTH_TEST);

    // only draw the collision sphere
    for (unsigned int x = 0; x < _kinematicShapes.size(); x++)
      drawKinematicShape(*_kinematicShapes[x]);
    
    //drawKinematicShape(*_kinematicShapes.back());

    //drawCollisionsOld(*_strandMesh);
    //drawCollisions(*_strandMesh);

    //drawKinematicConstraints(_strandMesh, _strandSolver);
    //drawPlaneConstraints(_strandMesh, _strandSolver);
    
  };
#endif

protected:
  //STRAND_MESH* _strandMesh;
  STRAND_MESH* _strandMesh;
  STRAND::TIMESTEPPER* _strandSolver;
  STRAND::STRETCHING* _strechingEnergy;

  unsigned int _substeps;

  string _filePath;
  bool _writeToFile;

  // GL drawing params
  bool _drawVertices;

  // strand generation parameters
  unsigned int _totalPoints;
  REAL _spacing;
  REAL _E;
  REAL _G;
  REAL _density;
  REAL _baseRadius;
  REAL _tipRadius;

  // how much to smash?
  REAL _smash;

  REAL _dt;
};

}

#endif

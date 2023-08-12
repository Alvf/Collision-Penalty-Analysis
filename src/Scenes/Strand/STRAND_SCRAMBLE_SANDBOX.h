#ifndef STRAND_SCRAMBLE_SANDBOX_H
#define STRAND_SCRAMBLE_SANDBOX_H

#include "Scenes/SIMULATION_SCENE.h"
#include "Geometry/TET_WISP_MESH.h"
#include "Timestepper/Strand/VOLUME_TIMESTEPPER.h"
#include "Timestepper/Strand/VOLUME_BDF_1.h"
#include "Hyperelastic/Strand/QUADRATIC_STRETCHING.h"
#include <random>
#include <float.h>

namespace HOBAK {

class STRAND_SCRAMBLE_SANDBOX : public SIMULATION_SCENE {
public:
  virtual void printSceneDescription() override
  {
    cout << "=====================================================================" << endl;
    cout << " Strand scramble scene" << endl;
    cout << "=====================================================================" << endl;
  }

  STRAND_SCRAMBLE_SANDBOX() 
  {
    //_E = 3.9e9;
    //_E = 3.9e7;
    _E = 3.9e4;
    REAL nu = 0.48;
    _G = _E / (2.0 * (1 + nu));
    _density = 1.3;
    //_baseRadius = 0.0037;
    //_tipRadius = 0.0037;
    _baseRadius = 0.2 * 0.5;
    _tipRadius = 0.2 * 0.5;
    
    //_jitter = 0.1; // stable
    //_jitter = 0.11; // stable at dt = 1 / 300
    //_jitter = 0.2; // stable at dt = 1 / 300
    //_jitter = 0.5; // stable at 1 / 3000
    _jitter = 1.0; // stable at 1 / 3000
    
    //_writeToFile = true;
    _writeToFile = false;
  };

  // for testing stretch
  virtual bool buildScene() override
  {
    _autoplay = false;
    
    // this will determine the MOV and JSON filenames
    char buffer[512];
    snprintf(buffer,256, "jitter=%4.2e", _jitter);
    string jitterString(buffer);
    _sceneName = "scramble_sandbox_" + jitterString;

    //const REAL rayleighAlpha = 0.0;
    //const REAL rayleighBeta = 0.0;
    //const REAL rayleighAlpha = 0.001;
    //const REAL rayleighBeta = 0.001;
    //const REAL rayleighAlpha = 0.0025;
    //const REAL rayleighBeta = 0.0025;

    //const int maxNewtonIterations = 10; // flickering stabilizes
    //const int maxNewtonIterations = 3;  // works with warm start

    REAL twistScale = 1.0;
    REAL bendScale = 1.0;
    REAL stretchScale = 1.0;

    // stable up to 2000
#if 0
    const int maxNewtonIterations = 1; 
    const REAL rayleighAlpha = 0.0;
    const REAL rayleighBeta = 0.0;
    _totalPoints = 10;
    _spacing = 0.1;   // stable 
    _coilHeight = 1.0;
    _gravity = VECTOR3(0,0,0);
#endif

    // stable up to 2000
#if 0
    const int maxNewtonIterations = 1; 
    const REAL rayleighAlpha = 0.0;
    const REAL rayleighBeta = 0.0;
    _totalPoints = 10;
    _spacing = 0.1;   // stable 
    _coilHeight = 1.0;
    _gravity = VECTOR3(0,-981,0);
#endif

    // stable up to 2000, long helix
#if 0
    const int maxNewtonIterations = 1; 
    const REAL rayleighAlpha = 0.0;
    const REAL rayleighBeta = 0.0;
    _totalPoints = 100;
    _spacing = 0.1;   // stable 
    _coilHeight = 1.0;
    _gravity = VECTOR3(0,-981,0);
#endif

    // stable up to 2000, long helix
#if 0
    const int maxNewtonIterations = 1; 
    const REAL rayleighAlpha = 0.0;
    const REAL rayleighBeta = 0.0;
    _totalPoints = 100;
    _spacing = 0.09; // NEW
    _coilHeight = 1.0;
    _gravity = VECTOR3(0,-981,0);
#endif

    // pops around frame 500, long helix
#if 0
    const int maxNewtonIterations = 1; 
    const REAL rayleighAlpha = 0.0;
    const REAL rayleighBeta = 0.0;
    _totalPoints = 100;
    _spacing = 0.085; // NEW
    _coilHeight = 1.0;
    _gravity = VECTOR3(0,-981,0);
#endif

    // unstable
#if 0
    const int maxNewtonIterations = 1; 
    const REAL rayleighAlpha = 0.0;
    const REAL rayleighBeta = 0.0;
    _totalPoints = 100;
    _spacing = 0.08; // NEW
    _coilHeight = 1.0;
    _gravity = VECTOR3(0,-981,0);
#endif

    // vibrates slightly
#if 0
    const int maxNewtonIterations = 1; 
    const REAL rayleighAlpha = 0.0;
    const REAL rayleighBeta = 0.0;
    _totalPoints = 10;
    _spacing = 0.1;   // stable 
    _coilHeight = 1e-15;
    _gravity = VECTOR3(0,0,0);
#endif

    // stable up to 2000
#if 0
    const int maxNewtonIterations = 1; 
    const REAL rayleighAlpha = 0.0025;
    const REAL rayleighBeta = 0.0025;
    _totalPoints = 7;
    _spacing = 0.075;
    _coilHeight = 1;
    _gravity = VECTOR3(0,0,0);
#endif

    // unstable from the beginning
#if 0
    const int maxNewtonIterations = 1; 
    const REAL rayleighAlpha = 0.0025;
    const REAL rayleighBeta = 0.0025;
    //_totalPoints = 7;
    _totalPoints = 15;
    _spacing = 0.075;
    _coilHeight = 2;
    _gravity = VECTOR3(0,-981,0);
#endif

    // unstable at ~900 steps
#if 0
    const int maxNewtonIterations = 1; 
    const REAL rayleighAlpha = 0.001;
    const REAL rayleighBeta = 0.001;
    _totalPoints = 7;
    _spacing = 0.075;
    _coilHeight = 1;
    _gravity = VECTOR3(0,0,0);
#endif

    // unstable around 200 steps
#if 0
    const int maxNewtonIterations = 1; 
    const REAL rayleighAlpha = 0.0;
    const REAL rayleighBeta = 0.0;
    _totalPoints = 7;
    _spacing = 0.075;
    _coilHeight = 1;
    _gravity = VECTOR3(0,0,0);
#endif

    // current best unstable test
#if 0
    const int maxNewtonIterations = 1; 
    //const REAL rayleighAlpha = 0.0;
    //const REAL rayleighBeta = 0.0;
    const REAL rayleighAlpha = 0.0;
    //const REAL rayleighBeta = 0.01; // stable
    //const REAL rayleighBeta = 0.001;  // unstable
    //const REAL rayleighBeta = 0.0025; // unstable
    //const REAL rayleighBeta = 0.005; // unstable
    //const REAL rayleighBeta = 0.0075; // unstable
    const REAL rayleighBeta = 0.009; // unstable
    _totalPoints = 7;
    _spacing = 0.05;
    _coilHeight = 1;
    _gravity = VECTOR3(0,0,0);

    //twistScale = 0.01;
    //twistScale = 1e-8;
    twistScale = 0;
    bendScale = 0;
    stretchScale = 1.0;
#endif

    // trying a stable test
#if 1
    const int maxNewtonIterations = 1; 
    //const REAL rayleighAlpha = 0.0;
    const REAL rayleighAlpha = 0.01;
    //const REAL rayleighBeta = 0.0; // unstable
    const REAL rayleighBeta = 0.001; // stable
    //const REAL rayleighBeta = 0.01; // stable
    _totalPoints = 7;
    _spacing = 0.1;
    //_coilHeight = 1;    // stable
    //_coilHeight = 0.1; // stable
    //_coilHeight = 0.01; // stable
    //_coilHeight = 0.001; // stable
    //_coilHeight = 1e-6; // stable
    _coilHeight = 1e-8; // stable (!!)
    //_gravity = VECTOR3(0,0,0);
    _gravity = VECTOR3(0,-981,0);

    //twistScale = 0.01;
    //twistScale = 1e-8;
    twistScale = 1.0;
    bendScale = 1.0;
    stretchScale = 1.0;
#endif

    _drawVertices = true;
    _totalStrands = 1;

    buildCoil();

    _strandMesh->scaleTwistingMu(twistScale);
    _strandMesh->scaleBendingMu(bendScale);
    _strandMesh->scaleStretchingMu(stretchScale);

    cout << " Strand length: " << _strandMesh->strandLength(0) << endl;

    const REAL k = 1.0;
    STRAND::QUADRATIC_STRETCHING* stretchingEnergy = new STRAND::QUADRATIC_STRETCHING(k);
    setUnsignedCollisions(1000.0,0.1); 
    _strandSolver = new STRAND::VOLUME_BDF_1(*_strandMesh, *stretchingEnergy, *_eeGeneral);
    
    _strandSolver->maxNewtonIterations() = maxNewtonIterations;

    _substeps = 1;
    const REAL dt = 1.0 / 30.0;
    //const REAL dt = 1.0 / 3000.0;
    _strandSolver->setDt(dt);

    // pin both ends of the strand
    //const int lastIndex = _strandMesh->vertices().size() - 1;
    //const REAL radius = 0.01;
    const REAL radius = 0.15;
    VECTOR3 sphereCenter;
    sphereCenter = _strandMesh->vertices()[0];
    _kinematicShapes.push_back(new SPHERE(sphereCenter, radius));
    _strandSolver->bindVertexToKinematic(_kinematicShapes.back(), 0);

    sphereCenter = _strandMesh->vertices()[1];
    _kinematicShapes.push_back(new SPHERE(sphereCenter, radius));
    _strandSolver->bindVertexToKinematic(_kinematicShapes.back(), 1);

#if 0
    sphereCenter = _strandMesh->vertices()[lastIndex - 1];
    _kinematicShapes.push_back(new SPHERE(sphereCenter, radius));
    _strandSolver->bindVertexToKinematic(_kinematicShapes.back(), lastIndex - 1);

    sphereCenter = _strandMesh->vertices()[lastIndex];
    _kinematicShapes.push_back(new SPHERE(sphereCenter, radius));
    _strandSolver->bindVertexToKinematic(_kinematicShapes.back(), lastIndex);
#endif

    sphereCenter = VECTOR3(0,-1.5,1.5);
    _kinematicShapes.push_back(new SPHERE(sphereCenter, 1.0));
    _strandSolver->addKinematicCollisionObject(_kinematicShapes.back());

    cout << " Total constrained nodes: " << _strandSolver->constrainedNodes().size() << endl;

    // apply the scramble
    //jitterVertices();

    // make the sphere a collision object
    //_strandSolver->addKinematicCollisionObject(_kinematicShapes.back());

    //_strandSolver->collisionsEnabled() = true;
    _strandSolver->collisionsEnabled() = false;
    _strandSolver->pcgEnabled() = false;
    //_strandSolver->pcgEnabled() = true;
    _strandSolver->hessianClampingEnabled() = true;
    //_strandSolver->hessianClampingEnabled() = false;
    //_strandSolver->collisionStiffness() = 1000;
    _strandSolver->collisionDampingBeta() = 0;
    _strandSolver->edgeEdgeSelfCollisionsOn() = false;

    _strandSolver->setRayeligh(rayleighAlpha, rayleighBeta);

    _strandMesh->bendingForceFilterEnabled() = false;
    //_strandMesh->bendingForceFilterEnabled() = true;

    _eye    = VECTOR3(0.321923, 2.62036, 5.06155);
    _lookAt = VECTOR3(0.257654, 2.17211, 4.16996);
    _up     = VECTOR3(-0.00502089, 0.893574, -0.448882);

    _worldCenter = VECTOR3(0,0,0);

    //_pauseFrame = 41;
    //_pauseFrame = 2000;
    _pauseFrame = 1000;
    _exitOnPause = false;
    //_pauseFrame = 1000;
    //_pauseFrame = -1;

    _filePath = string("../data/render/");

    snprintf(buffer,256, "E=%4.2e_G=%4.2e", _E, _G);
    string Estring(buffer);
    snprintf(buffer,256, "r0=%.3f_r1=%.3f", _baseRadius, _tipRadius);
    string Rstring(buffer);

    //_sceneName = _sceneName + string("_") + Sstring + string("_") + Estring + string("_") + Rstring;
    _sceneName = _sceneName + string("_") + Estring + string("_") + Rstring;
    cout << " Scene name: " << _sceneName.c_str() << endl;

    _filePath = string("../data/render/") + _sceneName + string("/");
    string mkdir("mkdir ");
    mkdir = mkdir + _filePath;
    system(mkdir.c_str());

    return true;
  }

  // jitter the vertices in space
  void jitterVertices()
  {
    std::mt19937 gen(123);
    std::uniform_real_distribution<REAL> rng(-_jitter, _jitter);

    vector<VECTOR3>& vertices = _strandMesh->vertices();
    for (unsigned int x = 0; x < vertices.size(); x++)
    {
      // skip constrained vertices
      if (_strandSolver->isKinematicallyConstrained(x)) continue;

      vertices[x][0] += rng(gen);
      vertices[x][1] += rng(gen);
      vertices[x][2] += rng(gen);
    }
    _strandMesh->updateProperties();
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
  }

  // build a single helix
  void buildCoil()
  {
    _drawVertices = false;

    vector<VECTOR3> restVertices;
    vector<vector<int> > strands;

    const REAL frequency = 1;
    const REAL amplitude = 0.5;
    const REAL radius = -0.25;

    //const REAL straighten = 1.0; // original
    //const REAL straighten = 2.0;
    //const REAL straighten = 4.0;  // stable
    //const REAL straighten = 8.0;  // stable
    //const REAL straighten = 16.0;  // stable
    //const REAL straighten = 32.0;  // stable
    //const REAL straighten = 128.0; // stable
    //const REAL straighten = 1024.0; // stable
    //const REAL straighten = 8 * 1024.0; // stable
    //const REAL straighten = 32 * 1024.0; // stable
    const REAL straighten = 256* 1024* 1024.0; // stable
    //const REAL straighten = 1024 * 1024* 1024.0; // dies!

    //const int stemPoints = 1;
    //const int stemPoints = 0;

    for (unsigned int x = 0; x < _totalPoints; x++)
    {
      const REAL yReal = _coilHeight * _spacing * (x + 1) + radius;
      const REAL theta = _spacing * (x + 1) * 4.0 * M_PI / straighten;
      const REAL xReal = -cos(frequency * theta) * amplitude * straighten;
      const REAL zReal = sin(frequency * theta) * amplitude * straighten;

      restVertices.push_back(VECTOR3(xReal, yReal, zReal));
    }

    // recenter everything at the origin
    const VECTOR3 original = restVertices[0];
    for (unsigned int x = 0; x < restVertices.size(); x++)
      restVertices[x] -= original;

    vector<int> strand0;
    for (unsigned int x = 0; x < restVertices.size(); x++)
      strand0.push_back(x);
    strands.push_back(strand0);
      
    _strandMesh = new TET_WISP_MESH(restVertices, strands,
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
      STRAND_MESH::writeSOBJFile(filename.c_str(), _strandMesh->vertices(), _strandMesh->strandIndices());
      cout << "done." << endl;

      TIMER gzipper("Gzipping");
      cout << " Gzipping ..." << flush;
      string gzip = string("gzip -f -9 ") + filename;
      system(gzip.c_str());
      cout << "done." << endl;
      gzipper.stop();
    }

    //const REAL time = _frameNumber * _strandSolver->dt();
    for (unsigned int x = 0; x < _substeps; x++)
    {
      _strandSolver->externalForces().setZero();
      _strandSolver->addGravity(_gravity);
      _strandSolver->solveDynamics(verbose);
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
      STRAND_MESH::writeSOBJFile(filename.c_str(), _strandMesh->vertices(), _strandMesh->strandIndices());
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
    glEnable(GL_DEPTH_TEST);
    drawAxes();
    TET_WISP_MESH* wisp = (TET_WISP_MESH*)_strandMesh;
    drawWispMesh(*wisp, -1, _drawVertices);

    // only draw the collision sphere
    for (unsigned int x = 0; x < _kinematicShapes.size(); x++)
      drawKinematicShape(*_kinematicShapes[x]);
    
    //drawKinematicShape(*_kinematicShapes.back());

    //drawCollisionsOld(*_strandMesh);
    //drawCollisions(*_strandMesh);

    //drawKinematicConstraints(_strandMesh, _strandSolver);
    //drawPlaneConstraints(_strandMesh, _strandSolver);
 
   /* 
    // draw some random points on the sphere
    glDisable(GL_DEPTH_TEST);
    glColor4f(1.0, 0,0,1);
    glPointSize(10);
    glBegin(GL_POINTS);
      for (unsigned int x = 0; x < _randomPoints.size(); x++)
        glVertex3f(_randomPoints[x][0], _randomPoints[x][1], _randomPoints[x][2]);
    glEnd(); 
    */
  };
#endif

protected:
  //STRAND_MESH* _strandMesh;
  TET_WISP_MESH* _strandMesh;
  STRAND::TIMESTEPPER* _strandSolver;
  STRAND::STRETCHING* _strechingEnergy;

  unsigned int _substeps;

  string _filePath;
  bool _writeToFile;

  // GL drawing params
  bool _drawVertices;

  // strand generation parameters
  unsigned int _totalPoints;
  unsigned int _totalStrands;
  REAL _spacing;
  REAL _density;
  REAL _baseRadius;
  REAL _tipRadius;
  REAL _coilHeight;

  vector<VECTOR3> _randomPoints;

  REAL _jitter;
};

}

#endif

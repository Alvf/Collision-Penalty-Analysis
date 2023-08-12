#ifndef STRAND_SCRAMBLE_H
#define STRAND_SCRAMBLE_H

#include "Scenes/SIMULATION_SCENE.h"
#include "Geometry/TET_WISP_MESH.h"
#include "Timestepper/Strand/VOLUME_TIMESTEPPER.h"
#include "Timestepper/Strand/VOLUME_BDF_1.h"
#include "Hyperelastic/Strand/QUADRATIC_STRETCHING.h"
#include <random>
#include <float.h>

namespace HOBAK {

class STRAND_SCRAMBLE : public SIMULATION_SCENE {
public:
  virtual void printSceneDescription() override
  {
    cout << "=====================================================================" << endl;
    cout << " Strand scramble scene" << endl;
    cout << "=====================================================================" << endl;
  }

  STRAND_SCRAMBLE() 
  {
    _E = 3.9e9;
    REAL nu = 0.48;
    _G = _E / (2.0 * (1 + nu));
    //_G = _E;
    _density = 1.3;
    //_baseRadius = 0.0037;
    //_tipRadius = 0.0037;
    _baseRadius = 0.2 * 0.5;
    _tipRadius = 0.2 * 0.5;
    
    //_jitter = 0.05; // stable
    //_jitter = 0.1; // stable
    //_jitter = 0.11; // stable at dt = 1 / 300
    //_jitter = 0.2; // stable at dt = 1 / 300
    _jitter = 0.5; // stable at 1 / 3000
    //_jitter = 1.0; // stable at 1 / 3000
    
    _writeToFile = true;
    //_writeToFile = false;
  };

  // for testing stretch
  virtual bool buildScene() override
  {
    _autoplay = false;

    
    // this will determine the MOV and JSON filenames
    char buffer[512];
    snprintf(buffer,256, "jitter=%4.2e", _jitter);
    string jitterString(buffer);
    _sceneName = "lifted_curls_" + jitterString;

    //_totalPoints = 50;
    _spacing = 1.0;

    //_gravity = VECTOR3(0,-981,0);
    _gravity = VECTOR3(0,0,0);

    _totalPoints = 100; // decent for type 4C

    _drawVertices = true;

    //_totalStrands = 15;
    //_totalStrands = 100;
    //_totalStrands = 1000;
    _totalStrands = 1;

    buildType4C();

    cout << " Strand length: " << _strandMesh->strandLength(0) << endl;

    const REAL k = 1.0;
    STRAND::QUADRATIC_STRETCHING* stretchingEnergy = new STRAND::QUADRATIC_STRETCHING(k);
    //_strandSolver = new STRAND::VOLUME_TIMESTEPPER(*_strandMesh, *stretchingEnergy);
    setUnsignedCollisions(1000.0,0.1); 
    _strandSolver = new STRAND::VOLUME_BDF_1(*_strandMesh, *stretchingEnergy, *_eeGeneral);
    
    //const int maxNewtonIterations = 10; // flickering stabilizes
    const int maxNewtonIterations = 3;  // works with warm start
    _strandSolver->maxNewtonIterations() = maxNewtonIterations;

    _substeps = 1;
    //const REAL dt = 1.0 / 3000.0;
    //const REAL dt = 1.0 / 300.0;
    const REAL dt = 1.0 / 30.0;
    _strandSolver->setDt(dt);

    const REAL radius = 15.25;
    VECTOR3 sphereCenter0(0,-radius,0); 
    _kinematicShapes.push_back(new SPHERE(sphereCenter0, radius));

    // attach hairs to the sphere
    for (unsigned int x = 0; x < _kinematicShapes.size(); x++)
      _strandSolver->attachKinematicSurfaceConstraints(_kinematicShapes[x]);

    // apply the scramble
    //scrambleVertices();
    //smashVertices();
    jitterVertices();

    // bump up the stiffness of attached edges to reflect the rigid attachment
    //((STRAND::VOLUME_TIMESTEPPER*)_strandSolver)->stiffenKinematicStrands();

    // make the sphere a collision object
    //_strandSolver->addKinematicCollisionObject(_kinematicShapes.back());

    //_strandSolver->collisionsEnabled() = true;
    _strandSolver->collisionsEnabled() = false;
    //_strandSolver->pcgEnabled() = false;
    _strandSolver->pcgEnabled() = true;
    _strandSolver->hessianClampingEnabled() = true;
    //_strandSolver->hessianClampingEnabled() = false;
    //_strandSolver->collisionStiffness() = 1000;
    _strandSolver->collisionDampingBeta() = 0;

    _strandMesh->bendingForceFilterEnabled() = false;
    //_strandMesh->bendingForceFilterEnabled() = true;

    _eye    = VECTOR3(16.6521, 13.8222, 42.1509);
    _lookAt = VECTOR3(16.3798, 13.7082, 41.1954);
    _up     = VECTOR3(0.0296061, 0.991492, -0.126738);

    _eye    = VECTOR3(6.351, 11.607, 59.9579);
    _lookAt = VECTOR3(6.22913, 11.4656, 58.9755);
    _up     = VECTOR3(0.0176662, 0.989338, -0.144556);

    _eye    = VECTOR3(3.245, 19.1677, 47.9845);
    _lookAt = VECTOR3(3.15827, 19.0252, 46.9985);
    _up     = VECTOR3(0.0176633, 0.989333, -0.144589);
    
    _eye    = VECTOR3(3.21116, 6.56164, 17.5514);
    _lookAt = VECTOR3(3.06594, 6.42116, 16.572);
    _up     = VECTOR3(0.017671, 0.989341, -0.144537);

    _worldCenter = VECTOR3(0,0,0);

    //_pauseFrame = 50;// okay for initial tests
    _pauseFrame = 200;
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

  // mix up the vertices in space
  void smashVertices()
  {
    vector<VECTOR3>& vertices = _strandMesh->vertices();
    for (unsigned int x = 0; x < vertices.size(); x++)
    {
      // skip constrained vertices
      if (_strandSolver->isKinematicallyConstrained(x)) continue;

      vertices[x][0] *= 1e-5;
      vertices[x][1] *= 1e-5;
      vertices[x][2] *= 1e-5;
    }
  }

  // Type O hair on LOIS scale
  void buildType4C()
  {
    _drawVertices = false;

    vector<VECTOR3> restVertices;
    vector<vector<int> > strands;

    //const REAL straightLength = _spacing * _totalPoints;

    const REAL frequency = 1;
    const REAL amplitude = 0.5;

    //_spacing = 0.714621;
    //_spacing = 1.0;
    _spacing = 0.1;

    const REAL radius = -0.25;

    const int stemPoints = 1;

    VECTOR3 scalpVertex;
    scalpVertex.setZero();

    //restVertices.push_back(VECTOR3(-1, radius - stemPoints, 0.0));
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
      //const REAL xReal = cos(frequency * theta) * amplitude;  // wrong ordering
      const REAL zReal = sin(frequency * theta) * amplitude;

      restVertices.push_back(VECTOR3(xReal, yReal, zReal));
    }

    const REAL curlLength = (scalpVertex - restVertices.back()).norm();
    cout << " Curl length: " << curlLength << endl;

    vector<int> strand0;
    for (unsigned int x = 0; x < restVertices.size(); x++)
      strand0.push_back(x);
    strands.push_back(strand0);
      
    _strandMesh = new TET_WISP_MESH(restVertices, strands,
                                    _E, _G, _density, _baseRadius, _tipRadius);
                                    //_E, _nu, _density, _baseRadius, _tipRadius);
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
#if 0
    drawStrandMesh(*_strandMesh, -1, true);
#else
    TET_WISP_MESH* wisp = (TET_WISP_MESH*)_strandMesh;
    drawWispMesh(*wisp, -1, _drawVertices);
#endif

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

  vector<VECTOR3> _randomPoints;

  REAL _jitter;
};

}

#endif

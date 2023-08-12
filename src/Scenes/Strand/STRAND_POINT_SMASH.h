#ifndef STRAND_POINT_SMASH_H
#define STRAND_POINT_SMASH_H

#include "Scenes/SIMULATION_SCENE.h"
#include "Geometry/TET_WISP_MESH.h"
#include "Timestepper/Strand/VOLUME_TIMESTEPPER.h"
#include "Timestepper/Strand/VOLUME_BDF_1.h"
#include "Hyperelastic/Strand/QUADRATIC_STRETCHING.h"
#include <random>
#include <float.h>

namespace HOBAK {

class STRAND_POINT_SMASH : public SIMULATION_SCENE {
public:
  virtual void printSceneDescription() override
  {
    cout << "=====================================================================" << endl;
    cout << " Strand scramble scene" << endl;
    cout << "=====================================================================" << endl;
  }

  STRAND_POINT_SMASH() {
    _E = 3.9e9;
    //_G = _E;
    REAL nu = 0.48;
    _G = _E / (2.0 * (1 + nu));
    _density = 1.3;
    //_baseRadius = 0.0037;
    //_tipRadius = 0.0037;
    //_baseRadius = 0.01;
    //_tipRadius = 0.01;
    _baseRadius = 0.1;
    _tipRadius = 0.1;
    
    //_smash = 0.75; // this recovers
    //_smash = 8.0e-1; // this recovers
    //_smash = 7.5e-1; // this recovers
    //_smash = 6.0e-1; // this diverges
    //_smash = 1e-5;  // diverges
    //_smash = 6e-1;  // diverges
    //_smash = 1e-1;  // diverges, Gauss-Newton diverges (DVT)
    
    _smash = 7e-1;  // recovers (DVT)
    _smash = 6e-1;  // diverges, Gauss-Newton recovers (DVT)
    _smash = 4e-1;  // diverges, Gauss-Newton diverges (DVT)
    
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
    _sceneName = "lifted_curls_" + smashString;

    _spacing = 1.0;

    //_gravity = VECTOR3(0,-981,0);
    _gravity = VECTOR3(0,0,0);

    _totalPoints = 100; // decent for type 4C
    //_totalPoints = 5; // DEBUG

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
    const REAL dt = 1.0 / 30.0 / _substeps;
    _strandSolver->setDt(dt);

    const REAL radius = 15.25;
    VECTOR3 sphereCenter0(0,-radius,0); 
    _kinematicShapes.push_back(new SPHERE(sphereCenter0, radius));

    // attach hairs to the sphere
    for (unsigned int x = 0; x < _kinematicShapes.size(); x++)
      _strandSolver->attachKinematicSurfaceConstraints(_kinematicShapes[x]);

    // apply the scramble
    //scrambleVertices();
    smashVertices();

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

    _eye    = VECTOR3(3.21116, 6.56164, 17.5514);
    _lookAt = VECTOR3(3.06594, 6.42116, 16.572);
    _up     = VECTOR3(0.017671, 0.989341, -0.144537);

    _worldCenter = VECTOR3(0,0,0);

    _pauseFrame = 100;
    //_pauseFrame = 200;
    //_pauseFrame = 1000;
    //_pauseFrame = -1;

    snprintf(buffer,256, "E=%4.2e", _E);
    string Estring(buffer);
    snprintf(buffer,256, "r0=%.3f_r1=%.3f", _baseRadius, _tipRadius);
    string Rstring(buffer);

    const REAL damping = _strandSolver->rayleighBeta();
    snprintf(buffer,256, "damping=%.3f", damping);
    string Dstring(buffer);

    int totalStrands = _strandMesh->totalStrands();
    cout << " total strands: " << totalStrands << endl;

    snprintf(buffer,256, "s=%d", _totalStrands);
    string Sstring(buffer);

    //_sceneName = _sceneName + string("_") + Sstring + string("_") + Estring + string("_") + Rstring;
    //_sceneName = _sceneName + string("_") + Sstring + string("_") + Estring + string("_") + Rstring + string("_") + Dstring;
    _sceneName = _sceneName + string("_") + Estring + string("_") + Rstring;
    cout << " Scene name: " << _sceneName.c_str() << endl;

    _filePath = string("../data/render/") + _sceneName + string("/");
    string mkdir("mkdir ");
    mkdir = mkdir + _filePath;
    system(mkdir.c_str());

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
  }

  // mix up the vertices in space
  void smashVertices()
  {
    // write out an unsmashed version
    if (_writeToFile)
    {
      string filename = _filePath + string("unsmashed.sobj");

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

    vector<VECTOR3>& vertices = _strandMesh->vertices();
    for (unsigned int x = 0; x < vertices.size(); x++)
    {
      // skip constrained vertices
      if (_strandSolver->isKinematicallyConstrained(x)) continue;

      vertices[x][0] *= _smash;
      //vertices[x][1] *= _smash;
      //vertices[x][2] *= _smash;
    }
  }

  // Type O hair on LOIS scale
  void buildType4C()
  {
    _drawVertices = false;

    vector<VECTOR3> restVertices;
    vector<vector<int> > strands;

    const REAL frequency = 1;
    const REAL amplitude = 0.5;

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
  
  // how much to smash?
  REAL _smash;
};

}

#endif

#ifndef DVT_SCRAMBLE_H
#define DVT_SCRAMBLE_H

#include "Scenes/SIMULATION_SCENE.h"
#include "Geometry/WISP_MESH.h"
#include "Geometry/STRAND_MESH_FASTER.h"
#include "Timestepper/Strand/TIMESTEPPER.h"
#include "Timestepper/Strand/BDF_1.h"
#include "Hyperelastic/Strand/QUADRATIC_STRETCHING.h"

namespace HOBAK {

class DVT_SCRAMBLE : public SIMULATION_SCENE {
public:
  virtual void printSceneDescription() override
  {
    cout << "=====================================================================" << endl;
    cout << " Hair wisp scene" << endl;
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

  DVT_SCRAMBLE() 
  {
    _E = 3.9e9;
    //_E = 1e3;
    //_G = _E;
    REAL nu = 0.48;
    _G = _E / (2.0 * (1 + nu));
    _density = 1.3; // 5% of the usual density
    
    _baseRadius = 0.2 * 0.5;
    _tipRadius = 0.2 * 0.5;
    
    _substeps = 1;
    
    _dt = 1.0 / 30.0;
    //_dt = 1.0 / 300.0; // GNF, stable for jitter = 0.2
    //_dt = 1.0 / 3000.0; // GNF, stable for jitter = 0.5
    //_dt = 1.0 / 30000.0;
    //_dt = 1.0 / 300000.0;// stable for full filter, jitter = 0.2
    //_dt = 1.0 / 3000000.0;// stable for full filter, jitter = 0.5

    //_jitter = 0.1; // stable
    //_jitter = 0.05; // stable
    //_jitter = 0.075; // unstable full, stable GNF
    //_jitter = 0.1; // unstable full, stable GNF
    //_jitter = 0.15; // stable GNF
    //_jitter = 0.2; // unstable
    
    //_jitter = 0.05; // stable
    //_jitter = 0.1; // unstable full, stable GNF
    _jitter = 0.2; // unstable
    //_jitter = 0.5; // unstable
    
    _writeToFile = true;
    //_writeToFile = false;
  };

  // for testing stretch
  virtual bool buildScene() override
  {
    _autoplay = false;

    // this will determine the MOV and JSON filenames
    char buffer[512];
    sprintf(buffer, "jitter=%4.2e", _jitter);
    string jitterString(buffer);
    _sceneName = "dvt_" + jitterString;

    //_totalPoints = 50;
    _spacing = 1.0;

    // based on photo of store-bought hair

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
    //smashVertices();
    //scrambleVertices();
    jitterVertices();
    //_strandSolver->position() = _strandMesh->getDisplacement();

    //_strandSolver->collisionsEnabled() = true;
    _strandSolver->collisionsEnabled() = false;
    //_strandSolver->pcgEnabled() = false;
    _strandSolver->pcgEnabled() = true;
    _strandSolver->hessianClampingEnabled() = true;
    //_strandSolver->hessianClampingEnabled() = false;
    _strandSolver->collisionStiffness() = 10;
    _strandSolver->collisionDampingBeta() = 0;
    _strandSolver->disablePreconditioner() = true;  // otherwise it makes a NaN

    _strandMesh->bendingForceFilterEnabled() = false;
    //_strandMesh->bendingForceFilterEnabled() = true;

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

    //_eye    = VECTOR3(16.6521, 13.8222, 42.1509);
    //_lookAt = VECTOR3(16.3798, 13.7082, 41.1954);
    //_up     = VECTOR3(0.0296061, 0.991492, -0.126738);

    _eye    = VECTOR3(3.21116, 6.56164, 17.5514);
    _lookAt = VECTOR3(3.06594, 6.42116, 16.572);
    _up     = VECTOR3(0.017671, 0.989341, -0.144537);
    
    _eye    = VECTOR3(3.245, 19.1677, 47.9845);
    _lookAt = VECTOR3(3.15827, 19.0252, 46.9985);
    _up     = VECTOR3(0.0176633, 0.989333, -0.144589);

    _eye    = VECTOR3(3.21116, 6.56164, 17.5514);
    _lookAt = VECTOR3(3.06594, 6.42116, 16.572);
    _up     = VECTOR3(0.017671, 0.989341, -0.144537);

    _worldCenter = VECTOR3(0,0,0);

    _pauseFrame = 50;
    //_pauseFrame = 200;
    //_pauseFrame = -1;

    _filePath = string("../data/render/");

    sprintf(buffer, "E=%4.2e", _E);
    string Estring(buffer);
    sprintf(buffer, "r0=%.3f_r1=%.3f", _baseRadius, _tipRadius);
    string Rstring(buffer);
    sprintf(buffer, "dt=%4.2e", _dt);
    string Tstring(buffer);

    _sceneName = _sceneName + string("_") + Estring + string("_") + Rstring + string("_") + Tstring;
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
    _strandMesh->updateProperties();
  }

  // mix up the vertices in space
  void smashVertices()
  {
    vector<VECTOR3>& vertices = _strandMesh->vertices();

    const REAL smash = 1e-5;
    //const REAL smash = 1e-10;

    for (unsigned int x = 0; x < vertices.size(); x++)
    {
      // skip constrained vertices
      if (_strandSolver->isKinematicallyConstrained(x)) continue;

      vertices[x][0] *= smash;
      vertices[x][1] *= smash;
      vertices[x][2] *= smash;
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
    //                              _E, _nu, _density, _radiusA, _radiusB);
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
      sprintf(buffer, "%04i", _frameNumber);
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
      sprintf(buffer, "%04i", _frameNumber);
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

  REAL _jitter;
  REAL _dt;
};

}

#endif

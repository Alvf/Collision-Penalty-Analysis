#ifndef WISP_SCENE_H
#define WISP_SCENE_H

#include "Scenes/SIMULATION_SCENE.h"
#include "Geometry/WISP_MESH.h"
#include "Geometry/STRAND_MESH_FASTER.h"
#include "Timestepper/Strand/TIMESTEPPER.h"
#include "Hyperelastic/Strand/QUADRATIC_STRETCHING.h"

namespace HOBAK {

class WISP_SCENE : public SIMULATION_SCENE {
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

  // for testing stretch
  virtual bool buildScene() override
  {
    // this will determine the MOV and JSON filenames
    _sceneName = "wisp";

    //_totalPoints = 50;
    _spacing = 1.0;
    //_E = 1e10; // default
    //_E = 5e4; // decent for type 1
    //_E = 5e4; // decent for type 4C
    _E = 2.5e4; // decent for sample match
    
    //_E = 1e4; // decent for type 4B
    //_nu = 0.36;
    _nu = 0.01;
    _density = 0.065; // 5% of the usual density

    // based on photo of A.M.'s hair, in centimeters
    //_baseRadius = 0.23;
    //_tipRadius = baseRadius / 5.0;

    // based on photo of store-bought hair
    _baseRadius = 1.0 * 0.5;
    //_baseRadius = 0.4 * 0.5;
    //_tipRadius = 0.4 * 0.5;
    _tipRadius = 0.2 * 0.5;
    //_tipRadius = 0.1 * 0.5;

    _gravity = VECTOR3(0,-981,0);
    //_gravity = VECTOR3(0,0,0);

    _totalPoints = 20;  // decent test for type 1 
    //_totalPoints = 50; // decent for type 4C
    _totalPoints = 100; // decent for type 4C
    //_totalPoints = 200; // decent for type 4C

    _drawVertices = true;

    // which type of hair?
    //buildType1();
    //buildType4B();
    
    //buildType4C();
    buildType4CAlternate();
    //buildSampleMatch();

    cout << " Strand length: " << _strandMesh->strandLength(0) << endl;

    const REAL k = 1.0;
    STRAND::QUADRATIC_STRETCHING* stretchingEnergy = new STRAND::QUADRATIC_STRETCHING(k);
    setUnsignedCollisions(10.0,0.1); 
    _strandSolver = new STRAND::TIMESTEPPER(*_strandMesh, *stretchingEnergy, *_eeGeneral);

    _substeps = 1;
    //const REAL dt = 1.0 / 3000.0;
    //const REAL dt = 1.0 / 300.0;
    const REAL dt = 1.0 / 30.0 / _substeps;
    _strandSolver->setDt(dt);

    VECTOR3 sphereCenter0(0,0,0); 
    const REAL radius = 15.0;
    _kinematicShapes.push_back(new SPHERE(sphereCenter0, radius));

    // attach hairs to the sphere
    for (unsigned int x = 0; x < _kinematicShapes.size(); x++)
      _strandSolver->attachKinematicSurfaceConstraints(_kinematicShapes[x]);

    // make the sphere a collision object
    _strandSolver->addKinematicCollisionObject(_kinematicShapes.back());

    //_strandSolver->collisionsEnabled() = true;
    _strandSolver->collisionsEnabled() = false;
    _strandSolver->pcgEnabled() = false;
    //_strandSolver->pcgEnabled() = true;
    _strandSolver->hessianClampingEnabled() = true;
    //_strandSolver->hessianClampingEnabled() = false;
    //_strandSolver->collisionStiffness() = 10;
    _strandSolver->collisionDampingBeta() = 0;

    _strandMesh->bendingForceFilterEnabled() = false;
    //_strandMesh->bendingForceFilterEnabled() = true;

    _eye    = VECTOR3(16.6521, 13.8222, 42.1509);
    _lookAt = VECTOR3(16.3798, 13.7082, 41.1954);
    _up     = VECTOR3(0.0296061, 0.991492, -0.126738);

    _worldCenter = VECTOR3(0,0,0);

    //_pauseFrame = 400;
    _pauseFrame = -1;

    //_writeToFile = true;
    _writeToFile = false;
    _filePath = string("../data/render/");

    char buffer[256];
    snprintf(buffer,256, "E=%4.2e", _E);
    string Estring(buffer);
    snprintf(buffer,256, "r0=%.3f_r1=%.3f", _baseRadius, _tipRadius);
    string Rstring(buffer);

    int totalStrands = _strandMesh->totalStrands();
    cout << " total strands: " << totalStrands << endl;

    snprintf(buffer,256, "s=%d", _totalPoints);
    string Sstring(buffer);

    _sceneName = _sceneName + string("_") + Sstring + string("_") + Estring + string("_") + Rstring;
    cout << " Scene name: " << _sceneName.c_str() << endl;

    return true;
  }

  // Type L hair on LOIS scale
  void buildType4B()
  {
    vector<VECTOR3> restVertices;
    vector<vector<int> > strands;

    const REAL theta = M_PI * 0.25; // unstable

    const REAL rCosTheta = _spacing * cos((M_PI - theta) * 0.5);
    const REAL rSinTheta = _spacing * sin((M_PI - theta) * 0.5);

    cout << " increments: " << rCosTheta << " " << rSinTheta << endl;
    for (unsigned int x = 0; x < _totalPoints; x++)
    {
      if (x <= 3)
      {
        restVertices.push_back(VECTOR3(0, x + 15.0 - 3, 0.5));
        continue;
      }
      const REAL sign = (x % 2) ? 0 : 1;
      restVertices.push_back(VECTOR3(sign * rSinTheta, rCosTheta * x + 15.0 - 3 * rCosTheta, 0.5));
    }

    // slight rotation
    for (unsigned int x = 0; x < _totalPoints; x++)
    {
      MATRIX3 R = rotationMatrix(VECTOR3(0,0,1), -M_PI * 0.25);
      restVertices[x] = R * restVertices[x];
    }
      
    vector<int> strand0;
    for (unsigned int x = 0; x < _totalPoints; x++)
      strand0.push_back(x);
    strands.push_back(strand0);
      
    _strandMesh = new WISP_MESH(restVertices, strands,
                                _E, _nu, _density, _baseRadius, _tipRadius);
  }

  // Type O hair on LOIS scale
  void buildType4CAlternate()
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

    const int stemPoints = 1;

    VECTOR3 scalpVertex;

    //cout << " increments: " << rCosTheta << " " << rSinTheta << endl;
    restVertices.push_back(VECTOR3(-1, 15.0 - stemPoints, 0.0));
    for (unsigned int x = 0; x < _totalPoints; x++)
    {
      if (x <= stemPoints)
      {
        restVertices.push_back(VECTOR3(0, x + 15.0 - stemPoints, 0.0));
        if (x == stemPoints)
          scalpVertex = restVertices.back();
        continue;
      }

      const REAL yReal = _spacing * (x - stemPoints + 1) + 15.0;
      const REAL theta = (_spacing * (x - stemPoints + 1) * 4.0 * M_PI);
      const REAL xReal = -cos(frequency * theta) * amplitude;
      //const REAL xReal = cos(frequency * theta) * amplitude;  // wrong ordering
      const REAL zReal = sin(frequency * theta) * amplitude;

      restVertices.push_back(VECTOR3(xReal, yReal, zReal));
    }

    const REAL curlLength = (scalpVertex - restVertices.back()).norm();
    cout << " Curl length: " << curlLength << endl;

    // slight rotation
    for (unsigned int x = 0; x < restVertices.size(); x++)
    {
      MATRIX3 R = rotationMatrix(VECTOR3(0,0,1), -M_PI * 0.125); // unstable after a while
      restVertices[x] = R * restVertices[x];
    }
      
    vector<int> strand0;
    for (unsigned int x = 0; x < restVertices.size(); x++)
      strand0.push_back(x);
    strands.push_back(strand0);
      
    _strandMesh = new WISP_MESH(restVertices, strands,
                                _E, _nu, _density, _baseRadius, _tipRadius);
  }

  // Type O hair on LOIS scale
  void buildType4C()
  {
    _drawVertices = false;

    vector<VECTOR3> restVertices;
    vector<vector<int> > strands;

    //const REAL theta = M_PI;

    //const REAL theta = M_PI * 0.5; // jumpy
    //const REAL theta = M_PI * 0.375; // jumpy
    //const REAL theta = M_PI * 0.275; // jumpy
    //const REAL theta = M_PI * 0.26; // jumpy
    //const REAL theta = M_PI * 0.25; // unstable
    //const REAL theta = M_PI * 0.175; // unstable
    //const REAL theta = M_PI * 0.125; // unstable

    //const REAL rCosTheta = _spacing * cos((M_PI - theta) * 0.5);
    //const REAL rSinTheta = _spacing * sin((M_PI - theta) * 0.5);

    //const REAL straightLength = _spacing * _totalPoints;

    //const REAL frequency = straightLength * 2;
    //const REAL frequency = straightLength * 2;
    const REAL frequency = 1;
    const REAL amplitude = 0.5;

    //_spacing = 0.714621;
    //_spacing = 1.0;
    _spacing = 0.1;

    const int stemPoints = 2;

    VECTOR3 scalpVertex;

    //cout << " increments: " << rCosTheta << " " << rSinTheta << endl;
    for (unsigned int x = 0; x < _totalPoints; x++)
    {
      if (x <= stemPoints)
      {
        restVertices.push_back(VECTOR3(0, x + 15.0 - stemPoints, 0.0));
        if (x == stemPoints)
          scalpVertex = restVertices.back();
        continue;
      }

      /*
      const REAL yReal = _spacing * x + 15.0;
      //const REAL yReal = 0.714621* x + 12.0;

      const REAL theta = (_spacing * x + 12.0) / (2.0 * M_PI);
      //const REAL theta = (0.714621* x + 12.0) / (2.0 * M_PI);

      const REAL xReal = cos(frequency * theta) * amplitude;
      const REAL zReal = sin(frequency * theta) * amplitude;
      */
      const REAL yReal = _spacing * (x - stemPoints + 1) + 15.0;
      const REAL theta = (_spacing * (x - stemPoints + 1) * 4.0 * M_PI);
      const REAL xReal = cos(frequency * theta) * amplitude;
      const REAL zReal = sin(frequency * theta) * amplitude;

      restVertices.push_back(VECTOR3(xReal, yReal, zReal));
      
    }

    const REAL curlLength = (scalpVertex - restVertices.back()).norm();
    cout << " Curl length: " << curlLength << endl;

    // slight rotation
    for (unsigned int x = 0; x < _totalPoints; x++)
    {
      MATRIX3 R = rotationMatrix(VECTOR3(0,0,1), -M_PI * 0.125); // unstable after a while
      restVertices[x] = R * restVertices[x];
    }
      
    vector<int> strand0;
    for (unsigned int x = 0; x < _totalPoints; x++)
      strand0.push_back(x);
    strands.push_back(strand0);
      
    //_strandMesh = new STRAND_MESH(restVertices, strands,
    //                              _E, _nu, _density, _radiusA, _radiusB);
    _strandMesh = new WISP_MESH(restVertices, strands,
                                _E, _nu, _density, _baseRadius, _tipRadius);
  }

  // match the hair sample we have
  void buildSampleMatch()
  {
    _drawVertices = false;

    vector<VECTOR3> restVertices;
    vector<vector<int> > strands;

    //const REAL straightLength = _spacing * _totalPoints;

    const REAL frequency = 0.5;
    const REAL amplitude = 1.25 * 0.5;

    //_spacing = 0.714621;
    //_spacing = 1.0;
    _totalPoints = 50;
    _spacing = 0.2;

    const int stemPoints = 2;

    VECTOR3 scalpVertex;

    //cout << " increments: " << rCosTheta << " " << rSinTheta << endl;
    for (unsigned int x = 0; x < _totalPoints; x++)
    {
      if (x <= stemPoints)
      {
        restVertices.push_back(VECTOR3(0, x + 15.0 - stemPoints, 0.0));
        if (x == stemPoints)
          scalpVertex = restVertices.back();
        continue;
      }

      const REAL yReal = _spacing * (x - stemPoints + 1) + 15.0;
      const REAL theta = (_spacing * (x - stemPoints + 1) * 4.0 * M_PI);
      const REAL xReal = cos(frequency * theta) * amplitude;
      const REAL zReal = sin(frequency * theta) * amplitude;

      restVertices.push_back(VECTOR3(xReal, yReal, zReal));
      
    }

    const REAL curlLength = (scalpVertex - restVertices.back()).norm();
    cout << " Curl length: " << curlLength << endl;

    // slight rotation
    for (unsigned int x = 0; x < _totalPoints; x++)
    {
      MATRIX3 R = rotationMatrix(VECTOR3(0,0,1), -M_PI * 0.125); // unstable after a while
      restVertices[x] = R * restVertices[x];
    }
      
    vector<int> strand0;
    for (unsigned int x = 0; x < _totalPoints; x++)
      strand0.push_back(x);
    strands.push_back(strand0);
      
    //_strandMesh = new STRAND_MESH(restVertices, strands,
    //                              _E, _nu, _density, _radiusA, _radiusB);
    _strandMesh = new WISP_MESH(restVertices, strands,
                                _E, _nu, _density, _baseRadius, _tipRadius);
  }

  // totally straight hair
  void buildType1()
  {
    vector<VECTOR3> restVertices;
    vector<vector<int> > strands;
    const REAL theta = M_PI; // stable

    const REAL rCosTheta = _spacing * cos((M_PI - theta) * 0.5);
    const REAL rSinTheta = _spacing * sin((M_PI - theta) * 0.5);

    const int rootVertices = 2;

    cout << " increments: " << rCosTheta << " " << rSinTheta << endl;
    for (unsigned int x = 0; x < _totalPoints; x++)
    {
      //const REAL sign = (x % 2) ? 0 : 1;
      //const REAL sign = (x % 2) ? 0 : 1;
      const REAL sign = (x % 2 || x < rootVertices) ? 0 : 1;
      restVertices.push_back(VECTOR3(sign * rSinTheta, rCosTheta * x + 15.0 - rootVertices * rCosTheta, 0.5));
    }

    // slight rotation
    for (unsigned int x = 0; x < _totalPoints; x++)
    {
      MATRIX3 R = rotationMatrix(VECTOR3(0,0,1), -M_PI * 0.25);
      restVertices[x] = R * restVertices[x];
    }
      
    vector<int> strand0;
    for (unsigned int x = 0; x < _totalPoints; x++)
      strand0.push_back(x);
    strands.push_back(strand0);
    _strandMesh = new WISP_MESH(restVertices, strands,
                                _E, _nu, _density, _baseRadius, _tipRadius);
  }

  // simulation loop
  virtual void stepSimulation(const bool verbose = true) override
  {
    const REAL time = _frameNumber * _strandSolver->dt();
    const VECTOR3 translation(0.2 * sin( time), 0,0);

    _kinematicShapes[0]->translation() += translation;

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
      STRAND_MESH::writeSOBJFile(filename.c_str(), _strandMesh->vertices(), _strandMesh->strandIndices());
      cout << "done." << endl;

      TIMER gzipper("Gzipping");
      cout << " Gzipping ..." << flush;
      string gzip = string("gzip -9 ") + filename;
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
  REAL _nu;
  REAL _density;
  REAL _baseRadius;
  REAL _tipRadius;
};

}

#endif

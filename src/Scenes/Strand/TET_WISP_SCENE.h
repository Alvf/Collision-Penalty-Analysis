#ifndef TET_WISP_SCENE_H
#define TET_WISP_SCENE_H

#include "Scenes/SIMULATION_SCENE.h"
#include "Geometry/TET_WISP_MESH.h"
#include "Timestepper/Strand/VOLUME_TIMESTEPPER.h"
#include "Hyperelastic/Strand/QUADRATIC_STRETCHING.h"

namespace HOBAK {

class TET_WISP_SCENE : public SIMULATION_SCENE {
public:
  virtual void printSceneDescription() override
  {
    cout << "=====================================================================" << endl;
    cout << " Hair wisp scene" << endl;
    cout << "=====================================================================" << endl;
  }

  TET_WISP_SCENE() {
    _E = 2.5e7; // looks okay for quadratic 4C
    //_E = 2.5e6; // looks okay for quadratic 4C
    //_E = 2.5e5; // looks okay for quadratic 4C
    //_E = 2.5e4; // decent for sample match
    //_nu = 0.01;
    _G = _E;
  };

  // for testing stretch
  virtual bool buildScene() override
  {
    // this will determine the MOV and JSON filenames
    _sceneName = "wisp";

    //_totalPoints = 50;
    _spacing = 1.0;
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

    _totalPoints = 100; // decent for type 4C
    //_totalPoints = 7; // DEBUG
    //_totalPoints = 8; // DEBUG

    _drawVertices = true;

    // which type of hair?
    //buildType1();
    //buildType4B();
    
    buildType4C();
    //buildSampleMatch();

    cout << " Strand length: " << _strandMesh->strandLength(0) << endl;

    const REAL k = 1.0;
    STRAND::QUADRATIC_STRETCHING* stretchingEnergy = new STRAND::QUADRATIC_STRETCHING(k);
    setUnsignedCollisions(10.0,0.1); 
    _strandSolver = new STRAND::VOLUME_TIMESTEPPER(*_strandMesh, *stretchingEnergy, *_eeGeneral);

    _substeps = 1;
    //const REAL dt = 1.0 / 3000.0;
    //const REAL dt = 1.0 / 300.0;
    const REAL dt = 1.0 / 30.0;
    //const REAL dt = 1.0 / 60.0;
    _strandSolver->setDt(dt);

    VECTOR3 sphereCenter0(0,0,0); 
    //const REAL radius = 15.0;
    const REAL radius = 15.25;
    //const REAL radius = 16.0;
    _kinematicShapes.push_back(new SPHERE(sphereCenter0, radius));

    // attach hairs to the sphere
    for (unsigned int x = 0; x < _kinematicShapes.size(); x++)
      _strandSolver->attachKinematicSurfaceConstraints(_kinematicShapes[x]);

    // bump up the stiffness of attached edges to reflect the rigid attachment
    ((STRAND::VOLUME_TIMESTEPPER*)_strandSolver)->stiffenKinematicStrands();

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

    /*
    _eye    = VECTOR3(32.8654, 20.1034, 7.7234);
    _lookAt = VECTOR3(32.0035, 19.7099, 7.40374);
    _up     = VECTOR3(-0.373167, 0.919277, -0.125191);
    */

    _eye    = VECTOR3(10.8074, 15.06, 5.70459);
    _lookAt = VECTOR3(10.3451, 14.8271, 4.84902);
    _up     = VECTOR3(-0.0261983, 0.968055, -0.249361);

    /*
    _eye    = VECTOR3(17.951, 15.5614, 11.1614);
    _lookAt = VECTOR3(17.5589, 15.3638, 10.2629);
    _up     = VECTOR3(-0.00628123, 0.977208, -0.212183);

    _eye    = VECTOR3(18.8276, 14.4529, 3.42433);
    _lookAt = VECTOR3(18.5438, 14.3279, 2.47363);
    _up     = VECTOR3(0.0250533, 0.99017, -0.137598);
    */

#if 1
    // see if pop and beginning and end
    _eye    = VECTOR3(22.7633, 17.265, 12.5873);
    _lookAt = VECTOR3(22.2279, 17.0408, 11.773);
    _up     = VECTOR3(-0.0262002, 0.968054, -0.249364);
#endif

    _worldCenter = VECTOR3(0,0,0);

    //_pauseFrame = 400;
    //_pauseFrame = -1;
    //_pauseFrame = 18;
    _pauseFrame = 94;

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
      
    _strandMesh = new TET_WISP_MESH(restVertices, strands,
                                    _E, _G, _density, _baseRadius, _tipRadius);
                                    //_E, _nu, _density, _baseRadius, _tipRadius);
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
      
    _strandMesh = new TET_WISP_MESH(restVertices, strands,
                                    _E, _G, _density, _baseRadius, _tipRadius);
                                    //_E, _nu, _density, _baseRadius, _tipRadius);
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
    _strandMesh = new TET_WISP_MESH(restVertices, strands,
                                    _E, _G, _density, _baseRadius, _tipRadius);
                                    //_E, _nu, _density, _baseRadius, _tipRadius);
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
    _strandMesh = new TET_WISP_MESH(restVertices, strands,
                                    _E, _G, _density, _baseRadius, _tipRadius);
                                    //_E, _nu, _density, _baseRadius, _tipRadius);
  }

  // simulation loop
  virtual void stepSimulation(const bool verbose = true) override
  {
    const REAL time = _frameNumber * _strandSolver->dt();
    const REAL amplitude = (0.2 * _strandSolver->dt()) / (1.0 / 30.0);

    //const VECTOR3 translation(0,0,0.2 * sin(time)); // back and forth
    const VECTOR3 translation(amplitude * sin(time), 0,0); // left and right
    //const VECTOR3 translation(0.5 * amplitude, 0,0); // left and right
    //const VECTOR3 translation(0, 0.1 * sin(time), 0);   // up and down
    //const VECTOR3 translation(0, 0.2 * sin(time), 0);   // up and down
    //const VECTOR3 translation(0.2 * sin(time * 10), 0,0);

    //_kinematicShapes[0]->translation() += translation;

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
      //TIMER::printTimingsPerFrame(_frameNumber);
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
#if 1
    drawStrandMesh(*_strandMesh, -1, true);
#else
    TET_WISP_MESH* wisp = (TET_WISP_MESH*)_strandMesh;
    drawWispMesh(*wisp, -1, _drawVertices);
#endif

    glEnable(GL_DEPTH_TEST);

    // only draw the collision sphere
    //for (unsigned int x = 0; x < _kinematicShapes.size(); x++)
    //  drawKinematicShape(*_kinematicShapes[x]);

    //drawInvertedTetStrandMesh(*_strandMesh);

    //drawTetForces(*_strandMesh, 2);
    //drawTetForces(*_strandMesh, 3);
    drawTetForces(*_strandMesh, 4);
    //drawTetForces(*_strandMesh, 5);
    
    //drawKinematicShape(*_kinematicShapes.back());

    //drawCollisionsOld(*_strandMesh);
    //drawCollisions(*_strandMesh);

    //drawKinematicConstraints(_strandMesh, _strandSolver);
    //drawPlaneConstraints(_strandMesh, _strandSolver);
    
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
  REAL _spacing;
  REAL _density;
  REAL _baseRadius;
  REAL _tipRadius;
};

}

#endif

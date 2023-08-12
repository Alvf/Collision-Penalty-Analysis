#ifndef STRAND_STRAIGHTEN_REDUX_H
#define STRAND_STRAIGHTEN_REDUX_H

#include "Scenes/SIMULATION_SCENE.h"
#include "Geometry/TET_WISP_MESH.h"
#include "Timestepper/Strand/VOLUME_TIMESTEPPER.h"
#include "Timestepper/Strand/VOLUME_BDF_1.h"
#include "Hyperelastic/Strand/QUADRATIC_STRETCHING.h"
#include <random>
#include <float.h>

namespace HOBAK {

class STRAND_STRAIGHTEN_REDUX : public SIMULATION_SCENE {
public:
  virtual void printSceneDescription() override
  {
    cout << "=====================================================================" << endl;
    cout << " Strand straighten scene" << endl;
    cout << "=====================================================================" << endl;
  }

  STRAND_STRAIGHTEN_REDUX() 
  {
    _E = 3.9e9;
    //_E = 1.0e8;
    //_nu = 0.48;
    //_G = _E;
    //_density = 0.065;
    REAL nu = 0.48;
    _G = _E / (2.0 * (1 + nu));
    _density = 1.3; // 5% of the usual density
    //_baseRadius = 0.0037;
    //_tipRadius = 0.0037;
    _baseRadius = 0.2 * 0.5;
    _tipRadius = 0.2 * 0.5;
    
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
    _sceneName = "strand_straighten_" + jitterString;

    //_totalPoints = 50;
    //_spacing = 1.0;

    _gravity = VECTOR3(0,-981,0);
    //_gravity = VECTOR3(-981,0,0);
    //_gravity = VECTOR3(0,-98.1,0);
    //_gravity = VECTOR3(0,0,0);

    _totalPoints = 20; // decent for type 4C

    _drawVertices = true;

    //_totalStrands = 15;
    //_totalStrands = 100;
    //_totalStrands = 1000;
    _totalStrands = 1;

    buildCoil();

    cout << " Strand length: " << _strandMesh->strandLength(0) << endl;

    const REAL k = 1.0;
    STRAND::QUADRATIC_STRETCHING* stretchingEnergy = new STRAND::QUADRATIC_STRETCHING(k);
    //_strandSolver = new STRAND::VOLUME_TIMESTEPPER(*_strandMesh, *stretchingEnergy);
    setUnsignedCollisions(1000.0,0.1); 
    _strandSolver = new STRAND::VOLUME_BDF_1(*_strandMesh, *stretchingEnergy, *_eeGeneral);
    
    //const int maxNewtonIterations = 10; // flickering stabilizes
    const int maxNewtonIterations = 3;  // works with warm start
    //const int maxNewtonIterations = 1;  // works with warm start
    _strandSolver->maxNewtonIterations() = maxNewtonIterations;

    _substeps = 1;
    //const REAL dt = 1.0 / 3000.0;
    //const REAL dt = 1.0 / 300.0;
    const REAL dt = 1.0 / 30.0;
    _strandSolver->setDt(dt);

    //const REAL radius = 15.25;
    //const REAL radius = 1.0;
    //VECTOR3 sphereCenter0(0,-radius,0); 
    const VECTOR3 sphereCenter0 = _strandMesh->restVertices()[0];
    cout << " Sphere center: " << sphereCenter0.transpose() << endl;
    _kinematicShapes.push_back(new SPHERE(sphereCenter0, 0.75));
    //_kinematicShapes.push_back(new SPHERE(sphereCenter0, 0.25));
    //const VECTOR3 sphereCenter1 = _strandMesh->restVertices().back();
    //_kinematicShapes.push_back(new SPHERE(sphereCenter1, 0.25));

    // attach hairs to the sphere
    for (unsigned int x = 0; x < _kinematicShapes.size(); x++)
      _strandSolver->attachKinematicSurfaceConstraints(_kinematicShapes[x]);

    // apply the straighten
    straightenVertices();

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

    _eye    = VECTOR3(7.29718, 3.0797, 17.5469);
    _lookAt = VECTOR3(7.17407, 2.90108, 16.5708);
    _up     = VECTOR3(0.0132358, 0.983286, -0.181586);

    /*
    _eye    = VECTOR3(10.3154, 32.1291, 174.712);
    _lookAt = VECTOR3(10.3056, 31.9476, 173.728);
    _up     = VECTOR3(0.0132357, 0.983283, -0.181605);

    _eye    = VECTOR3(11.6214, 46.0981, 250.485);
    _lookAt = VECTOR3(11.6034, 45.9169, 249.502);
    _up     = VECTOR3(0.0132376, 0.983301, -0.181503);
    */

    _eye    = VECTOR3(-3.03448, -1.46185, 22.7735);
    _lookAt = VECTOR3(-2.84879, -1.64754, 21.8086);
    _up     = VECTOR3(0.0350931, 0.982607, -0.182349);

    _worldCenter = VECTOR3(0,0,0);

    //_pauseFrame = 1000;
    //_pauseFrame = 200;
    _pauseFrame = 40;
    //_pauseFrame = -1;

    _filePath = string("../data/render/");

    snprintf(buffer,256, "E=%4.2e", _E);
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

  // Type O hair on LOIS scale
  void buildCoil()
  {
    _drawVertices = false;

    vector<VECTOR3> restVertices;
    vector<vector<int> > strands;

    //const REAL straightLength = _spacing * _totalPoints;

    const REAL frequency = 2.0;
    const REAL amplitude = 0.75;

    //_spacing = 0.1;
    //_spacing = 0.08;
    _spacing = 0.05;

    const REAL radius = -0.25;

    const int stemPoints = 0;

    VECTOR3 scalpVertex;
    scalpVertex.setZero();

    VECTOR3 first(-0.231763, -0.8, 0.713292);
    first[1] -= 0.5;
    restVertices.push_back(first);

    //restVertices.push_back(VECTOR3(-0.1, radius - stemPoints, 0.0));
    for (unsigned int x = 0; x < _totalPoints; x++)
    {
      /*
      if (x <= stemPoints)
      {
        restVertices.push_back(VECTOR3(0, x + radius - stemPoints, 0.0));
        if (x == stemPoints)
          scalpVertex = restVertices.back();
        continue;
      }
      */

      const REAL yReal = _spacing * (x - stemPoints + 1) + radius;
      const REAL theta = (_spacing * (x - stemPoints + 1) * 4.0 * M_PI);
      const REAL xReal = -cos(frequency * theta) * amplitude;
      //const REAL xReal = cos(frequency * theta) * amplitude;  // wrong ordering
      const REAL zReal = sin(frequency * theta) * amplitude;

      const VECTOR3 vertex = VECTOR3(xReal, 4.0 * yReal, zReal);

      cout << " vertex: " << vertex.transpose() << endl;
      //restVertices.push_back(VECTOR3(xReal, 4.0 * yReal, zReal));
      restVertices.push_back(vertex);
    }

    MATRIX3 R = rotationMatrix(VECTOR3(0,0,1), -M_PI);
    for (unsigned int x = 0; x < restVertices.size(); x++)
    {
      restVertices[x] = R * restVertices[x];
      //cout << " vertex " << x << ": " << restVertices[x].transpose() << endl;
    }

    const REAL curlLength = (scalpVertex - restVertices.back()).norm();
    cout << " Curl length: " << curlLength << endl;

    vector<int> strand0;
    for (unsigned int x = 0; x < restVertices.size(); x++)
      strand0.push_back(x);
    strands.push_back(strand0);
      
    _strandMesh = new TET_WISP_MESH(restVertices, strands,
                                    _E, _G, _density, _baseRadius, _tipRadius);
  }

  void straightenVertices()
  {
    vector<VECTOR3>& vertices = _strandMesh->vertices();

    // get length of first segment
    const REAL length = (vertices[0] - vertices[1]).norm();
    
    // get position of first vertex
    const VECTOR3 v0 = vertices[0];

    // displace the rest
    for (unsigned int x = 1; x < vertices.size(); x++)
      vertices[x] = v0 - VECTOR3(0,1,0) * x * length;
    _strandMesh->updateProperties();
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

    /*
    if (_frameNumber < 500)
    {
      _kinematicShapes[0]->translation() -= VECTOR3(0.3,0,0);
      _kinematicShapes[1]->translation() += VECTOR3(0.3,0,0);
    }
    if (_frameNumber == 500)
      _strandSolver->popBackKinematicConstraint();
      */

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
    glColor4f(1.0, 0.0, 0.0, 1.0);
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

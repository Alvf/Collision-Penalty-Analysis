#ifndef TET_MULTI_WISP_SCENE_H
#define TET_MULTI_WISP_SCENE_H

#include "Scenes/SIMULATION_SCENE.h"
#include "Geometry/TET_WISP_MESH.h"
#include "Timestepper/Strand/VOLUME_TIMESTEPPER.h"
#include "Timestepper/Strand/VOLUME_BDF_1.h"
#include "Hyperelastic/Strand/QUADRATIC_STRETCHING.h"
#include <random>
#include <float.h>

namespace HOBAK {

class TET_MULTI_WISP_SCENE : public SIMULATION_SCENE {
public:
  virtual void printSceneDescription() override
  {
    cout << "=====================================================================" << endl;
    cout << " Multi-hair wisp scene" << endl;
    cout << "=====================================================================" << endl;
  }

  // generate some random points on the sphere
  void buildRandomPoints()
  {
    std::mt19937 gen(123);
    std::uniform_real_distribution<REAL> spinRNG(0.0, 2 * M_PI);
    std::uniform_real_distribution<REAL> heightRNG(0.0, 1.0);

    unsigned int totalPoints = _totalStrands;
    //const REAL blueRadius = 0.05; // seems to work for 1000
    //const REAL blueRadius = 0.04; // seems to work for 2000
    REAL blueRadius = 0.05;
    if (_totalStrands == 2000)
      blueRadius = 0.04;
    else if (_totalStrands == 1000)
      blueRadius = 0.05;
    else if (_totalStrands == 100 || _totalStrands < 100)
      blueRadius = 0.05;
    else
    {
      std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
      cout <<" NEED TO SET BLUR NOISE RADIUS " << endl;
      exit(0);
    }

    while (_randomPoints.size() < totalPoints)
    {
      //cout << "Blue noise points: " << _randomPoints.size() << endl;
      const REAL theta = spinRNG(gen);
      const REAL u = heightRNG(gen);
      const REAL coeff = sqrt(1 - u * u);
      VECTOR3 point(coeff * cos(theta), u, coeff * sin(theta));

      // find the closest point
      REAL closest = FLT_MAX;
      for (unsigned int x = 0; x < _randomPoints.size(); x++)
      {
        const VECTOR3 diff = point - _randomPoints[x];
        if (diff.norm() < closest)
          closest = diff.norm();
      }

      // rejection sampling
      if (closest < blueRadius) continue;

      _randomPoints.push_back(point);
    }
    
    for (unsigned int x = 0; x < _randomPoints.size(); x++)  
      _randomPoints[x] *= 15.0;
  }

  TET_MULTI_WISP_SCENE() 
  {
    //_E = 2.5e5; // looks okay for quadratic 4C
    //_E = 2.5e7; // looks okay for quadratic 4C
    //_E = 2.5e7;
    //_E = 1.0e8;
    //_E = 7.5e7;
    _E = 5.0e7;
    //_E = 7.5e5; 
    //_E = 5e4; // decent for type 1
    //_E = 5e4; // decent for type 4C
    //_E = 2.5e4; // decent for sample match
    
    //_E = 1e4; // decent for type 4B
    //_nu = 0.36;
    //_nu = 0.01;
    _G = _E;
  };

  // for testing stretch
  virtual bool buildScene() override
  {
    _autoplay = true;

    // this will determine the MOV and JSON filenames
    //_sceneName = "multi_wisp";
    _sceneName = "multi_wisp_newton_3";

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

    _drawVertices = true;

    //_totalStrands = 15;
    //_totalStrands = 100;
    //_totalStrands = 1000;
    _totalStrands = 2000;
    buildRandomPoints();

    //buildType4C();
    buildMultiType4C();

    cout << " Strand length: " << _strandMesh->strandLength(0) << endl;

    const REAL k = 1.0;
    STRAND::QUADRATIC_STRETCHING* stretchingEnergy = new STRAND::QUADRATIC_STRETCHING(k);
    //_strandSolver = new STRAND::VOLUME_TIMESTEPPER(*_strandMesh, *stretchingEnergy);
    //_strandSolver = new STRAND::VOLUME_BDF_1(*_strandMesh, *stretchingEnergy);
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

    _strandSolver->collisionsEnabled() = true;
    //_strandSolver->collisionsEnabled() = false;
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

    _worldCenter = VECTOR3(0,0,0);

    _pauseFrame = 500;
    //_pauseFrame = 1000;
    //_pauseFrame = -1;

    //_writeToFile = true;
    _writeToFile = false;
    _filePath = string("../data/render/");

    char buffer[256];
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
    _sceneName = _sceneName + string("_") + Sstring + string("_") + Estring + string("_") + Rstring + string("_") + Dstring;
    cout << " Scene name: " << _sceneName.c_str() << endl;

    return true;
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

  // Type O hair on LOIS scale
  void buildMultiType4C()
  {
    _drawVertices = false;

    //const REAL straightLength = _spacing * _totalPoints;

    const REAL frequency = 1;
    const REAL amplitude = 0.5;

    //_spacing = 0.714621;
    //_spacing = 1.0;
    _spacing = 0.1;

    const int stemPoints = 1;

    vector<VECTOR3> wisp;
    //cout << " increments: " << rCosTheta << " " << rSinTheta << endl;
    wisp.push_back(VECTOR3(-1, 15.0 - stemPoints, 0.0));
    for (unsigned int x = 0; x < _totalPoints; x++)
    {
      if (x <= stemPoints)
      {
        wisp.push_back(VECTOR3(0, x + 15.0 - stemPoints, 0.0));
        continue;
      }

      const REAL yReal = _spacing * (x - stemPoints + 1) + 15.0;
      const REAL theta = (_spacing * (x - stemPoints + 1) * 4.0 * M_PI);
      const REAL xReal = -cos(frequency * theta) * amplitude;
      const REAL zReal = sin(frequency * theta) * amplitude;

      wisp.push_back(VECTOR3(xReal, yReal, zReal));
    }

    vector<VECTOR3> restVertices;
    vector<vector<int> > strands;

    // spin it along its axis too, for variety
    std::mt19937 gen(123);
    std::uniform_real_distribution<REAL> spinRNG(0.0, 2 * M_PI);

    // slight rotation
    for (unsigned int i = 0; i < _randomPoints.size(); i++)
    {
      vector<int> strand;

      VECTOR3 point = _randomPoints[i].normalized();
      VECTOR3 up(0,1,0);
      VECTOR3 axis = up.cross(point);
      axis.normalize();
      const REAL angle = acos(up.dot(point));

      const REAL theta = spinRNG(gen);
      MATRIX3 spin = rotationMatrix(up, theta);

      MATRIX3 R = rotationMatrix(axis, angle);
      for (unsigned int x = 0; x < wisp.size(); x++)
      {
        //restVertices.push_back(R * wisp[x]);
        restVertices.push_back(R * spin * wisp[x]);
        strand.push_back(restVertices.size() - 1);
      }
      strands.push_back(strand);
    }
      
    _strandMesh = new TET_WISP_MESH(restVertices, strands,
                                    _E, _G, _density, _baseRadius, _tipRadius);
                                    //_E, _nu, _density, _baseRadius, _tipRadius);
  }

  // simulation loop
  virtual void stepSimulation(const bool verbose = true) override
  {
    const REAL time = _frameNumber * _strandSolver->dt();
    //const VECTOR3 translation(0,0,0.2 * sin(time)); // back and forth
    //const VECTOR3 translation(0.2 * sin(time), 0,0); // left and right
    const VECTOR3 translation(0, 0.1 * sin(time), 0);   // up and down
    //const VECTOR3 translation(0.2 * sin(time * 10), 0,0);

    //_kinematicShapes[0]->translation() += translation;

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
};

}

#endif

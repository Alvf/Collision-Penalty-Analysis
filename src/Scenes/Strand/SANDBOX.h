#ifndef SANDBOX_H
#define SANDBOX_H

#include "Scenes/SIMULATION_SCENE.h"
#include "Geometry/TET_STRAND_MESH.h"
#include "Timestepper/Strand/TIMESTEPPER.h"
#include "Timestepper/Strand/VOLUME_TIMESTEPPER.h"
#include "Hyperelastic/Strand/QUADRATIC_STRETCHING.h"
#include "Hyperelastic/Strand/ISOTROPIC_BENDING.h"
#include "Hyperelastic/Strand/QUADRATIC_BENDING.h"
#include "Hyperelastic/Strand/QUADRATIC_UNIT_BENDING.h"
#include "Hyperelastic/Strand/ISOTROPIC_THETA.h"

namespace HOBAK {

class SANDBOX : public SIMULATION_SCENE {
public:
  virtual void printSceneDescription() override
  {
    cout << "=====================================================================" << endl;
    cout << " Sandbox strand scene" << endl;
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
    _sceneName = "sandbox";

    _totalPoints = 3;
    //_totalPoints = 50;
    //_totalPoints = 20;
    _spacing = 1.0;
    //_E = 1e8; // default
    _E = 1e9;
    _nu = 0.36;
    _density = 1.32;
    _radiusA = 0.02;
    _radiusB = 0.02;

    _gravity = VECTOR3(0,-981,0);
    //_gravity = VECTOR3(0,0,0);

#if 1
    // which type of hair?
    //_E = 1e5;
    _E = 1e10;
    _totalPoints = 5;
    //_totalPoints = 25;
    _totalPoints = 50;
    //_totalPoints = 500;
    //buildType1();
    buildType4CShortStem();
#elif 0
    _totalPoints = 11;
    buildType4B();
#else
    _totalPoints = 50;
    buildType4C();
#endif
    _strandMesh->twists()[0] = 1;
    std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;

    cout << " Strand length: " << _strandMesh->strandLength(0) << endl;

    const REAL k = _E * M_PI * _radiusA * _radiusB;
    STRAND::QUADRATIC_STRETCHING* stretchingEnergy = new STRAND::QUADRATIC_STRETCHING(k);
    setCollisions(10.0,0.001); 
    _strandSolver = new STRAND::VOLUME_TIMESTEPPER(*_strandMesh, *stretchingEnergy, *_eeGeneral);

    _substeps = 1;
    const REAL dt = 1.0 / 300.0;
    //const REAL dt = 1.0 / 30.0;
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

    // zoom out for Type 1
    _eye    = VECTOR3(21.8627, 15.9426, 59.9575);
    _lookAt = VECTOR3(21.6341, 15.8259, 58.991);
    _up     = VECTOR3(0.0296078, 0.991494, -0.126731);

    _eye    = VECTOR3(10.4504, 13.2706, 14.3831);
    _lookAt = VECTOR3(10.1183, 13.1608, 13.4462);
    _up     = VECTOR3(0.029611, 0.991494, -0.126722);

    _eye    = VECTOR3(11.2095, 15.5793, 5.67528);
    _lookAt = VECTOR3(10.5983, 15.2593, 4.95136);
    _up     = VECTOR3(-0.138005, 0.943713, -0.300598);

    _worldCenter = VECTOR3(0,0,0);

    _pauseFrame = -1;
    //_pauseFrame = 1000;
    //_pauseFrame = 57;
    //_pauseFrame = 58;
    _pauseFrame = 100;

    //_writeToFile = true;
    _writeToFile = false;
    _filePath = string("../data/render/");

    char buffer[256];
    sprintf(buffer, "E=%4.2e", _E);
    string Estring(buffer);
    sprintf(buffer, "r=%.3f", _radiusA);
    string Rstring(buffer);

    int totalStrands = _strandMesh->totalStrands();
    cout << " total strands: " << totalStrands << endl;

    if (totalStrands < 1000)
      sprintf(buffer, "s=%d", totalStrands);
    else
      sprintf(buffer, "s=%dK", totalStrands / 1000);
    string Sstring(buffer);

    _sceneName = _sceneName + string("_") + Sstring + string("_") + Estring + string("_") + Rstring;
    cout << " Scene name: " << _sceneName.c_str() << endl;

    return true;
  }

  // totally straight hair
  void buildType1()
  {
    vector<VECTOR3> restVertices;
    vector<vector<int> > strands;
    const REAL theta = M_PI; // stable

    const REAL rCosTheta = _spacing * cos((M_PI - theta) * 0.5);
    const REAL rSinTheta = _spacing * sin((M_PI - theta) * 0.5);

    const int stemPoints = 1;
    //const int stemPoints = 3;

    cout << " increments: " << rCosTheta << " " << rSinTheta << endl;
    for (unsigned int x = 0; x < _totalPoints; x++)
    {
      //const REAL sign = (x % 2) ? 0 : 1;
      //const REAL sign = (x % 2) ? 0 : 1;
      const REAL sign = (x % 2 || x < stemPoints) ? 0 : 1;
      restVertices.push_back(VECTOR3(sign * rSinTheta, rCosTheta * x + 15.0 - stemPoints* rCosTheta, 0.5));
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
    _strandMesh = new TET_STRAND_MESH(restVertices, strands,
                                  _E, _nu, _density, _radiusA, _radiusB);
  }

  // Type O hair on LOIS scale
  void buildType4CShortStem()
  {
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

    const int stemPoints = 1;

    //cout << " increments: " << rCosTheta << " " << rSinTheta << endl;
    for (unsigned int x = 0; x < _totalPoints; x++)
    {
      if (x <= stemPoints)
      {
        restVertices.push_back(VECTOR3(0, x + 15.0 - stemPoints - 0.01, 0.0));
        //restVertices.push_back(VECTOR3(0, x + 15.0 - stemPoints, 0.0));
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

      // winds clockwise
      const REAL theta = (_spacing * (x - stemPoints + 1) * 4.0 * M_PI);
      const REAL xReal = -cos(frequency * theta) * amplitude;
      const REAL zReal = sin(frequency * theta) * amplitude;

      restVertices.push_back(VECTOR3(xReal, yReal, zReal));
    }

    // slight rotation
    for (unsigned int x = 0; x < _totalPoints; x++)
    {
      MATRIX3 R = rotationMatrix(VECTOR3(0,0,1), -M_PI * 0.125); // unstable after a while
      //MATRIX3 R = rotationMatrix(VECTOR3(0,0,1), -M_PI * 0.25); // jitters
      //MATRIX3 R = rotationMatrix(VECTOR3(0,0,1), M_PI * 0.25); // jitters
      //MATRIX3 R = rotationMatrix(VECTOR3(0,0,1), -M_PI * 0.26); // jitters
      //MATRIX3 R = rotationMatrix(VECTOR3(0,0,1), -M_PI * 0.126); // unstable
      restVertices[x] = R * restVertices[x];
    }
      
    vector<int> strand0;
    for (unsigned int x = 0; x < _totalPoints; x++)
      strand0.push_back(x);
    strands.push_back(strand0);
      
    _strandMesh = new TET_STRAND_MESH(restVertices, strands,
                                  _E, _nu, _density, _radiusA, _radiusB);
  }

  // Type O hair on LOIS scale
  void buildType4C()
  {
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

    //const int stemPoints = 3;
    const int stemPoints = 0;

    //cout << " increments: " << rCosTheta << " " << rSinTheta << endl;
    for (unsigned int x = 0; x < _totalPoints; x++)
    {
      if (x <= stemPoints)
      {
        restVertices.push_back(VECTOR3(0, x + 15.0 - stemPoints, 0.0));
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

      // winds clockwise
      const REAL theta = (_spacing * (x - stemPoints + 1) * 4.0 * M_PI);
      const REAL xReal = cos(frequency * theta) * amplitude;
      const REAL zReal = sin(frequency * theta) * amplitude;

      restVertices.push_back(VECTOR3(xReal, yReal, zReal));
    }

    // slight rotation
    for (unsigned int x = 0; x < _totalPoints; x++)
    {
      MATRIX3 R = rotationMatrix(VECTOR3(0,0,1), -M_PI * 0.125); // unstable after a while
      //MATRIX3 R = rotationMatrix(VECTOR3(0,0,1), -M_PI * 0.25); // jitters
      //MATRIX3 R = rotationMatrix(VECTOR3(0,0,1), M_PI * 0.25); // jitters
      //MATRIX3 R = rotationMatrix(VECTOR3(0,0,1), -M_PI * 0.26); // jitters
      //MATRIX3 R = rotationMatrix(VECTOR3(0,0,1), -M_PI * 0.126); // unstable
      restVertices[x] = R * restVertices[x];
    }
      
    vector<int> strand0;
    for (unsigned int x = 0; x < _totalPoints; x++)
      strand0.push_back(x);
    strands.push_back(strand0);
      
    _strandMesh = new TET_STRAND_MESH(restVertices, strands,
                                  _E, _nu, _density, _radiusA, _radiusB);
  }

  // Type L hair on LOIS scale
  void buildType4B()
  {
    vector<VECTOR3> restVertices;
    vector<vector<int> > strands;

    //const REAL theta = M_PI * 0.5; // jumpy
    //const REAL theta = M_PI * 0.375; // jumpy
    //const REAL theta = M_PI * 0.275; // jumpy
    //const REAL theta = M_PI * 0.26; // jumpy
    const REAL theta = M_PI * 0.25; // unstable
    //const REAL theta = M_PI * 0.175; // unstable
    //const REAL theta = M_PI * 0.125; // unstable

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
      
    _strandMesh = new TET_STRAND_MESH(restVertices, strands,
                                  _E, _nu, _density, _radiusA, _radiusB);
  }

  // simulation loop
  virtual void stepSimulation(const bool verbose = true) override
  {
    //_kinematicShapes[1]->translation() -= VECTOR3(0,0, 0.01);
    //_kinematicShapes[1]->translation() -= VECTOR3(0,0, 0.075);
    //_kinematicShapes[1]->translation() -= VECTOR3(0,0, 0.1);
    //_kinematicShapes[1]->translation() -= VECTOR3(0,0, 0.2);

    for (unsigned int x = 0; x < _substeps; x++)
    {
      _strandSolver->externalForces().setZero();
      _strandSolver->addGravity(_gravity);
      _strandSolver->solveDynamics(verbose);
      //_strandSolver->solveNewton(verbose);

      //const int longest = _strandMesh->longestStrand();
      //cout << " Longest is strand " << longest << " which is " << _strandMesh->strandLength(longest) << endl;
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
    //drawStrandMeshOld(*_strandMesh);
    
    //drawTetStrandMesh(*_strandMesh, -1, true);
    drawStrandMesh(*_strandMesh, -1, true);
    _strandMesh->computeFs();

    //glColor4f(1.0, 0.0, 0.0, 0.1);
    //drawInvertedTetStrandMesh(*_strandMesh);
    
    const int highlight = 4;

    if (_strandMesh->isTetInverted()[highlight])
      glColor4f(1.0, 0.0, 0.0, 0.1);
    else
      glColor4f(0.0, 0.0, 1.0, 0.1);
    drawTet(*_strandMesh, 4);
    //drawTetForces(*_strandMesh, 4);
    //drawAnisotropy(*_strandMesh, 4);

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

  // accessors
  TET_STRAND_MESH* strandMesh() { return _strandMesh; };

protected:
  TET_STRAND_MESH* _strandMesh;
  STRAND::TIMESTEPPER* _strandSolver;
  STRAND::STRETCHING* _strechingEnergy;

  unsigned int _substeps;

  string _filePath;
  bool _writeToFile;

  // strand generation parameters
  unsigned int _totalPoints;
  REAL _spacing;
  REAL _E;
  REAL _nu;
  REAL _density;
  REAL _radiusA;
  REAL _radiusB;
};

}

#endif

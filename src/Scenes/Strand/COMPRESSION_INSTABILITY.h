#ifndef COMPRESSION_INSTABILITY_H
#define COMPRESSION_INSTABILITY_H

#include "Scenes/SIMULATION_SCENE.h"
#include "Geometry/STRAND_MESH.h"
#include "Geometry/STRAND_MESH_FASTER.h"
#include "Timestepper/Strand/TIMESTEPPER.h"
#include "Hyperelastic/Strand/QUADRATIC_STRETCHING.h"
#include "Hyperelastic/Strand/ISOTROPIC_BENDING.h"
#include "Hyperelastic/Strand/QUADRATIC_BENDING.h"
#include "Hyperelastic/Strand/QUADRATIC_UNIT_BENDING.h"
#include "Hyperelastic/Strand/ISOTROPIC_THETA.h"

namespace HOBAK {

class COMPRESSION_INSTABILITY : public SIMULATION_SCENE {
public:
  virtual void printSceneDescription() override
  {
    cout << "=====================================================================" << endl;
    cout << " Unstable strand compression scene" << endl;
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
    _sceneName = "compression_instability";

    vector<VECTOR3> restVertices;
    vector<vector<int> > strands;

    const unsigned int totalPoints = 20;
    //const unsigned int totalPoints = 50;
    const REAL spacing = 1.0;

    const REAL theta = M_PI; // stable
    //const REAL theta = M_PI * 0.25; // unstable
    //const REAL theta = M_PI * 0.5; // jumpy

    const REAL rCosTheta = spacing * cos((M_PI - theta) * 0.5);
    const REAL rSinTheta = spacing * sin((M_PI - theta) * 0.5);

    cout << " increments: " << rCosTheta << " " << rSinTheta << endl;
    for (unsigned int x = 0; x < totalPoints; x++)
    {
      const REAL sign = (x % 2 || x < 3) ? 0 : 1;
      restVertices.push_back(VECTOR3(sign * rSinTheta, rCosTheta * x + 15.0 - 3 * rCosTheta, 0.5));
    }

    // theta = M_PI * 0.25
    //const REAL compression = 0.5; // stable
    //const REAL compression = 0.25; // stable
    //const REAL compression = 0.15; // stable
    //const REAL compression = 0.13; // stable
    //const REAL compression = 0.1275; // unstable
    //const REAL compression = 0.125; // unstable
    //const REAL compression = 0.1; // unstable
   
    // theta = M_PI * 0.5
    //const REAL compression = 0.1275; // stable
    //const REAL compression = 0.1; // stable
    //const REAL compression = 0.075; //stable
    //const REAL compression = 0.05; // unstable
    //const REAL compression = 0.01; // unstable
    
    // theta = M_PI
    //const REAL compression = 1.0;
    //const REAL compression = 0.1;
    //const REAL compression = 0.01; // stable
    const REAL compression = 1e-6; // stable

    const int stemSegments = 3;

    vector<VECTOR3> vertices;
    for (unsigned int x = 0; x < totalPoints; x++)
    {
      bool atStem = (x < stemSegments);
      const REAL sign = (x % 2 || atStem) ? 0 : 1;

      if (atStem)
        vertices.push_back(VECTOR3(sign * rSinTheta, rCosTheta * x + 15.0 - stemSegments * rCosTheta, 0.5));
      else
        vertices.push_back(VECTOR3(sign * rSinTheta, compression * rCosTheta * (x - stemSegments) + 15.0, 0.5));
    }

    /*
    // slight rotation
    for (unsigned int x = 0; x < totalPoints; x++)
    {
      MATRIX3 R = rotationMatrix(VECTOR3(0,0,1), -M_PI * 0.25);
      restVertices[x] = R * restVertices[x];
    }
    */
      
    vector<int> strand0;
    for (unsigned int x = 0; x < totalPoints; x++)
      strand0.push_back(x);
    strands.push_back(strand0);
      
    //_gravity = VECTOR3(0,-981,0);
    _gravity = VECTOR3(0,0,0);

    //const REAL E = 1e7; // Taz settings
    //const REAL E = 1e5; // misbehaving is quite obvious
    //const REAL E = 1e3;
    //const REAL E = 3.727e10;
    //const REAL E = 1e4;
    //const REAL E = 3.727e13;
    //const REAL E = 3.727e12;
    //const REAL E = 3.727e11;  // looks good for curly hair
    //const REAL E = 3.727e10;
    const REAL E = 1e8;
    //const REAL E = 3.727e8;
    //const REAL E = 3.727e6;
    const REAL nu = 0.36;
    const REAL density = 1.32;
    //const REAL radiusA = 0.005;
    //const REAL radiusB = 0.005;
    const REAL radiusA = 0.02;
    const REAL radiusB = 0.02;

    char buffer[256];
    sprintf(buffer, "E=%4.2e", E);
    string Estring(buffer);
    sprintf(buffer, "r=%.3f", radiusA);
    string Rstring(buffer);

    int totalStrands = strands.size();
    cout << " total strands: " << totalStrands << endl;

    if (totalStrands < 1000)
      sprintf(buffer, "s=%d", totalStrands);
    else
      sprintf(buffer, "s=%dK", totalStrands / 1000);
    string Sstring(buffer);

    _sceneName = _sceneName + string("_") + Sstring + string("_") + Estring + string("_") + Rstring;
    cout << " Scene name: " << _sceneName.c_str() << endl;
    
#if 0
    _strandMesh = new STRAND_MESH_FASTER(restVertices, strands,
                                         E, nu, density, radiusA, radiusB);
#else
    _strandMesh = new STRAND_MESH(restVertices, vertices,strands,
                                  E, nu, density, radiusA, radiusB);
#endif
    setUnsignedCollisions(10.0,0.1);

    const REAL k = E * M_PI * radiusA * radiusB;
    cout <<" stretching stiffness: " << k << endl;
    STRAND::QUADRATIC_STRETCHING* stretchingEnergy = new STRAND::QUADRATIC_STRETCHING(k);
    _strandSolver = new STRAND::TIMESTEPPER(*_strandMesh, *stretchingEnergy, *_eeGeneral);

    _substeps = 1;
    //const REAL dt = 1.0 / 300.0;
    const REAL dt = 1.0 / 30.0;
    _strandSolver->setDt(dt);

    VECTOR3 sphereCenter0(0,0,0); 
    //VECTOR3 sphereCenter1 = restVertices[last];

    const REAL radius = 15.0;
    //const REAL radius = 0.1;
    _kinematicShapes.push_back(new SPHERE(sphereCenter0, radius));
    //_kinematicShapes.push_back(new SPHERE(sphereCenter1, radius));

    // attach hairs to the sphere
    for (unsigned int x = 0; x < _kinematicShapes.size(); x++)
      _strandSolver->attachKinematicSurfaceConstraints(_kinematicShapes[x]);

    // make the sphere a collision object
    //_strandSolver->addKinematicCollisionObject(_kinematicShapes.back());

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

    _eye    = VECTOR3(-6.30346, 25.7801, 29.4684);
    _lookAt = VECTOR3(-6.09867, 25.6976, 28.4931);
    _up     = VECTOR3(0,1,0);

    _worldCenter = VECTOR3(0,0,0);

    //_pauseFrame = 3;
    //_pauseFrame = 33;
    //_pauseFrame = 200;
    //_pauseFrame = 1000;
    _pauseFrame = 100;
    //_pauseFrame = 10;

    //_writeToFile = true;
    _writeToFile = false;
    _filePath = string("../data/render/");

    return true;
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
    drawStrandMesh(*_strandMesh, -1, true);
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
};

}

#endif

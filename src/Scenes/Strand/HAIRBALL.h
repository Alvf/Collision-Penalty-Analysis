#ifndef HAIRBALL_H
#define HAIRBALL_H

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

class HAIRBALL : public SIMULATION_SCENE {
public:
  virtual void printSceneDescription() override
  {
    cout << "=====================================================================" << endl;
    cout << " Strand grid collision scene" << endl;
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
    _sceneName = "hairball";

    vector<VECTOR3> restVertices;
    vector<vector<int> > strands;

    //readSOBJFile("../data/hairball/straight_1_sparse.sobj", restVertices, strands);
    //readSOBJFile("../data/hairball/straight_3_sparse.sobj", restVertices, strands);
    //readSOBJFile("../data/hairball/straight_103_sparse.sobj", restVertices, strands);
    STRAND_MESH::readSOBJFile("../data/hairball/straight_498_sparse.sobj", restVertices, strands);
    //readSOBJFile("../data/hairball/straight_1000.sobj", restVertices, strands);  // too big to handle
    //readSOBJFile("../data/hairball/strand_250.sobj", restVertices, strands);
    //readSOBJFile("../data/hairball/strand_218.sobj", restVertices, strands);

    // curly
    //readSOBJFile("../data/hairball/curly_125.sobj", restVertices, strands);

    _gravity = VECTOR3(0,-981,0);
    //_gravity = VECTOR3(0,0,0);

    const REAL E = 1e7; // Taz settings
    //const REAL E = 1e5; // misbehaving is quite obvious
    //const REAL E = 1e3;
    //const REAL E = 3.727e10;
    //const REAL E = 1e4;
    //const REAL E = 3.727e10;
    //const REAL E = 3.727e8;
    //const REAL E = 3.727e6;
    const REAL nu = 0.36;
    const REAL density = 1.32;
    const REAL radiusA = 0.005;
    const REAL radiusB = 0.005;

#if 0
    _strandMesh = new STRAND_MESH_FASTER(restVertices, strands,
                                         E, nu, density, radiusA, radiusB);
#else
    _strandMesh = new STRAND_MESH(restVertices, strands,
                                  E, nu, density, radiusA, radiusB);
#endif
    setUnsignedCollisions(10.0,0.1); 
    
    const REAL k = E * M_PI * radiusA * radiusB;
    cout <<" stretching stiffness: " << k << endl;
    STRAND::QUADRATIC_STRETCHING* stretchingEnergy = new STRAND::QUADRATIC_STRETCHING(k);
    _strandSolver = new STRAND::TIMESTEPPER(*_strandMesh, *stretchingEnergy, *_eeGeneral);

    // write out a file for a single strand
    //_strandMesh->writeStrand("strand_218.sobj", 218);

    //_substeps = 10;
    //const REAL dt = 1.0 / 30.0 / (REAL)_substeps;
    _substeps = 1;
    //const REAL dt = 1.0 / 300.0;
    const REAL dt = 1.0 / 30.0;
    _strandSolver->setDt(dt);

    VECTOR3 sphereCenter(27.317198, 32.682804, 27.317198);

    const REAL delta = 5.1;
    _kinematicShapes.push_back(new SPHERE(sphereCenter, 30 + delta));

    // attach hairs to the sphere
    for (unsigned int x = 0; x < _kinematicShapes.size(); x++)
      _strandSolver->attachKinematicSurfaceConstraints(_kinematicShapes[x]);

    // make the sphere a collision object
    _strandSolver->addKinematicCollisionObject(_kinematicShapes.back());

    _strandSolver->collisionsEnabled() = false;
    //_strandSolver->collisionsEnabled() = true;
    //_strandSolver->pcgEnabled() = false;
    _strandSolver->pcgEnabled() = true;
    _strandSolver->hessianClampingEnabled() = true;
    //_strandSolver->collisionStiffness() = 10;
    _strandSolver->collisionDampingBeta() = 0;

    _eye    = VECTOR3(83.4295, 66.3376, -124.49);
    _lookAt = VECTOR3(83.0904, 66.2669, -123.552);
    _up     = VECTOR3(-0.018703, 0.997481, 0.0684159);
    _worldCenter = sphereCenter;

    //_pauseFrame = 3;
    //_pauseFrame = 33;
    //_pauseFrame = 49;
    _pauseFrame = 100;

    _writeToFile = true;

    _filePath = string("../data/render/straight_line_3/");
    string mkdir("mkdir " + _filePath);
    system(mkdir.c_str());

    return true;
  }

  // simulation loop
  virtual void stepSimulation(const bool verbose = true) override
  {
    for (unsigned int x = 0; x < _substeps; x++)
    {
      _strandSolver->externalForces().setZero();
      _strandSolver->addGravity(_gravity);
      _strandSolver->solveDynamics(verbose);
      //_strandSolver->solveNewton(verbose);

      //const int longest = _strandMesh->longestStrand();
      //cout << " Longest is strand " << longest << " which is " << _strandMesh->strandLength(longest) << endl;
    }

    if (_writeToFile) {
      if (!writeFrameToFile()) {
        cerr << "Failed to save frame" << endl;
      }
    }

    _frameNumber++;
  };

  virtual bool writeFrameToFile() override {
    char buffer[256];
    sprintf(buffer, "%04i", _frameNumber);
    string frameName = string("_frame_") + string(buffer) + string(".sobj");
    string filename = _filePath + _sceneName + frameName;

    cout << " Writing file " << filename.c_str() << " ... " << flush;
    if (!STRAND_MESH::writeSOBJFile(filename.c_str(), _strandMesh->vertices(),
                                    _strandMesh->strandIndices())) {
      return false;
    }
    cout << "done." << endl;

    TIMER gzipper("Gzipping");
    cout << " Gzipping ..." << flush;
    string gzip = string("gzip -9 ") + filename + string(" &");
    system(gzip.c_str());
    cout << "done." << endl;
    gzipper.stop();
    return true;
  }

// drawScene MUST have this guard, so it can be deactivated for regression testing
#ifndef GL_DISABLED
  // what to draw?
  virtual void drawScene()
  {
    drawAxes();
    //drawStrand(*_strandMesh, 250);
    //drawStrandMeshOld(*_strandMesh, 250);
    drawStrandMesh(*_strandMesh);
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
};

}

#endif

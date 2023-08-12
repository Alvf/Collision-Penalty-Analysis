#ifndef HAIRBALL_UNSTABLE_H
#define HAIRBALL_UNSTABLE_H

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

// caching the unstable version of the hairball here for illustrative purposes later
class HAIRBALL_UNSTABLE : public SIMULATION_SCENE {
public:
  virtual void printSceneDescription() override
  {
    cout << "=====================================================================" << endl;
    cout << " Unstable hair ball scene" << endl;
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
    _sceneName = "hairball_unstable";

    vector<VECTOR3> restVertices;
    vector<vector<int> > strands;

    STRAND_MESH::readSOBJFile("../data/hairball/straight_498_sparse.sobj", restVertices, strands);

    _gravity = VECTOR3(0,-981,0);

    //const REAL E = 1e7; // Taz settings
    //const REAL E = 1e10;
    const REAL E = 3.727e10;
    const REAL nu = 0.36;
    const REAL density = 1.32;
    const REAL radiusA = 0.005;
    const REAL radiusB = 0.005;

#if 1
    _strandMesh = new STRAND_MESH_FASTER(restVertices, strands,
                                         E, nu, density, radiusA, radiusB);
#else
    _strandMesh = new STRAND_MESH(restVertices, strands,
                                  E, nu, density, radiusA, radiusB);
#endif
    setUnsignedCollisions(10.0,0.1);

    //STRAND::QUADRATIC_STRETCHING* stretchingEnergy = new STRAND::QUADRATIC_STRETCHING(100.0);
    const REAL k = E * M_PI * radiusA * radiusB;
    cout <<" stretching stiffness: " << k << endl;
    STRAND::QUADRATIC_STRETCHING* stretchingEnergy = new STRAND::QUADRATIC_STRETCHING(k);
    _strandSolver = new STRAND::TIMESTEPPER(*_strandMesh, *stretchingEnergy, *_eeGeneral);

    // write out a file for a single strand
    //_strandMesh->writeStrand("strand_218.sobj", 218);

    _substeps = 1;
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
    _strandSolver->pcgEnabled() = false;

    // stable
    //_strandSolver->hessianClampingEnabled() = true;
    //_strandMesh->bendingForceFilterEnabled() = true;

    // totally unstable
    //_strandSolver->hessianClampingEnabled() = false;
    //_strandMesh->bendingForceFilterEnabled() = false;
    
    // still not that stable
    _strandSolver->hessianClampingEnabled() = true;
    _strandMesh->bendingForceFilterEnabled() = false;

    _eye    = VECTOR3(83.4295, 66.3376, -124.49);
    _lookAt = VECTOR3(83.0904, 66.2669, -123.552);
    _up     = VECTOR3(-0.018703, 0.997481, 0.0684159);
    _worldCenter = sphereCenter;

    _pauseFrame = 60;

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

      const int longest = _strandMesh->longestStrand();
      cout << " Longest is strand " << longest << " which is " << _strandMesh->strandLength(longest) << endl;
    }

    _frameNumber++;
  };

// drawScene MUST have this guard, so it can be deactivated for regression testing
#ifndef GL_DISABLED
  // what to draw?
  virtual void drawScene()
  {
    drawAxes();
    //drawStrandMesh(*_strandMesh);
    drawStrand(*_strandMesh, 289);

    glEnable(GL_DEPTH_TEST);
    
    drawCollisionsOld(*_strandMesh);
  };
#endif

protected:
  STRAND_MESH* _strandMesh;
  STRAND::TIMESTEPPER* _strandSolver;
  STRAND::STRETCHING* _strechingEnergy;

  unsigned int _substeps;
};

}

#endif

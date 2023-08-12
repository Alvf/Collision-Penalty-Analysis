#ifndef STRAIGHT_STRAND_SCENE_H
#define STRAIGHT_STRAND_SCENE_H

#include "Scenes/SIMULATION_SCENE.h"
#include "Geometry/STRAND_MESH_DVT.h"
#include "Timestepper/Strand/TIMESTEPPER_DVT.h"
#include "Hyperelastic/Strand/QUADRATIC_STRETCHING.h"
#include "Hyperelastic/Strand/ISOTROPIC_BENDING.h"
#include "Hyperelastic/Strand/QUADRATIC_BENDING.h"
#include "Hyperelastic/Strand/QUADRATIC_UNIT_BENDING.h"
#include "Hyperelastic/Strand/ISOTROPIC_THETA.h"

namespace HOBAK {

class STRAIGHT_STRAND_SCENE : public SIMULATION_SCENE {
public:
  virtual void printSceneDescription() override
  {
    cout << "=====================================================================" << endl;
    cout << " Straight strand scene" << endl;
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
    _sceneName = "straight_strand";

    // procedurally generate a scene
    std::vector<VECTOR3> restVertices;

    restVertices.push_back(VECTOR3(-3.888,4.64077,-4.85428));
    restVertices.push_back(VECTOR3(-4.38894,5.06663,-5.05925));
    restVertices.push_back(VECTOR3(-4.89188,5.2493,-5.49283));
    restVertices.push_back(VECTOR3(-5.11457,5.50839,-6.09083));
    restVertices.push_back(VECTOR3(-5.24247,6.03597,-6.51464));
    restVertices.push_back(VECTOR3(-5.61924,6.57603,-6.71637));
    restVertices.push_back(VECTOR3(-6.1732,6.85173,-7.01876));
    restVertices.push_back(VECTOR3(-6.53695,7.02742,-7.57656));
    restVertices.push_back(VECTOR3(-6.66054,7.43436,-8.11826));
    restVertices.push_back(VECTOR3(-6.89482,8.01788,-8.39922));
    restVertices.push_back(VECTOR3(-7.40592,8.4286,-8.60995));
    restVertices.push_back(VECTOR3(-7.89796,8.60554,-9.05817));
    restVertices.push_back(VECTOR3(-8.1069,8.87787,-9.65524));
    restVertices.push_back(VECTOR3(-8.24114,9.41557,-10.0641));
    restVertices.push_back(VECTOR3(-8.63377,9.94573,-10.2618));
    restVertices.push_back(VECTOR3(-9.18663,10.208,-10.5778));
    restVertices.push_back(VECTOR3(-9.53417,10.3892,-11.1441));
    restVertices.push_back(VECTOR3(-9.65348,10.8113,-11.6751));
    restVertices.push_back(VECTOR3(-9.90256,11.3944,-11.9439));
    restVertices.push_back(VECTOR3(-10.4228,11.7898,-12.1614));
    restVertices.push_back(VECTOR3(-10.9031,11.9621,-12.624));
    restVertices.push_back(VECTOR3(-11.099,12.2482,-13.219));
    restVertices.push_back(VECTOR3(-11.2407,12.7951,-13.6129));
    restVertices.push_back(VECTOR3(-11.6488,13.3145,-13.8077));
    restVertices.push_back(VECTOR3(-12.1994,13.5641,-14.1377));
    restVertices.push_back(VECTOR3(-12.5306,13.7518,-14.7116));
    restVertices.push_back(VECTOR3(-12.6468,14.1888,-15.231));
    restVertices.push_back(VECTOR3(-12.9112,14.7703,-15.4885));
    restVertices.push_back(VECTOR3(-13.4396,15.1503,-15.7137));
    restVertices.push_back(VECTOR3(-13.9073,15.319,-16.1902));
    restVertices.push_back(VECTOR3(-14.0909,15.6194,-16.7822));
    restVertices.push_back(VECTOR3(-14.241,16.1746,-17.161));
    restVertices.push_back(VECTOR3(-14.6643,16.6825,-17.354));
    restVertices.push_back(VECTOR3(-15.2114,16.9199,-17.6984));
    restVertices.push_back(VECTOR3(-15.5264,17.1151,-18.2789));
    restVertices.push_back(VECTOR3(-15.6407,17.5668,-18.7861));
    restVertices.push_back(VECTOR3(-15.9207,18.1456,-19.0329));
    restVertices.push_back(VECTOR3(-16.4562,18.51,-19.2668));
    restVertices.push_back(VECTOR3(-16.9106,18.6764,-19.7569));

    // vector of all the strands
    vector<vector<int> > strandIndices;

    // store the indices of individual strands
    vector<int> strand0;
    for (unsigned int x = 0; x < restVertices.size(); x++)
      strand0.push_back(x);
    strandIndices.push_back(strand0);

    // center the strand
    VECTOR3 origin = restVertices[0];
    for (unsigned int x = 0; x < restVertices.size(); x++)
      restVertices[x] -= origin;

    _gravity = VECTOR3(0,981,0);
    //_gravity = VECTOR3(0,0,0);

    vector<VECTOR3> vertices = restVertices;
    //for (unsigned int x = 0; x < vertices.size(); x++)
    //  vertices[x][1] = 0;

    const REAL E = 3.727e10;
    const REAL nu = 0.36;
    const REAL density = 1.32;
    const REAL radiusA = 0.005;
    const REAL radiusB = 0.005;
    _strandMesh = new STRAND_MESH_DVT(restVertices, vertices,
                                      E, nu, density, radiusA, radiusB);

    // create the integrator
#if 0
    const REAL k = E * M_PI * radiusA * radiusB;
    cout <<" stretching stiffness: " << k << endl;
    STRAND::QUADRATIC_STRETCHING* stretchingEnergy = new STRAND::QUADRATIC_STRETCHING(k);
#else
    STRAND::QUADRATIC_STRETCHING* stretchingEnergy = new STRAND::QUADRATIC_STRETCHING(100.0);
#endif
    _strandSolver = new STRAND::TIMESTEPPER_DVT(*_strandMesh, *stretchingEnergy);
    _strandSolver->setDt(1.0 / 30.0);
    //_strandSolver->setDt(1.0 / 100.0);
    //_strandSolver->setDt(1.0 / 700.0);

    // constrain the first vertex
    VECTOR3 center0(0,0, 0.0);
    _kinematicShapes.push_back(new CUBE(center0, 1.25));
    _strandSolver->attachKinematicSurfaceConstraints(_kinematicShapes[0]);

    _eye    = VECTOR3(3.79907, 1.29043, 4.21273);
    _lookAt = VECTOR3(3.37669, 1.07395, 3.33254);
    _up     = VECTOR3(-0.19182, 0.970417, -0.146615);

    _eye    = VECTOR3(5.6, 0, 0.);
    _lookAt = VECTOR3(4.6, 0.0, 0.0);
    _up     = VECTOR3(0,1,0);
    
    _eye = VECTOR3( 2.2056111, -0.52623286, -1.5676666 );
    _up = VECTOR3( -0.23574428, -0.95587433, 0.17529663 );
    _lookAt = VECTOR3( 0.25, 0.22275163, -0.11349762 );

    _eye    = VECTOR3(7.53, -3.73121, -32.0716);
    _lookAt = VECTOR3(7.08389, -3.40145, -31.2396);
    _up     = VECTOR3(-0.134503, -0.943787, 0.301943);  

    _pauseFrame = 100;

    return true;
  }

  // simulation loop
  virtual void stepSimulation(const bool verbose = true) override 
  {
    _strandSolver->externalForces().setZero();
    _strandSolver->addGravity(_gravity);
    _strandSolver->solveDynamics(verbose);

    _frameNumber++;
  };

// drawScene MUST have this guard, so it can be deactivated for regression testing
#ifndef GL_DISABLED
  // what to draw?
  virtual void drawScene()
  {
    drawAxes();
    drawStrandMesh(*_strandMesh);

    glEnable(GL_DEPTH_TEST);
    for (unsigned int x = 0; x < _kinematicShapes.size(); x++)
      drawKinematicShape(*_kinematicShapes[x]);

    /*
    drawKinematicConstraints(_strandMesh, _strandSolver);
    //drawPlaneConstraints(_triangleMesh, _strandSolver);
    drawStrandTwistFreeFrames(*_strandMesh);
    drawStrandTwistFrames(*_strandMesh);
    drawStrandTwistForces(*_strandMesh);
    */
  };
#endif

protected:
  //STRAND_MESH* _strandMesh;
  STRAND_MESH_DVT* _strandMesh;
  STRAND::TIMESTEPPER_DVT* _strandSolver;
  STRAND::STRETCHING* _strechingEnergy;

  //string _strandMeshFilename;
};

}

#endif

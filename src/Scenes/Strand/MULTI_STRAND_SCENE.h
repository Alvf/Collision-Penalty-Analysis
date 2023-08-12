#ifndef MULTI_STRAND_SCENE_H
#define MULTI_STRAND_SCENE_H

#include "Scenes/SIMULATION_SCENE.h"
#include "Geometry/STRAND_MESH_DVT.h"
#include "Timestepper/Strand/TIMESTEPPER_DVT.h"
#include "Hyperelastic/Strand/QUADRATIC_STRETCHING.h"
#include "Hyperelastic/Strand/ISOTROPIC_BENDING.h"
#include "Hyperelastic/Strand/QUADRATIC_BENDING.h"
#include "Hyperelastic/Strand/QUADRATIC_UNIT_BENDING.h"
#include "Hyperelastic/Strand/ISOTROPIC_THETA.h"

namespace HOBAK {

class MULTI_STRAND_SCENE : public SIMULATION_SCENE {
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

#if 1
  // for testing stretch
  virtual bool buildScene() override
  {
    // this will determine the MOV and JSON filenames
    _sceneName = "multi_strand";

    // procedurally generate a scene
    std::vector<VECTOR3> restVertices;
    std::vector<VECTOR3> vertices;

    REAL xShift = 0;
    REAL zShift = 0;

    restVertices.push_back(VECTOR3(xShift,0,zShift));
    restVertices.push_back(VECTOR3(xShift,1,zShift));
    restVertices.push_back(VECTOR3(xShift,2,zShift));
    restVertices.push_back(VECTOR3(xShift,3,zShift));

    vertices.push_back(VECTOR3(0,xShift, zShift));
    vertices.push_back(VECTOR3(1,xShift, zShift));
    vertices.push_back(VECTOR3(2,xShift, zShift));
    vertices.push_back(VECTOR3(3,xShift, zShift));

    // vector of all the strands
    vector<vector<int> > strandIndices;

    // store the indices of individual strands
    vector<int> strand0;
    for (unsigned int x = 0; x < restVertices.size(); x++)
      strand0.push_back(x);
    strandIndices.push_back(strand0);

    xShift = 0;
    zShift = 0.2;
    const int nextStart = restVertices.size();
    restVertices.push_back(VECTOR3(xShift,0,zShift));
    restVertices.push_back(VECTOR3(xShift,1,zShift));
    restVertices.push_back(VECTOR3(xShift,2,zShift));
    restVertices.push_back(VECTOR3(xShift,3,zShift));
    
    vertices.push_back(VECTOR3(xShift,0.12,0 + zShift));
    vertices.push_back(VECTOR3(xShift,0.12,1 + zShift));
    vertices.push_back(VECTOR3(xShift,0.12,2 + zShift));
    vertices.push_back(VECTOR3(xShift,0.12,3 + zShift));

    // store the indices of individual strands
    vector<int> strand1;
    for (unsigned int x = nextStart ; x < restVertices.size(); x++)
      strand1.push_back(x);
    strandIndices.push_back(strand1);

    _gravity = VECTOR3(0,981,0);
    //_gravity = VECTOR3(0,0,0);

    const REAL E = 3.727e10;
    const REAL nu = 0.36;
    const REAL density = 1.32;
    const REAL radiusA = 0.005;
    const REAL radiusB = 0.005;
    _strandMesh = new STRAND_MESH_DVT(restVertices, vertices, strandIndices,
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

    // constrain the first vertex
    VECTOR3 center0(0,0, 0.0);
    _kinematicShapes.push_back(new CUBE(center0, 1.25));
    _strandSolver->attachKinematicSurfaceConstraints(_kinematicShapes[0]);

    _eye    = VECTOR3(3.72997, -0.0834072, -7.74737);
    _lookAt = VECTOR3(3.28849, 0.0942485, -6.86786);
    _up     = VECTOR3(-0.0702835, -0.984034, 0.16349);

    _pauseFrame = 100;

    return true;
  }
#else
  // for testing stretch
  virtual bool buildScene() override
  {
    // this will determine the MOV and JSON filenames
    _sceneName = "multi_strand";

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

    // center to origin
    const VECTOR3 v0 = restVertices[0];
    for (unsigned int x = 0; x < restVertices.size(); x++)
      restVertices[x] -= v0;

    // vector of all the strands
    vector<vector<int> > strandIndices;

    // store the indices of individual strands
    vector<int> strand0;
    for (unsigned int x = 0; x < restVertices.size(); x++)
      strand0.push_back(x);
    strandIndices.push_back(strand0);

    const int nextStart = restVertices.size();
   
    restVertices.push_back(VECTOR3(-2.83364,-6.23249,3.65291));
    restVertices.push_back(VECTOR3(-3.11762,-6.85697,3.71375));
    restVertices.push_back(VECTOR3(-3.12058,-7.50261,3.95344));
    restVertices.push_back(VECTOR3(-3.21184,-7.96367,4.45685));
    restVertices.push_back(VECTOR3(-3.61336,-8.34244,4.8687));
    restVertices.push_back(VECTOR3(-4.0293,-8.87857,4.98654));
    restVertices.push_back(VECTOR3(-4.13975,-9.54906,5.09858));
    restVertices.push_back(VECTOR3(-4.13464,-10.1067,5.50272));
    restVertices.push_back(VECTOR3(-4.39676,-10.4926,6.00936));
    restVertices.push_back(VECTOR3(-4.86291,-10.9345,6.25785));
    restVertices.push_back(VECTOR3(-5.1294,-11.5667,6.31835));
    restVertices.push_back(VECTOR3(-5.12573,-12.2058,6.57478));
    restVertices.push_back(VECTOR3(-5.2323,-12.6571,7.08405));
    restVertices.push_back(VECTOR3(-5.64551,-13.039,7.48112));
    restVertices.push_back(VECTOR3(-6.05001,-13.5861,7.58819));
    restVertices.push_back(VECTOR3(-6.14499,-14.257,7.71145));
    restVertices.push_back(VECTOR3(-6.1462,-14.8038,8.13007));
    restVertices.push_back(VECTOR3(-6.42588,-15.1858,8.63029));
    restVertices.push_back(VECTOR3(-6.89153,-15.6371,8.86217));
    restVertices.push_back(VECTOR3(-7.14029,-16.2764,8.9236));
    restVertices.push_back(VECTOR3(-7.13125,-16.9085,9.19694));
    restVertices.push_back(VECTOR3(-7.25377,-17.3502,9.71089));
    restVertices.push_back(VECTOR3(-7.67763,-17.7362,10.0926));
    restVertices.push_back(VECTOR3(-8.06971,-18.294,10.1899));
    restVertices.push_back(VECTOR3(-8.14992,-18.9645,10.3253));
    restVertices.push_back(VECTOR3(-8.15866,-19.5004,10.7577));
    restVertices.push_back(VECTOR3(-8.4556,-19.8792,11.2503));
    restVertices.push_back(VECTOR3(-8.91944,-20.3404,11.4659));
    restVertices.push_back(VECTOR3(-9.15033,-20.9861,11.5295));
    restVertices.push_back(VECTOR3(-9.1372,-21.6105,11.8199));
    restVertices.push_back(VECTOR3(-9.27622,-22.0432,12.3373));
    restVertices.push_back(VECTOR3(-9.70964,-22.4339,12.7031));
    restVertices.push_back(VECTOR3(-10.0884,-23.0023,12.7917));
    restVertices.push_back(VECTOR3(-10.1546,-23.6715,12.9401));
    restVertices.push_back(VECTOR3(-10.1721,-24.1965,13.3855));
    restVertices.push_back(VECTOR3(-10.4859,-24.5728,13.8695));
    restVertices.push_back(VECTOR3(-10.9466,-25.0442,14.0692));
    restVertices.push_back(VECTOR3(-11.1596,-25.6957,14.1362));
    restVertices.push_back(VECTOR3(-11.1436,-26.3118,14.4436));

    // center to origin
    const VECTOR3 v1 = restVertices[nextStart];
    for (unsigned int x = nextStart; x < restVertices.size(); x++)
      restVertices[x] -= v1;

    // spin into collision with the other strand
    const MATRIX3 R = rotationMatrix(VECTOR3(0,0,1), M_PI * 0.333);
    for (unsigned int x = nextStart; x < restVertices.size(); x++)
      restVertices[x] = R * restVertices[x];

    // store the indices of individual strands
    vector<int> strand1;
    for (unsigned int x = nextStart ; x < restVertices.size(); x++)
      strand1.push_back(x);
    strandIndices.push_back(strand1);

    _gravity = VECTOR3(0,981,0);
    //_gravity = VECTOR3(0,0,0);

    vector<VECTOR3> vertices = restVertices;
    for (unsigned int x = 0; x < vertices.size(); x++)
      vertices[x][1] = 0;

    const REAL E = 3.727e10;
    const REAL nu = 0.36;
    const REAL density = 1.32;
    const REAL radiusA = 0.005;
    const REAL radiusB = 0.005;
    _strandMesh = new STRAND_MESH_DVT(restVertices, vertices, strandIndices,
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
    _kinematicShapes.push_back(new CUBE(center0, 0.25));
    _strandSolver->attachKinematicSurfaceConstraints(_kinematicShapes[0]);

    _eye    = VECTOR3(3.72997, -0.0834072, -7.74737);
    _lookAt = VECTOR3(3.28849, 0.0942485, -6.86786);
    _up     = VECTOR3(-0.0702835, -0.984034, 0.16349);

    _pauseFrame = 100;

    return true;
  }
#endif

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

    drawCollisions(*_strandMesh);

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

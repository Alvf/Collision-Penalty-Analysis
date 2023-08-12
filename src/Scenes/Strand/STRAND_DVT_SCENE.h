#ifndef STRAND_DVT_SCENE_H
#define STRAND_DVT_SCENE_H

#include "Scenes/SIMULATION_SCENE.h"
#include "Geometry/STRAND_MESH_DVT.h"
#include "Timestepper/Strand/TIMESTEPPER_DVT.h"
#include "Hyperelastic/Strand/QUADRATIC_STRETCHING.h"
#include "Hyperelastic/Strand/ISOTROPIC_BENDING.h"
#include "Hyperelastic/Strand/QUADRATIC_BENDING.h"
#include "Hyperelastic/Strand/QUADRATIC_UNIT_BENDING.h"
#include "Hyperelastic/Strand/ISOTROPIC_THETA.h"

namespace HOBAK {

class STRAND_DVT_SCENE : public SIMULATION_SCENE {
public:
  virtual void printSceneDescription() override
  {
    cout << "=====================================================================" << endl;
    cout << " Debugging strand scene" << endl;
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
    _sceneName = "strand_dvt";

    // procedurally generate a scene
    std::vector<VECTOR3> restVertices;

    // L shape
    restVertices.push_back(VECTOR3(-1,0,0));
    restVertices.push_back(VECTOR3(0,0,0));
    restVertices.push_back(VECTOR3(1,0,0));
    restVertices.push_back(VECTOR3(1,0,-1));

    //_gravity = VECTOR3(0,100,0);
    _gravity = VECTOR3(0,0,0);

    vector<VECTOR3> vertices = restVertices;
    //const REAL angle = M_PI * 0.5;
    //vertices[3] = VECTOR3(1,sin(angle), -cos(angle)); // all the way

    vertices[2] = VECTOR3(0,1,0);
    vertices[3] = VECTOR3(1,1,0);

    const REAL E = 10;
    const REAL nu = 0.025 - 1.0;  // non-physical, but good for debugging
    const REAL density = 1e-3;
    const REAL radiusA = 1;
    const REAL radiusB = 1;
    _strandMesh = new STRAND_MESH_DVT(restVertices, vertices,
                                      E, nu, density, radiusA, radiusB);
    /*
    const int maxElements = restVertices.size();
    VECTOR restTwists(maxElements - 1);
    restTwists.setZero();

    // initialize the strand mesh
    _strandMesh = new STRAND_MESH(restVertices, restTwists);
   
    // give it a twist 
    std::vector<VECTOR3> vertices;
    vertices = restVertices;
    vertices[3] = VECTOR3(1,0,-1); // symmetric, preferred
    //vertices[3] = VECTOR3(2,-1,1); // not symmetric
    //vertices[3] = VECTOR3(1,-1,1); // symmetric
    //vertices[3] = VECTOR3(1,-1,-1); // symmetric
    //vertices[3] = VECTOR3(1.5,0,-1);  // not symmetric, debug with this one. Why is the angle so far off?
    //vertices[3] = VECTOR3(1.25,0,-1);  // not symmetric
    //vertices[3] = VECTOR3(1.75,0,-1);  // not symmetric, way off
    //vertices[3] = VECTOR3(2.0,0,-1);  // not symmetric, way off
    //_strandMesh->setPositions(vertices);

    // the material needs to see the twist too
    VECTOR twistAngles(_strandMesh->totalEdges());
    twistAngles.setZero();
    //twistAngles[1] = -M_PI * 0.5;
    //_strandMesh->setTwistAngles(twistAngles);
    _strandMesh->setPositionsAndTwistAngles(vertices, twistAngles);

    cout << " Twists: " << endl << _strandMesh->twistAngles() << endl;

    VECTOR& restBendAngles = _strandMesh->restBendAngles();
    for (int x = 0; x < restBendAngles.size(); x++)
      restBendAngles[x] = 0.5 * M_PI;
    restBendAngles[0] = M_PI;
    */

    // create the integrator
    STRAND::QUADRATIC_STRETCHING* stretchingEnergy = new STRAND::QUADRATIC_STRETCHING(100.0);
    _strandSolver = new STRAND::TIMESTEPPER_DVT(*_strandMesh, *stretchingEnergy);
    _strandSolver->setDt(1.0 / 100.0);

    // constrain the first vertex
    VECTOR3 center0(-1.0,-0.2, 0.0);
    _kinematicShapes.push_back(new CUBE(center0, 2.1));
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

    _pauseFrame = 9;

    /*
    _strandMesh->computeBendingForces();
    _strandMesh->computeBendingClampedHessian();
    _strandMesh->computeTwistingForces();
    _strandMesh->computeTwistingClampedHessian();
    */

    return true;
  }

  // simulation loop
  virtual void stepSimulation(const bool verbose = true) override 
  {
#if 1
    _strandSolver->solveQuasistatic(verbose);
#else
    _strandSolver->externalForces().setZero();
    _strandSolver->addGravity(_gravity);
    _strandSolver->solveDynamics(verbose);
#endif

#if 1
    if (_frameNumber == 9)
    {
      VECTOR positionsKnown(12);
      positionsKnown << -1, 0, 0, 0, 0, 0, 0.9994173123909491, 0.1243500915608951, -0.2291030282394728, 0.8139664087434525, 0.1410819371218592, -1.225833006238255;
      VECTOR thetasKnown(3);
      thetasKnown << 0, 0.7256850604191948, -1.252920931479943;

      VECTOR positions = _strandMesh->getPositions();
      VECTOR thetas = _strandMesh->thetas();

      const REAL positionsDiff = (positionsKnown - positions).norm();
      const REAL thetasDiff = (thetasKnown - thetas).norm();

      cout << " position diff: " << positionsDiff << endl;
      cout << " theta diff: " << thetasDiff << endl; 
      if (positionsDiff + thetasDiff < 1e-8)
        cout << " Regression test PASSED" << endl;
      else
        cout << " Regression test FAILED" << endl;

      assert(positionsDiff < 1e-8);
      assert(thetasDiff < 1e-8);
    }
#endif
    
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

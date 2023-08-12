#ifndef STRAND_LINE_SCENE_H
#define STRAND_LINE_SCENE_H

#include "Scenes/SIMULATION_SCENE.h"
#include "Geometry/STRAND_MESH_DVT.h"
#include "Timestepper/Strand/TIMESTEPPER_DVT.h"
#include "Hyperelastic/Strand/QUADRATIC_STRETCHING.h"
#include "Hyperelastic/Strand/ISOTROPIC_BENDING.h"
#include "Hyperelastic/Strand/QUADRATIC_BENDING.h"
#include "Hyperelastic/Strand/QUADRATIC_UNIT_BENDING.h"
#include "Hyperelastic/Strand/ISOTROPIC_THETA.h"

namespace HOBAK {

class STRAND_LINE_SCENE : public SIMULATION_SCENE {
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
    _sceneName = "strand_line";

    // procedurally generate a scene
    std::vector<VECTOR3> restVertices;

    // L shape
    restVertices.push_back(VECTOR3(-1,0,0));
    restVertices.push_back(VECTOR3(0,0,0));
    restVertices.push_back(VECTOR3(1,0,0));
    //restVertices.push_back(VECTOR3(1,1,0));
    //restVertices.push_back(VECTOR3(0,1,0));

    //_gravity = VECTOR3(0,-0.1,0);

    const int maxElements = restVertices.size();
    VECTOR restTwists(maxElements - 1);
    restTwists.setZero();

    // initialize the strand mesh
    _strandMesh = new STRAND_MESH_DVT(restVertices, restTwists);
   
    // give it a twist 
    std::vector<VECTOR3> vertices;
    vertices = restVertices;
    //vertices[3] = VECTOR3(2,-1,1); // not symmetric
    //vertices[3] = VECTOR3(1,-1,1); // symmetric
    //vertices[3] = VECTOR3(1,-1,-1); // symmetric
    //vertices[3] = VECTOR3(1,0,-1); // symmetric, easy
    //vertices[3] = VECTOR3(1.5,0,-1);  // not symmetric, debug with this one. Why is the angle so far off?
    //vertices[3] = VECTOR3(1.25,0,-1);  // not symmetric
    //vertices[3] = VECTOR3(1.75,0,-1);  // not symmetric, way off
    //vertices[3] = VECTOR3(2.0,0,-1);  // not symmetric, way off
    _strandMesh->setPositions(vertices);

    // the material needs to see the twist too
    VECTOR twistAngles(_strandMesh->totalEdges());
    twistAngles.setZero();
    twistAngles[1] = M_PI * 0.5;
    _strandMesh->setTwistAngles(twistAngles);

    cout << " Twists: " << endl << _strandMesh->twistAngles() << endl;

    VECTOR& restBendAngles = _strandMesh->restBendAngles();
    //for (int x = 0; x < restBendAngles.size(); x++)
    //  restBendAngles[x] = 0.5 * M_PI;
    restBendAngles[0] = M_PI;

    STRAND::QUADRATIC_STRETCHING* stretchingEnergy = new STRAND::QUADRATIC_STRETCHING(1000.0);
    const REAL currentStretch = _strandMesh->computeStretchingEnergy(*stretchingEnergy);
    cout << " Stretching energy: " << currentStretch << endl;

    STRAND::ISOTROPIC_BENDING* isotropicBendingEnergy = new STRAND::QUADRATIC_UNIT_BENDING(10000.0);
    
    // create the integrator
    _strandSolver = new STRAND::TIMESTEPPER_DVT(*_strandMesh, *stretchingEnergy, *isotropicBendingEnergy);

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

    _eye    = VECTOR3(1.59594, 1.46548, 2.27851);
    _lookAt = VECTOR3(1.1492, 0.958721, 1.54121);
    _up     = VECTOR3(-0.329399, 0.859388, -0.391086);    
    return true;
  }

  // simulation loop
  virtual void stepSimulation(const bool verbose = true) override {
    _strandSolver->externalForces().setZero();
    _strandSolver->addGravity(_gravity);
    //_strandSolver->solve(verbose);
    _strandSolver->solveWithTwist(verbose);

    _frameNumber++;
  };

// drawScene MUST have this guard, so it can be deactivated for regression testing
#ifndef GL_DISABLED
  // what to draw?
  virtual void drawScene()
  {
    //drawAxes();
    drawStrandMesh(*_strandMesh);

    glEnable(GL_DEPTH_TEST);
    for (unsigned int x = 0; x < _kinematicShapes.size(); x++)
      drawKinematicShape(*_kinematicShapes[x]);

    drawKinematicConstraints(_strandMesh, _strandSolver);
    //drawPlaneConstraints(_triangleMesh, _strandSolver);
    drawStrandTwistFreeFrames(*_strandMesh);
    drawStrandTwistFrames(*_strandMesh);
    //drawStrandTwistForces(*_strandMesh);
  };
#endif

protected:
  STRAND_MESH* _strandMesh;
  STRAND::TIMESTEPPER* _strandSolver;
  STRAND::STRETCHING* _strechingEnergy;

  //string _strandMeshFilename;
};

}

#endif

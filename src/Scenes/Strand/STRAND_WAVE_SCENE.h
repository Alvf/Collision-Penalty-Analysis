#ifndef STRAND_WAVE_SCENE_H
#define STRAND_WAVE_SCENE_H

#include "Scenes/SIMULATION_SCENE.h"
#include "Geometry/STRAND_MESH_DVT.h"
#include "Timestepper/Strand/TIMESTEPPER_DVT.h"
#include "Hyperelastic/Strand/QUADRATIC_STRETCHING.h"
#include "Hyperelastic/Strand/ISOTROPIC_BENDING.h"
#include "Hyperelastic/Strand/QUADRATIC_BENDING.h"
#include "Hyperelastic/Strand/QUADRATIC_UNIT_BENDING.h"
#include "Hyperelastic/Strand/ISOTROPIC_THETA.h"

namespace HOBAK {

class STRAND_WAVE_SCENE : public SIMULATION_SCENE {
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
    _sceneName = "strand";

    // procedurally generate a scene
    std::vector<VECTOR3> restVertices;
    std::vector<VECTOR3> vertices;

    // line
    for (int x = 0; x < 15; x++)
      restVertices.push_back(VECTOR3(-1 + x,0.01,0));

    vertices = restVertices;

    _gravity = VECTOR3(0,-0.1,0);

    const int maxElements = restVertices.size();
    VECTOR restThetas(maxElements - 1);
    restThetas.setZero();

    // initialize the strand mesh
    _strandMesh = new STRAND_MESH_DVT(restVertices, restThetas);
    _strandMesh->setPositions(vertices);

    VECTOR& restBendAngles = _strandMesh->restBendAngles();
    for (int x = 0; x < restBendAngles.size(); x++)
      restBendAngles[x] = 0.5 * M_PI;
    restBendAngles[0] = M_PI;

    STRAND::QUADRATIC_STRETCHING* stretchingEnergy = new STRAND::QUADRATIC_STRETCHING(1000.0);
    const REAL currentStretch = _strandMesh->computeStretchingEnergy(*stretchingEnergy);
    cout << " Stretching energy: " << currentStretch << endl;

    // create the integrator
    STRAND::ISOTROPIC_BENDING* isotropicBendingEnergy = new STRAND::QUADRATIC_UNIT_BENDING(10000.0);
    _strandSolver = new STRAND::TIMESTEPPER_DVT(*_strandMesh, *stretchingEnergy, *isotropicBendingEnergy);

    // constrain the first vertex
    VECTOR3 center0(-1.0,-0.2, 0.0);
    _kinematicShapes.push_back(new CUBE(center0, 2.1));
    _strandSolver->attachKinematicSurfaceConstraints(_kinematicShapes[0]);

    _eye    = VECTOR3(7.72177, 2.75295, 10.9793);
    _lookAt = VECTOR3(7.38745, 2.57138, 10.0545);
    _up     = VECTOR3(-0.185233, 0.974786, -0.124423);

    return true;
  }

  // simulation loop
  virtual void stepSimulation(const bool verbose = true) override {
    _strandSolver->externalForces().setZero();
    _strandSolver->addGravity(_gravity);
    _strandSolver->solve(verbose);

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

    /*
    drawVertex(*_triangleMesh, _arrowCounter);
    */
    /*
    SIMULATION_SCENE::drawScene();

    if (_drawFeature)
    {
      glDisable(GL_DEPTH_TEST);
      glPointSize(10.0);
      drawVertexFacePairs(*_tetMesh, _arrowCounter);
      //drawVertexFacePairs(*_tetMesh);
    }
    */
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

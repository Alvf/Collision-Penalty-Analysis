#ifndef STRAND_SCENE_H
#define STRAND_SCENE_H

#include "Scenes/SIMULATION_SCENE.h"
#include "Hyperelastic/Strand/QUADRATIC_STRETCHING.h"
#include "Hyperelastic/Strand/ISOTROPIC_BENDING.h"
#include "Hyperelastic/Strand/QUADRATIC_BENDING.h"
#include "Hyperelastic/Strand/QUADRATIC_UNIT_BENDING.h"
#include "Hyperelastic/Strand/ISOTROPIC_THETA.h"

namespace HOBAK {

class STRAND_SCENE : public SIMULATION_SCENE {
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

    // straight line
    restVertices.push_back(VECTOR3(0,0,0));
    restVertices.push_back(VECTOR3(1,0,0));
    restVertices.push_back(VECTOR3(2,0,0));

#if 0
    // C shape
    restVertices.push_back(VECTOR3(-1,0,0));
    restVertices.push_back(VECTOR3(0,0,0));
    restVertices.push_back(VECTOR3(1,0,0));
    restVertices.push_back(VECTOR3(1,1,0));
    restVertices.push_back(VECTOR3(0,1,0));
#endif

    vertices = restVertices;
    //vertices[1] *= 1.0;
    //vertices[2] *= 1.0;

    //_gravity = VECTOR3(0.1,0,0);
    _gravity = VECTOR3(0,-0.1,0);
    //_gravity = VECTOR3(0,0,0);

    const int maxElements = restVertices.size();
    VECTOR restTwists(maxElements - 1);
    restTwists.setZero();

    // initialize the strand mesh
    _strandMesh = new STRAND_MESH_DVT(restVertices, restTwists);
    _strandMesh->setPositions(vertices);

    VECTOR& restBendAngles = _strandMesh->restBendAngles();
    for (int x = 0; x < restBendAngles.size(); x++)
      restBendAngles[x] = 0.5 * M_PI;
    restBendAngles[0] = M_PI;

    //STRAND::QUADRATIC_STRETCHING* stretchingEnergy = new STRAND::QUADRATIC_STRETCHING(1.0);
    STRAND::QUADRATIC_STRETCHING* stretchingEnergy = new STRAND::QUADRATIC_STRETCHING(100.0);
    const REAL currentStretch = _strandMesh->computeStretchingEnergy(*stretchingEnergy);
    cout << " Stretching energy: " << currentStretch << endl;

    //STRAND::BERGOU_2010 bendingEnergy(1.0);
    //const REAL currentBending = _strandMesh->computeBendingEnergy(bendingEnergy);
    //cout << " Bending energy: " << currentBending << endl;

    cout << " Current positions: " << endl;
    cout << _strandMesh->getPositions() << endl;
    cout << " Current twistAngles: " << endl;
    cout << _strandMesh->twistAngles() << endl;

    // create the integrator
    //STRAND::ISOTROPIC_BENDING* isotropicBendingEnergy = new STRAND::ISOTROPIC_THETA(1.0);
    //STRAND::ISOTROPIC_BENDING* isotropicBendingEnergy = new STRAND::QUADRATIC_BENDING(10000.0);
    STRAND::ISOTROPIC_BENDING* isotropicBendingEnergy = new STRAND::QUADRATIC_UNIT_BENDING(10000.0);
    _strandSolver = new STRAND::TIMESTEPPER_DVT(*_strandMesh, *stretchingEnergy, *isotropicBendingEnergy);

    // constrain the first vertex
    VECTOR3 center0(-1.0,-0.2, 0.0);
    _kinematicShapes.push_back(new CUBE(center0, 2.1));
    _strandSolver->attachKinematicSurfaceConstraints(_kinematicShapes[0]);

    _eye    = VECTOR3(3.79907, 1.29043, 4.21273);
    _lookAt = VECTOR3(3.37669, 1.07395, 3.33254);
    _up     = VECTOR3(-0.19182, 0.970417, -0.146615);

    _pauseFrame = 172;

    return true;
  }

  /*
  // for testing twist
  virtual bool buildScene() override
  {
    // this will determine the MOV and JSON filenames
    _sceneName = "strand";

    // procedurally generate a scene
    std::vector<VECTOR3> restVertices;
    std::vector<VECTOR3> vertices;

    restVertices.push_back(VECTOR3(0,0,0));
    restVertices.push_back(VECTOR3(1,0,0));
    restVertices.push_back(VECTOR3(1,0,1));

    vertices = restVertices;

    const int maxElements = 3;
    VECTOR restThetas(maxElements - 1);
    VECTOR thetas(maxElements - 1);
    restThetas.setZero();
    thetas.setZero();

    // initialize the strand mesh
    _strandMesh = new STRAND_MESH(restVertices, restThetas);
    _strandMesh->setPositions(vertices);

    STRAND::QUADRATIC* stretchingEnergy = new STRAND::QUADRATIC(1.0);
    const REAL currentStretch = _strandMesh->computeStretchingEnergy(*stretchingEnergy);
    cout << " Stretching energy: " << currentStretch << endl;

    STRAND::BERGOU_2010 bendingEnergy(1.0);
    const REAL currentBending = _strandMesh->computeBendingEnergy(bendingEnergy);
    cout << " Bending energy: " << currentBending << endl;

    cout << " Current positions: " << endl;
    cout << _strandMesh->getPositions() << endl;
    cout << " Current thetas: " << endl;
    cout << _strandMesh->thetas() << endl;

    // create the integrator
    _strandSolver = new STRAND::TIMESTEPPER(*_strandMesh, *stretchingEnergy);

    // constrain the first vertex
    VECTOR3 center0(0.0,0.0, 0.0);
    _kinematicShapes.push_back(new CUBE(center0, 0.1));
    _strandSolver->attachKinematicSurfaceConstraints(_kinematicShapes[0]);

    _eye    = VECTOR3(2.63737, 1.1219, 2.54607);
    _lookAt = VECTOR3(1.99534, 0.743962, 1.87901);
    _up     = VECTOR3(-0.270845, 0.925762, -0.263826);

    return true;
  }
  */


  // simulation loop
  virtual void stepSimulation(const bool verbose = true) override {
    _strandSolver->externalForces().setZero();
    _strandSolver->addGravity(_gravity);
    _strandSolver->solve(verbose);

    //_kinematicShapes[1]->translation() += VECTOR3(0.0, 0.01, 0.0);
    //_kinematicShapes[1]->translation() += VECTOR3(0.0, 0.02, 0.0);
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

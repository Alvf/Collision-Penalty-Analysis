#ifndef SIMULATION_SCENE_H
#define SIMULATION_SCENE_H

#ifndef GL_DISABLED
#include "util/DRAW_GL.h"
#endif

#include "Geometry/CAPSULE.h"
#include "Geometry/CUBE.h"
#include "Geometry/SPHERE.h"

// volume classes
#include "Geometry/TET_MESH.h"
#include "Geometry/TET_MESH_FASTER.h"
#include "Hyperelastic/Volume/SNH.h"
#include "Hyperelastic/Volume/STVK.h"
#include "Hyperelastic/Volume/ARAP.h"
#include "Hyperelastic/Volume/LINEAR.h"
#include "Timestepper/Volume/BACKWARD_EULER_VELOCITY.h"
#include "Timestepper/Volume/BACKWARD_EULER_POSITION.h"
#include "Timestepper/Volume/NEWMARK.h"
#include "Timestepper/Volume/QUASISTATIC.h"

// shell classes
#include "Geometry/TRIANGLE_MESH.h"
#include "Geometry/TRIANGLE_MESH_FASTER.h"
#include "Timestepper/Shell/TIMESTEPPER.h"
#include "Hyperelastic/Shell/STVK.h"
#include "Hyperelastic/Shell/ARAP.h"
#include "Hyperelastic/Shell/BW_STRETCH.h"
#include "Hyperelastic/Shell/BW_SHEAR.h"
#include "Hyperelastic/Shell/BARAFF_WITKIN.h"
#include "Hyperelastic/Shell/BENDING_SPRING.h"
#include "Hyperelastic/Shell/QUADRATIC_F_BENDING.h"
#include "Hyperelastic/Shell/DIHEDRAL.h"

// strand classes
#include "Geometry/TET_WISP_MESH.h"
#include "Timestepper/Strand/TIMESTEPPER.h"

// collision classes
#include "Collision/C_PLANES_EE.h"
#include "Collision/C_SPECIAL_EE.h"
#include "Collision/UNSIGNED_SPRING_1D.h"
#include "Collision/UNSIGNED_LEN_PLANES.h"
#include "util/TIMER.h"

namespace HOBAK {

class SIMULATION_SCENE {
public:

  // initialize the scene
  SIMULATION_SCENE() {
    _pauseFrame = -2;
    _exitOnPause = true;
    _arrowCounter = -1;
    _leftArrow = false;
    _rightArrow = false;
    _drawFeature = false;
    _frameNumber = 0;
    _normalizedVertices = false;
    _sceneName = std::string("default");
    _initialA = MATRIX3::Identity();
    _initialTranslation = VECTOR3::Zero();

    // volume variables
    _volumeSolver = NULL;
    _tetMesh = NULL;
    _hyperelastic = NULL;
    _gravity.setZero();

    // shell variables
    _shellSolver = NULL;
    _triangleMesh = NULL;
    _strechingEnergy = NULL;
    _bendingEnergy = NULL;

    // strand variables
    _strandSolver = NULL;
    _strandMesh = NULL;
    _E = -1;
    _G = -1;
    _density = -1;
    _baseRadius = -1;
    _tipRadius = -1;

    _movieInterval = -1;
    _autoplay = true;
  
    _normalSpring = NULL;
    _lengthFunc = NULL;
    _cFunc = NULL;
    _vfGeneral = NULL;
  
    _lengthFuncEE = NULL;
    _normalSpringEE = NULL;
    _cFuncEE = NULL;
    _eeGeneral = NULL;
  };

  virtual ~SIMULATION_SCENE()
  {
    delete _tetMesh;
    delete _volumeSolver;
    delete _hyperelastic;

    delete _triangleMesh;
    delete _shellSolver;
    delete _strechingEnergy;
    delete _bendingEnergy;

    for (unsigned int x = 0; x < _kinematicShapes.size(); x++)
      delete _kinematicShapes[x];

    delete _normalSpring;
    delete _lengthFunc;
    delete _cFunc;
    delete _vfGeneral;
  
    delete _lengthFuncEE;
    delete _normalSpringEE;
    delete _cFuncEE;
    delete _eeGeneral;
  };

  // TODO: Build the actual scene. You have to implement this!
  virtual bool buildScene() = 0;

  // TODO: Describe the scene build built. You have to do this!
  virtual void printSceneDescription() = 0;

  const VECTOR3& eye() const         { return _eye; };
  const VECTOR3& lookAt() const      { return _lookAt; };
  const VECTOR3& up() const          { return _up; };
  const VECTOR3& worldCenter() const { return _worldCenter; };
  const VECTOR3& gravity() const  { return _gravity; };
  const int& pauseFrame() const   { return _pauseFrame; };
  const int& frameNumber() const  { return _frameNumber; };
  const bool& drawFeature() const { return _drawFeature; };
  const int& arrowCounter() const { return _arrowCounter; };
  const string& tetMeshFilename() const { return _tetMeshFilename; };
  const string& triangleMeshFilename() const { return _triangleMeshFilename; };
  const string& strandMeshFilename() const { return _strandMeshFilename; };
  const vector<KINEMATIC_SHAPE*>& kinematicShapes() const { return _kinematicShapes; };
  const bool& normalizedVertices() const { return _normalizedVertices; };
  const bool& exitOnPause() const { return _exitOnPause; };

  const MATRIX3& initialA() const           { return _initialA; };
  const VECTOR3& initialTranslation() const { return _initialTranslation; };

  REAL& E()              { return _E; };
  const REAL& E() const  { return _E; };
  //REAL& nu()             { return _nu; };
  //const REAL& nu() const { return _nu; };
  REAL& G()             { return _G; };
  const REAL& G() const { return _G; };
  REAL& strandDensity()             { return _density; };
  const REAL& strandDensity() const { return _density; };
  REAL& strandBaseRadius()             { return _baseRadius; };
  const REAL& strandBaseRadius() const { return _baseRadius; };
  REAL& strandTipRadius()             { return _tipRadius; };
  const REAL& strandTipRadius() const { return _tipRadius; };

  VECTOR3& eye()         { return _eye; };
  VECTOR3& lookAt()      { return _lookAt; };
  VECTOR3& up()          { return _up; };
  VECTOR3& worldCenter() { return _worldCenter; };
  VECTOR3& gravity()     { return _gravity; };
  int& pauseFrame()      { return _pauseFrame; };
  int& frameNumber()     { return _frameNumber; };
  bool& drawFeature()    { return _drawFeature; };
  int& arrowCounter()    { return _arrowCounter; };
  bool& leftArrow()      { return _leftArrow; };
  bool& rightArrow()     { return _rightArrow; };
  string& tetMeshFilename()  { return _tetMeshFilename; };
  string& triangleMeshFilename() { return _triangleMeshFilename; };
  string& strandMeshFilename() { return _strandMeshFilename; };
  bool& normalizedVertices() { return _normalizedVertices; };
  vector<KINEMATIC_SHAPE*>& kinematicShapes() { return _kinematicShapes; };

  const string sceneName() const { return _sceneName; };
  const string movieName() const { return _sceneName + std::string(".mov"); };
  const string jsonName() const  { return _sceneName + std::string(".json"); };
  const bool autoplay() const { return _autoplay; };

  const VOLUME::TIMESTEPPER* solver() const { return _volumeSolver; };
  VOLUME::TIMESTEPPER* solver()             { return _volumeSolver; };
  const SHELL::TIMESTEPPER* shellSolver() const { return _shellSolver; };
  SHELL::TIMESTEPPER* shellSolver()             { return _shellSolver; };
  const STRAND::TIMESTEPPER* strandSolver() const { return _strandSolver; };
  STRAND::TIMESTEPPER* strandSolver()             { return _strandSolver; };
  const TET_MESH* tetMesh() const { return _tetMesh; };
  TET_MESH* tetMesh()             { return _tetMesh; };
  const TRIANGLE_MESH* triangleMesh() const { return _triangleMesh; };
  TRIANGLE_MESH* triangleMesh()             { return _triangleMesh; };
  const STRAND_MESH* strandMesh() const { return _strandMesh; };
  STRAND_MESH* strandMesh()             { return _strandMesh; };

  // simulation loop
  virtual void stepSimulation(const bool verbose = true) {
    // _volumeSolver->externalForces().setZero();
    // _volumeSolver->addGravity(_gravity);
    // _volumeSolver->solve(verbose);

    // if (_writeToFile) {
    //   if (!writeFrameToFile()) {
    //     // Do something here that isn't too noisy.
    //   }
    // }

    // _frameNumber++;
  };

  // should we write a movie at this frame?
  const bool writeMovie() const { 
    if (_movieInterval == -1) return false;

    if ((_frameNumber % _movieInterval) == 0)
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " Frame number: " << _frameNumber << " movie interval: " << _movieInterval << endl;
      cout << " mod: " << _frameNumber % _movieInterval << endl;
      return true;
    }
    return false;
  };

  virtual bool writeFrameToFile() { return true; };

// drawScene MUST have this guard, so it can be deactivated for regression testing
#ifndef GL_DISABLED
  // what to draw?
  virtual void drawScene()
  {
    glEnable(GL_DEPTH_TEST);
    drawSurfaceTriangles(*_tetMesh, true);
    for (unsigned int x = 0; x < _kinematicShapes.size(); x++)
      drawKinematicShape(*_kinematicShapes[x]);
  };
#endif

  void setCollisions()
  {
    setCollisions(1000.0, 0.02);
  }

  void setCollisions(const REAL& mu, const REAL& eps)
  {
    _normalSpring = new ENERGY_1D(mu, eps);
    _lengthFunc = new SIGNED_LEN_PLANES();
    _cFunc = new C_PLANES();

    _normalSpringEE = new ENERGY_1D(mu, eps);
    _lengthFuncEE = new SIGNED_LEN_PLANES();
    _cFuncEE = new C_PLANES_EE();

    _vfGeneral = new ENERGY_12D(_normalSpring, _cFunc, _lengthFunc);
    _eeGeneral = new ENERGY_12D(_normalSpringEE, _cFuncEE, _lengthFuncEE);
  }
  
  void setUnsignedCollisions(const REAL& mu, const REAL& eps)
  {
    _normalSpring = new UNSIGNED_SPRING_1D(mu, eps);
    _lengthFunc = new UNSIGNED_LEN_PLANES();
    _cFunc = new C_PLANES();

    _normalSpringEE = new UNSIGNED_SPRING_1D(mu, eps);
    _lengthFuncEE = new UNSIGNED_LEN_PLANES();
    _cFuncEE = new C_PLANES_EE();

    _vfGeneral = new ENERGY_12D(_normalSpring, _cFunc, _lengthFunc);
    _eeGeneral = new ENERGY_12D(_normalSpringEE, _cFuncEE, _lengthFuncEE);
  }

protected:
  // set the positions to previous timestep, in case the user wants to
  // look at that instead of the current step
  void setToPreviousTimestep()
  {
    const VECTOR& old = _volumeSolver->positionOld();
    _tetMesh->setDisplacement(old);
  }

  // restore positions from previous timestep, in case the user just drew
  // the previous timestep, but now we want the state to be consistent
  // when drawing the next frame
  void restoreToCurrentTimestep()
  {
    const VECTOR& current = _volumeSolver->position();
    _tetMesh->setDisplacement(current);
  }

  // scene geometry
  TET_MESH* _tetMesh;
  TRIANGLE_MESH* _triangleMesh;
  TET_WISP_MESH* _strandMesh;
  vector<KINEMATIC_SHAPE*> _kinematicShapes;

  // solver and materials
  VOLUME::TIMESTEPPER* _volumeSolver;
  VOLUME::HYPERELASTIC* _hyperelastic;
  SHELL::TIMESTEPPER* _shellSolver;
  SHELL::STRETCHING* _strechingEnergy;
  SHELL::BENDING* _bendingEnergy;
  STRAND::TIMESTEPPER* _strandSolver;
  
  // simulation parameters
  VECTOR3 _gravity;

  // initial rotation-scale and translation of tet mesh
  MATRIX3 _initialA;
  VECTOR3 _initialTranslation;

  // drawing parameters
  int _pauseFrame;

  // counter that can be incremented and decremented by the arrow keys
  int _arrowCounter;

  // bools that can be toggled by the arrow keys
  bool _leftArrow;
  bool _rightArrow;

  // flag for whether or not to draw some user-specific feature
  bool _drawFeature;

  // what frame are we on?
  int _frameNumber;

  // should we save this frame?
  bool _writeToFile = false;

  // geometry filenames
  string _tetMeshFilename;
  string _triangleMeshFilename;
  string _strandMeshFilename;

  // did we normalize the vertices when we read them in?
  bool _normalizedVertices;

  // camera parameters
  VECTOR3 _eye;
  VECTOR3 _lookAt;
  VECTOR3 _up;
  VECTOR3 _worldCenter;

  // scene name, used to determine the JSON and MOV filenames
  std::string _sceneName;

  // The path to save the files to
  std::string _filePath;

  // should the scene automatically start simulating?
  bool _autoplay;

  // should we exit when we hit the pause frame?
  bool _exitOnPause;
 
  // what interval should we write out a movie at?
  int _movieInterval;
 
  // Young's modulus and shear modulus for STRANDS
  REAL _E;
  //REAL _nu;     // let's prefer shear modulus to Poisson's ratio
  REAL _G;

  // strand thickness and mass parameters
  REAL _density;
  REAL _baseRadius;
  REAL _tipRadius;

  // vertex-face collision support
  ENERGY_1D* _normalSpring;
  SIGNED_LEN_PLANES* _lengthFunc;
  C_PLANES* _cFunc;
  ENERGY_12D* _vfGeneral;
    
  // edge-edge collision support
  SIGNED_LEN_PLANES* _lengthFuncEE;
  ENERGY_1D* _normalSpringEE;
  C_PLANES* _cFuncEE;
  ENERGY_12D* _eeGeneral;
};

}

#endif

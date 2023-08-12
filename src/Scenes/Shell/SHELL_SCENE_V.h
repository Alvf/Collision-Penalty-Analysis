#ifndef SHELL_SCENE_V_H
#define SHELL_SCENE_V_H

#include "Scenes/SIMULATION_SCENE.h"
#include <filesystem>

namespace HOBAK {

#define USING_FILESYSTEM 0

#if USING_FILESYSTEM
#if __APPLE__
namespace fs = std::__fs::filesystem;
#else
namespace fs = std::filesystem;
#endif
#endif

class SHELL_SCENE_V : public SIMULATION_SCENE {
public:
  virtual void printSceneDescription() override
  {
    cout << "=====================================================================" << endl;
    cout << " Shell V bending scene" << endl;
    cout << "=====================================================================" << endl;
  }

  virtual bool buildScene() override
  {

    // V shape parameters
    n0 = 31; n1 = 21; vTheta = M_PI / 4.0;

    //_triangleMeshFilename = string("../data/objs/curtain_10x10.obj");
    // _triangleMeshFilename = string("../data/objs/curtain_25x25.obj");
    // //_triangleMeshFilename = string("../data/objs/curtain_50x50.obj");
    // vector<VECTOR3> vertices;
    // vector<VECTOR3I> triangles;
    // bool success = TRIANGLE_MESH::readObjFile(_triangleMeshFilename, vertices, triangles);
    // if (!success)
    // {
    //   cout << " Failed to open file " << _triangleMeshFilename.c_str() << endl;
    //   return false;
    // }
    // _triangleMesh = new TRIANGLE_MESH_FASTER(vertices, triangles);

    buildVCloth1();
    //_triangleMesh->setCollisionEps(spacing / 8.0);
    //_gravity = VECTOR3(0, 0, -1.0);
    _gravity = VECTOR3(0, -1.0, 0.0);

    const REAL stretchingStiffness = 8.0;
    const REAL bendingStiffness = 1e5;
    //_strechingEnergy = new SHELL::STVK(10.0, 1.0);
    _strechingEnergy = new SHELL::ARAP(stretchingStiffness, 0.0);
    //_strechingEnergy = new SHELL::BW_STRETCH(1.0, 0.0);
    //_strechingEnergy = new SHELL::BW_SHEAR(1.0, 0.0);
    //_strechingEnergy = new SHELL::BARAFF_WITKIN(1.0, 0.0);
    //_bendingEnergy = new SHELL::BENDING_SPRING(20);
    //_bendingEnergy = new SHELL::BENDING_SPRING(10);
    //_bendingEnergy = new SHELL::DIHEDRAL(0.1);
    // _bendingEnergy = new SHELL::DIHEDRAL(1);
    _bendingEnergy = new SHELL::QUADRATIC_F_BENDING(bendingStiffness);
    // _bendingEnergy = new SHELL::DIHEDRAL(50);
    //_bendingEnergy = new SHELL::DIHEDRAL(100);

    // this will determine the MOV and JSON filenames
    // string method = "ours";
    // string method = "bruteCholesky";
    string method = "rank4Cholesky";
    // string method = "GNCholesky";
    // string method = "unFIlCholesky";
    char buffer[100];
    sprintf(buffer, "_S%.1f_B%.1f", stretchingStiffness, bendingStiffness);
    _sceneName = "shell_bending_V_" + method + buffer;
    // _writeToFile = true;

#if USING_FILESYSTEM
    if(_writeToFile){
      if(!fs::exists("../data/renders"))
        fs::create_directory("../data/renders");
      string outDir = string("../data/renders/") + _sceneName + string("/");
      if(!fs::exists(outDir))
        fs::create_directory(outDir);
      if(!fs::exists(outDir + string("shell")))
        fs::create_directory(outDir + string("shell"));

      _shellFilePath = outDir + string("shell/");
    }
#endif

    setCollisions(2000.0, spacing / 8.0);
    _shellSolver = new SHELL::TIMESTEPPER(*_triangleMesh, *_strechingEnergy, *_bendingEnergy, *_vfGeneral, *_eeGeneral);

    /*
    //Build the collision penalties
    vector<REAL> springArgs;
    springArgs.push_back(2000.0);
    springArgs.push_back(spacing / 8.0);
    ENERGY_1D* normalSpring = new ENERGY_1D(springArgs);
    SIGNED_LEN_PLANES* lengthFunc = new SIGNED_LEN_PLANES();
    C_PLANES* cFunc = new C_PLANES();

    SIGNED_LEN_PLANES* lengthFuncEE = new SIGNED_LEN_PLANES();
    ENERGY_1D* normalSpringEE = new ENERGY_1D(springArgs);
    C_PLANES_EE* cFuncEE = new C_PLANES_EE();

    ENERGY_12D* vfGeneral = new ENERGY_12D(normalSpring, cFunc, lengthFunc);
    ENERGY_12D* eeGeneral = new ENERGY_12D(normalSpringEE, cFuncEE, lengthFuncEE);

    // build the time integrator
    _shellSolver = new SHELL::TIMESTEPPER(*_triangleMesh, *_strechingEnergy, *_bendingEnergy, *vfGeneral, *eeGeneral);
    */
    // _shellSolver -> edgeEdgeSelfCollisionsOn() = false;
    // _shellSolver -> vertexFaceSelfCollisionsOn() = false;
    // _shellSolver -> setRayeligh(0.001, 0.0);
    // _shellSolver->setDt(1.0 / 60.0);

    // cube on top
    // VECTOR3 center0(0.0, 0.0, 0.95);
    //_kinematicShapes.push_back(new CUBE(center0, 1.0));
    //_shellSolver->attachKinematicSurfaceConstraints(_kinematicShapes[0]);
    //_shellSolver->addKinematicCollisionObject(_kinematicShapes[0]);
    
    const REAL cubeScale = 8;
    VECTOR3 center0(-0.5 * cubeScale, -0.4 * cubeScale, 0.4 * cubeScale);
    // VECTOR3 center0(-0.5 * cubeScale, 0.4 * cubeScale, 0.4 * cubeScale);
    //VECTOR3 center1(0.0, -1.0, -0.01);
    _kinematicShapes.push_back(new CUBE(center0, cubeScale));
    // _kinematicShapes.push_back(new SPHERE(center1, 0.25));
    _shellSolver->attachKinematicSurfaceConstraints(_kinematicShapes[0]);
    _shellSolver->addKinematicCollisionObject(_kinematicShapes[0]);
    
    // _shellSolver->addKinematicCollisionObject(_kinematicShapes.back());

    //_eye    = VECTOR3(1.24999, 2.07096, 0.502227);
    //_lookAt = VECTOR3(0.777846, 1.20965, 0.314523);
    //_up     = VECTOR3(-0.0859995, -0.166908, 0.982213);

    // build
    // _eye    = VECTOR3(0.5, 0.8, 6.0);
    // _lookAt = VECTOR3(0.3, 0.0, 0.0);
    // _up     = VECTOR3(0.0, 1.0, 0.0);

    // build 1
    _eye    = VECTOR3(5.0, 3.0, 1.0);
    _lookAt = VECTOR3(0.3, 0.0, 0.0);
    _up     = VECTOR3(0.0, 1.0, 0.0);

    _worldCenter = _triangleMesh->getRestTranslation();
    //_worldCenter = VECTOR3(0, 0, 0);
    _pauseFrame = 400;

    // try scaling see if F follows
    //for (unsigned int x = 0; x < _triangleMesh->vertices().size(); x++)
    //  _triangleMesh->vertex(x)[2] *= 2.0;

    return true;
  }

  // simulation loop
  virtual void stepSimulation(const bool verbose = true) override {
    if (_writeToFile && _frameNumber == 0)
    {
      char buffer[256];
      sprintf(buffer, "%04i", _frameNumber);
      string filenameShell = _shellFilePath + _sceneName + string("_frame_") + string(buffer) + string(".obj");
      // cloth
      cout << " Writing file " << filenameShell.c_str() << " ... " << flush;
      writeOBJFile(filenameShell.c_str(), _triangleMesh->vertices(), _triangleMesh->triangles());
      cout << "done." << endl;

      // TIMER gzipper("Gzipping");
      // string gzip = string("gzip -f -9 ");
      // cout << " Gzipping cloth mesh ..." << flush;
      // system((gzip + filenameShell).c_str());
      // cout << "done." << endl;
      // gzipper.stop();
    }

    _shellSolver->externalForces().setZero();
    _shellSolver->addGravity(_gravity);
    _shellSolver->solve(verbose);

    if (_writeToFile)
    {
      char buffer[256];
      sprintf(buffer, "%04i", _frameNumber + 1);
      string filenameShell = _shellFilePath + _sceneName + string("_frame_") + string(buffer) + string(".obj");
      // cloth
      cout << " Writing file " << filenameShell.c_str() << " ... " << flush;
      writeOBJFile(filenameShell.c_str(), _triangleMesh->vertices(), _triangleMesh->triangles());
      cout << "done." << endl;

      // TIMER gzipper("Gzipping");
      // string gzip = string("gzip -f -9 ");
      // cout << " Gzipping cloth mesh ..." << flush;
      // system((gzip + filenameShell).c_str());
      // cout << "done." << endl;
      // gzipper.stop();
    }

    //_kinematicShapes[1]->translation() += VECTOR3(0.0, 0.01, 0.0);
    //_kinematicShapes[1]->translation() += VECTOR3(0.0, 0.02, 0.0);
    _frameNumber++;
  };

// drawScene MUST have this guard, so it can be deactivated for regression testing
#ifndef GL_DISABLED
  // what to draw?
  virtual void drawScene()
  {
    drawAxes();
    drawTriangleMesh(*_triangleMesh, true);

    glEnable(GL_DEPTH_TEST);
    for (unsigned int x = 0; x < _kinematicShapes.size(); x++)
      drawKinematicShape(*_kinematicShapes[x]);

    drawKinematicConstraints(_triangleMesh, _shellSolver);
    //drawPlaneConstraints(_triangleMesh, _shellSolver);

    drawVertex(*_triangleMesh, _arrowCounter);
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

  void buildVCloth(){
    // the origin is at one corner.
    vector<VECTOR3> vertices;
    vector<VECTOR3I> triangles;
    bool alter = true, alterRow = false;

    spacing = 1.0/(n1 - 1.0);
    // first part
    const REAL zStep0 = spacing;
    const REAL xStep0 = zStep0;
    VECTOR3 vert(0.0, 0.0, 0.0);
    int vId = 0; // next to add
    //first row
    for(int i = 0; i < n0; i++) {
      vertices.push_back(vert);
      vert(2) += zStep0;
    }
    for(int j = 0; j< (n1 - 1)/2; j++){
      vert(0) += xStep0; vert(2) = 0.0;
      vertices.push_back(vert);
      // add row
      alter = alterRow;
      for(int i = 0; i < n0 - 1; i++) {
        vert(2) += zStep0;
        vId = (int)vertices.size();
        vertices.push_back(vert);
        if(alter){
          triangles.push_back(VECTOR3I(vId, vId - 1, vId - n0));
          triangles.push_back(VECTOR3I(vId - 1, vId - 1 - n0, vId - n0));
        }
        else{
          triangles.push_back(VECTOR3I(vId, vId - 1, vId - n0 - 1));
          triangles.push_back(VECTOR3I(vId, vId - 1 - n0, vId - n0));
        }
        alter = !alter;
      }
      alterRow = !alterRow;
    }

    // bended part
    const REAL zStep1 = zStep0;
    const REAL xStep1 = -xStep0 * cos(vTheta);
    const REAL yStep1 = xStep0 * sin(vTheta);
    for(int j = 0; j< (n1 - 1)/2; j++){
      vert(0) += xStep1; vert(2) = 0.0; vert(1) += yStep1;
      vertices.push_back(vert);
      // add row
      alter = alterRow;
      for(int i = 0; i < n0 - 1; i++) {
        vert(2) += zStep1;
        vId = (int)vertices.size();
        vertices.push_back(vert);
        if(alter){
          triangles.push_back(VECTOR3I(vId, vId - 1, vId - n0));
          triangles.push_back(VECTOR3I(vId - 1, vId - 1 - n0, vId - n0));
        }
        else{
          triangles.push_back(VECTOR3I(vId, vId - 1, vId - n0 - 1));
          triangles.push_back(VECTOR3I(vId, vId - 1 - n0, vId - n0));
        }
        alter = !alter;
      }
      alterRow = !alterRow;
    }
    // cout<<"vertices: "<<endl;
    // for(int i = 0; i < vertices.size(); i++)
    //   cout<<vertices[i]<<endl;
    _triangleMesh = new TRIANGLE_MESH_FASTER(vertices, triangles);

  }

  void buildVCloth1(){
    // the origin is at one corner.
    vector<VECTOR3> vertices;
    vector<VECTOR3I> triangles;
    bool alter = true, alterRow = false;

    spacing = 1.0/(n1 - 1.0);
    // first part
    const REAL xStep0 = spacing;
    const REAL zStep0 = xStep0 * sin(0.5 * vTheta);
    const REAL yStep0 = - xStep0 * cos(0.5 * vTheta);
    VECTOR3 vert(0.0, 0.0, 0.0);
    int vId = 0; // next to add
    //first row
    for(int i = 0; i < n0; i++) {
      vertices.push_back(vert);
      vert(0) += xStep0;
    }
    for(int j = 0; j< (n1 - 1)/2; j++){
      vert(2) += zStep0; vert(1) += yStep0; vert(0) = 0.0;
      vertices.push_back(vert);
      // add row
      alter = alterRow;
      for(int i = 0; i < n0 - 1; i++) {
        vert(0) += xStep0;
        vId = (int)vertices.size();
        vertices.push_back(vert);
        if(alter){
          triangles.push_back(VECTOR3I(vId, vId - 1, vId - n0));
          triangles.push_back(VECTOR3I(vId - 1, vId - 1 - n0, vId - n0));
        }
        else{
          triangles.push_back(VECTOR3I(vId, vId - 1, vId - n0 - 1));
          triangles.push_back(VECTOR3I(vId, vId - 1 - n0, vId - n0));
        }
        alter = !alter;
      }
      alterRow = !alterRow;
    }

    // bended part
    const REAL xStep1 = xStep0;
    const REAL zStep1 = zStep0;
    const REAL yStep1 = -yStep0;
    for(int j = 0; j< (n1 - 1)/2; j++){
      vert(2) += zStep1; vert(0) = 0.0; vert(1) += yStep1;
      vertices.push_back(vert);
      // add row
      alter = alterRow;
      for(int i = 0; i < n0 - 1; i++) {
        vert(0) += xStep1;
        vId = (int)vertices.size();
        vertices.push_back(vert);
        if(alter){
          triangles.push_back(VECTOR3I(vId, vId - 1, vId - n0));
          triangles.push_back(VECTOR3I(vId - 1, vId - 1 - n0, vId - n0));
        }
        else{
          triangles.push_back(VECTOR3I(vId, vId - 1, vId - n0 - 1));
          triangles.push_back(VECTOR3I(vId, vId - 1 - n0, vId - n0));
        }
        alter = !alter;
      }
      alterRow = !alterRow;
    }
    // cout<<"vertices: "<<endl;
    // for(int i = 0; i < vertices.size(); i++)
    //   cout<<vertices[i]<<endl;
    _triangleMesh = new TRIANGLE_MESH_FASTER(vertices, triangles);

  }
#endif

protected:
  bool _writeToFile;
  string _shellFilePath;

  // V shape parameters
  int n0, n1; // number of vertices along the edge. n1 should be odd
  REAL vTheta; // angle of V
  REAL spacing;
};

}

#endif

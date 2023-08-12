/*
This file is part of HOBAK.

HOBAK is free software: you can redistribute it and/or modify it under the terms of
the GNU General Public License as published by the Free Software Foundation, either
version 3 of the License, or (at your option) any later version.

HOBAK is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with HOBAK.
If not, see <https://www.gnu.org/licenses/>.
*/
#ifndef FILE_IO_H
#define FILE_IO_H

#include "SETTINGS.h"
#include "Timestepper/Volume/TIMESTEPPER.h"
#include "Scenes/SIMULATION_SCENE.h"
#include "Scenes/Volume/JSON_VOLUME_SCENE.h"
#include "Scenes/Shell/JSON_SHELL_SCENE.h"
#include "Scenes/Strand/JSON_STRAND_SCENE.h"
#include <document.h>

namespace HOBAK {
  void initializeVolumeSceneJSON(const SIMULATION_SCENE& scene);
  void initializeShellSceneJSON(const SIMULATION_SCENE& scene);
  void initializeStrandSceneJSON(const SIMULATION_SCENE& scene);
  
  void recordVolumeFrameToJSON(const SIMULATION_SCENE& scene);
  void recordShellFrameToJSON(const SIMULATION_SCENE& scene);
  void recordStrandFrameToJSON(const SIMULATION_SCENE& scene);

  void writeVolumeSceneJSON(const char* filename);
  void writeShellSceneJSON(const char* filename);
  void writeStrandSceneJSON(const char* filename);
 
  // are we looking at a volume, shell, or strand file?
  string readSceneTypeJSON(const char* filename);

  bool readVolumeSceneJSON(const char* filename, JSON_VOLUME_SCENE& scene);
  bool readShellSceneJSON(const char* filename, JSON_SHELL_SCENE& scene);
  bool readStrandSceneJSON(const char* filename, JSON_STRAND_SCENE& scene);

  bool writeOBJFile(const char *filename, const vector<VECTOR3> &vertices,
                     const vector<VECTOR3I> &indices);
};

#endif

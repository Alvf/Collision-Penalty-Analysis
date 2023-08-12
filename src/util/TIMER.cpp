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
#include "TIMER.h"
#include <cassert>
#include <stdio.h>

using namespace std;

timeval TIMER::_tick;
timeval TIMER::_tock;
std::map<std::string, double> TIMER::_timings;
std::stack<std::string> TIMER::_callStack;

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
TIMER::TIMER(string blockName) 
{
  // introduces more than just a pragme dependency on OpenMP, so backing
  // off for now
  //const int inOMP = omp_in_parallel();
  //if (inOMP) return;

  // look at the back of the call stack,
  // if there's something there, then store its timing
  //
  // else, it's the first call, so set the global start
  if (_callStack.size() > 0)
  {
    string function = _callStack.top();
    gettimeofday(&_tock, 0);

    _timings[function] += timing();
  }
  
  _callStack.push(blockName);
  gettimeofday(&_tick, 0);

  _stopped = false;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
TIMER::~TIMER()
{
  // introduces more than just a pragme dependency on OpenMP, so backing
  // off for now
  //const int inOMP = omp_in_parallel();
  //if (inOMP) return;

  // don't stop it twice
  if (!_stopped)
  {
    stop();
  }
}

///////////////////////////////////////////////////////////////////////
// stop the timer manually
///////////////////////////////////////////////////////////////////////
void TIMER::stop()
{
  // introduces more than just a pragme dependency on OpenMP, so backing
  // off for now
  //const int inOMP = omp_in_parallel();
  //if (inOMP) return;

  // don't stop it twice
  if (_stopped) return;

  assert(_callStack.size() > 0);
  
  string function = _callStack.top();
  _callStack.pop();
  gettimeofday(&_tock, 0);

  _timings[function] += timing();
  gettimeofday(&_tick, 0);

  _stopped = true;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TIMER::printTimings()
{
  timeval now;
  gettimeofday(&now, 0);

  double currentTimer = timing(_tick, now);
  string currentName;
  if (_callStack.size() > 0) 
    currentName = _callStack.top();

  // create an inverse map so that it will sort by time
  map<double, string> inverseMap;
  map<string, double>::iterator forwardIter;
  double totalTime = 0.0;
  for (forwardIter = _timings.begin(); forwardIter != _timings.end(); forwardIter++)
  {
    string name = forwardIter->first;
    double time = forwardIter->second;

    // if we are dinside this function, add the timing
    if (name.compare(currentName) == 0)
      time += currentTimer;

    inverseMap[time] = name;
    totalTime += time;
  }

  // print the map out backwards since it sorts from least to greatest
  cout << "=====================================================================================" << endl;
  cout << " TIMING BREAKDOWN: " << endl;
  cout << "=====================================================================================" << endl;
  map<double,string>::reverse_iterator backwardIter;
  for (backwardIter = inverseMap.rbegin(); backwardIter != inverseMap.rend(); backwardIter++)
  {
    string name = (*backwardIter).second;
    double time = (*backwardIter).first;

    name = name + string("                                                            ");
    name = name.substr(0,55);

    char buffer[512];
    sprintf(buffer, "%lf", time / totalTime * 100.0);
    string timing(buffer);
    timing = timing + string("%                                                            ");
    timing = timing.substr(0,14);

    //cout << "[" << time / totalTime * 100.0 << "%\t]: "
    cout << "[" << timing.c_str() << "]: "
         << name.c_str() << " " << hours(time) << ":" << minutes(time) << ":" << seconds(time) << "s" << endl;
  }
  cout << "=====================================================================================" << endl;
  cout << " total time: " << totalTime << endl;
  cout << "===========================================================================" << endl;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TIMER::printTimingsPerFrame(const int frames)
{
  timeval now;
  gettimeofday(&now, 0);

  double currentTimer = timing(_tick, now);
  string currentName;
  if (_callStack.size() > 0) 
    currentName = _callStack.top();

  // create an inverse map so that it will sort by time
  map<double, string> inverseMap;
  map<string, double>::iterator forwardIter;
  double totalTime = 0.0;
  for (forwardIter = _timings.begin(); forwardIter != _timings.end(); forwardIter++)
  {
    string name = forwardIter->first;
    double time = forwardIter->second;

    // if we are dinside this function, add the timing
    if (name.compare(currentName) == 0)
      time += currentTimer;

    inverseMap[time] = name;
    totalTime += time;
  }

  // print the map out backwards since it sorts from least to greatest
  cout << "==============================================================================================" << endl;
  cout << " TIMING BREAKDOWN, FRAME " << frames << ": " << endl;
  cout << "==============================================================================================" << endl;
  map<double,string>::reverse_iterator backwardIter;
  char buffer[256];
  string timeString;
  string hoursString;
  string minutesString;
  string secondsString;
  for (backwardIter = inverseMap.rbegin(); backwardIter != inverseMap.rend(); backwardIter++)
  {
    string name = (*backwardIter).second;
    double time = (*backwardIter).first;

    string padding("                                                  ");

    name = name + padding;
    name = name.substr(0,55);

    sprintf(buffer, "%f", time / totalTime * 100.0);
    timeString = string(buffer);
    timeString = timeString + padding;
    timeString = timeString.substr(0,10);

    sprintf(buffer, "%02i", hours(time));
    hoursString = string(buffer);
    
    sprintf(buffer, "%02i", minutes(time));
    minutesString = string(buffer);

    sprintf(buffer, "%02i", seconds(time));
    secondsString = string(buffer);

    cout << "[" << timeString.c_str() << "%]: "
         << name.c_str() << " " << hoursString.c_str() << ":" 
                         << minutesString << ":"
                         << secondsString << " s total ";

    double perFrame = time / frames;
    sprintf(buffer, "%02i", hours(perFrame));
    hoursString = string(buffer);
    
    sprintf(buffer, "%02i", minutes(perFrame));
    minutesString = string(buffer);

    sprintf(buffer, "%02i", seconds(perFrame));
    secondsString = string(buffer);
    cout << " " << hoursString.c_str() << ":" << minutesString.c_str() << ":" << secondsString.c_str() << " s/frame" << endl;
  }
  sprintf(buffer, "%02i", hours(totalTime));
  hoursString = string(buffer);
    
  sprintf(buffer, "%02i", minutes(totalTime));
  minutesString = string(buffer);

  sprintf(buffer, "%02i", seconds(totalTime));
  secondsString = string(buffer);
  cout << "==============================================================================================" << endl;
  cout << " Total time: " << hoursString.c_str() << ":" << minutesString.c_str() << ":" << secondsString.c_str() << " (" << totalTime << "s) " << endl;
  cout << " Time per frame: " << totalTime / frames << endl;
  cout << "==============================================================================================" << endl;
}


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
// MIN_HEAP.h: interface for the MIN_HEAP class.
//
//////////////////////////////////////////////////////////////////////

#ifndef MIN_HEAP_H
#define MIN_HEAP_H

#include <iostream>
#include <vector>
#include <cmath>
#include <map>

#include "SETTINGS.h"

namespace HOBAK {

class HEAP_ENTRY {
public:
  HEAP_ENTRY() {
    distance = 0;
    index = 0;
    heapIndex = 0;
  };

  REAL distance;
  int index;
  int heapIndex;
};

// the actual heap
class MIN_HEAP  
{
public:
  
	MIN_HEAP();
	virtual ~MIN_HEAP();

  // heap ops
  void insert(HEAP_ENTRY& cell);
  void decreaseKey(int toChange, REAL newKey);
  HEAP_ENTRY popMin();
  HEAP_ENTRY heapMin();
  int size() { return _size; };

  // debugging
  void print();
 
  void clear() {
    _heap.clear();
    _heapIndex.clear();
    _size = 0;
  };

  bool empty() { return _size == 0; };
  
private:
  // tree traversal
  int parent(int i) { return i / 2; };
  int left(int i) { return i * 2; };
  int right(int i) { return i * 2 + 1; };
 
  // the heap
  std::vector<HEAP_ENTRY> _heap;
  
  // hash table mapping grid index to heap index
  std::map<int,int> _heapIndex;

  // size of current heap
  int _size;
  
  // enforce heap property 
  void heapify(int index);

  // swap two entries
  void swap(int index1, int index2);
};

} // HOBAK

#endif

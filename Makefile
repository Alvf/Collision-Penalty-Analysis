SHELL := /bin/bash -e

build_mac_omp :
	echo -e "\n==== Building for OSX with OpenMP ====\n";cd projects; cp -f include_top.mac_omp include_top.mk; cd simulateHammock; make depend; make clean; make -j4;cd ../../

build_mac_omp_strand :
	echo -e "\n==== Building for OSX with OpenMP ====\n";cd projects; cp -f include_top.mac_omp include_top.mk; cd simulateStrands; make depend; make clean; make -j4;cd ../../

mac_omp :
	make build_mac_omp; cd bin;./simulateScene

build_mac :
	echo -e "\n==== Building for OSX ====\n";cd projects; cp -f include_top.mac include_top.mk; cd simulateHammock; make depend; make clean; make -j4;cd ../../

build_mac_strand :
	echo -e "\n==== Building for OSX ====\n";cd projects; cp -f include_top.mac include_top.mk; cd simulateStrands; make depend; make clean; make -j3;cd ../../

mac :
	make build_mac;cd bin;./simulateHammock

build_linux :
	echo -e "\n==== Building for Linux ====\n";cd projects; cp -f include_top.linux include_top.mk; cd simulateHammock; make depend; make clean; make -j4;cd ../../

build_linux_strand :
	echo -e "\n==== Building for Linux ====\n";cd projects; cp -f include_top.linux include_top.mk; cd simulateStrands; make depend; make clean; make -j4;cd ../../

build_linux_strand_cli :
	echo -e "\n==== Building for Linux ====\n";cd projects; cp -f include_top.linux include_top.mk; cd simulateStrandsCLI; make depend; make clean; make -j10;cd ../../

linux :
	make build_linux;cd bin;./simulateHammock

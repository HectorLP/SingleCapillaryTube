
CLINKER = g++

LIBRARY_PATH = header/

OBJS = main.o DynamicSingleTube.o

InterfaceLocation: $(OBJS)
	-${CLINKER} -std=c++11 -g -o ~/ComputationCplus/SingleTube/InterfaceLocation $(OBJS)
	${RM} *.o
	
main.o: main.cpp
	-${CLINKER} -std=c++11 -c -I$(LIBRARY_PATH) main.cpp
	
DynamicSingleTube.o: src/DynamicSingleTube.cpp
	-${CLINKER} -std=c++11 -c -I$(LIBRARY_PATH) src/DynamicSingleTube.cpp
	
	
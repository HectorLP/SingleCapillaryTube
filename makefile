
CLINKER = g++

LIBRARY_PATH = include/

OBJS = main.o DynamicSingleTube.o

InterfaceLocation: $(OBJS)
	-${CLINKER} -std=c++11 -g -lboost_system -o ~/ComputationCplus/SingleTube/InterfaceLocation $(OBJS)
	${RM} *.o
	
main.o: main.cpp
	-${CLINKER} -std=c++11 -c -g -I$(LIBRARY_PATH) main.cpp
	
DynamicSingleTube.o: src/DynamicSingleTube.cpp
	-${CLINKER} -std=c++11 -c -g -lboost_system -I$(LIBRARY_PATH) src/DynamicSingleTube.cpp
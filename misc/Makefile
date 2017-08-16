CC = g++

include ../config.mk

TARGETS = 
SOURCES = $(TARGETS:=.cpp)
HEADERS = $(TARGETS:=.hpp)
OBJECTS = $(TARGETS:=.o)

all: $(OBJECTS)

%.o : %.cpp %.hpp
	$(CC) $(INCLUDE) $(CCFLAGS) -c $<

clean: 
	rm -rf *.o 


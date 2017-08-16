CC = g++
#CCFLAGS = -Wall -ansi -pedantic -D DEBUG -04 -funroll-loops
OBJECTS := $(patsubst %.cpp, %.o, $(wildcard *.cpp))
MISC_DIR = ./../misc
INCLUDE := -I$(MISC_DIR)

include ../config.mk

all: $(OBJECTS)

%.o : %.cpp
	$(CC) $(INCLUDE) $(CCFLAGS) -c $<

clean: 
	rm -rf *.o 

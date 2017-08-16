ROOT = .
SRC_DIR = $(ROOT)/src/
MISC_DIR = $(ROOT)/misc/
BIN_DIR = $(ROOT)/bin/
CC = g++
INCLUDE = -I$(MISC_DIR)
COPY = cp

include config.mk

LOCAL_TARGETS = convert_align strand_shift regions tags_in_regions

LOCAL_EXECS = $(LOCAL_TARGETS)
LOCAL_OBJECTS = $(LOCAL_TARGETS:=.o)
OTHER_DEPENDENCIES = data filterstream format kernel peakcall bamtools/BamReader bamtools/BamAlignment bamtools/BGZF bamtools/BamIndex
OTHER_OBJECTS = $(OTHER_DEPENDENCIES:=.o)

EXPANDED_OTHER_OBJECTS = $(addprefix $(MISC_DIR), $(OTHER_OBJECTS))

default: all

.PHONY: default all clean


all: other_dependencies local_objects local_execs place

place:
	mv $(LOCAL_TARGETS) $(BIN_DIR)

local_execs: $(LOCAL_EXECS)

local_objects:
	cd $(SRC_DIR) && $(MAKE) $(LOCAL_OBJECTS)

$(LOCAL_EXECS) : %: $(SRC_DIR)/%.o
	echo $(INCLUDE)
	$(CC) $(INCLUDE) $(CCFLAGS) $(SRC_DIR)$@.o $(EXPANDED_OTHER_OBJECTS) -o $@ $(LIBS)

other_dependencies: other_objects

other_objects:
	cd $(MISC_DIR) && $(MAKE) $(OTHER_OBJECTS)

clean:
	rm -rf $(SRC_DIR)*.o
	cd $(MISC_DIR) && $(MAKE) clean	

all_clean:
	rm -rf $(SRC_DIR)*.o
	cd $(MISC_DIR) && $(MAKE) clean

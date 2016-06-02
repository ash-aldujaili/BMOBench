VERSION = 2014.07.04
# adapted from http://coco.gforge.inria.fr/doku.php?id=mo-gecco2015 Dimo Brockhoff, Thanh-Do Tran, and Nikolaus Hansen
#----------------------------------#
CC     = gcc
CXX    = g++
LINKER = gcc
#----------------------------------#

CFLAGS   = -Wall -O3 -std=c99
CXXFLAGS = -Wall -O3

LDFLAGS =
INCLUDE_FLAGS =
LIBRARY_FLAGS =
MODULES =
LIBRARIES =

#--------------------------------------------------#
SRC_ROOT   =  .
INC_ROOT   =  $(basename $(shell cd $(SRC_ROOT);  find . -type d -print | sed "s|^./||" | sed "/build*/d";  cd ..))
LIB_ROOT   =  .
MODULES   +=  $(basename $(shell cd $(SRC_ROOT);  find . -type d -print | sed "s|^./||" | sed "/build*/d";  cd ..))

LDFLAGS +=
#--------------------------------------------------#

SRC_DIR   =  $(SRC_ROOT)  $(addprefix $(SRC_ROOT)/, $(MODULES))
INC_DIR   =  $(INC_ROOT)  $(SRC_ROOT)
LIB_DIR   =  $(addprefix $(LIB_ROOT)/, $(LIBRARIES))
BUILD_DIR =  build
OBJ_DIR   =  $(BUILD_DIR)  $(addprefix $(BUILD_DIR)/, $(MODULES))

HEADERS_C    = $(foreach sdir, $(SRC_DIR), $(wildcard $(sdir)/*.h))
HEADERS_CPP    = $(foreach sdir, $(SRC_DIR), $(wildcard $(sdir)/*.hpp))
SOURCES_C    = $(foreach sdir, $(SRC_DIR), $(wildcard $(sdir)/*.c))
SOURCES_CPP  = $(foreach sdir, $(SRC_DIR), $(wildcard $(sdir)/*.cpp))
OBJECTS_C    = $(patsubst $(SRC_ROOT)/%.c, $(BUILD_DIR)/%.o, $(SOURCES_C))
OBJECTS_CPP  = $(patsubst $(SRC_ROOT)/%.cpp, $(BUILD_DIR)/%.o, $(SOURCES_CPP))
SOURCES      = $(SOURCES_C)  $(SOURCES_CPP)
OBJECTS      = $(OBJECTS_C)  $(OBJECTS_CPP)
INCLUDES     = $(foreach idir, $(sort $(SRC_DIR) $(INC_DIR)), $(wildcard $(idir)/*.h))

INC_ATTACH  = $(addprefix -I./, $(INC_DIR))
LIB_ATTACH  = $(addprefix -L./, $(LIB_DIR))

EXECUTABLE = $(BUILD_DIR)/runme

vpath %.c $(SRC_DIR)
vpath %.cpp $(SRC_DIR)

define COMPILE_RULES
$1/%.o: %.c $(INCLUDES)
	$(CC) -c $(CFLAGS) $$< -I./$$(<D) $(INC_ATTACH) $(INCLUDE_FLAGS) -o $$@
$1/%.o: %.cpp $(INCLUDES)
	$(CXX) -c $(CXXFLAGS) $$< -I./$$(<D) $(INC_ATTACH) $(INCLUDE_FLAGS) -o $$@
endef

.PHONY: all checkdirs clean backup

all: checkdirs $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(LINKER) $^ $(LIB_ATTACH) $(LIBRARY_FLAGS) -o $@ $(LDFLAGS) -lm
	@cp $@ ./ #copy the binary to the folder of Makefile to facilitate testing

checkdirs: $(OBJ_DIR)

$(OBJ_DIR):
	@mkdir -p $@

clean:
	-rm -rf $(BUILD_DIR)   # $(OBJ_DIR) $(EXECUTABLE)

backup:
	@read -p "Enter a backup name: "  MYFILENAME; \
	tar -czf $$MYFILENAME.tar.gz  src inc Makefile; \
	echo "'$$MYFILENAME.tar.gz' created to backup current version."

$(foreach odir, $(OBJ_DIR), $(eval $(call COMPILE_RULES, $(odir))))

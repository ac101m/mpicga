# Output files
RELEASE_EXEC ?= bin/mpicga
DEBUG_EXEC ?= bin/mpicga-debug
PATTERN_EXEC ?= bin/pattern
TEST_EXEC ?= bin/test

# Directory controls
OBJ_DIR_BASE ?= obj
OBJ_DIR_RELEASE ?= $(OBJ_DIR_BASE)/release
OBJ_DIR_DEBUG ?= $(OBJ_DIR_BASE)/debug
SRC_DIRS ?= src
INC_DIRS ?= include
MAIN_SRC_DIR ?= src/main

# Compiler configuration
CXX = mpicxx
INC_FLAGS := $(addprefix -I,$(INC_DIRS))
BASE_FLAGS ?= -MMD -MP -m64 -std=c++11 -fopenmp
DEBUG_FLAGS ?= $(INC_FLAGS) $(BASE_FLAGS) -g
RELEASE_FLAGS ?= $(INC_FLAGS) $(BASE_FLAGS) -O3

# Sources which define main functions
MAIN_SRCS := $(shell find $(SRC_DIRS) -name *.cpp | grep $(MAIN_SRC_DIR))
MAIN_OBJS_RELEASE := $(MAIN_SRCS:%=$(OBJ_DIR_RELEASE)/%.o)
MAIN_OBJS_DEBUG := $(MAIN_SRCS:%=$(OBJ_DIR_DEBUG)/%.o)
MAIN_DEPS := $(MAIN_OBJS:.o=.d)

# "Subordinate" sources which do not define mains
SUB_SRCS := $(shell find $(SRC_DIRS) -name *.cpp | grep -v $(MAIN_SRC_DIR))
SUB_OBJS_RELEASE := $(SUB_SRCS:%=$(OBJ_DIR_RELEASE)/%.o)
SUB_OBJS_DEBUG := $(SUB_SRCS:%=$(OBJ_DIR_DEBUG)/%.o)
SUB_DEPS := $(SUB_OBJS:.o=.d)

# C++ object compilation - debug symbols - no optimisations
$(OBJ_DIR_DEBUG)/%.cpp.o: %.cpp
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(DEBUG_FLAGS) -c $< -o $@

# C++ release compilation - release - optimisations etc
$(OBJ_DIR_RELEASE)/%.cpp.o: %.cpp
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(RELEASE_FLAGS) -c $< -o $@

# Release build target
RELEASE_OBJS := $(SUB_OBJS_RELEASE) obj/release/src/main/mpicga.cpp.o
release: $(RELEASE_OBJS)
	@$(MKDIR_P) $(dir $(RELEASE_EXEC))
	$(CXX) $(RELEASE_OBJS) -o $(RELEASE_EXEC) $(LDFLAGS) -fopenmp

# Debug build target
DEBUG_OBJS := $(SUB_OBJS_DEBUG) obj/debug/src/main/mpicga.cpp.o
debug: $(DEBUG_OBJS)
	@$(MKDIR_P) $(dir $(DEBUG_EXEC))
	$(CXX) $(DEBUG_OBJS) -o $(DEBUG_EXEC) $(LDFLAGS) -fopenmp

# Build target for pattern generator
PATTERN_OBJS := $(SUB_OBJS_DEBUG) obj/debug/src/main/pattern.cpp.o
pattern: $(PATTERN_OBJS)
	@$(MKDIR_P) $(dir $(PATTERN_EXEC))
	$(CXX) $(PATTERN_OBJS) -o $(PATTERN_EXEC) $(LDFLAGS) -fopenmp

# Build target for unit test program
TEST_OBJS := $(SUB_OBJS_DEBUG) obj/debug/src/main/test.cpp.o
test: $(TEST_OBJS)
	@$(MKDIR_P) $(dir $(TEST_EXEC))
	$(CXX) $(TEST_OBJS) -o $(TEST_EXEC) $(LDFLAGS) -fopenmp

# Clean, be careful with this
.PHONY: clean
clean:
	@$(RM) -rv $(OBJ_DIR_BASE)

# Include dependencies
-include $(MAIN_DEPS) $(SUB_DEPS)

# Make directory
MKDIR_P ?= mkdir -p

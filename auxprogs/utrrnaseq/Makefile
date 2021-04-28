# Makefile for utrrnaseq

CXX ?= g++
CXXFLAGS := -O0 -g3 -pedantic -Wall -fmessage-length=0 ${CXXFLAGS}

SRC_DIR := src
OBJ_DIR := Debug/src
SRC_FILES := $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES))
DEP_FILES := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.d,$(SRC_FILES))

.PHONY: all clean test

# All Target
all: utrrnaseq

# Tool invocations
utrrnaseq: $(OBJ_FILES)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) $(LIBS) -o $@ $(OBJ_FILES)
	@echo 'Finished building target: $@'
	@echo ' '
	
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler: $(CXX)'
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) -c  -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

test: CPPFLAGS += -DTEST_MODE # Target-specific Variable
test: utrrnaseq
	@./utrrnaseq

clean:
	rm -rf $(OBJ_FILES) $(DEP_FILES) utrrnaseq
	@echo ' '

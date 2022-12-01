# Makefile for utrrnaseq

include ../../common.mk

CXX ?= g++
CXXFLAGS := -Wall -O2 -pedantic -fmessage-length=0 -std=c++11 ${CXXFLAGS}
INCLS   += $(INCLUDE_PATH_BOOST) # set boost include path in INCLUDE_PATH_BOOST, if boost is not installed system-wide
LDFLAGS += $(LIBRARY_PATH_BOOST) # set boost library path in LIBRARY_PATH_BOOST, if boost is not installed system-wide

SRC_DIR  := src
OBJ_DIR  := build
TEST_DIR := build/test
SOURCES  := $(wildcard $(SRC_DIR)/*.cpp)
OBJS     := $(SOURCES:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
DEPS     := $(SOURCES:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.d)

.PHONY: all prep clean prep_test test clean_test

# All Target
all: prep utrrnaseq

prep:
	mkdir -p $(OBJ_DIR)

# Tool invocations
utrrnaseq: $(OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

$(OBJS): $(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler: $(CXX)'
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLS) -c -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

test: CPPFLAGS += -DTEST_MODE # see Target-specific Variable in make
test: CXXFLAGS += -g3
test: LIBS += -lboost_unit_test_framework
test: clean prep utrrnaseq prep_test
	mv $(OBJ_DIR)/*.o $(TEST_DIR)/
	mv $(OBJ_DIR)/*.d $(TEST_DIR)/
	mv utrrnaseq $(TEST_DIR)/utrrnaseq_selftest
	@$(TEST_DIR)/utrrnaseq_selftest

prep_test:
	mkdir -p $(TEST_DIR)

clean_test:
	rm -rf $(TEST_DIR) *.gff

clean: clean_test
	rm -rf $(OBJS) $(DEPS) $(OBJ_DIR) utrrnaseq
	@echo ' '

-include $(DEPS)

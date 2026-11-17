# Compiler and linker flags for the optional build switches MIMALLOC, LTO and PGO, which are
# set in common.mk or on the command line, e.g. make LTO=true
# Included by src/Makefile and src/unittests/Makefile, after CXXFLAGS has been set with :=

ifneq (,$(findstring $(MIMALLOC),1 true True TRUE))
	LDLIBS   += $(MIMALLOC_LIB)
endif
ifneq (,$(findstring $(LTO),1 true True TRUE))
	CXXFLAGS += -flto=auto
endif
# AUGUSTUS_PGO_GENERATE: the workers of --cpus write their profile counters before _exit
ifeq ($(PGO),gen)
	CXXFLAGS += -fprofile-generate -fprofile-dir=$(PGO_DIR) -fprofile-update=single -DAUGUSTUS_PGO_GENERATE
else ifeq ($(PGO),use)
	CXXFLAGS += -fprofile-use -fprofile-dir=$(PGO_DIR) -fprofile-partial-training -Wno-missing-profile
endif

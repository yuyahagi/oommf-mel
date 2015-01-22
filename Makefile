TARGET = ../../darwin/oxs
SRCS = $(shell ls *.cc)
HEADS = $(shell ls *.h)
OBJS = ../../darwin/$(subst .cc,.o,$(SRCS))

LOCALDIR = ../../local
OOMMFDIR = ../../../..
TESTMIF = testdata/Hyz_02_run.mif
TESTMIF_STEP = testdata/Hyz_02_run_step.mif
TESTMIF_STEP_FILELIST = testdata/Hyz_02_run_step_filelist.mif

CXXFLAGS =

RM = rm
CP = cp
MV = mv
PUSHD = pushd
CXX = g++
CC = gcc
SED = sed
OOMMF = tclsh $(OOMMFDIR)/oommf.tcl

#all: $(TARGET)
all: pimake

#$(TARGET): $(OBJS) $(HEADS)
pimake:
	$(CP) $(SRCS) $(LOCALDIR)
	$(CP) $(HEADS) $(LOCALDIR)
	$(OOMMF) pimake -cwd $(OOMMFDIR)

#depend:
#	$(CXX) -MM -MG $(SRCS) > Makefile.depend
#	cat Makefile.depend

clean:
	$(RM) -f $(OBJS) $(TARGET) *~ \#*\#

test:
	$(OOMMF) oxsii $(TESTMIF)

test_step:
	$(OOMMF) oxsii $(TESTMIF_STEP)

test_step_filelist:
	$(OOMMF) oxsii $(TESTMIF_STEP_FILELIST)

#-include Makefile.depend


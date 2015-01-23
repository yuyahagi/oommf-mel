TARGET = ../../darwin/oxs
SRCS = $(shell ls *.cc)
HEADS = $(shell ls *.h)
COMMONSRCS = yy_mel_util.cc
COMMONHEADS = yy_mel_util.h
OBJSDIR = ../../darwin
OBJS = $(addprefix $(OBJSDIR)/,$(subst .cc,.o,$(SRCS)))

LOCALDIR = ../../local
OOMMFDIR = ../../../..

TESTMIF = testdata/Hyz_02_run.mif
TESTMIF_STEP = testdata/Hyz_02_run_step.mif
TESTMIF_STEP_FILELIST = testdata/Hyz_02_run_step_filelist.mif

.SUFFIXES: .cc .h .o

RM = rm
CP = cp
OOMMF = tclsh $(OOMMFDIR)/oommf.tcl

all: $(OBJS)
	$(OOMMF) pimake -cwd $(OOMMFDIR)

$(OBJSDIR)/%.o: %.cc %.h $(COMMONSRCS) $(COMMONHEADS)
	$(CP) $*.cc $(LOCALDIR)
	$(CP) $*.h $(LOCALDIR)

clean:
	$(RM) -f $(OBJS) $(TARGET) *~ \#*\#

test:
	$(OOMMF) oxsii $(TESTMIF)

test_step:
	$(OOMMF) oxsii $(TESTMIF_STEP)

test_step_filelist:
	$(OOMMF) oxsii $(TESTMIF_STEP_FILELIST)

# Dependency rule
#-include Makefile.depend

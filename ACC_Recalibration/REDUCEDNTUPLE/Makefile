VERSION := $(shell $(ROOTSYS)/bin/root-config --version | cut -b1-4)
ifeq ($(VERSION),5.27)
VERSIONP=
else
VERSIONP = $(VERSION)
endif

MARCH:=`root-config --arch`
CXX:=`root-config --cxx`

SRC_U= analyze.C 

OBJS_U=$(SRC_U:%.C=%.o)

CPPFLAGS=-D_PGTRACK_
CXXFLAGS= -O3 -Wno-write-strings $(CPPFLAGS)

ROOTDYNAMICLIBS := $(shell root-config --glibs)
ROOTDYNAMICLIBS += -lMinuit

EXTRALIBS= -ldl -lssl -lMinuit -lTreePlayer -lProof -lProofPlayer -lgfortran -lGeom -lEG

ifndef ONMACOSX
        EXTRALIBS += -lTMVA -lXMLIO -lMLP
        SOFLAGS = -shared
else
        EXTRALIBS += -lTMVA -lXMLIO -lMLP
        SOFLAGS = --dynamic --shared --dynamiclib -undefined dynamic_lookup
endif

GBATCH:= $(AMSSRC)/
GBATCHLIBS:= $(GBATCH)/lib/$(MARCH)$(VERSIONP)/
GBATCHINC:= $(GBATCH)/include/

AMSDYNAMICLIBS = $(GBATCHLIBS)/ntuple_slc6_PG.so
AMSSTATICLIBS = -L$(GBATCHLIBS) -lntuple_slc6_PG_dynamic

ACQTINC:= -I$(AMSSRC)/ACsoft/interface/
ACQTINC2:= -I$(AMSSRC)/ACsoft/acroot/

ifdef AMS_ACQT_INTERFACE
ifndef ONMACOSX
        AC_LIBS=-L$(QTDIR)/lib -lQtCore
else
        AC_LIBS = -F$(QTDIR) -framework QtCore
endif
endif

INCLUDES=-I./ -I$(ROOTSYS)/include -I$(GBATCHINC) $(ACQTINC) $(ACQTINC2)

#default: libAntiPG.so libAntiPGa.a static 
default: static 

static: analyze
shared: analyze_shared

analyze: $(SRC_U)
	$(CXX) -o $@ $(CXXFLAGS) $(INCLUDES) $^ $(AMSSTATICLIBS) `root-config --libs` $(ROOTDYNAMICLIBS) $(EXTRALIBS) $(AC_LIBS)

#analyze: $(SRC_U) AntiPGDict.o libAntiPGa.a
#	$(CXX) -o $@ $(CXXFLAGS) $(INCLUDES) $^ $(AMSSTATICLIBS) `root-config --libs` $(ROOTDYNAMICLIBS) $(EXTRALIBS) $(AC_LIBS) -L./ -lAntiPGa

analyze_shared: $(SRC_U) AntiPGDict.o libAntiPG.so
	$(CXX) -o $@ $(CXXFLAGS) $(INCLUDES) $^ $(AMSDYNAMICSLIBS) `root-config --libs` $(ROOTDYNAMICLIBS) $(EXTRALIBS) $(AC_LIBS) -L./ -lAntiPG

#------------------------------------------------------------------------------------------------------------------------------------------

TARGETLIB := AntiPG
LIBOBJS = AntiPG.o

lib$(TARGETLIB): shared-$(TARGETLIB) static-$(TARGETLIB)

shared-$(TARGETLIB): lib$(TARGETLIB).so

lib$(TARGETLIB).so: $(LIBOBJS) AntiPGDict.o
	@echo "Creating shared library \"$@\" ..."
	$(CXX) $(SOFLAGS) $(AMSDYNAMICLIBS) $(ROOTDYNAMICLIBS) $(EXTRALIBS) $(AC_LIBS) -o $@ $^

static-$(TARGETLIB): lib$(TARGETLIB)a.a

lib$(TARGETLIB)a.a:  $(LIBOBJS) AntiPGDict.o
	@echo "Creating static library \"$@\" ..."
	$(AR) rv $@ $^

AntiPGDict.cxx: AntiPG.o
	@$(ROOTSYS)/bin/rootcint -f $@ -c -p $(INCLUDES) $(CXXFLAGS) AntiPG.h LinkDef.h

%.o: %.cxx
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

%.o: %.C
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -fv *.o
	rm -fv analyze analyze_shared
	rm -fv AntiPGDict*
	rm -fv libAntiPG*


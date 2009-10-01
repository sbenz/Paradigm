LIBDAIDIR ?= /projects/sysbio/apps/${MACHTYPE}
BOOST_DIR ?= /projects/sysbio/apps/${MACHTYPE}

# Standard include directories
CCINC=-I${LIBDAIDIR}/include -I$(BOOST_DIR)/include/boost-1_38
CPPFLAGS=-O3 -W -Wall -Wextra -fPIC ${CCINC}

LIBDAIFLAGS=-DDAI_WITH_BP -DDAI_WITH_MF -DDAI_WITH_HAK -DDAI_WITH_LC -DDAI_WITH_TREEEP -DDAI_WITH_JTREE -DDAI_WITH_MR -DDAI_WITH_GIBBS
LIB_DIR=-L${LIBDAIDIR}/lib -L$(BOOST_DIR)/lib
LIBS=-ldai
LIBFLAGS=${LIBDAIFLAGS} ${LIB_DIR} ${LIBS}

SOURCES=configuration.cpp \
	evidencesource.cpp \
	pathwaytab.cpp \
	externVars.cpp


OBJECTS=$(SOURCES:.cpp=.o)

ALLSOURCES=$(SOURCES) pathwaytab2daifg.cpp main.cpp
ALLOBJECTS=$(ALLSOURCES:.cpp=.o)

EXECUTABLES=hgFactorGraph pathwaytab2daifg

all: $(EXECUTABLES)

include $(ALLSOURCES:.cpp=.d)

hgFactorGraph: main.o ${OBJECTS} 
	${CXX} ${CPPFLAGS} -o $@ $< ${OBJECTS} ${LIBFLAGS} 

pathwaytab2daifg: pathwaytab2daifg.o ${OBJECTS} 
	${CXX} ${CPPFLAGS} -o $@ $< ${OBJECTS} ${LIBFLAGS} 

clean:
	rm -f ${EXECUTABLES} ${ALLOBJECTS} $(ALLSOURCES:.cpp=.d)

%.d: %.cpp
	@set -e; rm -f $@; \
	$(CXX) -M $(CPPFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

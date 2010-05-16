## System configuration for library locations
## Defining any of these as environment variables will override these locations
LIBDAI_DIR ?= /projects/sysbio/apps/${MACHTYPE}
BOOST_DIR ?= /projects/sysbio/apps/${MACHTYPE}
BOOST_VER ?= boost-1_38

LIBDAI_INC ?= ${LIBDAI_DIR}/include
LIBDAI_LIB ?= ${LIBDAI_DIR}/lib

BOOST_INC ?= ${BOOST_DIR}/include/${BOOST_VER}
BOOST_LIB ?= ${BOOST_DIR}/lib

## Compilation Configuration
CCINC=-I${LIBDAI_INC} -I${BOOST_INC}
VERSION:=$(shell git describe --always) $(shell git diff --shortstat)
CPPFLAGS=-O3 -W -Wall -Wextra -fPIC ${CCINC} -D'VERSION="${VERSION}"'
LIBDAIFLAGS=-DDAI_WITH_BP -DDAI_WITH_MF -DDAI_WITH_HAK -DDAI_WITH_LC -DDAI_WITH_TREEEP -DDAI_WITH_JTREE -DDAI_WITH_MR -DDAI_WITH_GIBBS
LIB_DIR=-L${LIBDAI_LIB}
LIBS=-ldai
LIBFLAGS=${LIBDAIFLAGS} ${LIB_DIR} ${LIBS}
CPPFLAGS +=${LIBDAIFLAGS}
DEPDIR=.deps
DF=$(DEPDIR)/$(*).d

## Source files and executables
SOURCES=configuration.cpp \
	evidencesource.cpp \
	pathwaytab.cpp \
	externVars.cpp

OBJECTS=$(SOURCES:.cpp=.o)

ALLSOURCES=$(SOURCES) pathwaytab2daifg.cpp main.cpp
ALLOBJECTS=$(ALLSOURCES:.cpp=.o)

EXECUTABLES=paradigm pathwaytab2daifg

all: $(EXECUTABLES)

-include $(addprefix $(DEPDIR)/,$(ALLSOURCES:.cpp=.d))

paradigm: main.o ${OBJECTS} 
	${CXX} ${CPPFLAGS} -o $@ $< ${OBJECTS} ${LIBFLAGS} 

pathwaytab2daifg: pathwaytab2daifg.o ${OBJECTS} 
	${CXX} ${CPPFLAGS} -o $@ $< ${OBJECTS} ${LIBFLAGS} 

clean:
	rm -f ${EXECUTABLES} ${ALLOBJECTS}
	rm -Rf $(DEPDIR)

tests: $(EXECUTABLES)
	cd testdata && sh runtests.sh

%.o: %.cpp
	@mkdir -p $(DEPDIR)
	$(COMPILE.cpp) -MMD -MF $(DF) $(OUTPUT_OPTION) $<

TAGS:
	etags *.h *.cpp

CPPFLAGS = -c -O3 -Wall -static
LDFLAGS = 
SOURCES = src/GlyRotHelper.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = bin/GlyRotHelper
INCLUDE = /path/to/gromacs-5.1/src

ifeq "$(origin GMXLDLIB)" "undefined"
  $(error "GMXLDLIB not found, please source GMXRC")
else
  export PKG_CONFIG_PATH:=${PKG_CONFIG_PATH}:${GMXLDLIB}/pkgconfig
endif

CPPFLAGS += `pkg-config --cflags libgromacs`
LDFLAGS += `pkg-config --libs libgromacs`

CPPFLAGS += -I$(INCLUDE)

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CXX) $(CPPFLAGS) $< -o $@

clean:
	rm $(OBJECTS) $(EXECUTABLE)

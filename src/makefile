CCX = clang++
#CFLAGS = -c -Wall -g -Os
CFLAGS = -std=c++11 $(shell python-config --includes) $(shell python -c "import numpy; print '-I' + numpy.get_include()") $(shell python-config --ldflags)
LD = $(CCX)
LDFLAGS = -lfoo
INCLUDE = -I.

TARGET = Simulator.exe

OBJECTS = $(patsubst %.cpp, %.o, $(shell find . -name "*.cpp"), $(shell find . -name "*.h"))

$(info $(OBJECTS))

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(LD) -o $@ $^ $(LDFLAGS)

#You don't even need to be explicit here,
#compiling C files is handled automagically by Make.

%.o: %.cpp
	$(CCX) $(CFLAGS) $(INCLUDE) $^ -o $@

clean:
	rm $(TARGET) $(OBJECTS)

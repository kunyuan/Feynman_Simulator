CXX=icc
#CXX=clang++
CXXFLAGS=-O3

ODIR=.

_OBJ = utility.o complex.o timer.o benchmark.o fft.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.c
	$(CXX) -c -o $@ $< $(CXXFLAGS)

benchmark.exe: $(OBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core

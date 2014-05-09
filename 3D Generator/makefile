SRC = FiniteAmortizedNoise3D.cpp InfiniteAmortizedNoise3D.cpp CPUtime.cpp Main.cpp MurmurHash3.cpp
EXE = generator

all: $(SRC) $(EXE)

$(EXE): $(SRC)
	g++ -O3 -ffast-math -o $(EXE) $(SRC)

cleanup: 
	rm -f  $(EXE) *.tga

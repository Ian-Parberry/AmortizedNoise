SRC = MurmurHash3.cpp FiniteAmortizedNoise2D.cpp CPUtime.cpp InfiniteAmortizedNoise2D.cpp infiniteperlin.cpp infinitesmoothperlin.cpp Main.cpp perlin.cpp
EXE = evaluator

all: $(SRC) $(EXE)

$(EXE): $(SRC)
	g++ -O3 -ffast-math -o $(EXE) $(SRC)

cleanup: 
	rm -f  $(EXE) 

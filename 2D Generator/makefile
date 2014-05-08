SRC = FiniteAmortizedNoise2D.cpp InfiniteAmortizedNoise2D.cpp CPUtime.cpp main.cpp MurmurHash3.cpp SaveImage.cpp
EXE = generator

all: $(SRC) $(EXE)

$(EXE): $(SRC)
	g++ -O3 -ffast-math -o $(EXE) $(SRC)

cleanup: 
	rm -f  $(EXE) *.tga

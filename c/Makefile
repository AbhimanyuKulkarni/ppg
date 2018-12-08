SRC		:= src
BIN		:= bin
OBJS	:= $(BIN)/main.o \
					$(BIN)/ppg_algorithm.o \
					$(BIN)/datagen.o \
					$(BIN)/setup.o \
					$(BIN)/fixed_point.o \
					$(BIN)/fixed_point_experiments.o

all: $(BIN)/c_reconstruction $(BIN)/fixed_point_expt

$(BIN)/c_reconstruction: $(BIN)/main.o $(BIN)/ppg_algorithm.o \
													$(BIN)/datagen.o $(BIN)/setup.o
	gcc $^ -o $@ -lm

$(BIN)/fixed_point_expt: $(BIN)/fixed_point.o \
													$(BIN)/fixed_point_experiments.o
	gcc $^ -o $@ -lm

$(BIN)/%.o: $(SRC)/%.c
	gcc -c -O3 -Wall $^ -o $@

.PHONY: clean

clean:
	rm $(BIN)/*

INCLUDE	:= include
SRC			:= src
BIN			:= bin

FRAC_TYPE ?= -DPPG_USE_FLOAT64
CFLAGS	:= -c -g -Wall -I$(INCLUDE) $(FRAC_TYPE)

OBJS	:= $(BIN)/main.o \
					$(BIN)/ppg_algorithm.o \
					$(BIN)/datagen.o \
					$(BIN)/setup.o \
					$(BIN)/fxp_ops.o \
					$(BIN)/fixed_point_experiments.o

all: $(BIN)/c_reconstruction $(BIN)/fixed_point_expt

$(BIN)/c_reconstruction: $(BIN)/main.o $(BIN)/ppg_algorithm.o \
													$(BIN)/datagen.o $(BIN)/setup.o \
													$(BIN)/fxp_ops.o
	gcc $^ -o $@ -lm

$(BIN)/fixed_point_expt: $(BIN)/fxp_ops.o \
													$(BIN)/fixed_point_experiments.o
	gcc $^ -o $@ -lm

$(BIN)/%.o: $(SRC)/%.c
	gcc $(CFLAGS) $^ -o $@

.PHONY: clean

clean:
	rm $(BIN)/*

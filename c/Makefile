SRC		:= src
BIN		:= bin
OBJS	:= $(BIN)/main.o \
					$(BIN)/ppg_algorithm.o \
					$(BIN)/datagen.o \
					$(BIN)/setup.o

all: $(BIN)/c_reconstruction

$(BIN)/c_reconstruction: $(OBJS)
	gcc $^ -o $@ -lm

$(BIN)/%.o: $(SRC)/%.c
	gcc -c -O3 -Wall $^ -o $@

.PHONY: clean

clean:
	rm $(BIN)/*

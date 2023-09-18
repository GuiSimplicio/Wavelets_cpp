#adding rule to Makefile

ROOTLIB := $(shell root-config --libs)
ROOTINC := $(shell root-config --incdir)

#Makefile complete and modified

# Makefile

BINDIR := bin
LIBDIR := lib

CCFLAGS := -pedantic

CC := g++ -std=c++14 

# src/ (declarcaoes de funcoes, de calsses + codigo)
# main/ (programas principais)
# bin/ (temporarios, .o, .exe)
# lib/ (bibliotecas) biblioteca FC

# making library
# - static: .a
# - shareable: .so

VPATH = main:src

ROOTLIB := $(shell root-config --libs)
ROOTINC := $(shell root-config --incdir)

SRC := $(wildcard src/*.C)
OBJ := $(patsubst %.C, $(BINDIR)/%.o, $(notdir $(SRC)))
INC := $(wildcard src/*.h)

lib: $(LIBDIR)/libFC.a

$(LIBDIR)/libFC.a: $(OBJ) 
	@echo make lib...
	ar ruv $@ $^
	ranlib $@

%.exe: $(BINDIR)/%.o $(LIBDIR)/libFC.a 
	@echo compiling and linking... 
	$(CC) -I src $< -o $(BINDIR)/$@ -L lib -l FC $(ROOTLIB) -lMathMore

#Integral: $(BINDIR)/Integral.o $(LIBDIR)/libFC.a 
#	@echo Compiling and linking... 
#	$(CC) -I src $(BINDIR)/Integral.o -o $(BINDIR)/Integral.exe -L lib -l FC $(ROOTLIB) -lMathMore
#	@echo
# SE QUISERES MAIS QUE UMA SÓ COM UM MAKE.
#trab03: $(BINDIR)/main_trab03_1.o $(BINDIR)/main_trab03_2.o $(LIBDIR)/libFC.a 
#	@echo Compiling and linking... 
#	$(CC) -I src $(BINDIR)/main_trab03_2.o -o $(BINDIR)/main_$@_2.exe -L lib -l FC $(ROOTLIB)
#	@echo
#	$(CC) -I src $(BINDIR)/main_trab03_1.o -o $(BINDIR)/main_$@_1.exe -L lib -l FC $(ROOTLIB)
#	@echo

$(BINDIR)/%.o: %.C | $(INC)
	@echo compiling... $<
	$(CC) -I src -I $(ROOTINC) -c $< -o $@


######### clean

tilde := $(wildcard */*~) $(wildcard *~)
exe := $(wildcard */*.exe) $(wildcard *.exe)
obj := $(wildcard */*.o) $(wildcard *.o) $(wildcard */*.so) $(wildcard */*.pcm) $(wildcard */*.d)

clean:
	@echo cleaning dir...
	rm -f $(exe) $(obj) $(tilde) $(LIBDIR)/libFC.a



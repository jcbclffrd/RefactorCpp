CC = g++

GSL_DIR = /usr

INC = -I$(GSL_DIR)
# Debug build: use -g for debug symbols, -O0 for no optimization
# Release build: use -O3 for optimization
DEBUG_FLAGS = -g -O0
RELEASE_FLAGS = -O3

# Default to debug build - change to $(RELEASE_FLAGS) for optimized build
CFLAGS = $(DEBUG_FLAGS) $(INC) 

LIB_DIRS = -L$(GSL_DIR)/lib
LIBS = -lgsl -lgslcblas -lm -ljsoncpp
LFLAGS = $(DEBUG_FLAGS) $(LIB_DIRS) $(LIBS)

BIN =  seq2exp mcp_demo
all: $(BIN)

clean:
	rm -f $(BIN) *.o

Tools.o : Tools.h Tools.cpp
	$(CC) $(CFLAGS) -c Tools.cpp
ExprPredictor.o : Tools.h ExprPredictor.h  ExprPredictor.cpp
	$(CC) $(CFLAGS) -c ExprPredictor.cpp
seq2exp.o : Tools.h ExprPredictor.h seq2exp.cpp
	$(CC) $(CFLAGS) -c seq2exp.cpp
mcp_tools.o : mcp_tools.h mcp_tools.cpp ExprPredictor.h Tools.h
	$(CC) $(CFLAGS) -c mcp_tools.cpp
mcp_demo.o : mcp_tools.h mcp_demo.cpp
	$(CC) $(CFLAGS) -c mcp_demo.cpp

seq2exp : Tools.o ExprPredictor.o seq2exp.o 
	$(CC) -o $@ Tools.o ExprPredictor.o seq2exp.o $(LFLAGS) 

mcp_demo : Tools.o ExprPredictor.o mcp_tools.o mcp_demo.o
	$(CC) -o $@ Tools.o ExprPredictor.o mcp_tools.o mcp_demo.o $(LFLAGS) 



# Compliler
CXX := g++

# Compiler flags
CXX_FLAGS := -std=c++11
DBG_FLAGS := DEBUG
NAIVE_IMPL := NAIVE
CSR_IMPL := CSR
CSC_IMPL := CSC
COLUMN_BLOCK_CSR_IMPL := COLUMN_BLOCK_CSR
ROW_BLOCK_CSR_IMPL := ROW_BLOCK_CSR
CACHE_BLOCK_IMPL := CACHE_BLOCK

# Targets
NAIVE := naive
CSR := csr
CSC := csc
COL_BLK_CSR := cloumn_csr
ROW_BLK_CSR := row_csr
CACHE_BLOCK := cache-block

# Include directories
INC := .

# Rules
naive: CXX_FLAGS += -D$(DBG_FLAGS) -D$(NAIVE_IMPL)
naive: naive.o main.o
	$(CXX) -o $(NAIVE) $+

csr: CXX_FLAGS += -D$(DBG_FLAGS) -D$(CSR_IMPL)
csr: csr.o main.o naive.o
	$(CXX) -o $(CSR) $+

csc: CXX_FLAGS += -D$(DBG_FLAGS) -D$(CSC_IMPL)
csc: csc.o main.o naive.o
	$(CXX) -o $(CSC) $+

column-csr: CXX_FLAGS += -D$(DBG_FLAGS) -D$(COLUMN_BLOCK_CSR_IMPL)
column-csr: csr_column_block.o main.o naive.o
	$(CXX) -o $(COL_BLK_CSR) $+

naive.o: naive.cpp
	$(CXX) -I $(INC) $(CXX_FLAGS) -o $@ -c $<

csr.o: csr.cpp
	$(CXX) -I $(INC) $(CXX_FLAGS) -o $@ -c $<

csc.o: csc.cpp
	$(CXX) -I $(INC) $(CXX_FLAGS) -o $@ -c $<

csr_column_block.o: csr_column_block.cpp
	$(CXX) -I $(INC) $(CXX_FLAGS) -o $@ -c $<

main.o: main.cpp utils.h
	$(CXX) -I $(INC) $(CXX_FLAGS) -o $@ -c $<
	#$(CXX) -fverbose-asm -O3 -g -S main.cpp

clean:
	@rm -f $(NAIVE) $(CSR) $(CSC) *.o *.s


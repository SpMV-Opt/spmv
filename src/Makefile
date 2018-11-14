# Compliler
CXX := g++

# Compiler flags
CXX_FLAGS := -std=c++11
DBG_FLAGS := -DDEBUG
NAIVE_IMPL := -DNAIVE
CSR_IMPL := -DCSR
CSC_IMPL := -DCSC
CSR_OMP_IMPL := -DCSR_OMP
CSC_OMP_IMPL := -DCSC_OMP
COLUMN_BLOCK_CSR_IMPL := -DCOLUMN_BLOCK_CSR
ROW_BLOCK_CSR_IMPL := -DROW_BLOCK_CSR
CACHE_BLOCK_IMPL := -DCACHE_BLOCK

# Targets
NAIVE := naive
CSR := csr
CSC := csc
CSR_OMP := csr_omp
CSC_OMP := csc_omp
COL_BLK_CSR := cloumn_csr
ROW_BLK_CSR := row_csr
CACHE_BLOCK := cache-block

# Include directories
INC := .

# Rules
naive: CXX_FLAGS += $(DBG_FLAGS) $(NAIVE_IMPL)
naive: naive.o main.o
	$(CXX) -o $(NAIVE) $+

csr: CXX_FLAGS += $(DBG_FLAGS) $(CSR_IMPL)
csr: csr.o main.o naive.o
	$(CXX) -o $(CSR) $+

csr-omp: CXX_FLAGS += $(DBG_FLAGS) $(CSR_OMP_IMPL)
csr-omp: csr.o main.o naive.o
	$(CXX) -o $(CSR_OMP) $+

csc: CXX_FLAGS += $(DBG_FLAGS) $(CSC_IMPL)
csc: csc.o main.o naive.o
	$(CXX) -o $(CSC) $+

csc-omp: CXX_FLAGS += $(DBG_FLAGS) $(CSC_OMP_IMPL)
csc-omp: csc.o main.o naive.o
	$(CXX) -o $(CSC_OMP) $+

column-csr: CXX_FLAGS += $(DBG_FLAGS) $(COLUMN_BLOCK_CSR_IMPL)
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


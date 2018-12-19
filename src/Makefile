# Compliler
CXX := g++
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	OMP_CXX := g++
endif
ifeq ($(UNAME_S),Darwin)
	OMP_CXX := g++-8
endif

# Compiler flags
CXX_FLAGS := -std=c++11
OMP_FLAGS := -fopenmp
NAIVE_IMPL := -DNAIVE
CSR_IMPL := -DCSR
CSC_IMPL := -DCSC
CSR_OMP_IMPL := -DCSR_OMP
CSC_OMP_IMPL := -DCSC_OMP
COLUMN_BLOCK_CSR_IMPL := -DCOLUMN_BLOCK_CSR
ROW_BLOCK_CSR_IMPL := -DROW_BLOCK_CSR
CACHE_BLOCK_IMPL := -DCACHE_BLOCK

ifeq ($(VER), TRUE)
	DBG_FLAGS := -DDEBUG -DRESULT_VERIFY
else
	DBG_FLAGS :=
endif
PAPI_LIB := -lpapi

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
	$(CXX) $(PAPI_LIB) -o $(NAIVE) $+

csr: CXX_FLAGS += $(DBG_FLAGS) $(CSR_IMPL)
csr: csr.o main.o naive.o
	$(CXX) $(PAPI_LIB) -o $(CSR) $+

csr-omp: CXX_FLAGS += $(DBG_FLAGS) $(CSR_OMP_IMPL)
csr-omp: csr_omp.o main.o naive.o
	$(OMP_CXX) $(PAPI_LIB) $(OMP_FLAGS) -o $(CSR_OMP) $+

csc: CXX_FLAGS += $(DBG_FLAGS) $(CSC_IMPL)
csc: csc.o main.o naive.o
	$(CXX) $(PAPI_LIB) -o $(CSC) $+

csc-omp: CXX_FLAGS += $(DBG_FLAGS) $(CSC_OMP_IMPL)
csc-omp: csc_omp.o main.o naive.o
	$(OMP_CXX) $(PAPI_LIB) $(OMP_FLAGS) -o $(CSC_OMP) $+

naive.o: naive.cpp
	$(CXX) -I $(INC) $(CXX_FLAGS) -o $@ -c $<

csr.o: csr.cpp
	$(CXX) -I $(INC) $(CXX_FLAGS) -o $@ -c $<

csc.o: csc.cpp
	$(CXX) -I $(INC) $(CXX_FLAGS) -o $@ -c $<

csr_omp.o: csr_omp.cpp
	$(CXX) -I $(INC) $(OMP_FLAGS) $(CXX_FLAGS) -o $@ -c $<

csc_omp.o: csc_omp.cpp
	$(CXX) -I $(INC) $(OMP_FLAGS) $(CXX_FLAGS) -o $@ -c $<

csr_column_block.o: csr_column_block.cpp
	$(CXX) -I $(INC) $(CXX_FLAGS) -o $@ -c $<

main.o: main.cpp utils.h
	$(CXX) -I $(INC) $(CXX_FLAGS) -o $@ -c $<
	#$(CXX) -fverbose-asm -O3 -g -S main.cpp

clean:
	@rm -f *.o *.s


# Compliler
CXX := g++
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    OMP_FLAGS := -fopenmp -DOMP
	PAPI_LIB := -lpapi
PAPI_FLAG := -DPAPI
endif
ifeq ($(UNAME_S),Darwin)
    OMP_FLAGS :=
	PAPI_LIB :=
	PAPI_FLAG :=
endif

# Compiler flags
CXX_FLAGS := -std=c++11
NAIVE_IMPL := -DNAIVE
CSR_IMPL := -DCSR
BCSR_IMPL := -DBCSR
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

# Targets
NAIVE := naive
CSR := csr
BCSR := bcsr
CSC := csc
CSR_OMP := csr_omp
CSC_OMP := csc_omp
COL_BLK_CSR := cloumn_csr
ROW_BLK_CSR := row_csr
CACHE_BLOCK := cache-block

# Include directories
INC := .

# Rules
naive: CXX_FLAGS += $(DBG_FLAGS) $(NAIVE_IMPL) $(PAPI_FLAG)
naive: naive.o main.o
	$(CXX) -o $(NAIVE) $+ $(PAPI_LIB)

csr: CXX_FLAGS += $(DBG_FLAGS) $(CSR_IMPL) $(PAPI_FLAG)
csr: csr.o main.o naive.o
	$(CXX) -o $(CSR) $+ $(PAPI_LIB)

bcsr: CXX_FLAGS += $(DBG_FLAGS) $(BCSR_IMPL) $(PAPI_FLAG)
bcsr: bcsr.o naive.o bcsr_main.o
	$(CXX) -o $(BCSR) $+ $(PAPI_LIB) 

csr-omp: CXX_FLAGS += $(DBG_FLAGS) $(CSR_OMP_IMPL) $(PAPI_FLAG)
csr-omp: csr_omp.o main.o naive.o
	$(CXX) $(OMP_FLAGS) -o $(CSR_OMP) $+ $(PAPI_LIB)

csc: CXX_FLAGS += $(DBG_FLAGS) $(CSC_IMPL) $(PAPI_FLAG)
csc: csc.o main.o naive.o
	$(CXX) -o $(CSC) $+ $(PAPI_LIB)

csc-omp: CXX_FLAGS += $(DBG_FLAGS) $(CSC_OMP_IMPL) $(PAPI_FLAG)
csc-omp: csc_omp.o main.o naive.o
	$(CXX) $(OMP_FLAGS) -o $(CSC_OMP) $+ $(PAPI_LIB)

naive.o: naive.cpp
	$(CXX) -I $(INC) $(CXX_FLAGS) -o $@ -c $<

csr.o: csr.cpp
	$(CXX) -I $(INC) $(CXX_FLAGS) -o $@ -c $<

bcsr.o: bcsr.cpp
	$(CXX) -I $(INC) $(CXX_FLAGS) -o $@ -c $<

bcsr_main.o: bcsr_main.cpp utils.h
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


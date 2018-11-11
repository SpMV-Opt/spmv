# Compliler
CXX := g++

# Compiler flags
CXX_FLAGS := -std=c++11
DBG_FLAGS := DEBUG
NAIVE_IMPL := NAIVE
CSR_IMPL := CSR
CACHE_BLOCK_IMPL := CACHE_BLOCK

# Targets
NAIVE := naive
CSR := csr
CACHE_BLOCK := cache-block

# Include directories
INC := .

# Rules
naive: CXX_FLAGS += -D$(DBG_FLAGS) -D$(NAIVE_IMPL)
naive: naive.o main.o
	$(CXX) -o $(NAIVE) $+

csr: CXX_FLAGS += -D$(DBG_FLAGS) -D$(CSR_IMPL)
csr: csr.o main.o
	$(CXX) -o $(NAIVE) $+

naive.o: naive.cpp
	$(CXX) -I $(INC) $(CXX_FLAGS) -o $@ -c $<

csr.o: csr.cpp
	$(CXX) -I $(INC) $(CXX_FLAGS) -o $@ -c $<

main.o: main.cpp
	$(CXX) -I $(INC) $(CXX_FLAGS) -o $@ -c $<
	#$(CXX) -fverbose-asm -O3 -g -S main.cpp

clean:
	@rm -f $(NAIVE) $(CSR) *.o *.s


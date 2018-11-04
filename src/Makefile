all:
	g++ main.cpp -o foo
	g++ -fverbose-asm -O3 -g -S main.cpp
clean:
	@rm -f foo *.s

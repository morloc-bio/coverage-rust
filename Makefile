.PHONY: all test clean

all: test

test:
	morloc make -o test test.loc
	./test | tee test.out
	tail -n 1 test.out | grep -qx true

clean:
	rm -rf test test-build test.out

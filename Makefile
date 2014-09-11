

all: install

%.o: %.f
gfortran -c -fpic -Wall $<

%.o: %.c
gcc -c -fpic -Wall $<

libcsnopt.so: snopt.o snerror.o snset.o snget.o sn_open.o
gcc -shared -o $@ $^ -lsnopt7 -lsnblas -lsnprint7 -lgfortran

install: libcsnopt.so
sudo cp $< /usr/lib

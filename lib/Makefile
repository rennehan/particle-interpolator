# If " final link failed: nonrepresentable section on output" happens
# uncomment line below and comment CC = gcc.
CC = gcc -fPIC
#CC = gcc  

default: libinterpolate.a

libinterpolate.a: interpolate.o
	ar rcs $@ $^
    
interpolate.o: interpolate.c interpolate.h 
	$(CC) -c $<

clean:
	rm *.o *.a

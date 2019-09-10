LIB_DIR = lib

default: particleinterpolator

particleinterpolator: setup.py particleinterpolator.pyx $(LIB_DIR)/libinterpolate.a
	python setup.py build && python setup.py install && rm -f particleinterpolator.c && rm -rf build

$(LIB_DIR)/libinterpolate.a:
	make -C $(LIB_DIR) libinterpolate.a

clean:
	rm *.so

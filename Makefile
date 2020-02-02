CPPFLAGS = -Wall -Wextra

wfc: wfc.cpp
	g++ $(CPPFLAGS) -o wfc.o wfc.cpp

.PHONY: clean

clean:
	rm -f wfc.o

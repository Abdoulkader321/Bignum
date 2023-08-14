CFLAGS   =-Wall -Wextra -std=c++11 -g -O2
CPPFLAGS =-I../include -DDEBUG
OBJFILES = Main.o Bignum.o
TARGET   = Main

all:$(TARGET)

$(TARGET): $(OBJFILES)
	g++ $(CFLAGS) -o $@ $^ $(LDFLAGS)

Rsa: Rsa.o Bignum.o
	g++ $(CFLAGS) -o $@ $^ $(LDFLAGS)

Rsa.o: Rsa.cpp Bignum.cpp Bignum.hpp
	g++ $(CFLAGS) $(CPPFLAGS) -c $<

Main.o: Main.cpp Bignum.cpp Bignum.hpp
	g++ $(CFLAGS) $(CPPFLAGS) -c $<

Bignum.o: Bignum.cpp Bignum.hpp
	g++ $(CFLAGS) $(CPPFLAGS) -c $<

clean:
	@rm -f $(OBJFILES) $(TARGET) *~

.PHONY: all clean

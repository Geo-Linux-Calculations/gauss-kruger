CC        = gcc
CFLAGS    = -Wall
LDFLAGS   = -lm

default: cli lib
lib: libgausskruger.a libgausskruger.so
cli: gausskruger
all: default

.PHONY: clean
clean:
	rm -f src/*.o gausskruger libgausskruger.a libgausskruger.so

# CLI tool
gausskruger: src/gausskruger_cli.o src/gausskruger.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) -s

# Static library
libgausskruger.a: src/gausskruger.o
	ar rcs $@ $^

# Shared library
libgausskruger.so: src/gausskruger.o
	$(CC) $(CFLAGS) -shared $^ -o $@ $(LDFLAGS) -s

src/gausskruger.o: src/gausskruger.c src/gausskruger.h
src/gausskruger_cli.o: src/gausskruger_cli.c

CC=gcc
INCPATH=./inc/
OUTPATH=./out/
SRCPATH=./src/
CFLAGS=-std=c99 --coverage -I$(INCPATH) -Wall -Wpedantic -Werror -Wextra
LFLAGS=--coverage
SOURSES=main.c
OBJECTS=$(OUTPATH)$(SOURSES:.c=.o)
EXECUTABLE=app.exe
.PHONY: clean
.NOTPARALLEL: release debug

release: CFLAGS=-std=c99 -I$(INCPATH) -Wall -Wpedantic -Werror -Wextra
release: LFLAGS= -o2
release: $(OUTPATH)release.lastbuildstate app.exe

debug: CFLAGS=-std=c99 --coverage -I$(INCPATH) -Wall -Wpedantic -Werror -Wextra -ggdb -g3
debug:$(OUTPATH)debug.lastbuildstate app.exe

$(OUTPATH)debug.lastbuildstate:
	mkdir -p $(OUTPATH)
	rm -fv *.exe
	rm -fv $(OUTPATH)*.o
	rm -fv $(OUTPATH)*.gcno
	rm -fv $(OUTPATH)*.gcda
	rm -fv $(OUTPATH)release.lastbuildstate
	touch $(OUTPATH)debug.lastbuildstate

$(OUTPATH)release.lastbuildstate:
	mkdir -p $(OUTPATH)
	rm -fv *.exe
	rm -fv $(OUTPATH)*.o
	rm -fv $(OUTPATH)*.gcno
	rm -fv $(OUTPATH)*.gcda
	rm -fv $(OUTPATH)debug.lastbuildstate
	touch $(OUTPATH)release.lastbuildstate
	
$(EXECUTABLE): $(OBJECTS)
	mkdir -p $(OUTPATH)
	$(CC) $(LFLAGS) $(OBJECTS) -o $@

$(OBJECTS): $(SOURSES)
	mkdir -p $(OUTPATH)
	$(CC) $(CFLAGS) -c $< -o $(OBJECTS)

clean:
	mkdir -p $(OUTPATH)
	rm -fv *.exe
	rm -fv $(OUTPATH)*.o
	rm -fv $(OUTPATH)*.gcno
	rm -fv $(OUTPATH)*.gcda
	rm -fv $(OUTPATH)*.lastbuildstate

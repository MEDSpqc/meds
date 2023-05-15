SHELL := /bin/bash

PARAM = toy

LIBS = 
CC := gcc

ifdef DEBUG
CFLAGS := -g -Wall -D$(PARAM) -DDEBUG
OBJDIR := debug/$(PARAM)
else
CFLAGS := -O3 -Wall -D$(PARAM)
OBJDIR := build/$(PARAM)
endif

EXES = test bench

TARGETS := ${EXES:%=$(OBJDIR)/%}

.PHONY: default all clean

default: $(EXES)

OBJECTS = meds.o util.o seed.o osfreq.o fips202.o matrixmod.o bitstream.o
HEADERS = $(wildcard *.h) params.h

BUILDOBJ := ${OBJECTS:%=$(OBJDIR)/%}


params.h: gen_param.py params.py
	python $< > $@


$(EXES) : % :
	@make $(OBJDIR)/$(@F)

RUN: test
	@echo ""
	@$(OBJDIR)/test

RUN_ALL:
	for p in `./params.py -l`; do \
	echo "\nRunning par set" $$p; \
	make RUN PARAM=$$p; \
	echo ""; \
	done;

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(BUILDOBJ) : $(OBJDIR)/%.o: %.c $(HEADERS) | $(OBJDIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(TARGETS) : $(OBJDIR)/%: %.c $(BUILDOBJ)
	$(CC) $(@F).c $(BUILDOBJ) $(CFLAGS) $(LIBS) -o $@

BENCH:
	make bench PARAM=$(PARAM)
	$(OBJDIR)/bench | ./proc_bench.py >> bench.txt
	cat bench.txt

BENCH_ALL:
	rm -f bench.txt
	for p in `./params.py -l`; do \
	echo "Running par set" $$p; \
	make PARAM=$$p BENCH; \
	echo ""; \
	done

all:
	for p in `./params.py -l`; do \
	make PARAM=$$p; \
	make PARAM=$$p DEBUG=true; \
	done

clean:
	rm -rf params.h bench.txt build/ debug/ __pycache__/


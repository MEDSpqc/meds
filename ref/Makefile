SHELL := /bin/bash

PARAM = toy

LIBS = -lssl -lcrypto
CC := gcc

INC := ../lib

ifdef DEBUG
CFLAGS := -g -Wall -D$(PARAM) -DDEBUG
OBJDIR := debug/$(PARAM)
else
CFLAGS := -O3 -Wall -D$(PARAM)
OBJDIR := build/$(PARAM)
endif

ifdef LOG
CFLAGS += -DDEBUG
endif

CFLAGS += -I$(OBJDIR) -INIST -I$(INC)

EXES = test bench KAT_test

TARGETS := ${EXES:%=$(OBJDIR)/%}

.PHONY: default all clean

default: $(EXES)

OBJECTS = meds.o util.o seed.o osfreq.o fips202.o matrixmod.o bitstream.o randombytes.o
HEADERS = $(wildcard *.h) params.h $(OBJDIR)/api.h

BUILDOBJ := ${OBJECTS:%=$(OBJDIR)/%}


params.h: params.py
	python $< -p > $@

$(OBJDIR)/api.h: params.py $(OBJDIR)
	python $< -a $(PARAM) > $@


$(EXES) : % :
	@make $(OBJDIR)/$(@F)

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(BUILDOBJ) : $(OBJDIR)/%.o: %.c $(HEADERS) | $(OBJDIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(TARGETS) : $(OBJDIR)/%: %.c $(BUILDOBJ)
	$(CC) $(@F).c $(BUILDOBJ) $(CFLAGS) $(LIBS) -o $@


RUN: test
	@echo ""
	@$(OBJDIR)/test

RUN_ALL:
	for p in `./params.py -l`; do \
	echo "\nRunning par set" $$p; \
	make RUN PARAM=$$p; \
	echo ""; \
	done;

.PHONY: debug KAT

KAT: 
	$(MAKE) KAT_test PARAM=$(PARAM)
	$(OBJDIR)/KAT_test

KAT_ALL:
	@for p in `./params.py -l`; do \
	echo "Running par set" $$p; \
	make PARAM=$$p KAT; \
	echo ""; \
	done

BENCH:
	make bench PARAM=$(PARAM)
	$(OBJDIR)/bench | ./proc_bench.py >> bench.txt

BENCH_ALL:
	rm -f bench.txt
	for p in `./params.py -l`; do \
	echo "Running par set" $$p; \
	make PARAM=$$p BENCH; \
	echo ""; \
	cat bench.txt; \
	done

all:
	for p in `./params.py -l`; do \
	make PARAM=$$p; \
	make PARAM=$$p DEBUG=true; \
	done


clean:
	rm -rf build/ debug/ PQCsignKAT_*


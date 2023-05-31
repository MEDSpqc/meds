PARAM=toy

SK_SIZE := $(shell ./get_sk_size.py $(PARAM))

all: PQCsignKAT_$(SK_SIZE).rsp

MEDS.py: MEDS.sage.py
	sage -preparse MEDS.sage
	cp MEDS.sage.py MEDS.py

PQCsignKAT_$(SK_SIZE).rsp: KAT_test.sage MEDS.py SeedTree.py
	sage KAT_test.sage -c 5 $(PARAM) 3> PQCsignKAT_$(SK_SIZE).log

clean:
	@for p in `./params.py -l`; do \
	rm -rf IO-$$p.*; \
	done


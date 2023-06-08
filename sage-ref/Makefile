PARAM=toy

MEDS.py: MEDS.sage
	sage -preparse MEDS.sage
	cp MEDS.sage.py MEDS.py

debug:
	make -C ../ref DEBUG=true RUN 2> /tmp/c.out
	./MEDS.sage $(PARAM) 3> /tmp/sage.out
	diff /tmp/c.out /tmp/sage.out | head -n 10

KAT: MEDS.py
	make -C ../ref KAT PARAM=$(PARAM)
	./KAT_check.sage ../ref/PQCsignKAT_$(PARAM).rsp $(PARAM)

KAT_ALL:
	@for p in `./params.py -l`; do \
	echo "\nRunning par set" $$p "\n"; \
	make PARAM=$$p KAT; \
	echo ""; \
	done

clean:
	@for p in `./params.py -l`; do \
	rm -rf IO-$$p.*; \
	done


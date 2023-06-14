PARAM=toy

meds.py: meds.sage
	sage -preparse meds.sage
	cp meds.sage.py meds.py

debug:
	make -C ../ref DEBUG=true RUN 2> /tmp/c.out
	./meds.sage $(PARAM) 3> /tmp/sage.out
	diff /tmp/c.out /tmp/sage.out | head -n 10

KAT: meds.py
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


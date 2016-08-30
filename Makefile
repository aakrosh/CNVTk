all:
	cd src && $(MAKE)
	
install:
	mkdir -p bin
	cp src/count_in_bins bin/
	cp src/normalize_counts bin/
	cp src/segment_bins_using_CBS bin/
	cp src/genotype_segments bin/

.PHONY: clean

clean:
	-rm -rf bin
	cd src && $(MAKE) clean

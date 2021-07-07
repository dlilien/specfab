#
# Makefile
#
TOPTARGETS := all clean lib demo

SUBDIRS := specfab/src specfab/demo

$(TOPTARGETS): $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

.PHONY: $(TOPTARGETS) $(SUBDIRS) python

python:
	python setup.py install

clean:
	rm -rf dist wheelhouse build

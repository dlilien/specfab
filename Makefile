#
# Makefile
#
TOPTARGETS := all clean lib

SUBDIRS := specfab/src

$(TOPTARGETS): $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

.PHONY: $(TOPTARGETS) $(SUBDIRS)


# vim:ft=make
#

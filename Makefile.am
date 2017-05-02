## regenerate the ChangeLog on distribution
dist-hook: ChangeLog
ChangeLog: configure.ac
	git log --date-order --date=short | \
            awk '/^Author/ {sub(/\\$$/,""); getline t; print $$0 t; next}; 1' | \
            sed -e 's/^Author: //g' | \
            sed -e 's/>Date:   \([0-9]*-[0-9]*-[0-9]*\)/>\t\1/g' | \
            sed -e 's/^\(.*\) \(\)\t\(.*\)/\3    \1    \2/g' > ChangeLog && \
	cat ChangeLog.old >> ChangeLog 

## point to places where stuff can be found
if GUI
  MAYBE_GUI=src-gui
endif
SUBDIRS=dat src tools doc tests $(MAYBE_GUI)

EXTRA_DIST = AUTHORS COPYING INSTALL NEWS README.md THANKS ChangeLog

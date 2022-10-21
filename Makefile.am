## point to places where stuff can be found
SUBDIRS=dat src tools

EXTRA_DIST = AUTHORS COPYING INSTALL NEWS README.md THANKS ChangeLog

.PHONY: syntax
syntax:
	wget -qO - 'raw.githubusercontent.com/aoterodelaroza/aoterodelaroza.github.io/master/_critic2/12_syntax.md' | awk '/~~~/{a=!a;next}a{print}' > syntax.txt


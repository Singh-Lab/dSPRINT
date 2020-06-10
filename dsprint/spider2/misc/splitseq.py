#!/usr/bin/env python
import re, sys

if len(sys.argv) != 3:
	print("usage: RUN file ipos")
	exit(1)
ipos = int(argv[2])
fpo = ""
for x in open(argv[1]):
	if x.startswith('>'):
		ss = re.split(r'[>|:\s]+', x)
		if ss[0] == '': del ss[0]
		if fpo != "":
			fpo.close()
		elif argv[1] != "--":
			print 'Are you sure to use name "%s" from line:\n%s' %(ss[ipos], x.rstrip())
			if argv[1] != 'STDIN': raw_input()
		fn = ss[ipos]
#		fn = fn[:-1].lower() + fn[-1]
		fpo = file(fn + ".seq", "w")
		print >>fpo, x,
	elif fpo != "":
		print >>fpo, x,

#!/usr/bin/env python
from sys import *

dict_ASA0 = dict(zip("ACDEFGHIKLMNPQRSTVWY",
					(115, 135, 150, 190, 210, 75, 195, 175, 200, 170,
		185, 160, 145, 180, 225, 115, 140, 155, 255, 230)))
if __name__ == '__main__':
	if len(argv)<2:
		print >>stderr, 'usage: RUN *.spd3'
		exit(1)

	for x in open(argv[1]):
		if x[0] == '#': continue
		ss = x.split()
		rASA = float(ss[3]) / dict_ASA0[ss[1]] * 100.
		print ss[1], '%.2f' % rASA


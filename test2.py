#!/bin/python
import os



for x in xrange(1,11):
  d = {'sigma':x} 
  filename = "sigma" + str(x) + ".txt"
  d['filename'] = filename
  cmd = './main2.native abcd --window 10 --sigma %(sigma)f --iter 1000 --no-hbonds --prefix ../bpti_db/features/features_ --suffix .txt --files 5 > %(filename)s' % d
  print cmd
  os.system(cmd)


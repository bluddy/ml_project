#!/bin/python
import os



for x in xrange(5,21):
  d = {'window':x} 
  filename = "window" + str(x) + ".txt"
  d['filename'] = filename
  cmd = './main2.native abcd --window %(window)d --sigma 1 --iter 1000 --no-hbonds --prefix ../bpti_db/features/features_ --suffix .txt --files 5 > %(filename)s' % d
  print cmd
  os.system(cmd)


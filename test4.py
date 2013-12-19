#!/bin/python
import os



for x in xrange(1,11):
  d = {'file':x*50} 
  filename = "files" + str(d['file']) + ".txt"
  d['filename'] = filename
  cmd = './main2.native abcd --no-hbonds --window 10 --sigma 1 --iter 1000  --prefix ../bpti_db/features/features_ --suffix .txt --files %(file)s > %(filename)s' % d
  print cmd
  os.system(cmd)


#!/usr/bin/python
#Modified to accept file with multiple reads -BZ
import os
from os.path import join, abspath, isfile, isdir, exists, basename
import time

import sys
import fmindex
import tarfile

def diff_time(start, end):
    return int((end - start) * 1000)

def main():
    if not len(sys.argv) in [4]:
        print 'Usage: '
        print '  %s index reads_file output_file' % sys.argv[0]
        os.abort()
    else:
        if not isfile(sys.argv[1]):
            print "Index file doesn't exist"
            os.abort()
        tim = time.clock
        
        t_start = tim()
        
        idx = fmindex.load(sys.argv[1])
        t_load = tim()
        outFile=open(sys.argv[3], 'w')
        for line in open(sys.argv[2], 'r'):#expects a file with 1 read/line
            seq=line.strip()
            c = idx.count(seq)        
            m = idx.search(seq)
            outFile.write('%s\tcount:%d\tmatches:%s\n' % (seq, c, str(m)))
        t_search = time.clock()
        print "load: %sms" % diff_time(t_start, t_load)
        print "search: %sms" % diff_time(t_load, t_search)
if __name__ == '__main__':
    main()
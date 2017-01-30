# -*- coding: utf-8 -*-

import pickle
#import cPickle as pickle
import pickle
import bwt

# burrow wheeler transform
bw  = bwt.BurrowsWheeler()
# burrow wheeler inverse
bwi = bwt.BurrowsWheeler()

def save(filename, idx):
    f = open(filename, 'w')
    pickle.dump(idx,f)

def load(filename):
    f = open(filename)
    idx = pickle.load(f)
    return idx

def index(data):
    #return FMSimpleIndex(data)
    #return FMFullIndex(data)
    return FMCheckpointing(data)

class FMSimpleIndex(object):   
    def __init__(self, data):
        self.data = bw.transform(data)
        self.offset = {}
        self._build(data)
    
    def _build(self, data):
        """ build the index """
        self.occ = bwt.calc_first_occ(self.data)
    
    def _occ(self, qc):
        """ get the first occurance of letter qc in left-column"""
        c = self.occ.get(qc)
        if c == None:
            return 0
        return c
    
    def _count(self, idx, qc):
        """ count the occurances of letter qc (rank of qc) upto position idx """
        if not qc in self.occ.keys(): return 0
        c = 0
        for i in xrange(idx):
            if self.data[i] == qc:
                c += 1
        return c
    
    def _lf(self, idx, qc):
        #-------Your code should start here---------
        """ return the nearset lf mapping for letter qc at position idx """
        return self._occ(qc) + self._count(idx, qc)
        #--------Your code should end here----------
    
    def _walk(self, idx):
        """ find the offset in position idx of transformed string
            from the beginning """
        # walk to the beginning using lf mapping
        # this is same as inverse of burrow wheeler transformation
        # from arbitrary location
        r = 0
        i = idx 
        while self.data[i] != bw.EOS:
            if self.offset.get(i):
                # we have cached the location and can use it
                r += self.offset[i]
                break
            r += 1
            i = self._lf(i, self.data[i])
        
        # save the offset of some idx for faster searches
        if not self.offset.get(idx):
            self.offset[i] = r
        return r
    
    def bounds(self, q):
        """ find the first and last suffix positions for query q """
        """These are positions in the BWT string"""
        """This is the meat of the FM search algorithm"""
        top = 0
        bot = len(self.data)

        #-------Your code should start here---------
        #iterate over letters in query string q in reverse
        for qc in q[::-1]:
            # calculate top (which maps the position in the last column to the position in the first column)
            top = self._lf(top, qc) 
            # calculate bottom
            bot = self._lf(bot, qc)
            if top == bot: return (-1,-1)#since bottom is non-inclusive, top==bot implies that the string was not found.
        #--------Your code should end here----------
        return (top,bot)
    
    def search(self, q):
        """ search the positions of query q """
        
        # find the suffixes for the query
        top, bot = self.bounds(q)
        matches = []
        # find the location of the suffixes
        # by walking the reverse text from that position
        # with lf mapping
        #basically does an inverse BWT for each match to get the positions
        for i in range(top, bot):
            pos = self._walk(i)
            matches.append(pos)
        return sorted(matches)
    
    def count(self, q):
        """ count occurances of q in the index """
        """ does so by running the FM search algorithm"""
        top, bot = self.bounds(q)
        return bot - top
    
    def getOriginal(self):
        return bwi.inverse(self.data)
    
    def RLE(self):
        output = []
        last = ''
        k = 0
        for i in range(len(self.data)):
            ch = self.data[i]
            if ch == last:
                k += 1
            else:
                if k > 0:
                    output.append((last, k))
                last = ch
                k = 1
        output.append((last, k))
        return output
'''
class FMFullIndex(FMSimpleIndex):
    """ creates full LF index for each letter, space inefficient """
    
    def __init__(self, data):
        self.data = bw.transform(data)
        self.offset = {}
        self._build()
    
    def _build(self):
        """ build the index """
        occ = bwt.calc_first_occ(self.data)
        
        # FM Index
        FM = {}
        for i, c in enumerate(self.data):
            # we'll store the nearest LF mapping for each letter
            # space inefficient
            for x, v in occ.items():
                FM[(i,x)] = v
            occ[c] += 1
        i = len(self.data)
        for x, v in occ.items():
            FM[(i,x)] = v
        del occ
        
        self.FM = FM
    
    def _lf(self, idx, qc):
        return self.FM[(idx,qc)]
'''
class FMCheckpointing(FMSimpleIndex):
    """ creates LF index with checkpoints """
    
    def __init__(self, data, step = 50):
        self.data = bw.transform(data)#get BWT of input string
        self.offset = {}
        self.step = step
        self._build()
    
    def _build(self):
        """ build the index """
        #It seems that occ and C are backwards, as defined by the wikipedia page?
        self.occ = bwt.calc_first_occ(self.data)#returns C[c], a dictionary for each letter c, of the number of lexically smaller letters in the string
        self.C = bwt.calc_checkpoints(self.data, self.step) #returns the list of checkpoints, seperated by distance "step", each of which is a dict of the
                                                            #number of occurences of each character up to that point (noninclusive).
    
    def _count(self, idx, qc):
        """ count the occurances of letter qc (rank of qc) upto position idx """
        count = bwt.count_letter_with_checkpoints(self.C, self.step, self.data, idx, qc)
        return count
    
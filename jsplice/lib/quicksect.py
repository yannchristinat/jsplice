'''This file is part of jSplice.

    jSplice is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    jSplice is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with jSplice.  If not, see <http://www.gnu.org/licenses/>
    
    (C) copyright ETH Zurich; 2014 Institute for Molecular Health Sciences, Krek group
    Developer: Yann Christinat
    
    This file is an adaption of the quicksect.py script in the bx-python library. 
    No licensing or copyright information was found for this file.
    https://bitbucket.org/james_taylor/bx-python/raw/ebf9a4b352d3/lib/bx/intervals/operations/quicksect.py
'''


import math, random

class IntervalTree( object ):
    def __init__( self ):
        self.chroms = {}
        
    def insert( self, region):
        chrom = region.chrom
        strand = region.strand
        
        if chrom in self.chroms:
            if strand in self.chroms[chrom]:
                self.chroms[chrom][strand] = self.chroms[chrom][strand].insert( region )
            else:
                self.chroms[chrom][strand] = IntervalNode( region )
        else:
            self.chroms[chrom] = {strand: IntervalNode( region )}
            
    def overlaps( self, region):
        chrom = region.chrom
        strand = region.strand
        if chrom in self.chroms and strand in self.chroms[chrom]:
            return self.chroms[chrom][strand].overlaps( region )
        else:
            return None
        
    def contains(self, region):
        chrom = region.chrom
        strand = region.strand
        if chrom in self.chroms and strand in self.chroms[chrom]:
            return self.chroms[chrom][strand].contains( region )
        else:
            return False
            
#WARNING: assumes that the chromosome and strand have already been separated!
class IntervalNode( object ):
    def __init__( self, region ):
        # Python lacks the binomial distribution, so we convert a
        # uniform into a binomial because it naturally scales with
        # tree size.  Also, python's uniform is perfect since the
        # upper limit is not inclusive, which gives us undefined here.
        self.priority = math.ceil( (-1.0 / math.log(.5)) * math.log( -1.0 / (random.uniform(0,1) - 1)))
        self.start = region.start
        self.end = region.end
        self.maxend = self.end
        self.minend = self.end
        self.left = None
        self.right = None
        self.region = region
    def insert( self, region):
        root = self
        if region.start > self.start:
            # insert to right tree
            if self.right:
                self.right = self.right.insert( region )
            else:
                self.right = IntervalNode(region)
            # rebalance tree
            if self.priority < self.right.priority:
                root = self.rotateleft()
        else:
            # insert to left tree
            if self.left:
                self.left = self.left.insert( region )
            else:
                self.left = IntervalNode(region)
            # rebalance tree
            if self.priority < self.left.priority:
                root = self.rotateright()
        if root.right and root.left: 
            root.maxend = max( root.end, root.right.maxend, root.left.maxend )
            root.minend = min( root.end, root.right.minend, root.left.minend )
        elif root.right: 
            root.maxend = max( root.end, root.right.maxend )
            root.minend = min( root.end, root.right.minend )
        elif root.left:
            root.maxend = max( root.end, root.left.maxend )
            root.minend = min( root.end, root.left.minend )
        return root

    def rotateright( self ):
        root = self.left
        self.left = self.left.right
        root.right = self
        if self.right and self.left: 
            self.maxend = max(self.end, self.right.maxend, self.left.maxend)
            self.minend = min(self.end, self.right.minend, self.left.minend )
        elif self.right:
            self.maxend = max(self.end, self.right.maxend)
            self.minend = min(self.end, self.right.minend)
        elif self.left:
            self.maxend = max(self.end, self.left.maxend)
            self.minend = min(self.end, self.left.minend )
        return root
        
    def rotateleft( self ):
        root = self.right
        self.right = self.right.left
        root.left = self
        if self.right and self.left: 
            self.maxend = max(self.end, self.right.maxend, self.left.maxend)
            self.minend = min(self.end, self.right.minend, self.left.minend )
        elif self.right:
            self.maxend = max(self.end, self.right.maxend)
            self.minend = min(self.end, self.right.minend)
        elif self.left:
            self.maxend = max(self.end, self.left.maxend)
            self.minend = min(self.end, self.left.minend )
        return root
            
    def overlaps(self, region):
        rlist = list()
        if region.end <= self.end and region.start >= self.start: 
            rlist.append( self.region )
        if self.left and region.end <= self.left.maxend:
            rlist += self.left.overlaps( region )
        if self.right and region.start >= self.start:
            rlist += self.right.overlaps( region )
        
        return rlist
    
    def contains(self, region):
        if self.region==region:
            return True
        else:
            lfound = False
            if self.left and region.end <= self.left.maxend:
                lfound = self.left.contains( region )
            rfound = False
            if self.right and region.start >= self.start:
                rfound = self.right.contains( region )
            return rfound or lfound

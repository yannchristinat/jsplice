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
    
'''

#Experimental design

import os.path, math

# Describe a single sample or experiment
class Experiment:
    def __init__(self, expName, condName, bamFile, junctionFile, bedCountFile):
        self.experimentName = expName
        self.conditionName = condName
        self.bamFile = bamFile
        self.junctionFile = junctionFile
        self.bedCountFile = bedCountFile
        
    def __str__(self):
        return self.experimentName+'\t'+self.conditionName+'\t'+self.bamFile+'\t'+self.junctionFile+'\t'+str(self.bedCountFile)
        
# Abstract experimental design class. Provide a structure to store elements related to several experiment and conditions.
class AbstractDesign(dict):
    def addElement(self, expName, condition, element):
        if expName not in self:
            self[expName] = {condition: element}
        else:
            self[expName][condition] = element
            
    def getConditions(self):
        conds = list()
        for exp in self:
            for cond in self[exp]:
                if cond not in conds:
                    conds.append(cond)
        return conds
    
    def __len__(self):
        l=0
        for exp in self:
            l += len(self[exp])
        return l
    
    def equals(self,design, threshold):
        try:
            if len(self)!=len(design):
                return False
            for k in self:
                if len(self[k])!=len(design[k]):
                    return False
                for c in self[k]:
                    if math.fabs(self[k][c]-design[k][c])>threshold:
                        return False
        except:
            return False
        return True
        
# Experimental design class as used by jSplice
class ExperimentalDesign(AbstractDesign):
    NOBAM = 1
    JXNONLY = 2
    
    # Reads in the text file and instanciate the object
    def __init__(self, designFilename, flag):
        self.valid = True
        
        mainDir = os.path.dirname(designFilename)
                
        nobam = flag==ExperimentalDesign.NOBAM
        jxnonly = flag==ExperimentalDesign.JXNONLY
                 
        f = open(designFilename,'r')
        line = f.readline()
        lineNb=0
        
            
        while line:
            lineNb+=1
            if line.strip()=="" or line.strip()[0]=='#':
                line = f.readline() #skip comments and empty lines
                continue
            
            
            if line.strip()=='':
                line = f.readline()
                lineNb+=1
                continue
            toks = line.strip().split('\t')
            if len(toks)<3:
                print 'ERROR: Missing item on line '+str(lineNb)+'. Each line should have at least 3 tab-separated fields: exp_name, cond_name, junc_bed_file, [bam_file|covbed_file]'
            else:
                expName = toks[0]
                condName = toks[1]
                
                juncFile = os.path.expanduser(toks[2])
                if not os.path.isabs(juncFile): 
                    juncFile = os.path.join(mainDir,juncFile)
                juncFile = os.path.abspath(juncFile)
                
                if jxnonly:
                    bFile = None
                else:
                    bFile = os.path.expanduser(toks[3])
                    if not os.path.isabs(bFile):
                        bFile = os.path.join(mainDir,bFile)
                    bFile = os.path.abspath(bFile)
                
                if not jxnonly:
                    if not os.path.isfile(bFile):
                        self.valid = False
                        if nobam:
                            print 'ERROR: BED file "'+bFile+'" does not exists.'
                        else:
                            print 'ERROR: BAM file "'+bFile+'" does not exists.'
                if not os.path.isfile(juncFile):
                    self.valid = False
                    print 'ERROR: junction file "'+juncFile+'" does not exists.'
                
                if self.valid:
                    if jxnonly:
                        exp = Experiment(expName, condName, None, juncFile, None)
                    else: 
                        if nobam:
                            exp = Experiment(expName, condName, None, juncFile, bFile) #bamFile is actually a BED file
                        else:
                            exp = Experiment(expName, condName, bFile, juncFile, None)
                    self.addElement(expName, condName, exp)
                
                        
            line = f.readline()
            
        f.close()
            
        #Check that all experiments have the same two conditions
        conds = self.getConditions()
        if len(conds)!=2:
            print 'ERROR: Your design file should contain only two conditions!'
            self.cond0 = None
            self.cond1 = None
        else:
            self.cond0 = conds[0]
            self.cond1 = conds[1]
        
            for exp in self:
                if self.cond0 not in self[exp]:
                    print 'ERROR: Experiment '+exp+' is lacking condition '+self.cond0
                    self.valid = False
                if self.cond1 not in self[exp]:
                    print 'ERROR: Experiment '+exp+' is lacking condition '+self.cond1
                    self.valid = False
                
        if len(self)==0:
            print 'ERROR: Your design file is empty!'
            self.valid = False
        
        
    def getBEDCountFiles(self):
        files = []
        for exp in self:
            for cond in self[exp]:
                files.append(self[exp][cond].bedCountFile)
        
        return files
    
    def getJunctionFiles(self):
        files = []
        for exp in self:
            for cond in self[exp]:
                files.append(self[exp][cond].junctionFile)
        
        return files
    
    def getBAMFiles(self):
        files = []
        for exp in self:
            for cond in self[exp]:
                files.append(self[exp][cond].bamFile)
        
        return files
    
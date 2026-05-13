#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from rdkit import Chem
import sys,os

class PDBQtScoreParser:
    def __init__(self, score_column=-1, score_filter=None, name_filter="Name:"):
        self.score_data=dict()
        self.ligand_index=0
        self.score_column=score_column
        self.score_filter=score_filter
        self.name_filter=name_filter
    
    def loadFile(self, filename):
        infile=open(filename,"r")
        remark_start=False
        my_name=None
        my_score=None
        for l in infile:
            l=l.strip()
            if not len(l): continue
            if not l.startswith("REMARK"):
                if remark_start and my_score is not None:
                    if (my_name is None): my_name="ligand_"+str(self.ligand_index)
                    self.score_data[my_name]=my_score
                    self.ligand_index+=1
                remark_start=False
                my_score=None
                my_name=None
                continue
            else: remark_start=True
            l=l.replace("REMARK","").strip()
            
            is_name=True
            is_score=True
            if self.name_filter is not None and l.find(self.name_filter)==-1: is_name=False
            if self.score_filter is not None and l.find(self.score_filter)==-1: is_score=False
            if not (is_name or is_score):  continue
            
            if is_name:
                my_name=l.strip().split()[-1].strip()
                continue
            if my_score is None:
                l=l.strip().split()
                for k in l:
                    try:
                        score_num=float(k)
                        my_score=score_num
                        break
                    except: pass
    
    def getScoreDict(self): return self.score_data
    def dumpLigandsToFile(self,input_filename,outfile_obj,namelist,logfile=None):
        infile=open(input_filename,"r")
        remark_start=False
        my_name=None
        my_score=None
        lines=[]
        if logfile is not None: logfile.write("[INFO]: Searching for poses in "+input_filename+"\n")
        else: print("[INFO]: No logfile provided for dumpLigandsToFile()")
        for l in infile:
            l=l.strip()
            if not len(l): continue
            if not l.startswith("REMARK"):
                if remark_start and my_score is not None:
                    if (my_name is None): my_name="ligand_"+str(self.ligand_index)
                    self.ligand_index+=1
                lines.append(l)
                remark_start=False
                continue
            else:
                if not remark_start:
                    if my_name in namelist:
                        if logfile is not None: logfile.write("\tFound ligand: "+my_name+"\n")
                        for m in lines: outfile_obj.write(m+"\n")
                    my_score=None
                    my_name=None
                    lines=[]
                remark_start=True
            l=l.replace("REMARK","").strip()
            
            is_name=True
            is_score=True
            if self.name_filter is not None and l.find(self.name_filter)==-1: is_name=False
            if self.score_filter is not None and l.find(self.score_filter)==-1: is_score=False
            if is_name: my_name=l.strip().split()[-1].strip()
            
            lines.append("REMARK "+l)
        if len(lines) and my_name is not None:
            if my_name in namelist:
                if logfile is not None: logfile.write("\tFound ligand: "+my_name+"\n")
                for m in lines: outfile_obj.write(m+"\n")

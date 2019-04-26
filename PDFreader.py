#!/usr/bin/env python
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# PyWAD is open-source software and consists of a set of modules written in python for the WAD-Software medical physics quality control software. 
# The WAD Software can be found on https://github.com/wadqc
# 
# The pywad package includes modules for the automated analysis of QC images for various imaging modalities. 
# PyWAD has been originaly initiated by Dennis Dickerscheid (AZN), Arnold Schilham (UMCU), Rob van Rooij (UMCU) and Tim de Wit (AMC) 
#
#
# Changelog:
#   20180720: first complete version

from __future__ import print_function

__version__ = '201807020'
__author__ = 'jmgroen'

import os
# this will fail unless wad_qc is already installed
from wad_qc.module import pyWADinput
from wad_qc.modulelibs import wadwrapper_lib

import numpy as np
import scipy
if not 'MPLCONFIGDIR' in os.environ:
    import pkg_resources
    try:
        #only for matplotlib < 3 should we use the tmp work around, but it should be applied before importing matplotlib
        matplotlib_version = [int(v) for v in pkg_resources.get_distribution("matplotlib").version.split('.')]
        if matplotlib_version[0]<3:
            os.environ['MPLCONFIGDIR'] = "/tmp/.matplotlib" # if this folder already exists it must be accessible by the owner of WAD_Processor 
    except:
        os.environ['MPLCONFIGDIR'] = "/tmp/.matplotlib" # if this folder already exists it must be accessible by the owner of WAD_Processor 

        import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.

# we need pydicom to read out dicom tags
try:
    import pydicom as dicom
except ImportError:
    import dicom

def logTag():
    return "[PDFreader] "

def acqdatetime_series(data, results, action):
    """
    Read acqdatetime from dicomheaders and write to IQC database

    Workflow:
        1. Read only headers
    """
    try:
        import pydicom as dicom
    except ImportError:
        import dicom
    try:
        params = action['params']
    except KeyError:
        params = {}

    ## 1. read only headers
    dcmInfile = dicom.read_file(data.series_filelist[0][0], stop_before_pixels=True)

    dt = wadwrapper_lib.acqdatetime_series(dcmInfile)
    
    results.addDateTime('AcquisitionDateTime', dt)
     

def header_series(data, results, action):
    
    # function based on pyWAD function by A. Schilham
    
    # get the first (and only) file
    instances = data.getAllInstances()
    
    if len(instances) != 1:
        print('%s Error! Number of instances not equal to 1 (%d). Exit.'%(logTag(),len(instances)))
    instance=instances[0]
    

    
    # look in the config file for tags and write them as results, nested tags are supported 2 levels
    for key in action['tags']:
        varname=key
        tag=action['tags'][key]
        if tag.count('/')==0:
            value=instance[dicom.tag.Tag(tag.split(',')[0],tag.split(',')[1])].value
        elif tag.count('/')==1:
            tag1=tag.split('/')[0]
            tag2=tag.split('/')[1]
            value=instance[dicom.tag.Tag(tag1.split(',')[0],tag1.split(',')[1])][0]\
            [dicom.tag.Tag(tag2.split(',')[0],tag2.split(',')[1])].value
        elif tag.count('/')==2:
            tag1=tag.split('/')[0]
            tag2=tag.split('/')[1]
            tag3=tag.split('/')[2]
            value=instance[dicom.tag.Tag(tag1.split(',')[0],tag1.split(',')[1])][0]\
            [dicom.tag.Tag(tag2.split(',')[0],tag2.split(',')[1])][0]\
            [dicom.tag.Tag(tag3.split(',')[0],tag3.split(',')[1])].value
        else:
            # not more then 2 levels...
            value='too many levels'

        # write results
        results.addString(varname, str(value)[:min(len(str(value)),100)])    
    
        

        
def PDF(data, results, action):

    import pylab as plt 

    try:
        params = action['params']
    except KeyError:
        params = {}
    
    # assume that there is 1 file with multiple images
    instances = data.getAllInstances()
    instance=instances[0]
    
    # print(instance)
    
    # pixel_data=instance.pixel_array

    # check MIME Type Tag
    try:
        mimetype=instance[dicom.tag.Tag('0042','0012')].value
        print('MIME Type: ',mimetype)
    except:
        print('MIME Type Tag not present!')
        quit()
    
    if 'pdf' in mimetype:
        print('Encapsulated document should be a pdf.')
    else:
        print('Not a pdf, quitting.')
        quit()
    
    
    stream=instance[dicom.tag.Tag('0042','0011')].value
        
    stream=stream.decode("utf-8")
    # print(stream)
   
    # finale check: a pdf file should start with '%PDF'
    if stream[0:4]=='%PDF':
        print("it's a PDF!")
    else:
        print('Not a pdf, quitting.')
        quit()
    
    
    # import string
    # translator = str.maketrans('', '', string.punctuation)
    
    
    
    for key in action['texts']:
        pre_index=stream.find(key['pre'])
        if pre_index==-1:
            print(key['name'], ' - pre not found')
       
        newstream=stream[pre_index:len(stream)]
        
        post_index=newstream.find(key['post'])
        if post_index==-1:
            print(key['name'], ' - post not found')
 
        substream=stream[pre_index+len(key['pre']):post_index+pre_index]

        #remove all whites
        substream = ''.join(substream.split())
        # print(substream)
        
        if key['type']=='string':
            results.addString(key['name'], str(substream)[:min(len(str(substream)),100)]) 
        elif key['type']=='float':
            results.addFloat(key['name'], float(substream))
        else:
            print('other type?')
        
        # print(key['name'],' = ',substream)

            
if __name__ == "__main__":
    #import the pyWAD framework and get some objects
    data, results, config = pyWADinput()

    # look in the config for actions and run them
    for name,action in config['actions'].items():
        if name=='ignore':
            s='s'
        
        # save acquisition time and date as result        
        elif name == 'acqdatetime':
           acqdatetime_series(data, results, action)

        # save whatever tag is requested as result
        elif name == 'header_series':
           header_series(data, results, action)

        # run the PDF analysis
        elif name == 'pdf_series':
            PDF(data, results, action)

    results.write()

    # all done

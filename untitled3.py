# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 14:35:42 2015

@author: windows
"""

from sendNew2_registerdatabase import *
from newFeatures import *
import processDicoms

file_ids = open("viewStage1T2updatedFeatures.txt",'r')
line = file_ids.readline()
print line
fileline = line.split()
lesion_id = int(fileline[0] )
StudyID = fileline[1] 
PatientID = fileline[2]  
dateID = fileline[3]
Diagnosis = fileline[4:]

SendNew2DB = SendregisterNew()
[img_folder, cond, BenignNMaligNAnt, Diagnosis, casesFrame, MorNMcase, T2info] = SendNew2DB.querylocalDatabase(lesion_id)        
        
AccessionN = casesFrame['exam_a_number_txt']
DicomExamNumber = casesFrame['exam_img_dicom_txt']
## for old DicomExamNumber         
#AccessionN = DicomExamNumber
dateID = casesFrame['exam_dt_datetime']
if 'proc_proc_side_int' in casesFrame.keys():
    finding_side = casesFrame['proc_proc_side_int']
else:
    finding_side = casesFrame['exam_find_side_int ']

pathSegment = 'C:\Users\windows\Documents'+os.sep+'repoCode-local'+os.sep+'addnewLesion'+os.sep+'segmentations'
nameSegment = casesFrame['lesionfile'] 

# account for different filenames
if not '.vtk' in nameSegment:
    anameSegment = str(int(StudyID))+'_'+DicomExamNumber+'_'+nameSegment+'.vtk'
    nameSegment = anameSegment    
    
DynSeries_id = MorNMcase['DynSeries_id']
T2Series_id = MorNMcase['T2Series_id']


DicomDirectory = img_folder+os.sep+str(int(StudyID))+os.sep+AccessionN
if not os.path.exists(DicomDirectory):
    DicomDirectory = img_folder+os.sep+str(int(StudyID))+os.sep+DicomExamNumber
 
DynSeries_id = MorNMcase['DynSeries_id']
T2Series_id = MorNMcase['T2Series_id']

mhaDirectory = DicomDirectory
fixed_path =  DicomDirectory+os.sep+DynSeries_id+'_'+finding_side+'.mha'
moving_path =  DicomDirectory+os.sep+T2Series_id+'.mha'

path_rootFolder = os.path.dirname(os.path.abspath(__file__))
print path_rootFolder
#[SendNew2DB.T2_pos_pat, SendNew2DB.T2_pos_pat]  = SendNew2DB.DICOM2mha(path_rootFolder, img_folder, StudyID, AccessionN, DicomExamNumber, MorNMcase, finding_side)
[T2_pos_pat, T2_ori_pat] = processDicoms.get_T2_pos_ori(DicomDirectory, T2Series_id)  

SendNew2DB.T2_pos_pat=T2_pos_pat#[-22.39, -150.9, 115]
SendNew2DB.T2_ori_pat=T2_ori_pat#[-0, 1, 0, -0, -0, -1]
ouput_warped_image_elastix = mhaDirectory+os.sep+T2Series_id+'_warped_elastix.mha'

#############################                  
###### 2) Check segmentation accuracy with annotations
#############################       
print "\n Preload volumes and segmentation..."
[series_path, phases_series] = SendNew2DB.preloadSegment(img_folder, lesion_id, StudyID, AccessionN, DynSeries_id, pathSegment, nameSegment)

print "\n Visualize and load..."
SendNew2DB.loadSegment(pathSegment, nameSegment)

 
if(finding_side=='Right'):
    SendNew2DB.loadDisplay.visualizemha(fixed_path, ouput_warped_image_elastix, SendNew2DB.T2_pos_pat, SendNew2DB.T2_ori_pat, interact=True)

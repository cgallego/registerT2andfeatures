# -*- coding: utf-8 -*-

"""
This converst DICOM to mhd, registers T2 to DCE-MRI precontrast then
creates record in local database, segment and extract features

Created on Fri Jul 17 13:01:58 2015

@ author (C) Cristina Gallego, University of Toronto
"""

import os, os.path
import sys
from sys import argv, stderr, exit
import numpy as np
import dicom
import psycopg2
import pandas as pd

from sqlalchemy import Column, Integer, String
import datetime

from sqlalchemy.orm import sessionmaker
from sendNew2_registerdatabase import *
from newFeatures import *
   
  
if __name__ == '__main__':    
    # Get Root folder ( the directory of the script being run)
    path_rootFolder = os.path.dirname(os.path.abspath(__file__))
    print path_rootFolder
       
    # Open filename list
    print sys.argv[1]
    file_ids = open(sys.argv[1],"r")
    file_ids.seek(0)
    
    line = file_ids.readline()
    print line
       
    while ( line ) : 
        # Get the line: id  	id	CAD study #	Accession #	Date of Exam	Pathology Diagnosis 
        fileline = line.split()
        lesion_id = int(fileline[0] )
        StudyID = fileline[1] 
        AccessionN = fileline[2]  
        dateID = fileline[3]
        Diagnosis = fileline[4:]
  
        #############################
        print "\n Adding record radiology"
        SendNew2DB = SendregisterNew()
        radioinfo = SendNew2DB.queryRadioData(StudyID, dateID)
        radioinfo = radioinfo.iloc[0]
        ## send to new database    
        SendNew2DB.addRecordDB_radiology(lesion_id, radioinfo)
        
        #############################
        ###### 1) Querying Research database for clinical, pathology, radiology data
        #############################
        [img_folder, cond, BenignNMaligNAnt, Diagnosis, casesFrame, MorNMcase, T2info] = SendNew2DB.querylocalDatabase(lesion_id)        
        
        AccessionN = casesFrame['exam_a_number_txt']
        DicomExamNumber = casesFrame['exam_img_dicom_txt']
        ## for old DicomExamNumber         
        #AccessionN = DicomExamNumber
        dateID = casesFrame['exam_dt_datetime']
        finding_side = casesFrame['exam_find_side_int']
        if(finding_side=='L'):
            finding_side='Left'
        if(finding_side=='R'):
            finding_side='Right'

        pathSegment = 'C:\Users\windows\Documents'+os.sep+'repoCode-local'+os.sep+'addnewLesion'+os.sep+'segmentations'
        nameSegment = casesFrame['lesionfile'] 
        DynSeries_id = MorNMcase['DynSeries_id']
        T2Series_id = MorNMcase['T2Series_id']
        
        # account for different filenames
        if not '.vtk' in nameSegment:
            anameSegment = str(int(StudyID))+'_'+DicomExamNumber+'_'+nameSegment+'.vtk'
            nameSegment = anameSegment
            
        #############################
        ###### 2) Registration of T1-DCE pre contrast with T2
        #############################
        # Get right images in mha
        [mhaDirectory, fixed_path, moving_path]  = SendNew2DB.DICOM2mha(path_rootFolder, img_folder, StudyID, AccessionN, DicomExamNumber, MorNMcase, finding_side)
        print "\n fixed_path: %s" % fixed_path
        print " moving_path: %s" % moving_path
        
        ### do elastix registration
        regis_path = path_rootFolder+ os.sep +'registration'+ os.sep +'elastix'
        FNULL = open(os.devnull, 'w')    #use this if you want to suppress output to stdout from the subprocess
        elastix_exe = regis_path + os.sep +'elastix.exe'
        elastix_affine_pars = regis_path+ os.sep + 'elastix_pars_affine.txt'
        elastix_bspline_pars = regis_path + os.sep + 'elastix_pars_bspline.txt'
        ouput_warped_image_elastix = moving_path[:-4]+'_id'+str(lesion_id)+'_warped_elastix.mha'
        
        output_path = mhaDirectory+os.sep+'elastix_id'+str(lesion_id)
        TransformParameters_f = output_path+os.sep+'TransformParameters.0.txt'
        if not os.path.exists(output_path):
            os.mkdir(output_path)    
        
        #Running multiple registrations in succession: affine first, the output of affine registration as input to the Bspline:
        #Bspline registered result: result.1.mha
        #affine registered result: result.0.mha
#        print "\n Running registrations in succession: affine first..."
#        #elastix_cmd = elastix_exe + ' -f ' + fixed_path+ ' -m ' + moving_path +' -out ' + output_path + ' -p ' + elastix_affine_pars 
#        elastix_cmd_result = subprocess.call(elastix_cmd, stdout=FNULL, stderr=FNULL, shell=False)
        elastix_cmd_result=0
        if elastix_cmd_result==0: #if success
            print" \n Successfull registration!!"
            #rename result.1.mha, result.1.mha is the warped image name fixed in elastix executable
#            os.rename(output_path+os.sep+"result.0.mha",ouput_warped_image_elastix)
            TransformParameters = open(TransformParameters_f,"r")
            pars = TransformParameters.readline()
            # An affine transformation is defined as: Tμ(x) = A(x − c) + t + c, (taken from elastix manual page9)
            # where the matrix A has no restrictions. This means that the image can be translated, rotated, scaled,
            # and sheared. The parameter vector μ is formed by the matrix elements aij and the translation vector.
            # In 3D, this gives a vector of length 12: μ = (a11, a12, a13, a21, a22, a23, a31, a32, a33, tx, ty, tz)^T 
            while 'TransformParameters' not in TransformParameters.readline():
                print pars
                pars = TransformParameters.readline()
                print pars
            affine_pars = pars[pars.find('(')+1:pars.find(')')].split()[1:13]
        else:
            print" \n Unsuccessfull registration check for errors in %s" % output_path
                
        
        #############################                  
        ###### 3) Visualize and load, check with annotations
        #############################       
        print "\n Preload volumes and segmentation..."
        [series_path, phases_series, annotationsfound] = SendNew2DB.preloadSegment(img_folder, lesion_id, StudyID, AccessionN, DynSeries_id, pathSegment, nameSegment)
        path_T2Series = series_path+os.sep+T2Series_id 
        
        print "\n Visualize and load DCE lesion..."
        newnameSegment = str(int(StudyID))+'_'+DicomExamNumber+'_'+str(lesion_id)+'.vtk'
        SendNew2DB.loadSegment(pathSegment, newnameSegment)
        
        ############################# IF REPEATING SEGMENTATION
        #LesionZslice = SendNew2DB.loadDisplay.zImagePlaneWidget.GetSliceIndex()
        #SendNew2DB.segmentLesion(LesionZslice, newnameSegment)
        #############################
        
        print "\n Load registered T2 scan.."
        SendNew2DB.loadDisplay.visualizemha(fixed_path, SendNew2DB.T2_pos_pat, SendNew2DB.T2_ori_pat, interact=False)
        #SendNew2DB.T2_pos_pat[0] = 71 
            
        print "\n Visualize and load DCE lesion..."
        SendNew2DB.loadT2transSegment(pathSegment, newnameSegment, affine_pars)
        
        if annotationsfound:
           print "\n Check with annotations..."
           SendNew2DB.loadAnnotations(annotationsfound)
        
        #############################
        # 4) Extract Lesion and Muscle Major pectoralies signal                                   
        ############################# 
        line_muscleVOI = T2info['bounds_muscleSI']
        line_muscleVOI = line_muscleVOI.rstrip()
        l = line_muscleVOI[line_muscleVOI.find('[')+1:line_muscleVOI.find(']')].split(",")
        bounds_muscleSI = [float(l[0]), float(l[1]), float(l[2]), float(l[3]), float(l[4]), float(l[5]) ]
        print "\n bounds_muscleSI from file:"
        print bounds_muscleSI
        
        #############################
        # 5) Extract T2 features                            
        #############################
        [T2_muscleSI, muscle_scalar_range, bounds_muscleSI, 
         T2_lesionSI, lesion_scalar_range, LMSIR,
         morphoT2features, textureT2features] = SendNew2DB.T2_extract(T2Series_id, SendNew2DB.lesionT2trans, bounds_muscleSI)    
                 
        #############################
        ###### 6) Extract Dynamic features
        #############################
        [dyn_inside, dyn_contour] = SendNew2DB.extract_dyn(series_path, phases_series, SendNew2DB.lesion3D)
        
        #############################
        ###### 7) Extract Morphology features
        ############################# 
        morphofeatures = SendNew2DB.extract_morph(series_path, phases_series, SendNew2DB.lesion3D)
        
        #############################        
        ###### 8) Extract Texture features
        #############################
        texturefeatures = SendNew2DB.extract_text(series_path, phases_series, SendNew2DB.lesion3D)       
        
        #############################
        ###### 9) Extract new features from each DCE-T1 and from T2 using segmented lesion
        #############################
        newfeatures = newFeatures(SendNew2DB.load, SendNew2DB.loadDisplay)
        [deltaS, t_delta, centerijk] = newfeatures.extract_MRIsamp(series_path, phases_series, SendNew2DB.lesion3D, SendNew2DB.lesionT2trans, T2Series_id, SendNew2DB.T2_pos_pat, SendNew2DB.T2_ori_pat)
        
        # generate nodes from segmantation 
        [nnodes, curveT, earlySE, dce2SE, dce3SE, lateSE, ave_T2, prop] = newfeatures.generateNodesfromKmeans(deltaS['i0'], deltaS['j0'], deltaS['k0'], deltaS, centerijk, T2Series_id)    
        [kmeansnodes, d_euclideanNodes] = prop
        
        # pass nodes to lesion graph
        G = newfeatures.createGraph(nnodes, curveT, prop)                   

        [degreeC, closenessC, betweennessC, no_triangles, no_con_comp] = newfeatures.analyzeGraph(G)        
        network_measures = [degreeC, closenessC, betweennessC, no_triangles, no_con_comp]
                
        #############################
        ###### 10) End and Send record to DB
        #############################
        SendNew2DB.addRecordDB_lesion(newnameSegment, StudyID, DicomExamNumber, dateID, casesFrame, finding_side, MorNMcase, cond, Diagnosis, lesion_id, BenignNMaligNAnt,  DynSeries_id, T2Series_id)
                           
        SendNew2DB.addRecordDB_features(lesion_id, dyn_inside, dyn_contour, morphofeatures, texturefeatures)   
        
        if annotationsfound:
            SendNew2DB.addRecordDB_annot(lesion_id, SendNew2DB.annot_attrib, SendNew2DB.eu_dist_mkers, SendNew2DB.eu_dist_seg)

        SendNew2DB.addRecordDB_T2(lesion_id, T2Series_id, T2info, morphoT2features, textureT2features, T2_muscleSI, muscle_scalar_range, bounds_muscleSI, T2_lesionSI, lesion_scalar_range, LMSIR)
    
        #############################
        print "\n Adding record case to stage1"
        SendNew2DB.addRecordDB_stage1(lesion_id, d_euclideanNodes, earlySE, dce2SE, dce3SE, lateSE, ave_T2, network_measures)

        ## continue to next case
        line = file_ids.readline()
        print line
       
    file_ids.close()

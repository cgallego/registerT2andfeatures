# -*- coding: utf-8 -*-
"""
Process all pipeline for a record and
Send a record to database

Created on Fri Jul 25 15:46:19 2014

@ author (C) Cristina Gallego, University of Toronto
----------------------------------------------------------------------
 """
import os, os.path
import sys
import shutil
from sys import argv, stderr, exit
import numpy as np
import dicom
import psycopg2
import pandas as pd

import processDicoms
from inputs_init import *
from display import *
from features_dynamic import *
from features_morphology import *
from features_texture import *
from features_T2 import *
from segment import *
import pylab      
import annot
 
from sqlalchemy import Column, Integer, String
import datetime

from sqlalchemy.orm import sessionmaker
from add_newrecords import *

# to query biomatrix only needed
import database
from base import Base, engine
from query_database import *
from query_localdatabase import *

class SendregisterNew(object):
    """
    USAGE:
    =============
    Send2DB = SendregisterNew()
    """
    def __init__(self): 
        self.dataInfo = []
        self.queryData = Query() 
        self.load = Inputs_init()
        self.newrecords = AddNewRecords()
        self.queryBio = Query()
        self.querylocal = Querylocal()
        
        # Create only 1 display
        self.loadDisplay = Display()
        self.createSegment = Segment()
        self.loadDynamic = Dynamic()
        self.loadMorphology = Morphology()
        self.loadTexture = Texture()
        self.T2 = features_T2()  
        
        
    def queryRadioData(self, StudyID, dateID):
        """ Querying without known condition (e.g: mass, non-mass) if benign by assumption query only findings"""
        #############################
        ###### 1) Querying Research database for clinical, pathology, radiology data
        #############################
        print "Executing SQL connection..."

        try:
            ############# Query biomatrix            
            redateID = datetime.date(int(dateID[0:4]), int(dateID[5:7]), int(dateID[8:10]))
            radiologyinfo = self.queryBio.queryBiomatrix(StudyID, redateID)
                                    
        except Exception:
            print "Not able to query biomatrix"
            pass
            
        return radiologyinfo  
        
        
    def querylocalDatabase(self, lesion_id):
        """ Querying without known condition (e.g: mass, non-mass) if benign by assumption query only findings"""
        #############################
        ###### 1) Querying Research database for clinical, pathology, radiology data
        #############################
        print "Executing local connection..."
           
        # perform query
        try:
            ############# by lesion id
            l = self.querylocal.queryby_lesionid(lesion_id)
            self.casesFrame = pd.Series(l.__dict__)
            self.casesFrame.to_csv(sys.stdout)
            
            is_T2 = list(l.f_T2)
            T2 = pd.Series(is_T2[0])
            T2info =  pd.Series(T2.__dict__)

        except Exception:
            return -1
               
        # if correctly proccess
        #slice data, get only 1 record        
        is_mass = list(l.mass_lesion)
        if(is_mass):
            print "\n MASS:"
            cond = 'mass'
            mass = pd.Series(is_mass[0])
            mass_Case =  pd.Series(mass.__dict__)
            mass_Case.to_csv(sys.stdout) 
            # decide if it's a mass or nonmass
            MorNMcase = mass_Case
            
        is_nonmass = list(l.nonmass_lesion)
        if(is_nonmass):
            print "\n NON-MASS:"
            cond = 'nonmass'
            nonmass = pd.Series(is_nonmass[0])
            nonmass_Case =  pd.Series(nonmass.__dict__)
            nonmass_Case.to_csv(sys.stdout) 
            # decide if it's a mass or nonmass
            MorNMcase = nonmass_Case
            
        is_foci = list(l.foci_lesion)
        if(is_foci):
            print "\n FOCI:"
            cond = 'foci'
            foci = pd.Series(is_foci[0])
            foci_Case =  pd.Series(foci.__dict__)
            foci_Case.to_csv(sys.stdout) 
            # decide if it's a mass or nonmass
            MorNMcase = foci_Case
            
        print "\n----------------------------------------------------------"          
        img_folder = 'Z:\Breast\DICOMS'
        BenignNMaligNAnt = self.casesFrame['lesion_label'][-1]
        Diagnosis = self.casesFrame['lesion_diagnosis']
                        
        return img_folder, cond, BenignNMaligNAnt, Diagnosis, self.casesFrame, MorNMcase, T2info
        
        
    def DICOM2mha(self, path_rootFolder, img_folder, StudyID, AccessionN, DicomExamNumber, MorNMcase, finding_side):
        # perform conversion
        DicomDirectory = img_folder+os.sep+str(int(StudyID))+os.sep+AccessionN
        if not os.path.exists(DicomDirectory):
            DicomDirectory = img_folder+os.sep+str(int(StudyID))+os.sep+DicomExamNumber
         
        DynSeries_id = MorNMcase['DynSeries_id']
        T2Series_id = MorNMcase['T2Series_id']

        ### split bilateral series into left and right
        [self.Left_pos_pat, self.Left_ori_pat, self.Right_pos_pat, self.Right_ori_pat] = processDicoms.get_LorR_from_bilateral(DicomDirectory, DynSeries_id, T2Series_id)
        
        # get T2 position and orientation
        [self.T2_pos_pat, self.T2_ori_pat, self.T2fatsat] = processDicoms.get_T2_pos_ori(DicomDirectory, T2Series_id)  

        # convert DynSeries_id
        # Left
#        cmd = path_rootFolder+os.sep+'dcm23d'+os.sep+'bin'+os.sep+'dcm23d.exe -i '+DicomDirectory+os.sep+DynSeries_id+os.sep+'Left'+' -o '+DicomDirectory+os.sep+'temp'
#        print '\n---- Begin conversion of ' + DynSeries_id + 'Left to mha...' ;
#        p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
#        p1.wait()
#        
#        os.chdir(DicomDirectory+os.sep+'temp')
#        proc = subprocess.Popen('ls -a', stdout=subprocess.PIPE)
#        output = proc.stdout.read()
#        filetomove = output.split('\n')[2]
#        os.chdir(DicomDirectory)
#        proc = subprocess.Popen(['mv', 'temp'+os.sep+filetomove, DynSeries_id+'_Left.mha'], stdout=subprocess.PIPE)
#        proc.wait()
#        
#        # Right
#        cmd = path_rootFolder+os.sep+'dcm23d'+os.sep+'bin'+os.sep+'dcm23d.exe -i '+DicomDirectory+os.sep+DynSeries_id+os.sep+'Right'+' -o '+DicomDirectory+os.sep+'temp'
#        print '\n---- Begin conversion of ' + DynSeries_id + 'Right to mha...' ;
#        p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
#        p1.wait()
#        
#        os.chdir(DicomDirectory+os.sep+'temp')
#        proc = subprocess.Popen('ls -a', stdout=subprocess.PIPE)
#        output = proc.stdout.read()
#        filetomove = output.split('\n')[2]
#        os.chdir(DicomDirectory)
#        proc = subprocess.Popen(['mv', 'temp'+os.sep+filetomove, DynSeries_id+'_Right.mha'], stdout=subprocess.PIPE)
#        proc.wait()
#                
#        # convert T2Series_id
#        cmd = path_rootFolder+os.sep+'dcm23d'+os.sep+'bin'+os.sep+'dcm23d.exe -i '+DicomDirectory+os.sep+T2Series_id+' -o '+DicomDirectory+os.sep+'temp'
#        print '\n---- Begin conversion of ' + T2Series_id + ' to mha...' ;
#        p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
#        p1.wait()
#        
#        os.chdir(DicomDirectory+os.sep+'temp')
#        proc = subprocess.Popen('ls -a', stdout=subprocess.PIPE)
#        output = proc.stdout.read()
#        filetomove = output.split('\n')[2]
#        os.chdir(DicomDirectory)
#        proc = subprocess.Popen(['mv', 'temp'+os.sep+filetomove, T2Series_id+'.mha'], stdout=subprocess.PIPE)
#        proc.wait()
#
#        # remove temporal dirs, always when using  processDicoms.get_LorR_from_bilateral(abspath_PhaseID) create auxiliary Left and Right dirs
#        shutil.rmtree('temp')
        shutil.rmtree(DicomDirectory+os.sep+DynSeries_id+os.sep+'Left')
        shutil.rmtree(DicomDirectory+os.sep+DynSeries_id+os.sep+'Right')
        
        mhaDirectory = DicomDirectory
        moving_path =  DicomDirectory+os.sep+DynSeries_id+'_'+finding_side+'.mha'
        fixed_path =  DicomDirectory+os.sep+T2Series_id+'.mha'
                
        return mhaDirectory, fixed_path, moving_path
                        
        
    def preloadSegment(self, img_folder, lesion_id, StudyID, AccessionN, DynSeries_id, pathSegment, nameSegment):
        """ extract Series info and Annotations, then load"""
        AccessionN_loc = img_folder+os.sep+str(int(StudyID))+os.sep+AccessionN
        SeriesIDall = processDicoms.get_immediate_subdirectories(AccessionN_loc)
        
        #############################
        ###### Start by Loading 
        print "Start by loading volumes..."
        os.chdir(AccessionN_loc)
        if not os.path.exists('DynPhases'):
            os.mkdir('DynPhases') 
            
        #############################
        # Reveal annotations                                      
        #############################
        annotflag = False 
        annotationsfound = []
#        for iSer in SeriesIDall:
#            exam_loc = img_folder+os.sep+str(int(StudyID))+os.sep+AccessionN+os.sep+iSer
#            print "Path Series annotation inspection: %s" % iSer
#            os.chdir(img_folder+os.sep+str(int(StudyID))+os.sep+AccessionN+os.sep+'DynPhases')
#            annotationsfound, annotflag = annot.list_ann(exam_loc, annotflag, annotationsfound) 
#            
        
        [series_path, phases_series, self.lesionID_path] = self.load.readVolumes(img_folder, int(StudyID), AccessionN, DynSeries_id, lesion_id)
        print "Path to series location: %s" % series_path 
        print "List of pre and post contrast volume names: %s" % phases_series
        print "Path to new lesion segmentation: %s" % self.lesionID_path
        
        """ Load a previously existing segmentation"""
        print "\n Load Segmentation..."
        print "Path to lesion segmentation: %s" % pathSegment+os.sep+nameSegment
        self.lesion3D = self.load.loadSegmentation(pathSegment, nameSegment)
        print "Data Structure: %s" % self.lesion3D.GetClassName()
        print "Number of points: %d" % int(self.lesion3D.GetNumberOfPoints())
        print "Number of cells: %d" % int(self.lesion3D.GetNumberOfCells())
        
        return series_path, phases_series, annotationsfound
                
        
    def loadSegment(self, pathSegment, newnameSegment):
        #############################
        ###### Load Segmentation and visualize
        #############################           
        print "Path to lesion segmentation: %s" % pathSegment+os.sep+newnameSegment        
        self.loadDisplay.addSegment(self.lesion3D, (0,1,0), interact=False)                             
        self.createSegment.saveSegmentation(self.lesionID_path, self.lesion3D, lesionfilename=newnameSegment) 
        
        print "\n Visualize volumes..."
        self.loadDisplay.visualize(self.load.DICOMImages, self.load.image_pos_pat, self.load.image_ori_pat, sub=True, postS=1, interact=False)
                    
        # finally transform centroid world coords to ijk indexes
        im_pt = [0,0,0]
        ijk = [0,0,0]
        pco = [0,0,0]
        pixId_sliceloc = self.loadDisplay.transformed_image.FindPoint(self.loadDisplay.lesion_centroid)
        self.loadDisplay.transformed_image.GetPoint(pixId_sliceloc, im_pt) 
        io = self.loadDisplay.transformed_image.ComputeStructuredCoordinates( im_pt, ijk, pco)
        if io:
            self.lesion_centroid_ijk = ijk
            print "\n Lesion centroid"
            print self.lesion_centroid_ijk
        
        return              

    def segmentLesion(self, LesionZslice, newnameSegment):  
        """ Create a segmentation and check with annotations, if any"""
        
        print "\n Displaying picker for lesion segmentation"
        segorNot = 0 #int(raw_input('type 1 to segment: '))
        if segorNot == 1:
            seeds = self.loadDisplay.display_pick(self.load.DICOMImages, self.load.image_pos_pat, self.load.image_ori_pat, 4, LesionZslice)
            self.lesion3D = self.createSegment.segmentFromSeeds(self.load.DICOMImages, self.load.image_pos_pat, self.load.image_ori_pat, seeds, self.loadDisplay.iren1, self.loadDisplay.xImagePlaneWidget, self.loadDisplay.yImagePlaneWidget,  self.loadDisplay.zImagePlaneWidget)
            self.loadDisplay.addSegment(self.lesion3D, (0,1,1), interact=True)
            self.loadDisplay.picker.RemoveAllObservers()        

        axis_lengths = self.loadDisplay.extract_segment_dims(self.lesion3D)
        print axis_lengths
        self.eu_dist_seg = float( sqrt( axis_lengths[0] + axis_lengths[1])) # only measure x-y euclidian distance betweeen extreme points
        print "eu_dist_seg : " 
        print self.eu_dist_seg 
                    
        ###### loadSegmentation
        print "Data Structure: %s" % self.lesion3D.GetClassName()
        print "Number of points: %d" % int(self.lesion3D.GetNumberOfPoints())
        print "Number of cells: %d" % int(self.lesion3D.GetNumberOfCells())
        
        # finally transform centroid world coords to ijk indexes
        im_pt = [0,0,0]
        ijk = [0,0,0]
        pco = [0,0,0]
        pixId_sliceloc = self.loadDisplay.transformed_image.FindPoint(self.loadDisplay.lesion_centroid)
        self.loadDisplay.transformed_image.GetPoint(pixId_sliceloc, im_pt) 
        io = self.loadDisplay.transformed_image.ComputeStructuredCoordinates( im_pt, ijk, pco)
        if io:
            self.lesion_centroid_ijk = ijk
            print "\n Lesion centroid"
            print self.lesion_centroid_ijk
        
        #############################
        # 4) Parse annotations (display and pick corresponding to lesion)
        #############################                             
        self.createSegment.saveSegmentation(self.lesionID_path, self.lesion3D, lesionfilename=newnameSegment) 
        
        return 
        
        
    def loadT2transSegment(self, pathSegment, newnameSegment, ap):
        #############################
        ###### Load Segmentation affine_pars visualize
        ############################# 
        # Create matrix 4x4    
        affineM = vtk.vtkMatrix4x4()      
        affineM.DeepCopy((float(ap[0]), float(ap[1]), float(ap[2]), float(ap[9]),
                float(ap[3]), float(ap[4]), float(ap[5]), float(ap[10]),
                float(ap[6]), float(ap[7]), float(ap[8]), float(ap[11]),
                0, 0, 0, 1))
        
        # Set up the axes   
        #I_affineM = vtk.vtkMatrix4x4() 
        #affineM.Invert(affineM, I_affineM)
        
        affineT = vtk.vtkTransform()
        affineT.Concatenate(affineM)
        affineT.Update()      
        
        T2trans = vtk.vtkTransformPolyDataFilter()
        T2trans.SetInput(self.lesion3D)
        T2trans.SetTransform(affineT)
        T2trans.Update()
        
        # get lesionT2trans        
        self.lesionT2trans = T2trans.GetOutput()
        T2transnameSegment = newnameSegment[:-4]+'_T2trans.vtk'
    
        print "\n Visualize lesionT2trans based on affine pars of DCE-pre to T2 space..."
        self.loadDisplay.addSegment(self.lesionT2trans, (1,0,0), interact=False)                             
        transorNot = 1 #int(raw_input('type 1 to use trans T2 segmentation?: '))
        if transorNot == 1:
            self.createSegment.saveSegmentation(self.lesionID_path, self.lesionT2trans, lesionfilename=T2transnameSegment)
        else:
            self.lesionT2trans = self.lesion3D
        
        return 
        
    def loadAnnotations(self, annotationsfound):
        """ Explore annotations, pick one that overlays the segmentation if any"""      
        print annotationsfound
        
        if annotationsfound:
            print "\nLoading annotations..." 
            annots_dict_list = self.loadDisplay.extract_annot(annotationsfound)
            print "\nDisplay annotations:" 
            self.loadDisplay.display_annot(self.load.DICOMImages, self.load.image_pos_pat, self.load.image_ori_pat, annots_dict_list, interact=True)
        
        else:
            print "\n####################"
            print "No Annotations"
            print "####################"
                
        chooseAnnot = int(raw_input('\n Enter # corresponding to Lesion Annotation or 0 to skip: ') )
        if chooseAnnot != 0:
            #self.casesFrame['LesionAnnot']= str(annots_dict_list[chooseAnnot-1])
            self.annot_attrib = annots_dict_list[chooseAnnot-1]
            pi = self.annot_attrib['pi_2display']
            pf = self.annot_attrib['pf_2display'] 
            
            #############################
            ###### Compare manual marker distance with auto segmentation length for validation
            #############################
            self.eu_dist_mkers = float( sqrt( (pi[0]-pf[0])**2 + (pi[1]-pf[1])**2 + (pi[2]-pf[2])**2 ) )           
            print "eu_dist_mkers: " 
            print self.eu_dist_mkers
            
            axis_lengths = self.loadDisplay.extract_segment_dims(self.lesion3D)
            self.eu_dist_seg =  float(sqrt( axis_lengths[0] + axis_lengths[1]))  # only measure x-y euclidian distance betweeen extreme points
            print "eu_dist_seg : " 
            print self.eu_dist_seg 
       
        else:
            self.annot_attrib=[]
            self.eu_dist_mkers = []
            self.eu_dist_seg = []

        return 
        
    def extract_dyn(self, series_path, phases_series, lesion3D):            
        #############################
        ###### Extract Dynamic features
        #############################
        print "\n Extract Dynamic contour features..."
        dyn_contour = self.loadDynamic.extractfeatures_contour(self.load.DICOMImages, self.load.image_pos_pat, self.load.image_ori_pat, series_path, phases_series, lesion3D)
        print "\n=========================================="
        print dyn_contour
                
        print "\n Extract Dynamic inside features..."
        dyn_inside = self.loadDynamic.extractfeatures_inside(self.load.DICOMImages, self.load.image_pos_pat, self.load.image_ori_pat, series_path, phases_series, lesion3D)
        print dyn_inside
        print "\n=========================================="
 
        pylab.close('all') 
        
        return dyn_inside, dyn_contour
        
    def extract_morph(self, series_path, phases_series, lesion3D):      
        #############################
        ###### Extract Morphology features
        #############################
        print "\n Extract Morphology features..."
        morphofeatures = self.loadMorphology.extractfeatures(self.load.DICOMImages, self.load.image_pos_pat, self.load.image_ori_pat, series_path, phases_series, lesion3D)
        print "\n=========================================="
        print morphofeatures
        print "\n=========================================="

        pylab.close('all') 
        
        return morphofeatures
        
    def extract_text(self, series_path, phases_series, lesion3D):  
        #############################        
        ###### Extract Texture features
        #############################
        print "\n Extract Texture features..."
        texturefeatures = self.loadTexture.extractfeatures(self.load.DICOMImages, self.load.image_pos_pat, self.load.image_ori_pat, series_path, phases_series, lesion3D, self.loadMorphology.VOI_efect_diameter, self.loadMorphology.lesion_centroid )
        print "\n=========================================="
        print texturefeatures
        print "\n=========================================="

        pylab.close('all')  
        
        return texturefeatures
        
        
    def T2_extract(self, T2SeriesID, lesion3D, m_bounds):            
        #############################        
        ###### Extract T2 features, Process T2 and visualize
        #############################
        if T2SeriesID != 'NONE':  
            # Do extract_muscleSI 
            m_bounds[4] = m_bounds[4]+self.load.readjust_origin
            m_bounds[5] = m_bounds[5]+self.load.readjust_origin
            [T2_muscleSI, muscle_scalar_range, bounds_muscleSI]  = self.T2.load_muscleSI_mha(self.loadDisplay.warpT2_mha, self.T2_pos_pat, self.T2_ori_pat, m_bounds, self.loadDisplay.iren1)
            print "ave. T2_muscleSI: %d" % mean(T2_muscleSI)
            
            ############################# Do extract_muscleSI 
            #[T2_muscleSI, muscle_scalar_range, bounds_muscleSI]  = self.T2.extract_muscleS_mhaI(self.loadDisplay.warpT2_mha, self.T2_pos_pat, self.T2_ori_pat, self.loadDisplay.iren1, self.loadDisplay.renderer1, self.loadDisplay.picker, self.loadDisplay.xImagePlaneWidget, self.loadDisplay.yImagePlaneWidget, self.loadDisplay.zImagePlaneWidget)
            #print "ave. T2_muscleSI: %d" % mean(T2_muscleSI)
            #self.loadDisplay.iren1.Start()
            #############################
            
            # Do extract_lesionSI          
            [T2_lesionSI, lesion_scalar_range]  = self.T2.extract_lesionSI_mha(self.loadDisplay.warpT2_mha, lesion3D, self.T2_pos_pat, self.T2_ori_pat, self.loadDisplay)
            print "ave. T2_lesionSI: %d" % mean(T2_lesionSI)
            
            LMSIR = mean(T2_lesionSI)/mean(T2_muscleSI)
            print "LMSIR: %d" % LMSIR
                    
            #############################
            # Extract morphological and margin features from T2                                   
            #############################
            print "\n Extract T2 Morphology features..."
            morphoT2features = self.T2.extractT2morphology(self.loadDisplay.warpT2_mha, lesion3D, self.T2_pos_pat, self.T2_ori_pat)
            print "\n=========================================="
            print morphoT2features
            print "\n Extract T2 Texture features..."
            textureT2features = self.T2.extractT2texture(self.loadDisplay.warpT2_mha, lesion3D, self.T2_pos_pat, self.T2_ori_pat)
            print textureT2features
            print "\n=========================================="
            
            pylab.close('all')
        else:
            T2_muscleSI=[]; muscle_scalar_range=[]; bounds_muscleSI=[]; T2_lesionSI=[]; lesion_scalar_range=[]; LMSIR=[]; morphoT2features=[]; textureT2features=[];
        
        return T2_muscleSI, muscle_scalar_range, bounds_muscleSI, T2_lesionSI, lesion_scalar_range, LMSIR, morphoT2features, textureT2features
        
        
        
    def addRecordDB_lesion(self, Lesionfile, fStudyID, DicomExamNumber, dateID, casesFrame, finding_side, MorNMcase, cond, Diagnosis, 
                           lesion_id, BenignNMaligNAnt,  SeriesID, T2SeriesID):
                                       
        #############################
        ###### Send record to DB
        ## append collection of cases
        #############################  
        print "\n Adding record case to DB..."
        if 'proc_pt_procedure_id' in casesFrame.keys():
            self.newrecords.lesion_2DB(Lesionfile, fStudyID, casesFrame['anony_dob_datetime'], DicomExamNumber, str(casesFrame['exam_a_number_txt']), dateID, str(casesFrame['exam_mri_cad_status_txt']), 
                           str(casesFrame['cad_latest_mutation_status_int']), casesFrame['exam_find_mri_mass_yn'], casesFrame['exam_find_mri_nonmass_yn'], casesFrame['exam_find_mri_foci_yn'], finding_side, str(casesFrame['proc_pt_procedure_id']), 
                            casesFrame['proc_proc_dt_datetime'], str(casesFrame['proc_proc_side_int']), str(casesFrame['proc_proc_source_int']),  str(casesFrame['proc_proc_guid_int']), 
                            str(casesFrame['proc_proc_tp_int']), str(casesFrame['proc_lesion_comments_txt']), str(casesFrame['proc_original_report_txt']), str(casesFrame['find_curve_int']), 
                            str(casesFrame['find_mri_dce_init_enh_int']), str(casesFrame['find_mri_dce_delay_enh_int']), casesFrame['BIRADS'],
                            str(cond)+str(BenignNMaligNAnt),  Diagnosis)
        
        if not 'proc_pt_procedure_id' in casesFrame.keys():
            self.newrecords.lesion_2DB(Lesionfile, fStudyID, DicomExamNumber, str(casesFrame['exam_a_number_txt']), dateID, str(casesFrame['exam_mri_cad_status_txt']), 
                           str(casesFrame['cad_latest_mutation']), casesFrame['exam_find_mri_mass_yn'], casesFrame['exam_find_mri_nonmass_yn'], casesFrame['exam_find_mri_foci_yn'], finding_side, 'NA', 
                            datetime.date(9999, 12, 31), 'NA', 'NA', 'NA', 'NA', str(casesFrame['proc_lesion_comment_txt']), 'NA', str(casesFrame['find_curve_int']), str(casesFrame['find_mri_dce_init_enh_int']), str(casesFrame['find_mri_dce_delay_enh_int']), casesFrame['BIRADS'], cond+BenignNMaligNAnt,  Diagnosis)
                            
        if "mass" == cond:
            self.newrecords.mass_2DB(lesion_id, str(BenignNMaligNAnt), SeriesID, T2SeriesID, MorNMcase['find_mammo_n_mri_mass_shape_int'], MorNMcase['find_mri_mass_margin_int'] )

        if "nonmass" == cond: 
            self.newrecords.nonmass_2DB(lesion_id, str(BenignNMaligNAnt), SeriesID, T2SeriesID, MorNMcase['find_mri_nonmass_dist_int'], MorNMcase['find_mri_nonmass_int_enh_int'])
        
        if "foci" == cond: 
            self.newrecords.foci_2DB(lesion_id, str(BenignNMaligNAnt), SeriesID, T2SeriesID, MorNMcase['mri_foci_distr_int'])
       
       
        return
        
        
    def addRecordDB_features(self, lesion_id, dyn_inside, dyn_contour, morphofeatures, texturefeatures): 
        
        # SEgmentation details
        self.newrecords.segment_records_2DB(lesion_id, self.loadDisplay.lesion_bounds[0], self.loadDisplay.lesion_bounds[1], self.loadDisplay.lesion_bounds[2], self.loadDisplay.lesion_bounds[3], self.loadDisplay.lesion_bounds[4], self.loadDisplay.lesion_bounds[5],
                                                    self.loadDisplay.no_pts_segm, self.loadDisplay.VOI_vol, self.loadDisplay.VOI_surface, self.loadDisplay.VOI_efect_diameter, str(list(self.loadDisplay.lesion_centroid)), str(self.lesion_centroid_ijk))
                                                    
        # send features
        # Dynamic
        self.newrecords.dyn_records_2DB(lesion_id, dyn_inside['A.inside'], dyn_inside['alpha.inside'], dyn_inside['beta.inside'], dyn_inside['iAUC1.inside'], dyn_inside['Slope_ini.inside'], dyn_inside['Tpeak.inside'], dyn_inside['Kpeak.inside'], dyn_inside['SER.inside'], dyn_inside['maxCr.inside'], dyn_inside['peakCr.inside'], dyn_inside['UptakeRate.inside'], dyn_inside['washoutRate.inside'], dyn_inside['maxVr.inside'], dyn_inside['peakVr.inside'], dyn_inside['Vr_increasingRate.inside'], dyn_inside['Vr_decreasingRate.inside'], dyn_inside['Vr_post_1.inside'],
                               dyn_contour['A.contour'], dyn_contour['alpha.contour'], dyn_contour['beta.contour'], dyn_contour['iAUC1.contour'], dyn_contour['Slope_ini.contour'], dyn_contour['Tpeak.contour'], dyn_contour['Kpeak.contour'], dyn_contour['SER.contour'], dyn_contour['maxCr.contour'], dyn_contour['peakCr.contour'], dyn_contour['UptakeRate.contour'], dyn_contour['washoutRate.contour'], dyn_contour['maxVr.contour'], dyn_contour['peakVr.contour'], dyn_contour['Vr_increasingRate.contour'], dyn_contour['Vr_decreasingRate.contour'], dyn_contour['Vr_post_1.contour'] )
        
        # Morphology
        self.newrecords.morpho_records_2DB(lesion_id, morphofeatures['min_F_r_i'], morphofeatures['max_F_r_i'], morphofeatures['mean_F_r_i'], morphofeatures['var_F_r_i'], morphofeatures['skew_F_r_i'], morphofeatures['kurt_F_r_i'], morphofeatures['iMax_Variance_uptake'], 
                                                  morphofeatures['iiMin_change_Variance_uptake'], morphofeatures['iiiMax_Margin_Gradient'], morphofeatures['k_Max_Margin_Grad'], morphofeatures['ivVariance'], morphofeatures['circularity'], morphofeatures['irregularity'], morphofeatures['edge_sharp_mean'],
                                                  morphofeatures['edge_sharp_std'], morphofeatures['max_RGH_mean'], morphofeatures['max_RGH_mean_k'], morphofeatures['max_RGH_var'], morphofeatures['max_RGH_var_k'] )
        # Texture
        self.newrecords.texture_records_2DB(lesion_id, texturefeatures['texture_contrast_zero'], texturefeatures['texture_contrast_quarterRad'], texturefeatures['texture_contrast_halfRad'], texturefeatures['texture_contrast_threeQuaRad'], 
                                                  texturefeatures['texture_homogeneity_zero'], texturefeatures['texture_homogeneity_quarterRad'], texturefeatures['texture_homogeneity_halfRad'], texturefeatures['texture_homogeneity_threeQuaRad'], 
                                                  texturefeatures['texture_dissimilarity_zero'], texturefeatures['texture_dissimilarity_quarterRad'], texturefeatures['texture_dissimilarity_halfRad'], texturefeatures['texture_dissimilarity_threeQuaRad'], 
                                                  texturefeatures['texture_correlation_zero'], texturefeatures['texture_correlation_quarterRad'], texturefeatures['texture_correlation_halfRad'], texturefeatures['texture_correlation_threeQuaRad'], 
                                                  texturefeatures['texture_ASM_zero'], texturefeatures['texture_ASM_quarterRad'], texturefeatures['texture_ASM_halfRad'], texturefeatures['texture_ASM_threeQuaRad'], 
                                                  texturefeatures['texture_energy_zero'], texturefeatures['texture_energy_quarterRad'], texturefeatures['texture_energy_halfRad'], texturefeatures['texture_energy_threeQuaRad'] )
        return


    def addRecordDB_annot(self, lesion_id, annot_attrib, eu_dist_mkers, eu_dist_seg):
        # Send annotation if any
        if annot_attrib:
            self.newrecords.annot_records_2DB(lesion_id, annot_attrib['AccessionNumber'], annot_attrib['SeriesDate'], annot_attrib['SeriesNumber'], annot_attrib['SliceLocation'], annot_attrib['SeriesDescription'], annot_attrib['PatientID'], annot_attrib['StudyID'], annot_attrib['SeriesInstanceUID'], annot_attrib['note'], annot_attrib['xi'], annot_attrib['yi'], annot_attrib['xf'], annot_attrib['yf'], 
                                                    str(annot_attrib['pi_ijk']), str(annot_attrib['pi_2display']), str(annot_attrib['pf_ijk']), str(annot_attrib['pf_2display']),
                                                    eu_dist_mkers, eu_dist_seg)
        
        return
      
      
    def addRecordDB_T2(self, lesion_id, T2SeriesID, T2info, morphoT2features, textureT2features, T2_muscleSI, muscle_scalar_range, bounds_muscleSI, T2_lesionSI, lesion_scalar_range, LMSIR):
                                                              
        # T2 relative signal, morphology and texture
        if T2SeriesID != 'NONE':                                                       
            self.newrecords.t2_records_2DB(lesion_id, T2info['find_t2_signal_int'], str(list(self.loadDisplay.T2dims)), str(list(self.loadDisplay.T2spacing)), str(self.T2fatsat), mean(T2_muscleSI), std(T2_muscleSI), str(muscle_scalar_range), str(bounds_muscleSI), mean(T2_lesionSI), std(T2_lesionSI), str(lesion_scalar_range), LMSIR, 
                                            morphoT2features['T2min_F_r_i'], morphoT2features['T2max_F_r_i'], morphoT2features['T2mean_F_r_i'], morphoT2features['T2var_F_r_i'], morphoT2features['T2skew_F_r_i'], morphoT2features['T2kurt_F_r_i'], morphoT2features['T2grad_margin'], morphoT2features['T2grad_margin_var'], morphoT2features['T2RGH_mean'], morphoT2features['T2RGH_var'], 
                                            textureT2features['T2texture_contrast_zero'], textureT2features['T2texture_contrast_quarterRad'], textureT2features['T2texture_contrast_halfRad'], textureT2features['T2texture_contrast_threeQuaRad'], 
                                            textureT2features['T2texture_homogeneity_zero'], textureT2features['T2texture_homogeneity_quarterRad'], textureT2features['T2texture_homogeneity_halfRad'], textureT2features['T2texture_homogeneity_threeQuaRad'], 
                                            textureT2features['T2texture_dissimilarity_zero'], textureT2features['T2texture_dissimilarity_quarterRad'], textureT2features['T2texture_dissimilarity_halfRad'], textureT2features['T2texture_dissimilarity_threeQuaRad'], 
                                            textureT2features['T2texture_correlation_zero'], textureT2features['T2texture_correlation_quarterRad'], textureT2features['T2texture_correlation_halfRad'], textureT2features['T2texture_correlation_threeQuaRad'], 
                                            textureT2features['T2texture_ASM_zero'], textureT2features['T2texture_ASM_quarterRad'], textureT2features['T2texture_ASM_halfRad'], textureT2features['T2texture_ASM_threeQuaRad'], 
                                            textureT2features['T2texture_energy_zero'], textureT2features['T2texture_energy_quarterRad'], textureT2features['T2texture_energy_halfRad'], textureT2features['T2texture_energy_threeQuaRad'])
    
    def addRecordDB_stage1(self, lesion_id, d_euclidean, earlySE, dce2SE, dce3SE, lateSE, ave_T2, network_meas):        
        
        # Send to database lesion info
        self.newrecords.stage1_2DB(lesion_id, d_euclidean, earlySE, dce2SE, dce3SE, lateSE, ave_T2, network_meas)
              
        return


    
    def addRecordDB_radiology(self, lesion_id, radioinfo):        
        
        # Send to database lesion info
        self.newrecords.radiology_2DB(lesion_id, radioinfo['cad.cad_pt_no_txt'], radioinfo['cad.latest_mutation'], radioinfo['exam.exam_dt_datetime'],
                        radioinfo['exam.mri_cad_status_txt'], radioinfo['exam.comment_txt'], 
                        str(radioinfo['exam.original_report_txt']),
                        radioinfo['exam.sty_indicator_rout_screening_obsp_yn'], 
                        radioinfo['exam.sty_indicator_high_risk_yn'], radioinfo['exam.sty_indicator_high_risk_brca_1_yn'], radioinfo['exam.sty_indicator_high_risk_brca_2_yn'], radioinfo['exam.sty_indicator_high_risk_brca_1_or_2_yn'], 
                        radioinfo['exam.sty_indicator_high_risk_at_yn'], radioinfo['exam.sty_indicator_high_risk_other_gene_yn'],
                        radioinfo['exam.sty_indicator_high_risk_prior_high_risk_marker_yn'], radioinfo['exam.sty_indicator_high_risk_prior_personal_can_hist_yn'], radioinfo['exam.sty_indicator_high_risk_hist_of_mantle_rad_yn'],
                        radioinfo['exam.sty_indicator_high_risk_fam_hist_yn'], radioinfo['exam.sty_indicator_add_eval_as_folup_yn'], radioinfo['exam.sty_indicator_folup_after_pre_exam_yn'], 
                        radioinfo['exam.sty_indicator_pre_operative_extent_of_dis_yn'], radioinfo['exam.sty_indicator_post_operative_margin_yn'], radioinfo['exam.sty_indicator_pre_neoadj_trtmnt_yn'],
                        radioinfo['exam.sty_indicator_prob_solv_diff_img_yn'], radioinfo['exam.sty_indicator_scar_vs_recurr_yn'], radioinfo['exam.sty_indicator_folup_recommend_yn'], 
                        radioinfo['exam.sty_indicator_prior_2_prophy_mast_yn'])
              
        return
    
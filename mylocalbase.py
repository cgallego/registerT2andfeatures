# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 16:10:52 2014

@ author (C) Cristina Gallego, University of Toronto
"""

from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker


# configure Session class with desired options
Session = sessionmaker()
#myengine = create_engine('sqlite:///C:\\Users\\windows\\Documents\\repoCode-local\\stage1features\\stage1localData.db', echo=False) # now on, when adding new cases
newengine = create_engine('sqlite:///C:\\Users\\windows\\Documents\\repoCode-local\\registerT2andfeatures\\newregisterStage1T2updatedFeatures.db', echo=False) # now on, when adding new cases
localengine = create_engine('sqlite:///C:\\Users\\windows\\Documents\\repoCode-local\\addnewLesion\\stage1T2updatedFeatures.db', echo=False) # when querying local data
#myradengine = create_engine('sqlite:///C:\\Users\\windows\\Documents\\MyPubli and Abstracts\\Radiology2014-revision\\revision-revise\\addressing-revisions\\stage1localData_forRadiology-revisiondb.db',echo=False)


# later, we create the engine
Base = declarative_base(newengine)

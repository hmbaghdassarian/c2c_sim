#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import shutil
import sys
import io
import zipfile
import gzip
import requests
import subprocess
from multiprocessing import Pool
import warnings

from googleapiclient.http import MediaIoBaseDownload
from googleapiclient.discovery import build
from apiclient import discovery
from httplib2 import Http
import oauth2client
from oauth2client import file, client, tools

import pandas as pd

sys.path.insert(1, '../scripts/') # comment out in python script
from load_environmental_variables import *


# In[ ]:


print('Generate relevant data directories')
if not os.path.isdir(local_data_path):
    os.mkdir(local_data_path)
if not os.path.isdir(local_data_path + 'figures/'):
    os.mkdir(local_data_path + 'figures/')
if not os.path.isdir(local_data_path + 'raw/'):
    os.mkdir(local_data_path + 'raw/')
if not os.path.isdir(local_data_path + 'interim/'):
    os.mkdir(local_data_path + 'interim/')
if not os.path.isdir(local_data_path + 'processed/'):
    os.mkdir(local_data_path + 'processed/')


# # Google Drive

# In[ ]:


# print('Download google drive data')
# # this will only work if you have a client_id.json in root_path

# try:
#     obj = lambda: None
#     lmao = {"auth_host_name":'localhost', 'noauth_local_webserver':'store_true', 'auth_host_port':[8080, 8090], 'logging_level':'ERROR'}
#     for k, v in lmao.items():
#         setattr(obj, k, v)

#     # authorization boilerplate code
#     SCOPES = 'https://www.googleapis.com/auth/drive.readonly'
#     store = file.Storage(root_path + 'token.json')
#     creds = store.get()
#     # The following will give you a link if token.json does not exist, the link allows the user to give this app permission
#     if not creds or creds.invalid:
#         flow = client.flow_from_clientsecrets(root_path + 'client_id.json', SCOPES)
#         creds = tools.run_flow(flow, store, obj)
# except:
#     raise ValueError('Have you created a client_id.json?')
#     raise ValueError ('Have you created a .env in root directory with appropriate variables?')
#     raise ValueError('Current cell may first need to be run in jupyter notebook (see projects/notebooks/download_dataset.ipynb) to create a token.json')
    
    
# service = discovery.build('drive', 'v2', http=creds.authorize(Http()))
# # get ID and names of all files in a specified folder (by the folder_id in .env file)
# children = service.children().list(folderId=folder_id).execute()
# children_ids = [child['id'] for child in children.get('items', [])]
# file_names = [service.files().get(fileId=file_id).execute()['title'] for file_id in children_ids]

# # parse
# files = list(zip(children_ids, file_names))
# # files = [file for file in files if 'fastq.gz' in file[1]] # COMMENT OUT IF REUSING IN OTHER PROJECTS

# # download files from drive

# DRIVE = discovery.build('drive', 'v3', http=creds.authorize(Http()))
# count = 1
# for file in files:
#     file_id = file[0]
#     file_name = file[1]
#     print(str(count) + ' of ' + str(len(files)) + ' files: ' + file_name)
#     # if you get the shareable link, the link contains this id, replace the file_id below
#     request = DRIVE.files().get_media(fileId=file_id)

#     # replace the filename and extension in the first field below
#     fh = io.FileIO(file_name, mode='w')
#     downloader = MediaIoBaseDownload(fh, request)
#     done = False
#     while done is False:
#         status, done = downloader.next_chunk()
#         print("Download %d%%." % int(status.progress() * 100))

#     # move files to raw data path
#     str_ = 'mv ' + os.getcwd() + '/' + file_name + ' ' + local_data_path + 'raw/' + file_name
#     os.system(str_)
    
#     # untar
#     if file_name[-4:] == '.tar':
#         os.mkdir(local_data_path + 'raw/' + file_name[:-4])
#         untar = 'tar xf ' + local_data_path + 'raw/' + file_name 
#         untar += ' -C ' + local_data_path + 'raw/' + file_name[:-4]
#         os.system(untar)
    
#     count += 1


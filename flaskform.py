#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 12:20:49 2017

@author: jerrywzy
"""
import os
import Selenzy
from flask import Flask, flash, render_template, request, redirect, url_for, send_from_directory
from werkzeug import secure_filename
import pandas as pd
import numpy as np
import argparse
import uuid

app = Flask(__name__)



def arguments():
    parser = argparse.ArgumentParser(description='Options for the webserver')
    parser.add_argument('upload_folder', 
                        help='Upload folder')
    parser.add_argument('p2env', 
                        help='Specify path to python 2 environment directory')
    parser.add_argument('datadir',
                        help='specify data directory for required databases files, please end with slash')
    arg = parser.parse_args()
    return arg


def allowed_file(filename):
    return filename 

def file_path(uniqueid, filename):
    uniquefolder = os.path.join(app.config['UPLOAD_FOLDER'], uniqueid)
    uniquename = os.path.join(uniquefolder, filename)
    return uniquename


# @ signifies a decorator - wraps a function and modifies its behaviour
@app.route('/')
def upload_form():
    return render_template("my_form.html")
    

@app.route('/uploader', methods=['GET', 'POST'])
def upload_file():
    if request.method == 'POST':
        # check if post request has file part
        if 'file' not in request.files:
            flash('No file part')
            return redirect(request.url)
        file = request.files['file']   
        # if user does not select file, browser also submit empty part without fielname
        if file.filename == '':
            flash("No file selected")
            return redirect (request.url)
        if file and allowed_file(file.filename):
            direction = 0
            targets = request.form['targets']
            if request.form.get('direction'):
                direction = 1
            filename = secure_filename(file.filename)
            uniqueid = str(uuid.uuid4())
            uniquename = file_path(uniqueid, filename)
            uniquefolder = os.path.dirname(uniquename)
            if not os.path.exists(uniquefolder):
                os.mkdir(uniquefolder)
            file.save(uniquename)
            return redirect(url_for('analysed_file', session=uniqueid, filename=filename, targets=targets, direction=direction))
        
    return upload_form
    
@app.route('/results/files/<filename>')
def results_file(filename):
    return send_from_directory(app.config['UPLOAD_FOLDER'], filename)

@app.route('/results/<session>/<filename>?<targets>?<direction>')    
def analysed_file(session, filename, targets, direction):  
    filename = file_path(session, filename)
    filenameshort = os.path.splitext(filename)[0]
    realfile = (''.join(list(filter(str.isdigit, filenameshort))))
    csvfile = os.path.basename(filenameshort)+".csv"
    Selenzy.analyse(filename, 
                    app.config['PYTHON2'],
                    targets,
                    app.config['DATA_FOLDER'],  
                    os.path.dirname(filename),
                    csvfile,
                    int(direction)) # this creates CSV file in Uploads directory
    data = pd.read_csv(file_path(session, csvfile))
    data.index = data.index + 1
    return render_template('results.html', tables=data.to_html(), query=realfile, csvfile=csvfile)

    
if __name__== "__main__":  #only run server if file is called directly

    arg = arguments()
    ALLOWED_EXTENSIONS = set(['rxn', ' '])

    app.config['UPLOAD_FOLDER'] = os.path.abspath(arg.upload_folder)
    app.config['DATA_FOLDER'] = os.path.abspath(arg.datadir)
    app.config['PYTHON2'] = os.path.abspath(arg.p2env)

    app.run(debug=True)

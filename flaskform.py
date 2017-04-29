#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 12:20:49 2017

@author: jerrywzy
"""
import os
import SeqFind
from flask import Flask, flash, render_template, request, redirect, url_for, send_from_directory
from werkzeug import secure_filename
import pandas as pd
import numpy as np

UPLOAD_FOLDER = '/home/jerrywzy/Python Scripts/uploads'
ALLOWED_EXTENSIONS = set(['rxn', ' '])

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

def allowed_file(filename):
    return filename 


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
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            return redirect(url_for('analysed_file', filename=filename, targets=targets, direction=direction))
        
    return upload_form
    
@app.route('/results/files/<filename>')
def results_file(filename):
    return send_from_directory(app.config['UPLOAD_FOLDER'], filename)

@app.route('/results/<filename>?<targets>?<direction>')    
def analysed_file(filename, targets, direction):  
#    file = os.path.join(app.config['UPLOAD_FOLDER'], filename)
    filenameshort = os.path.splitext(filename)[0]
    realfile = (''.join(list(filter(str.isdigit, filenameshort))))
    if "." in filename:
        filenamepure = filenameshort.rsplit('/', 1)[-1]
        csvfile = "results_"+filenamepure+".csv"
    else:
        csvfile = "results_"+filename+".csv"
    SeqFind.analyse("uploads/"+filename, targets, direction) # this creates CSV file in uploads folder
    data = pd.read_csv(os.path.join(app.config['UPLOAD_FOLDER'], csvfile))
    data.index = data.index + 1
    return render_template('results.html', tables=data.to_html(), query=realfile, csvfile=csvfile)

    
if __name__== "__main__":  #only run server if file is called directly
    app.run(debug=True)

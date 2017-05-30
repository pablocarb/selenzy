#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 12:20:49 2017

@author: jerrywzy
"""
import os
import Selenzy
from flask import Flask, flash, render_template, request, redirect, url_for, send_from_directory, jsonify
from flask_restful import Resource, Api
from werkzeug import secure_filename
import pandas as pd
import numpy as np
import argparse
import uuid
import json

app = Flask(__name__)
api = Api(app)
app.config['SECRET_KEY'] = str(uuid.uuid4())


class RestGate(Resource):
    """ REST interface, returns api info """
    def get(self):
        return {'app': 'Selenzy', 'version': '1.0', 'author': 'Synbiochem'}

class RestQuery(Resource):
    """ REST interface to Selenzy, by default it does not run the MSA to be faster """
    def post(self):
        args = request.json
        if 'smarts' in args:
            rxntype = 'smarts'
            rxninfo = args['smarts']
            if 'targets' in args:
                targets = args['targets']
            else:
                targets = '20'
            if 'direction' in args:
                direction = int(args['direction'])
            else:
                direction = 0
            if 'noMSA' in args:
                noMSA = args['noMSA']
            else:
                noMSA = True
            data, csvfile, session = run_session(rxntype, rxninfo, targets, direction, noMSA)
            return jsonify({'app': 'Selenzy', 'version': '1.0', 'author': 'Synbiochem', 'data': data.to_json()})


api.add_resource(RestGate, '/REST')

api.add_resource(RestQuery, '/REST/Query')

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
    

def run_session(rxntype, rxninfo, targets, direction, noMSA):
    uniqueid = str(uuid.uuid4())
    uniquefolder = os.path.join(app.config['UPLOAD_FOLDER'], uniqueid)
    if not os.path.exists(uniquefolder):
        os.mkdir(uniquefolder)
    if rxntype == 'rxn':
        filename = secure_filename(rxninfo.filename)
        uniquename = file_path(uniqueid, filename)
        rxninfo.save(uniquename)
        rxninfo = uniquename
    csvfile = "selenzy_results.csv"
    Selenzy.analyse(['-'+rxntype, rxninfo], 
                    app.config['PYTHON2'],
                    targets,
                    app.config['DATA_FOLDER'],  
                    uniquefolder,
                    csvfile,
                    pdir = int(direction),
                    NoMSA = noMSA
                ) # this creates CSV file in Uploads directory
    data = pd.read_csv(file_path(uniqueid, csvfile))
    data.index = data.index + 1
    return data, csvfile, uniqueid

@app.route('/uploader', methods=['GET', 'POST'])
def upload_file():
    if request.method == 'POST':

        # check if post request has smarts part
        smarts = request.form['smarts']
        if len(smarts) > 0:
            rxntype = 'smarts'
            rxninfo = smarts
        # check if post request has file part
        elif 'file' in request.files:
            rxntype = 'rxn'
            file = request.files['file']   
            if file.filename == '' or not allowed_file(file.filename):
                flash("No file selected")
                return redirect (request.url)
            rxninfo = file
        else:
            flash('No reaction')
            return redirect(request.url)
        direction = 0
        noMSA = False
        targets = request.form['targets']
        if request.form.get('direction'):
            direction = 1
        if request.form.get('noMSA'):
            noMSA = True
        data, csvfile, sessionid = run_session(rxntype, rxninfo, targets, direction, noMSA)
        return render_template('results.html', tables=data.to_html(), csvfile=csvfile, sessionid=sessionid)

    return render_template("my_form.html")
    
@app.route('/results/<session>/files/<filename>')
def results_file(session,filename):
    return send_from_directory(os.path.join(app.config['UPLOAD_FOLDER'], session), filename)
    
if __name__== "__main__":  #only run server if file is called directly

    arg = arguments()
    ALLOWED_EXTENSIONS = set(['rxn', ' '])

    app.config['UPLOAD_FOLDER'] = os.path.abspath(arg.upload_folder)
    app.config['DATA_FOLDER'] = os.path.abspath(arg.datadir)
    app.config['PYTHON2'] = os.path.abspath(arg.p2env)

    app.run(host="0.0.0.0",port=5000, debug=True)
#    app.run(port=5000, debug=True)

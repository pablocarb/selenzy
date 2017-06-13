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
from flask import session
from werkzeug import secure_filename
import pandas as pd
import numpy as np
import argparse
import uuid
import json
import csv



app = Flask(__name__)
api = Api(app)
app.config['SECRET_KEY'] = str(uuid.uuid4())

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

def save_rxn(rxninfo):
    filename = secure_filename(rxninfo.filename)
    uniquename = file_path(session['uniqueid'], filename)
    rxninfo.save(uniquename)
    session['rxninfo'] = uniquename
    return uniquename

def init_session():
    reset_session()
    uniqueid = session['uniqueid']
    uniquefolder = os.path.join(app.config['UPLOAD_FOLDER'], uniqueid)
    if not os.path.exists(uniquefolder):
        os.mkdir(uniquefolder)
    session['uniquefolder'] = uniquefolder
    session['rxnifo'] = None
    session['status'] = False

def reset_session():
    uniqueid = str(uuid.uuid4())
    session['uniqueid'] = uniqueid

def run_session(rxntype, rxninfo, targets, direction, noMSA):
    uniqueid = session['uniqueid']
    uniquefolder = session['uniquefolder']
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

    
def retrieve_session(csvinfo):
    uniqueid = session['uniqueid']
    uniquefolder = os.path.join(app.config['UPLOAD_FOLDER'], uniqueid)
    if not os.path.exists(uniquefolder):
        os.mkdir(uniquefolder)
    filename = secure_filename(csvinfo.filename)
    uniquename = file_path(uniqueid, filename)
    csvinfo.save(uniquename)
    data = pd.read_csv(uniquename)
    data.index = data.index + 1
    csvfile = os.path.basename(uniquename)
    return data, csvfile, uniqueid



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
            init_session()
            data, csvfile, session = run_session(rxntype, rxninfo, targets, direction, noMSA)
            return jsonify({'app': 'Selenzy', 'version': '1.0', 'author': 'Synbiochem', 'data': data.to_json()})


api.add_resource(RestGate, '/REST')

api.add_resource(RestQuery, '/REST/Query')


# @ signifies a decorator - wraps a function and modifies its behaviour
@app.route('/')
def upload_form():
    if 'username' not in session:
        return redirect(url_for('login'))
    return render_template("my_form.html", username=session['username'])


@app.route('/login', methods=['GET', 'POST'])
def login():
    if app.debug == True:
      session['username'] = 'debug'  
      init_session()
      return redirect(url_for('upload_form'))
    if request.method == 'POST':
        session['username'] = request.form['username']
        init_session()
        return redirect(url_for('upload_form'))
    return '''
    <form method="post">
    <p><input type=text name=username>
    <p><input type=submit value=Login>
    </form>
    '''

@app.route('/logout')
def logout():
    session.pop('username', None)
    return redirect(url_for('login'))
        

@app.route('/display', methods=['POST'])
def display_reaction():
    """ Display the reaction """
    if request.method == 'POST':
        if 'file' in request.files and len(request.files['file'].filename) > 0:
            rxntype = 'rxn'
            fileinfo = request.files['file']   
            if fileinfo.filename == '' or not allowed_file(fileinfo.filename):
                flash("No file selected")
                return redirect (request.url)
            rxninfo = save_rxn(fileinfo)
            svgstream = Selenzy.display_reaction(rxninfo)
            svg = svgstream.decode('utf-8')
            success = False
            if len(svg) > 0:
                session['rxninfo'] = rxninfo
                session['rxntype'] = 'rxn'
                session['status'] = True
                success = True
                return json.dumps( {'data': svg, 'status': session['status'], 'success': success} )
        elif len(request.form['smarts']) > 0:
            rxninfo = request.form['smarts']
            svgstream = Selenzy.display_reaction(rxninfo)
            svg = svgstream.decode('utf-8')
            success = False
            if len(svg) > 0:
                session['rxninfo'] = rxninfo
                session['rxntype'] = 'smarts'
                session['status'] = True
                success = True
                return json.dumps( {'data': svg, 'status': session['status'], 'success': success} )


@app.route('/sorter', methods=['POST'])
def sort_table():
    """ Sorts table """
    if request.method == 'POST':
        jfilter = json.loads(request.values.get('filter'))
        print(jfilter)
        try:
            filt = [int(x) for x in jfilter]
        except:
            return
        session = json.loads(request.values.get('session'))
        csvname = os.path.basename(json.loads(request.values.get('csv')))
        csvfile = os.path.join(app.config['UPLOAD_FOLDER'], session, csvname)
        head, rows = Selenzy.read_csv(csvfile)
        sortrows = Selenzy.sort_rows(rows, filt)
        Selenzy.write_csv(csvfile, head, sortrows)
        data = pd.read_csv(csvfile)
        data.index = data.index + 1
        return json.dumps( {'data': {'csv':  data.to_html()}} )
                          


@app.route('/uploader', methods=['GET', 'POST'])
def upload_file():
    if request.method == 'POST':
        # check if post request has smarts part
        if 'csv' in request.files  and len(request.files['csv'].filename) > 0:
            fileinfo = request.files['csv']   
            if fileinfo.filename == '' or not allowed_file(fileinfo.filename):
                flash("No file selected")
                return redirect (request.url)
            data, csvfile, sessionid = retrieve_session(fileinfo)
            return render_template('results.html', tables=data.to_html(), csvfile=csvfile, sessionid=sessionid)
        else:
            try:
                rxninfo = session['rxninfo']
                rxntype = session['rxntype']
            except:
                return redirect(url_for('login'))

        direction = 0
        noMSA = False
        targets = request.form['targets']
        if request.form.get('direction'):
            direction = 1
        if request.form.get('noMSA'):
            noMSA = True
        
        data, csvfile, sessionid = run_session(rxntype, rxninfo, targets, direction, noMSA)
        return render_template('results.html', tables=data.to_html(), csvfile=csvfile, sessionid=sessionid)
    elif request.method == 'GET':
        smarts = request.args.get('smarts')
        if smarts is None:
            return redirect(request.url)
        rxntype = 'smarts'
        rxninfo = smarts
        direction = 0
        noMSA = False
        targets = 20
        init_session()
        session['rxninfo'] = rxninfo
        session['rxntype'] = rxntype
        data, csvfile, sessionid = run_session(rxntype, rxninfo, targets, direction, noMSA)
        return render_template('results.html', tables=data.to_html(), csvfile=csvfile, sessionid=sessionid)
    return render_template("my_form.html")
    
@app.route('/results/<sessionid>/files/<filename>')
def results_file(sessionid,filename):
    return send_from_directory(os.path.join(app.config['UPLOAD_FOLDER'], sessionid), filename)
    
if __name__== "__main__":  #only run server if file is called directly

    arg = arguments()
    ALLOWED_EXTENSIONS = set(['rxn', ' '])

    app.config['UPLOAD_FOLDER'] = os.path.abspath(arg.upload_folder)
    app.config['DATA_FOLDER'] = os.path.abspath(arg.datadir)
    app.config['PYTHON2'] = os.path.abspath(arg.p2env)

    app.run(host="0.0.0.0",port=5000, debug=True)
#    app.run(port=5000, debug=True)

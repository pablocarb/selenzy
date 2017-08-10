#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 12:20:49 2017

@author: jerrywzy
"""
import os, subprocess, glob, time, shutil
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
app.config['MARVIN'] = False

def arguments():
    parser = argparse.ArgumentParser(description='Options for the webserver')
    parser.add_argument('upload_folder', 
                        help='Upload folder')
    parser.add_argument('datadir',
                        help='Data directory for required databases files')
    parser.add_argument('-d', action='store_true',
                        help='Run in debug mode (no preload)')
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
    try:
        uniquename = file_path(session['uniqueid'], filename)
    except:
        init_session()
        uniquename = file_path(session['uniqueid'], filename)
    rxninfo.save(uniquename)
    outname = file_path(session['uniqueid'], session['uniqueid'])
    rxninfo = Selenzy.sanitizeRxn(uniquename, outname)
    session['rxninfo'] = rxninfo
    return rxninfo

def init_session():
    maintenance()
    reset_session()
    uniqueid = session['uniqueid']
    uniquefolder = os.path.join(app.config['UPLOAD_FOLDER'], uniqueid)
    if not os.path.exists(uniquefolder):
        os.mkdir(uniquefolder)
    session['uniquefolder'] = uniquefolder
    session['rxnifo'] = None
    session['status'] = False
    session['username'] = session['uniqueid']

def reset_session():
    uniqueid = str(uuid.uuid4())
    session['uniqueid'] = uniqueid

def run_session(rxntype, rxninfo, targets, direction, host, noMSA):
    uniqueid = session['uniqueid']
    uniquefolder = session['uniquefolder']
    csvfile = "selenzy_results.csv"
    success, app.config['TABLES'] = Selenzy.analyse(['-'+rxntype, rxninfo], 
                                                    targets,
                                                    app.config['DATA_FOLDER'],  
                                                    uniquefolder,
                                                    csvfile,
                                                    pdir = int(direction),
                                                    host = host,
                                                    NoMSA = noMSA,
                                                    pc = app.config['TABLES']
    ) # this creates CSV file in Uploads directory
    if success:
        data = pd.read_csv(file_path(uniqueid, csvfile))
        data.index = data.index + 1
        data.rename_axis('Select', axis="columns")
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
    data.rename_axis('Select', axis="columns")
    return data, csvfile, uniqueid


def maintenance(expDay=10):
    secs = expDay*24*60*60
    for folder in glob.glob(os.path.join(app.config['UPLOAD_FOLDER'], '*')):
        name = os.path.basename(folder)
        if name.startswith('debug'):
            continue
        modiftime = os.path.getmtime(folder)
        lapse = time.time() - modiftime
        if lapse > secs:
            # Double check that this an upload folder
            rxnfile = os.path.join(folder, name+'.rxn')
            if os.path.exists(rxnfile):
                try:
                    for x in glob.glob(os.path.join(folder, '*')):
                        os.unlink(x)
                except:
                    pass
                try:
                    os.rmdir(folder)
                except:
                    pass

        
class RestGate(Resource):
    """ REST interface, returns api info """
    def get(self):
        return {'app': 'Selenzy', 'version': '1.0', 'author': 'Synbiochem'}

class RestQuery(Resource):
    """ REST interface to Selenzy, by default it does not run the MSA to be faster. 
    We init an independent session for the REST request."""
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
            if 'host' in args:
                host = args['host']
            else:
                host = '83333'
            init_session()
            try:
                data, csvfile, session = run_session(rxntype, rxninfo, targets, direction, host, noMSA)
                return jsonify({'app': 'Selenzy', 'version': '1.0', 'author': 'Synbiochem', 'data': data.to_json()})
            except:
                return jsonify({'app': 'Selenzy', 'version': '1.0', 'author': 'Synbiochem', 'data': None})

class RestSource(Resource):
    """ REST interface, returns api info """
    def get(self):
        orgs = {}
        for seq in app.config['ORG']:
            orgs[app.config['ORG'][seq][1]] = app.config['ORG'][seq][0]
        return jsonify({'app': 'Selenzy', 'version': '1.0', 'author': 'Synbiochem', 'data': orgs})


api.add_resource(RestGate, '/REST')

api.add_resource(RestQuery, '/REST/Query')

api.add_resource(RestSource, '/REST/Source')


@app.errorhandler(404)
def page_not_found(e):
    return redirect(url_for('upload_form'))

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
    else:
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


@app.route('/msa', methods=['POST'])
def post_msa():
    """ Post safely the MSA """
    if request.method == 'POST':
        sessionid = json.loads(request.values['sessionid'])
        msafile = os.path.join(app.config['UPLOAD_FOLDER'], sessionid, 'sequences_aln.fasta')
        treefile = os.path.join(app.config['UPLOAD_FOLDER'], sessionid, 'sequences.dnd')
        if os.path.exists(msafile) and os.path.exists(treefile):
            msa = open(msafile).readlines()
            tree = open(treefile).readlines()
            return json.dumps({'msa': ''.join(msa), 'tree': ' '.join(tree)})
    return redirect ( url_for('upload_form') )

@app.route('/msaview', methods=['GET'])
def display_msa():
    """ Display the MSA """
    if request.method == 'GET':
        if 'id' in request.values:
            sessionid = request.values['id']
            msafile = os.path.join(app.config['UPLOAD_FOLDER'], sessionid, 'sequences_aln.fasta')
            if os.path.exists(msafile):
                return render_template('viewmsa.html', sessionid=sessionid)
    return redirect ( url_for('upload_form') )

@app.route('/display', methods=['POST'])
def display_reaction(marvin=app.config['MARVIN']):
    """ Display the reaction """
    if request.method == 'POST':
        size = (600,400)
        if 'file' in request.files and len(request.files['file'].filename) > 0:
            fileinfo = request.files['file']   
            if fileinfo.filename == '' or not allowed_file(fileinfo.filename):
                flash("No file selected")
                return redirect (request.url)
            rxninfo = save_rxn(fileinfo)
            success = True
            if len(rxninfo) == 0:
                success = False
                data = ''
            else:
                if marvin:
                    svgstream = Selenzy.display_reaction(rxninfo, outfolder=session['uniquefolder'], outname = str(uuid.uuid4()), marvin=True)
                    data = svgstream.decode('utf-8')
                    if len(data) == 0:
                        success = False
                else:
                    outfile, size = Selenzy.display_reaction(rxninfo, outfolder=session['uniquefolder'], outname = str(uuid.uuid4()), marvin=False)
                    if len(outfile) == 0:
                        success = False
                    data = os.path.join('/results', session['uniqueid'], 'files', os.path.basename(outfile))
                session['rxninfo'] = rxninfo
                session['rxntype'] = 'smarts'
                session['status'] = True
                success = True
            return json.dumps( {'data': data, 'status': session['status'], 'success': success, 'svg': marvin, 'size': size} )
        elif len(request.form['smarts']) > 0:
            outname = file_path(session['uniqueid'], session['uniqueid'])
            rxninfo = Selenzy.sanitizeSmarts(request.form['smarts'], outname)
            success = True
            if marvin:
                svgstream = Selenzy.display_reaction(rxninfo, outfolder=session['uniquefolder'], outname = str(uuid.uuid4()), marvin=True)
                data = svgstream.decode('utf-8')
                if len(data) == 0:
                    success = False
            else:
                outfile, size = Selenzy.display_reaction(rxninfo, outfolder=session['uniquefolder'], outname = str(uuid.uuid4()), marvin=False)
                if len(outfile) == 0:
                    success = False
                data = os.path.join('/results', session['uniqueid'], 'files', os.path.basename(outfile))
            session['rxninfo'] = rxninfo
            session['rxntype'] = 'smarts'
            session['status'] = True
            return json.dumps( {'data': data, 'status': session['status'], 'success': success, 'svg': marvin, 'size': size} )


@app.route('/sorter', methods=['POST'])
def sort_table():
    """ Sorts table """
    if request.method == 'POST':
        jfilter = json.loads(request.values.get('filter'))
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
        data.rename_axis('Select', axis="columns")
        return json.dumps( {'data': {'csv':  data.to_html(), 'filter': filt}} )

@app.route('/remover', methods=['POST'])
def delete_rows():
    """ Sorts table """
    if request.method == 'POST':
        selrows = json.loads(request.values.get('filter'))
        session = json.loads(request.values.get('session'))
        csvname = os.path.basename(json.loads(request.values.get('csv')))
        outdir = os.path.join(app.config['UPLOAD_FOLDER'], session)
        csvfile = os.path.join(outdir, csvname)
        head, rows = Selenzy.read_csv(csvfile)
        filt = []
        for i in selrows:
            try:
                index = int(i) - 1
                filt.append(index)
            except:
                continue
        newtargets = []
        newrows = []
        for j in range(0, len(rows)):
            if j not in filt:
                newtargets.append(rows[j][0])
                newrows.append(rows[j])
        fastaFile = os.path.join(outdir, "sequences.fasta")
        Selenzy.write_fasta(fastaFile, newtargets, app.config['TABLES'])
        # Avoid issues with sequence ids
        fastaShortNameFile = os.path.join(outdir, "seqids.fasta")
        Selenzy.write_fasta(fastaShortNameFile, newtargets, app.config['TABLES'], short=True)
        # Recompute MSA if exists
        dndFile = os.path.join(outdir, 'sequences.dnd')
        if os.path.exists(dndFile):
            cons = Selenzy.doMSA(fastaShortNameFile, outdir)
            for i in range(0, len(newrows)):
                try:
                    newrows[i][head.index('Consv. Score')] = cons[newrows[i][head.index('Seq. ID')]]
                except:
                    pass
        Selenzy.write_csv(csvfile, head, newrows)
        data = pd.read_csv(csvfile)
        data.index = data.index + 1
        data.rename_axis('Select', axis="columns")
        return json.dumps( {'data': {'csv':  data.to_html()}} )
                          


@app.route('/debug', methods=['GET'])
def show_table():
    if app.debug == True:
        csvfile = os.path.join('uploads', 'debug', 'selenzy_results.csv')
        data = pd.read_csv(csvfile)
        data.index = data.index + 1
        sessionid = 'debug'
        data.rename_axis('Select', axis="columns")
        return render_template('results.html', tables=data.to_html(), csvfile=csvfile, sessionid=sessionid, flags={'fasta': False, 'msa': False})
    else:
        return redirect ( url_for('upload_form') )

@app.route('/results', methods=['GET', 'POST'])
def upload_file():
    if request.method == 'POST':
        """ The POST request should come from an already initalised session """
        if 'uniqueid' not in session:
            return redirect ( url_for('upload_form') )
        # check if post request has smarts part
        if 'csv' in request.files  and len(request.files['csv'].filename) > 0:
            fileinfo = request.files['csv']   
            if fileinfo.filename == '' or not allowed_file(fileinfo.filename):
                flash("No file selected")
                return redirect (request.url)
            data, csvfile, sessionid = retrieve_session(fileinfo)
            return render_template('results.html', tables=data.to_html(), csvfile=csvfile, sessionid=sessionid, flags={'fasta': False, 'msa': False})
        else:
            try:
                rxninfo = session['rxninfo']
                rxntype = session['rxntype']
            except:
                return redirect(url_for('login'))

        direction = 0
        noMSA = False
        targets = request.form['targets']
        host = request.form['host']
        if request.form.get('direction'):
            direction = 1
        if request.form.get('noMSA'):
            noMSA = True
        try:
            data, csvfile, sessionid = run_session(rxntype, rxninfo, targets, direction, host, noMSA)
            return render_template('results.html', tables=data.to_html(), csvfile=csvfile, sessionid=sessionid, flags={'fasta': True, 'msa': not noMSA})
        except:
            return redirect( url_for("upload_form") )
    elif request.method == 'GET':
        """ A GET request would require an independently initialised session """
        init_session()
        smarts = request.args.get('smarts')
        if smarts is None:
            return redirect( url_for("upload_form") )
        host = request.args.get('host')
        if smarts is None:
            host = '83333'
        rxntype = 'smarts'
        rxninfo = smarts
        direction = 0
        noMSA = False
        targets = 20
        session['rxninfo'] = rxninfo
        session['rxntype'] = rxntype
        try:
            data, csvfile, sessionid = run_session(rxntype, rxninfo, targets, direction, host, noMSA)
            return render_template('results.html', tables=data.to_html(), csvfile=csvfile, sessionid=sessionid, flags={'fasta': True, 'msa': not noMSA})
        except:
            return redirect( url_for("upload_form") )
    return redirect( url_for("upload_form") )
    
@app.route('/results/<sessionid>/files/<filename>')
def results_file(sessionid,filename):
    return send_from_directory(os.path.join(app.config['UPLOAD_FOLDER'], sessionid), filename)

# Reconfigure for gunicorn
if __name__== "__main__":  #only run server if file is called directly

    arg = arguments()
    ALLOWED_EXTENSIONS = set(['rxn', 'smi', ' '])

    app.config['UPLOAD_FOLDER'] = os.path.abspath(arg.upload_folder)
    app.config['DATA_FOLDER'] = os.path.abspath(arg.datadir)

    if arg.d:
        app.config['DEBUG'] = True
        app.config['PRELOAD'] = False
    else:
        app.config['DEBUG'] = False
        app.config['PRELOAD'] = True        

    app.config['ORG'] = Selenzy.seqOrganism(arg.datadir, "seq_org.tsv")

    if app.config['PRELOAD']:
        app.config['TABLES'] = Selenzy.readData(arg.datadir)
    else:
        app.config['TABLES'] = None

    app.run(host="0.0.0.0",port=5000, debug=app.config['DEBUG'], threaded=True)
#    app.run(port=5000, debug=True)

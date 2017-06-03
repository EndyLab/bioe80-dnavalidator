"""
Flask Documentation:     http://flask.pocoo.org/docs/
Jinja2 Documentation:    http://jinja.pocoo.org/2/documentation/
Werkzeug Documentation:  http://werkzeug.pocoo.org/documentation/

This file creates your application.
"""

import os
from flask import Flask, render_template, request, redirect, url_for
from flask_sqlalchemy import SQLAlchemy

import pandas as pd

import re

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

app = Flask(__name__)

app.config['SECRET_KEY'] = os.environ.get('SECRET_KEY', 'this_should_be_configured')
app.config['SQLALCHEMY_DATABASE_URI'] = os.environ.get('DATABASE_URL', 'postgres:///anton')

db = SQLAlchemy(app)

class Submission(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    submission_time = db.Column(db.DateTime, server_default=db.func.now())

    suid = db.Column(db.String(100))

    jcvi_g1 = db.Column(db.String)
    jcvi_g2 = db.Column(db.String)
    jcvi_gu = db.Column(db.String)

    ecoli_g1 = db.Column(db.String)
    ecoli_g2 = db.Column(db.String)
    ecoli_gu = db.Column(db.String)

    ec_g1 = db.Column(db.String)
    ec_g2 = db.Column(db.String)

    def __repr__(self):
        return "{} {} {}".format(self.suid, self.submission_time, self.jcvi_g1)

    def genes(self):
        return {
            'JCVI-G1': self.jcvi_g1,
            'JCVI-G2': self.jcvi_g2,
            'JCVI-GU': self.jcvi_gu,
            'ECOLI-G1': self.ecoli_g1,
            'ECOLI-G2': self.ecoli_g2,
            'ECOLI-GU': self.ecoli_gu,
            'EC-G1': self.ec_g1,
            'EC-G2': self.ec_g2
        }

### Load sequence data

jcvi_genes = pd.read_csv('./data/aad6253-Hutchison-SM-database-S1.csv', skiprows=[0])
syn3_genes = pd.concat([jcvi_genes[jcvi_genes['KD_S3'] == 'k'], jcvi_genes[jcvi_genes['KD_S3'] == 'r']])
assignments = pd.read_csv('./data/assignments.csv')

data = assignments.merge(syn3_genes[['Locus', 'ProteinSequence']], left_on='Known Gene 1', right_on='Locus', suffixes=['', '-G1'])
data = data.merge(syn3_genes[['Locus', 'ProteinSequence']], left_on='Known Gene 2', right_on='Locus', suffixes=['', '-G2'])
data = data.merge(syn3_genes[['Locus', 'ProteinSequence']], left_on='Unknown Gene', right_on='Locus', suffixes=['', '-GU'])

def get_student(suid):
    students = data['Email'].str.split('@').str[0].str.match(suid)
    if not students.any() or suid == '':
        return False

    student_genes = data[students]
    name = student_genes['First Name'].iloc[0] + " " + student_genes['Last Name'].iloc[0]
    genes = {
        'JCVI-G1': student_genes['Locus'].iloc[0],
        'JCVI-G2': student_genes['Locus-G2'].iloc[0],
        'JCVI-GU': student_genes['Locus-GU'].iloc[0]
    }

    aas = {
            'JCVI-G1': student_genes['ProteinSequence'].iloc[0],
            'JCVI-G2': student_genes['ProteinSequence-G2'].iloc[0],
            'JCVI-GU': student_genes['ProteinSequence-GU'].iloc[0]
    }

    return name, genes, aas

def validate_gene(dna, aa):
    notes = list()

    # Clean up DNA sequence input
    dna = re.sub('[^A-Z]+', '', dna.strip().upper())

    if dna is None or dna == "":
        return True, []

    if len(dna) % 3 != 0:
        notes.append("Sequence doesn't obey triplet code (length not a multiple of three).")

    if not dna[:3] == "ATG":
        notes.append("No Start")

    try:
        prot = Seq(dna, generic_dna).translate(table=11)
        prot_myco = Seq(dna, generic_dna).translate(table=4)
    except:
        notes.append("Error translating DNA")
        return False, notes

    if len(prot) == 0:
        notes.append("Sequence does not translate to protein.")

    if not (len(prot) > 0 and prot[-1] == "*"):
        notes.append("No stop codon.")
        prot += "*"

    if aa is not None:
        if not prot[:-1] == aa:
            if prot_myco[:-1] == aa:
                notes.append("Not codon-optimized for E. coli.")
            else:
                notes.append("Translated DNA does not match target amino acid sequence.")

    if len(notes) == 0:
        return True, ['DNA meets validation criteria.']

    return False, notes

def validate_genes(genes, aas=dict()):
    dnas = dict()
    okays = dict()
    gene_feedback = dict()

    for gene in ['JCVI-G1', 'JCVI-G2', 'JCVI-GU', 'ECOLI-G1', 'ECOLI-G2', 'ECOLI-GU', 'EC-G1', 'EC-G2']:
        dna = genes[gene]
        if dna is None or dna.strip() == "":
            continue

        aa = None
        if gene in aas:
            aa = aas[gene]

        okay, notes = validate_gene(dna, aa)

        dnas[gene] = dna
        okays[gene] = okay
        if notes:
            gene_feedback[gene] = notes

    return  okays, dnas, gene_feedback

###
# Routing for your application.
###

@app.route('/')
def home():
    """Render website's home page."""
    return render_template('home.html')

@app.route('/dna', methods=['GET'])
def dna():
    """Render the main DNA validation page."""

    suid = request.args.get('suid', '')
    student = get_student(suid)

    if not student:
        return redirect(url_for('home'))

    name, genes, aas = student

    okay, dna, feedback = dict(), dict(), dict()

    last_sub = Submission.query.filter_by(suid=suid).order_by(db.desc(Submission.submission_time)).first()
    if last_sub:
        okay, dna, feedback = validate_genes(last_sub.genes(), aas)

    return render_template('dna.html', suid=suid, name=name, genes=genes, okay=okay, dna=dna, feedback=feedback)

@app.route('/dna', methods=['POST'])
def dna_submit():
    suid = request.form['suid']
    student = get_student(suid)

    if not student:
        return redirect(url_for('home'))

    name, genes, aas = student

    sub = Submission()
    sub.suid = suid

    okays, dnas, gene_feedback = validate_genes(request.form, aas)

    for gene, dna in dnas.items():
        setattr(sub, gene.lower().replace("-", "_"), dna)

    db.session.add(sub)
    db.session.commit()

    return render_template('dna.html', suid=suid, name=name, genes=genes, okay=okays, dna=dnas, feedback=gene_feedback)

@app.route('/about/')
def about():
    """Render the website's about page."""
    return render_template('about.html')


###
# The functions below should be applicable to all Flask apps.
###

@app.route('/<file_name>.txt')
def send_text_file(file_name):
    """Send your static text file."""
    file_dot_text = file_name + '.txt'
    return app.send_static_file(file_dot_text)


@app.after_request
def add_header(response):
    """
    Add headers to both force latest IE rendering engine or Chrome Frame,
    and also to cache the rendered page for 10 minutes.
    """
    response.headers['X-UA-Compatible'] = 'IE=Edge,chrome=1'
    # response.headers['Cache-Control'] = 'public, max-age=600'
    return response


@app.errorhandler(404)
def page_not_found(error):
    """Custom 404 page."""
    return render_template('404.html'), 404


if __name__ == '__main__':
    app.run(debug=True)

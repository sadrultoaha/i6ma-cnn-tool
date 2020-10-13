from flask import Flask, request, jsonify, render_template, redirect, send_file,flash
import os
from werkzeug.utils import secure_filename
from Bio import SeqIO
from model import prediction
import pandas as pd

UPLOAD_FOLDER = "static/uploads/"
ALLOWED_EXTENSIONS = ["FASTA","fasta"]

app = Flask(__name__)

app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config["ALLOWED_EXTENSIONS"] = ALLOWED_EXTENSIONS
app.config['SECRET_KEY'] = b'\xf1)\xfbn\xf95\x1a,\xad\xe5/\x1b\x87\nzJ'
#app.config["IMAGE_UPLOADS"] = "/uploads"
#app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024

ENV = 'prod'
if ENV == 'dev':
    app.debug = True
else:
    app.debug = False

def allowed_files(filename):
    # We only want files with a . in the filename
    if not "." in filename:
        return False
    # Split the extension from the filename
    
    ext = filename.rsplit(".", 1)[1]
    #print(ext.upper())
    # Check if the extension is in ALLOWED_EXTENSIONS
    if ext.upper() in app.config["ALLOWED_EXTENSIONS"]:
        return True
    else:
        return False

@app.route('/')
def home():
	return render_template('index.html')

@app.route('/info')
def info():
	return render_template('info.html')

@app.route("/", methods=["GET", "POST"])
def upload_image():
    if request.method == "POST":
        if request.files:
            fasta = request.files["fasta"]
            if fasta.filename == "":
                #print("No filename")
                flash('Warning! Select a fasta file before upload and run')
                return redirect(request.url)
            if allowed_files(fasta.filename):
                filename = secure_filename(fasta.filename)
                path = os.path.join(app.config["UPLOAD_FOLDER"], filename)
                #print(path)
                fasta.save(path)
                #print("fasta saved")
                #data = handle_fa(path)
                prediction(path)
                if os.path.exists(path):
                    if filename != 'sample.fasta':
                        os.remove(path)
                filename = 'prediction.csv'
                return redirect('/downloadfile/'+ filename)
            else:
                #print("That file extension is not allowed")
                flash('Error! That file extension is not allowed, upload only fasta file')
                return redirect(request.url)
    return render_template('index.html')


# Download API

@app.route('/fasta_file')
def return_fasta():
    file_path =  UPLOAD_FOLDER + 'sample.fasta'
    return send_file(file_path, as_attachment=True)

@app.route("/downloadfile/<filename>", methods = ['GET'])
def download_file(filename):
    flash('Prediction is successfully compleated !')
    return render_template('download.html',value=filename)

@app.route('/return-files/<filename>')
def return_files_tut(filename):
    file_path =  filename
    return send_file(file_path, as_attachment=True)

if __name__ == "__main__":
    app.run()
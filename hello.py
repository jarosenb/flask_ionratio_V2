# Flask wrapper for the isotope calculator.
__author__ = 'jakerosenberg'

from flask import Flask, make_response, request
import os
from output_formatter import result_table



app = Flask(__name__)
uploadfolder = './uploads'

@app.route('/')
def form():
    return """
        <html>
        <body>
        <h1>a-ion ratio calculator V2</h1>
        <form action="/results" method="post" enctype="multipart/form-data">

        Sequence (add modifications like in Prospector, IAM = +57.021464):<br>
        <input type="text" name="seq" size = "100%"> <br><br>
        Maximum charge state: <br>
        <input type="maxc" name="maxcharge"> <br><br>
        Upload MzML file: <br>
        <input type="file" name="data_file" /> <br><br>
        <input type="submit" />

        </form>


        </body>
        </html>
        """


@app.route('/results', methods=["POST"])
def transform_view():
    file = request.files['data_file']
    seq = request.form['seq']
    maxc = int(request.form['maxcharge'])
    if not file:
        return "No file"
    file.save('./uploads/' + file.filename)
    file_contents = './uploads/' + file.filename

    the_result = result_table(file_contents, seq, maxc)

    os.system('rm {0}'.format(file_contents))

    return the_result



app.run(debug = True)


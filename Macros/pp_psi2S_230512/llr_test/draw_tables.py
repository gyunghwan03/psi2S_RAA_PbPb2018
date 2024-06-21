import os
import glob
import subprocess

def compile_to_jpeg(tex_file):
    subprocess.run(['pdflatex', tex_file])
    pdf_file = os.path.splitext(tex_file)[0] + '.pdf'

    jpeg_file = os.path.splitext(tex_file)[0] + '.jpeg'
    subprocess.run(['convert', '-density', '300', pdf_file, '-quality', '90', jpeg_file])


txt_files = glob.glob('*.txt')
for txt_file in txt_files:
    compile_to_jpeg(txt_file)

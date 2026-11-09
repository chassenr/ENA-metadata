import sys
import csv
import os
import re
import requests
import xml.etree.ElementTree as ET
from xml.dom import minidom
from os import listdir
import sys
from tqdm import tqdm
import gzip

file = sys.argv[1] # input file name (xml)
outdir = sys.argv[2] # folder where csv-files should be saved

print(file)

tree = ET.parse(file)
root = tree.getroot()

pid = root.findall('.//SAMPLE/IDENTIFIERS')[0].find('PRIMARY_ID').text

ltag = [] # create a list with all values marked as "TAG"
lvalue = [] # create a list with all values marked as "VALUE"

def clean_text(text):
    if text is None:
        return ""
    return re.sub(r"\s+", " ", text).strip()

for x in root.findall('.//SAMPLE_ATTRIBUTE'):  # parse all from SAMPLE_ATTRIBUTES
    tag = clean_text(x.find('TAG').text)
    value = clean_text(x.find('VALUE').text)
    ltag.append(tag)
    lvalue.append(value)

filename = f"{pid}.csv"
path = os.path.join(outdir, filename)

csvfile = open(path,'w',encoding='utf-8')
csvfile_writer = csv.writer(csvfile, delimiter='\t', lineterminator='\n')

csvfile_writer.writerow(["PRIMARY ID", *ltag])
csvfile_writer.writerow([pid, *lvalue])

csvfile.close()

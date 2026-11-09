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

pid = root.findall('.//EXPERIMENT/IDENTIFIERS')[0].find('PRIMARY_ID').text

def clean_text(text):
    if text is None:
        return ""
    return re.sub(r"\s+", " ", text).strip()

ddesc = clean_text(root.findtext('.//DESIGN/DESIGN_DESCRIPTION'))

ltag = [] # create a list with all values marked as "TAG"
lvalue = [] # create a list with all values marked as "VALUE"

for x in root.findall('.//EXPERIMENT_ATTRIBUTE'): # parse all from EXPERIMENT_ATTRIBUTES
    # add rule for if either TAG or VALUE is not present.
    tag_elem = x.find('TAG')
    value_elem = x.find('VALUE')    

    tag = clean_text(tag_elem.text if tag_elem is not None else None)
    value = clean_text(value_elem.text if value_elem is not None else None)

    # Only keep pairs where both are non-empty
    if tag and value:
        ltag.append(tag)
        lvalue.append(value)


# LIBRARY_DESCRIPTOR
lib = root.find('.//LIBRARY_DESCRIPTOR')
lib_dict = {}
locus = lib.find("TARGETED_LOCI/LOCUS")
if locus is not None:
    lib_dict["TARGETED_LOCUS_NAME"] = locus.get("locus_name", "")
    lib_dict["TARGETED_LOCUS_DESC"] = locus.get("description", "")

filename = f"{pid}.csv"
path = os.path.join(outdir, filename)

csvfile = open(path,'w',encoding='utf-8')
csvfile_writer = csv.writer(csvfile, delimiter='\t', lineterminator='\n')

csvfile_writer.writerow(["PRIMARY ID", "DESIGN DESCRIPTION", *ltag, *lib_dict.keys()])
csvfile_writer.writerow([pid, ddesc, *lvalue, *lib_dict.values()])

csvfile.close()

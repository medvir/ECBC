#! /opt/python2.7/bin/python2.7

import json
jsonFile = open('/data/AbX/ECBC/Test_S1_ab-analysis/Test_S1_IgG1_ECBC_panda_blastout.json', 'r')
values = json.load(jsonFile)
jsonFile.close()

for criteria in values['criteria']:
    for key, value in criteria.iteritems():
        print key, 'is:', value
    print ''
	
	


import json
from pprint import pprint

with open('/data/AbX/ECBC/Test_S1_ab-analysis/Test_S1_IgG1_ECBC_panda_blastout.json') as data_file:
	data = json.load(data_file)

pprint(data)

data["maps"][0]["id"]
data["masks"]["id"]
data["om_points"]



	
	
	
#su -c 'service mongod start'
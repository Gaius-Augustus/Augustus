#!/usr/bin/env python3

import json

config = {
    'dbname':   'aug_vertebrates',
    'dbhost':   '127.0.0.1',
    'dbuser':   'augustus',
    'dbpasswd': 'aug_passwd',
    'cpuno':    '2' 
    }


with open('testconfig.json', 'w') as file:
    json.dump(config, file, indent=4, sort_keys=True)

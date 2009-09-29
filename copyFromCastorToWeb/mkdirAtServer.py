#!/usr/bin/env python
import urllib
import urllib2
import sys, os, string, fileinput
from sys import argv

url = 'http://higgs.skku.ac.kr/test/zip.php'
folder = sys.argv[1] 
values = {'dirn' : folder, 'modez' : 'md' }

data = urllib.urlencode(values)
req = urllib2.Request(url, data)
response = urllib2.urlopen(req)
the_page = response.read()

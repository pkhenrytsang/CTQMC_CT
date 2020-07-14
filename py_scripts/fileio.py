# File IO operations
# Created by Pak Ki Henry Tsang, Jun 15 2020

def CreateFolder(path,quiet=True):
	import os
	try: 
		os.makedirs(path) 
		if not(quiet):
			print("Directory '%s' created successfully" % path) 
	except OSError as error: 
		if not(quiet):
			print("Directory '%s' can not be created" % path) 

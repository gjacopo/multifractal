## Home directory
HOME = /home/grazzja

## Variables for Makefile locations
# You may have to change these lines
MAKEDIR=make
REL=..
CUR=.

include $(REL)/$(MAKEDIR)/makefile.defs

# List of source files to compile
# In each papckage directory, the FILELIST file gives the list
# of files that need to be compiled.
FILELIST=FileList
DIR=$(CUR)/

# Location of the sources
include $(CUR)/$(FILELIST)
# provides the variables 

include $(REL)/$(MAKEDIR)/makefile.flags

include $(REL)/$(MAKEDIR)/makefile.rules

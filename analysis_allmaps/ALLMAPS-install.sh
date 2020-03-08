#!/bin/bash

# Python dependencies
#wget https://bootstrap.pypa.io/ez_setup.py -O - | python - --user
#~/.local/bin/easy_install jcvi

# Concorde
wget http://www.math.uwaterloo.ca/tsp/concorde/downloads/codes/linux24/concorde.gz
gunzip concorde.gz
chmod u+x concorde

# Kent utilities
# Please download liftOver directly from UCSC
# Either from: http://hgdownload.cse.ucsc.edu/downloads.html for personal and non-profit
#          Or: https://genome-store.ucsc.edu/ for commercial access
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faSize
chmod u+x faSize liftOver
cp concorde faSize liftOver ~/.local/bin/

# Setup python
export PYTHONPATH=$HOME/.local/lib/python2.7/site-packages:$PYTHONPATH
export PATH=$HOME/.local/bin:$PATH

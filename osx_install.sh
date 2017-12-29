#!/bin/bash
#-----------------------------------------------------------------------------------------------
# osx_install v 1.0.3
# Michael G. Campana, 2017
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------------
curl -sSL https://get.rvm.io | bash -s stable
source ~/.rvm/scripts/rvm
rvm install 2.4.1
rvm --default use 2.4.1
gem install tk
mkdir $HOME/baitstools
chmod +x *.rb
mv *.rb $HOME/baitstools/
echo 'export PATH="$PATH:$HOME/baitstools"' >> $HOME/.bash_profile
exec bash

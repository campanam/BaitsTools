#!/bin/bash
#-----------------------------------------------------------------------------------------------
# osx_install v 1.6.5
# Michael G. Campana, 2017-2019
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------------
curl -sSL https://get.rvm.io | bash -s stable
source ~/.rvm/scripts/rvm
rvm install 2.6.5
rvm --default use 2.6.5
gem install tk
mkdir $HOME/baitstools
chmod +x *.rb
mv *.rb $HOME/baitstools/
echo 'export PATH="$PATH:$HOME/baitstools"' >> $HOME/.bash_profile
exec bash

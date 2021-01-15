#!/bin/bash
#-----------------------------------------------------------------------------------------------
# osx_install v 1.6.8
# Michael G. Campana, 2017-2020
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------------
curl -sSL https://get.rvm.io | bash -s stable
source ~/.rvm/scripts/rvm
rvm install 2.7.2
rvm --default use 2.7.2
gem install tk
gem install shell
mkdir $HOME/baitstools
chmod +x *.rb
mv *.rb $HOME/baitstools/
echo 'export PATH="$PATH:$HOME/baitstools"' >> $HOME/.bash_profile
exec bash

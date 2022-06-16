#!/bin/bash
#-----------------------------------------------------------------------------------------------
# osx_install v 1.7.5
# Michael G. Campana, 2017-2022
# Smithsonian's National Zoo and Conservation Biology Institute
#-----------------------------------------------------------------------------------------------

curl -sSL https://get.rvm.io | bash -s stable
source ~/.rvm/scripts/rvm
rvm install 3.0
rvm --default use 3.0
gem build baitstools.gemspec
gem install ./baitstools-1.7.5.gem

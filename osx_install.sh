#!/bin/bash
#-----------------------------------------------------------------------------------------------
# osx_install v 1.7.8
# Michael G. Campana, 2017-2023
# Smithsonian's National Zoo and Conservation Biology Institute
#-----------------------------------------------------------------------------------------------

curl -sSL https://get.rvm.io | bash -s stable
source ~/.rvm/scripts/rvm
rvm install 3.1.2
rvm --default use 3.1.2
gem build baitstools.gemspec
gem install ./baitstools-1.7.8.gem

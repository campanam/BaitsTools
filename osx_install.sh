#!/bin/bash
curl -sSL https://get.rvm.io | bash -s stable --ruby=2.4.1
gem install tk
mkdir $HOME/baitstools
chmod +x *.rb
mv *.rb $HOME/baitstools/
echo 'export PATH="$PATH:$HOME/baitstools"' >> $HOME/.bash_profile
exec bash

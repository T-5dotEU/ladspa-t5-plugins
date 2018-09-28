#!/bin/bash

cd ..
VERSION=`cat .previous_version`
echo -n "Version (previous version was $VERSION: "
read VERSION

sudo mkdir -p debian/DEBIAN
sudo mkdir -p debian/usr/lib/ladspa
sudo chown -R jh debian
cat control.in | sed "s#Version: _VERSION_#Version: ${VERSION}#" > debian/DEBIAN/control

sudo cp plugins/*.so debian/usr/lib/ladspa
sudo chmod 644 debian/usr/lib/ladspa/*.so
sudo chown -R root.root debian/
dpkg --build debian && \
    mv debian.deb builds/ladspa-t5-plugins_${VERSION}_amd64.deb && \
    aptly repo add t-5 builds && \
    echo -n "$VERSION" > .previous_version

sudo rm -rdf debian

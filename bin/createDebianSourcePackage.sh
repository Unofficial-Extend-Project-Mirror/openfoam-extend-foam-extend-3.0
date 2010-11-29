#!/bin/sh

# TODO: checkout correct tag
GPGKey=${1}
ubuntuVersion="10.04"
myTag="1.6-extRelease"
pushOrigSources="-sa" ## yes
# pushOrigSources="-sd" ## no

## create debian changelog for a release
# git-dch --release
## or for a snapshot
# git-dch --snapshot


## create debian source package with *.orig.tar.gz and patches
git-buildpackage \
	--git-builder="debuild -i\.git -I.git -S ${pushOrigSources} -k${GPGKey}" \
	--git-upstream-branch=master \
	--git-debian-branch=packaging/ubuntu/${ubuntuVersion} \
	--git-ignore-new

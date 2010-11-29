#!/bin/sh
#set -x

for name in `ls openfoam-1.6-ext/usr/lib/OpenFOAM-1.6-ext/applications/bin/`; do

manname=`echo -n ${name}.1`
cat > $manname <<EOF
.TH NAME 1
.\" NAME $name, SECTION 1
.SH NAME
$name \-
.SH SYNOPSIS
.B $name
.br
.SH DESCRIPTION
This manual page documents briefly the
.BR $name
command.

You can find information in the openfoam manual page and on the
openfoam web page http://www.openfoam.org.

.SH AUTHOR
The Debian Scientific Computing Team
EOF
echo debian/$manname /usr/share/man/man1
done

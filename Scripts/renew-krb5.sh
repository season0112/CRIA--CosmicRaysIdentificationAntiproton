#!/bin/bash
cd $HOME
mv .cern.keytab .cern.keytab.old
read -s -p "Password: " PASSWORD
ktutil << EOF
addent -password -p sili@CERN.CH -k 1 -e arcfour-hmac-md5
$PASSWORD
wkt .cern.keytab
quit
EOF

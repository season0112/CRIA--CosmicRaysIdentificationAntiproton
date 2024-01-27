cd $1
find . -name "$2" | xargs ac_merge --resultdir $3 --resultprefix $4 --input

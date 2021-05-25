#!/bin/bash
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
cd $DIR

VERSION=`cat README.md | grep VERSION | cut -d " " -f 2`

DIST=$DIR/dist
echo $DIST
if [ ! -d "$DIST" ]; then
    mkdir $DIST
else
  echo 'already exists'
  rm -rf $DIST
  mkdir $DIST
fi

PY_DIR=`ls -d PyLink_*`
for py in $PY_DIR
do
   cd $DIR
   echo $py
   cp -r $py ${DIST}/
   cp PyLink/* ${DIST}/${py}/
   rm -rf ${DIST}/${py}/.idea/
   cp README.md ${DIST}/${py}/
   cp requirements* ${DIST}/${py}/
   cd $DIST

   if [[ $py == *"Windows" ]]; then
      zip -r ${py}-${VERSION}.zip ${py}
      rm -rf ${py}
   else
      tar cfpz ${py}-${VERSION}.tar.gz ${py}
      rm -rf ${py}
   fi
done

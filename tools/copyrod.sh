#!/bin/bash
# Create a copy for production server
SOURCE=../../selenzy
TARGET=../../selenzyPro
rsync -a --delete --exclude='.git/' --exclude='uploads/*' --exclude='notes/' --exclude='tools/' --exclude='log/*' --exclude='old/' --exclude='*~' $SOURCE/ $TARGET
cd ../../
tar -czvf selenzy.tar.gz selenzyPro
mv selenzy.tar.gz /var/www/html/selenzy

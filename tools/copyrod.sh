#!/bin/bash
# Create a copy for production server
SOURCE=../../selenzy
TARGET=../../selenzyPro
rsync -a --delete --exclude='.git/' --exclude='uploads/*' --exclude='notes/' --exclude='tools/' $SOURCE/ $TARGET

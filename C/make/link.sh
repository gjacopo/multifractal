#!/bin/bash

for d in *; do \
if [ -d  ~/Codes/$d ]; then \
cd $d; \
for f in *; do \
if [ -d "$f" ]; then \
if ! [ -d ~/Codes/$d/"$f-M" ]; then \
# echo ln -s `pwd`/$f ~/Codes/$d/"$f-M" ; \
echo "linking $f"; \
 ln -s `pwd`/$f ~/Codes/$d/"$f-M" ; \
fi; \
fi; \
done; \
cd ..; \
fi; \
done;

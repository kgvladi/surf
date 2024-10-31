#!/bin/bash
path=/home/scr/kgvladi/10000_700
user=kgvladi
let m=12
for l in `seq 0 $m`
 do
   ./gs $l > loggs$l &
 done

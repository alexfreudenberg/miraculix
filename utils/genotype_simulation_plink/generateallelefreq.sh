#!/bin/bash
awk 'NR>1 {print NR-1,1-$5}' $1 > $2


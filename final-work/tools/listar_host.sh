#!/bin/bash
nmap -sP $1 > .nmap_out.txt
./host_up.py .nmap_out.txt $2
rm .nmap_out.txt

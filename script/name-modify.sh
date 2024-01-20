
mutation=$1

echo "$mutation" |tr '[a-z]' '[A-Z]'|sed 's/DELINS/delins/g'|sed 's/DUP/dup/g'|sed 's/INS/ins/g'|sed 's/DEL/del/g'|sed 's/ //g'|sed 's/GLY/G/g'|sed 's/ALA/A/g'|sed 's/VAL/V/g'|sed 's/LEU/L/g'|sed 's/ILE/I/g'|sed 's/PRO/P/g'|sed 's/PHE/F/g'|sed 's/TYR/Y/g'|sed 's/TRP/W/g'|sed 's/SER/S/g'|sed 's/THR/T/g'|sed 's/CYS/C/g'|sed 's/MET/M/g'|sed 's/ASN/N/g'|sed 's/GLN/Q/g'|sed 's/ASP/D/g'|sed 's/GLU/E/g'|sed 's/LYS/K/g'|sed 's/ARG/R/g'|sed 's/HIS/H/g'



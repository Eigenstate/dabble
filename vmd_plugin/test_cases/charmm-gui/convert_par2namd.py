#!/usr/bin/env python

import sys

def main(argv):
    charmm_toppar_file = sys.argv[1]
    charmm_toppar_dir = sys.argv[2]
    input_file = open(charmm_toppar_file, 'r')
    prm_num = 0
    str_num = 0

    # count paramter files
    for line in input_file:
        if ( ".prm" in line ) or ( ".str" in line ) :
            for i in range(len(line.split())):
                string = line.split()[i]
                if ".prm" in string:
                    if "toppar" not in string:
                        continue
                    prm_num = prm_num + 1
                if ".str" in string:
                    str_num = str_num + 1

    # list parameter files
    prm_files = [""] * prm_num
    str_files = [""] * str_num
    prm_num = 0
    str_num = 0

    input_file = open(charmm_toppar_file, 'r')
    for line in input_file:
        if ( ".prm" in line ) or ( ".str" in line ) :
            for i in range(len(line.split())):
                string = line.split()[i]
                if ".prm" in string:
                    if "toppar" not in string:
                        continue
                    prm_files[prm_num] = string
                    prm_num = prm_num + 1
                if ".str" in string:
                    str_files[str_num] = string
                    str_num = str_num + 1
                   
    # convert parameter files into NAMD readable format    
    for i in range(prm_num):
        input_file = open(prm_files[i], 'r')
        output_file = open("namd/"+prm_files[i], 'w')
        for line in input_file:
            if line.startswith("ATOM") or line.startswith("MASS") or ( "BOM" in line ) or ( "WRN" in line ) or line.startswith("set") or line.startswith("if"):
                line = "!" + line
            output_file.write(line)

    for i in range(str_num):
        input_file = open(str_files[i], 'r')
        output_file = open("namd/"+str_files[i], 'w')
        flag4par = 0
        flag4nbfix = 0
        for line in input_file:
            if ( "read para" in line):                
                flag4par = 1
            if ( flag4par == 1 ):
                if ( "read para" in line ):
                    continue
                if line.startswith("ATOM") or line.startswith("MASS") or ( "BOM" in line ) or ( "WRN" in line ) or line.startswith("set") or line.startswith("if"):
                   line = "!" + line
                if ( "NBFix between carboxylate and sodium" in line ):
                    flag4nbfix = 1
                    continue
                if ( flag4nbfix == 1 ) and ( "*" in line ):
                    flag4nbfix = 0
                    continue

                output_file.write(line)                

if __name__ == "__main__":
   main(sys.argv[1:])

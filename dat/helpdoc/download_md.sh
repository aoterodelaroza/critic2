#! /bin/bash

list="00_main 01_inputoutput 02_arithmetics 03_crystal 04_molecule 05_write 06_fields 07_graphics 08_cpsearch 09_gradientpath 10_basinplot 11_integrate 12_nciplot 13_structure 14_stm 15_misc"

rm -f *.md
for i in $list ; do
    wget -q --output-document - https://raw.githubusercontent.com/aoterodelaroza/aoterodelaroza.github.io/master/_critic2/98_manual/${i}.md | awk '/---/{a = !a;next}!a{print}' > ${i}.md
done

echo "--------- 95% vs. Median"
for i in $(seq 1 1 9); do echo "Area $i $(awk '{if (NR>1) {x=$5/$4; y=$10/$9; z=$15/$14; print NR-1,x,y,z}}' Area_Results_$i.txt | sort -k1 -g | tail -n1)"; done
echo "---------"
echo "--------- Median vs. 5%"
for i in $(seq 1 1 9); do echo "Area $i $(awk '{if (NR>1) {x=$4/$3; y=$9/$8; z=$14/$13; print NR-1,x,y,z}}' Area_Results_$i.txt | sort -k1 -g | tail -n1)"; done

echo "---------"
echo "---------"
echo "---------"

echo "--------- 95% vs. Median"
for i in $(seq 1 1 9); do echo "Area $i $(awk '{if (NR>1) {x=$5/$4; y=$10/$9; z=$15/$14; print NR-1,x,y,z}}' Area_Results_$i.txt | sort -k2 -g | tail -n1)"; done
echo "---------"
echo "--------- Median vs. 5%"
for i in $(seq 1 1 9); do echo "Area $i $(awk '{if (NR>1) {x=$4/$3; y=$9/$8; z=$14/$13; print NR-1,x,y,z}}' Area_Results_$i.txt | sort -k2 -g | tail -n1)"; done

echo "---------"
echo "---------"
echo "---------"
echo "--------- 95% vs. Median"
for i in $(seq 1 1 9); do echo "Area $i $(awk '{if (NR>1) {x=$5/$4; y=$10/$9; z=$15/$14; print NR-1,x,y,z}}' Area_Results_$i.txt | sort -k3 -g | tail -n1)"; done
echo "---------"
echo "--------- Median vs. 5%"
for i in $(seq 1 1 9); do echo "Area $i $(awk '{if (NR>1) {x=$4/$3; y=$9/$8; z=$14/$13; print NR-1,x,y,z}}' Area_Results_$i.txt | sort -k3 -g | tail -n1)"; done
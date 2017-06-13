col_num=$(head -n1 $1 | grep -o "	" | wc -l)

echo "Number of columns in file: "$(($col_num+1))

row_num=$(wc -l $1 | awk '{print $1}')

echo "Number of lines in file: "$row_num

hom_num=$(grep "Homo sapiens" $1 | wc -l)

echo "Number of lines with Homo sapiens:"$hom_num

shortest_length=$(cut -f7 $1 | sort -n | head -n2 | tail -n1)

echo "Shortest sequence length: "$shortest_length

unique_species_num=$(cut -f6 $1 | uniq | wc -l)

echo "Number of unique spiecies: "$unique_species_num

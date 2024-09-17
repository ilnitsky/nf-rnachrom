[ -z "$1" ] && echo "Usage: $0 description_file path_to_input_data" >&2 && exit 1
cat "$1"
d1=$(cat "$1"  |  sed 's/[+?\-s][^\)]*\]//g')
echo ./alpha1 -s -e -t -i $2 -d $d1 -l 14 -m 14
./alpha1 -s -e -t -i $2 -d $d1 -l 14 -m 14
d2=$(cat "$1"  |  sed 's/[?!<][^\)]*)//g; s/b[^\)]*)/\n/g')
echo ./alpha2 "$2"_DNA "$2"_RNA $d2
./alpha2 "$2"_DNA "$2"_RNA $d2

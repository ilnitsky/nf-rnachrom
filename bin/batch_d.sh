[ -z "$1" ] && echo "Usage: $0 description_file path_to_input_data1 path_to_input_data2" >&2 && exit 1
cat "$1"
d2=$(cat "$1"  |  sed 's/[?!<][^\)]*)//g; s/b[^\)]*)/\n/g')
echo ./alpha2 "$2" "$3" $d2
./alpha2 "$2" "$3" $d2

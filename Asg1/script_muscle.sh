find $1 -name "*" -exec sh -c 'echo '$2'/$(basename "{}")_conv ;muscle -in "{}" -out '$2'"$(basename "{}")_conv"' \;

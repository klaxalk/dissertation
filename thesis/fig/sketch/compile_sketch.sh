# get the path to this script
APP_PATH=`dirname "$0"`
APP_PATH=`( cd "$APP_PATH" && pwd )`

cd "$APP_PATH"

NAME=${1%.sk}

sketch -T "$1" > "$NAME.tex"

ls

if [ ! "$?" -eq 0 ]; then

  mv "$NAME.tex" > "error.txt"
  
else

  latex "$NAME.tex"
  rm "$NAME.tex"

  dvips -Ppdf -G0 "$NAME.dvi"
  rm "$NAME.dvi"

  ps2pdf "$NAME.ps"
  pdfcrop "$NAME.pdf" "$NAME.pdf"

  rm "$NAME.ps"
  rm "$NAME.log"
  rm "$NAME.aux"
  rm "texput.log"
  
fi

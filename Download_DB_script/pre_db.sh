wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1qhc17ZSRop0hp5lrM2sAQiffv0PXuQ6r' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1qhc17ZSRop0hp5lrM2sAQiffv0PXuQ6r" -O DB_Pre.zip && rm -rf /tmp/cookies.txt  &&\


unzip DB_Pre.zip &&\
rm  DB_Pre.zip


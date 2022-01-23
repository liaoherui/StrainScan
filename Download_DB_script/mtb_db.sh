wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1jjpDVVN7yD4f7r8yp_30cIgfTJ5KoGUt' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1jjpDVVN7yD4f7r8yp_30cIgfTJ5KoGUt" -O DB_Mtb.zip && rm -rf /tmp/cookies.txt  &&\


unzip DB_Mtb.zip &&\
rm  DB_Mtb.zip


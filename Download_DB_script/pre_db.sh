wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1yclHCQYpatDChOO2lVjQorQ43K4H1ltz' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1yclHCQYpatDChOO2lVjQorQ43K4H1ltz" -O DB_Pre.zip && rm -rf /tmp/cookies.txt  &&\


unzip DB_Pre.zip &&\
rm  DB_Pre.zip


wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1otmTt98xu8YUiTh4kQDKo3e6wMaOpWSV' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1otmTt98xu8YUiTh4kQDKo3e6wMaOpWSV" -O DB_Ecoli.zip && rm -rf /tmp/cookies.txt  &&\


unzip DB_Ecoli.zip &&\
rm  DB_Ecoli.zip


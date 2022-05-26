wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=18pPGyHODRwYV_d-l-xOKiLkbj7WcHCR_' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=18pPGyHODRwYV_d-l-xOKiLkbj7WcHCR_" -O DB_Mtb.zip && rm -rf /tmp/cookies.txt  &&\


unzip DB_Mtb.zip &&\
rm  DB_Mtb.zip


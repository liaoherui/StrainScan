wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1ENEfaKN2Mjrt1RHGIylhLIcP4NYogIFb' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1ENEfaKN2Mjrt1RHGIylhLIcP4NYogIFb" -O DB_Lcr.zip && rm -rf /tmp/cookies.txt  &&\


unzip DB_Lcr.zip &&\
rm  DB_Lcr.zip

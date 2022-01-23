wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1wwGubVO9F4r0pwHqQWQ2_NC151WDn96I' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1wwGubVO9F4r0pwHqQWQ2_NC151WDn96I" -O DB_Sep.zip && rm -rf /tmp/cookies.txt  &&\


unzip DB_Sep.zip &&\
rm  DB_Sep.zip


wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1BAoi5u4JuPTapULjbZRgaBBW5JrNYd7P' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1BAoi5u4JuPTapULjbZRgaBBW5JrNYd7P" -O DB_Akk.zip && rm -rf /tmp/cookies.txt  &&\


unzip DB_Akk.zip &&\
rm  DB_Akk.zip


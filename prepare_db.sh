wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1q_1n8DyCU1GjDIOESuEqWbRV10R6fO5y' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1q_1n8DyCU1GjDIOESuEqWbRV10R6fO5y" -O database.tgz && rm -rf /tmp/cookies.txt
tar -zxvf database.tgz
rm database.tgz

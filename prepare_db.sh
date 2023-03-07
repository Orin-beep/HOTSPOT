wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1VIbfp35X5JMiA7BfOS3lNBycT73uFcTN' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1VIbfp35X5JMiA7BfOS3lNBycT73uFcTN" -O database.tgz && rm -rf /tmp/cookies.txt
tar -zxvf database.tgz
rm database.tgz

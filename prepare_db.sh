wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1pvFZ6HWPmdAA0ZSOXS8xK4IFgkJXU7sA' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1pvFZ6HWPmdAA0ZSOXS8xK4IFgkJXU7sA" -O database.tgz && rm -rf /tmp/cookies.txt
tar -zxvf database.tgz
rm database.tgz

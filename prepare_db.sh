wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1ZSTz3kotwF8Zugz_aBGDtmly8BVo9G4T' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1ZSTz3kotwF8Zugz_aBGDtmly8BVo9G4T" -O database.tar.gz && rm -rf /tmp/cookies.txt
tar -zxvf database.tar.gz
rm database.tar.gz

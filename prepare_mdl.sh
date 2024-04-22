wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1bnA1osvYDgYBi-DRFkP-HrvcnvBvbipF' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1bnA1osvYDgYBi-DRFkP-HrvcnvBvbipF" -O models.tar.gz && rm -rf /tmp/cookies.txt
tar -zxvf models.tar.gz
rm models.tar.gz

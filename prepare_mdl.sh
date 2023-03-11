wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1L6ogZhdAWJ7Ns8Hz59m7W2kFPMcMTC6u' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1L6ogZhdAWJ7Ns8Hz59m7W2kFPMcMTC6u" -O models.tgz && rm -rf /tmp/cookies.txt
tar -zxvf models.tgz
rm models.tgz

mkdir models

gdown https://drive.google.com/drive/folders/1DR2W_FzqbslPykQ9T_Q4r3aoOD0TxBtb -O models/mdl1 --folder

python uncompresse.py models/mdl1

gdown https://drive.google.com/drive/folders/1HjZqFS0EoTapbo4K24aYd-e3bbT0EVk9 -O models/mdl2 --folder

python uncompresse.py models/mdl2

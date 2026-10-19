pkill -f f02e_baixar_ortologos.sh; pkill -f f02e3_lotes_reverso.sh; pkill -f "datasets download gene"
sleep 2; ps -eo pid,args | grep -E "datasets|f02e" | grep -v grep
cd ~/rapinantes/locos/ortologos && rm -rf *.tmp *.rtmp *.zip; ls -d *.d

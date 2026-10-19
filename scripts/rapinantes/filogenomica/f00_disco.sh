df -h / | tail -1
du -sh ~/nf-work/filogenia/work ~/rapinantes ~/genomas ~/envs ~/.pixi /tmp 2>/dev/null
du -sh ~/* 2>/dev/null | sort -h | tail -8
docker system df 2>/dev/null | head -5

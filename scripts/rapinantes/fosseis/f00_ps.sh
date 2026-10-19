#!/usr/bin/env bash
ps -eo pid,pcpu,etime,args --sort=-pcpu | grep -E "iqtree3|astral" | grep -v grep | cut -c1-220 | head -5
uptime

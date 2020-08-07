#!/bin/sh
excite -p 240 f0 | mglsadf -m 49 -c 0 -a 0.55 -p 240 mgc > data.mgcep.syn
gwave data.mgcep.syn | xgr

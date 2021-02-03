#!/bin/bash
# proxy use for GIT --JB
(echo "CONNECT $1:$2 HTTP/1.0"; echo; cat ) | netcat proxyout.lanl.gov 8080| (read a; read a; cat )


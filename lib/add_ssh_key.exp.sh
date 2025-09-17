#!/usr/bin/expect
spawn ssh-add /Users/mingda/.ssh/id_rsa_jasmin
expect "Enter passphrase"
send "19930517\r"
expect eof
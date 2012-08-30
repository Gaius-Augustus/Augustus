CFLAGS := -Wall -ansi -pedantic -O2 ${CFLAGS}
CC      = g++
INCLS	= -I/usr/include/mysql -I/usr/local/include/mysql++
LIBS 	=  -lmysqlclient -lmysqlpp

faload2db: faload2db.cc
	$(CC) $(CFLAGS) $(INCLS) $(LIBS) -o faload2db faload2db.cc


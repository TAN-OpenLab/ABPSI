# MUPSI
MUPSI: Malicious and Updatable Multi-party Delegated Private Set Intersection Protocol

This project is based on Feather:  https://github.com/AydinAbadi/Feather/tree/master/Feather-implementation

# How to run MUPSI?

1. Download the above files, and put them in a folder "MUPSI".
2. cd MUPSI
3. Install
  * GMP: https://gmplib.org/
  * Cryptopp: https://www.cryptopp.com
  * NTL: https://www.shoup.net/ntl
  * Bloom filter: http://www.partow.net/programming/bloomfilter/index.html
  * Boost: https://www.boost.org/users/history/
4. run
  * g++ -c Rand.cpp -c Hashtable.cpp -c Polynomial.cpp -c Server.cpp -c Client.cpp
  * g++  -I$home/homeDirectory/include -I/homeDirectory/bloom_filter/bloom_filter.hpp Rand.o Hashtable.o Polynomial.o Server.o Client.o test.cpp  -o test  -lntl -lgmpxx -lgmp -lcryptopp -lpthread
  * ./test
    
  Note that in the above "homeDirectory" should be replaced with the name of your machine home directory.
  

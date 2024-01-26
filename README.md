# mersenne prime finder
This is a program to find Mersenne primes. It is broadly divided into the following parts:

- Determine that the Mersenne number is not prime by using 2p + 1.
 Euler and Lagrange proved the following about Sophie Germain primes: if p ≡ 3 (mod 4) and p > 3, then the prime 2p+1 divides the Mersenne number Mp. (https://t5k.org/glossary/page.php?sort=SophieGermainPrime)
 The basic algorithm was implemented using Fermat's little theorem.
If p is a prime and a is any integer not divisible by p, then a p − 1 − 1 is divisible by p.(https://en.wikipedia.org/wiki/Fermat%27s_little_theorem)


- Mersenne primes are determined using the Lucas-Lehmer primality test.
In mathematics, the Lucas–Lehmer test (LLT) is a primality test for Mersenne numbers. The test was originally developed by Édouard Lucas in 1878[1] and subsequently proved by Derrick Henry Lehmer in 1930. (https://en.wikipedia.org/wiki/Lucas%E2%80%93Lehmer_primality_test)

- The currently known Mersenne primes are:
2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217, 4253, 4423, 9689, 9941, 11213, 19937, 21701, 23209, 44497, 86243, 110503, 132049, 216091, 756839, 859433, 1257787, 1398269, 2976221, 3021377, 6972593, 13466917, 20996011, 24036583, 25964951, 30402457, 32582657, 37156667, 42643801, 43112609, 57885161, 74207281, 77232917, 82589933. (sequence A000043 in the OEIS)
# build
 gcc mersenne.c -o mersenne
# run
  ./mersenne option.config
  
# mersenne database
  The mersenne.dat file contains unresolved Mersenne numbers.
# question
 Kyojun.Kim ( kkjms2@gmail.com )
# reference
- https://t5k.org/
- https://en.wikipedia.org/wiki/Lucas%E2%80%93Lehmer_primality_test
- https://en.wikipedia.org/wiki/Mersenne_prime
- https://en.wikipedia.org/wiki/Fermat%27s_little_theorem

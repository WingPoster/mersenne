# Mersenne Prime Finder (MPF)

&nbsp;&nbsp;&nbsp;This program finds unknown Mersenne primes. A Mersenne number(Mp) is a number that is 1 less than the power of 2. The Mersenne number for exponent n can be expressed as $M_n = 2^n - 1$. The reason why finding Mersenne primes is difficult is because it requires calculation of very large numbers. The largest Mersenne prime known so far is $2^{82589933} - 1$. These Mersenne primes were officially discovered on December 7, 2018. It is not yet known whether there are infinite number of Mersenne primes.<br>
The currently known Mersenne primes are:

```2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217, 4253, 4423, 9689, 9941, 11213, 19937, 21701, 23209, 44497, 86243, 110503, 132049, 216091, 756839, 859433, 1257787, 1398269, 2976221, 3021377, 6972593, 13466917, 20996011, 24036583, 25964951, 30402457, 32582657, 37156667, 42643801, 43112609, 57885161, 74207281, 77232917, 82589933. (sequence A000043 in the OEIS)```

&nbsp;&nbsp;&nbsp;This program includes two algorithms. First, an algorithm is used to divide a Mersenne number by a prime number p and determine that it is not a Mersenne prime if the remainder is 0. Second, Lucas–Lehmer primality test (LLT) is used to determine whether the Mersenne number is prime. The reason the algorithm is divided into two areas is that LLT takes calculation time from several minutes to several months as Mersenne's exponent increases, so if the LLT algorithm is used for all Mersenne numbers, search may not be possible. In addition, if the remainder is 0 as a result of dividing a Mersenne prime by a prime p using the general trial factoring method (TFM), it can be determined that it is not a prime. This can be determined to be not a Mersenne number by dividing Mersenne numbers by a relatively small prime number p (p < 100,000,000). However, although it is possible to determine that this is not a prime number, it is difficult to determine whether the relevant Mersenne number is a Mersenne prime number. In order to determine a Mersenne prime, it is necessary to confirm that there are prime factors, but if the smallest number among the prime factors is a large number greater than 100 million, the general TFM takes a lot of time and may be impossible to calculate.<br>

| Process | Full Name |Description |
| --- |  --- | --- |
| LLT | Lucas-Lehmer Primality Testing |If $p$ is an odd prime, then $M_p = 2^p -1$ is prime if and only if $M_p$ evenly divides $S_{p-2}$, where $S_0 = 4$ and for $k > 0$, $S_k = S^{2}_{k-1} -2$.|
| TFM | Trial Factoring Method |Trial factoring is a factorization method that checks every prime less than the square root of the number in question to see if any of them divide that number, trying to finding a factor. If none are found, the number in question is prime. Otherwise, it is a composite number.|

&nbsp;&nbsp;&nbsp;Therefore, this program uses an efficient strategy to find Mersenne primes. First, we use the TFM algorithm to find the prime factor p that is the prime factor of Mersenne, and then we find p that divides the Mersenne number. This procedure continues until p is 10 billion or more and excludes numbers with prime factors among Mersenne numbers from the list of possible primes because they cannot be Mersenne primes. However, Mersenne numbers with large prime factors are not easy to find using TFM. Therefore, the deterministic LLT method is used to determine whether a Mersenne number for which no prime factors have been found is prime. Although this is a time-consuming task, the prime number determination process is performed only on the Mersenne number list, excluding Mersenne numbers with primarily small prime factors, so the operation can be efficiently optimized. Running the two algorithms in parallel finds whether the Mersenne number is a Mersenne prime, starting from the smallest number, and Mersenne primes are found sequentially from the smallest number to the largest number.<br>

The mathematical background knowledge used in this program is as follows:
| Mathematical Background | Reference |
| --- | --- |
| Euler and Lagrange proved the following about Sophie Germain primes: if p ≡ 3 (mod 4) and p > 3, then the prime 2p+1 divides the Mersenne number Mp.  | <a href="https://t5k.org/glossary/page.php?sort=SophieGermainPrime">Sophie Germain Prime</a> |
| The Lucas–Lehmer test (LLT) is a primality test for Mersenne numbers. The test was originally developed by Édouard Lucas in 1878 and subsequently proved by Derrick Henry Lehmer in 1930. | <a href=https://en.wikipedia.org/wiki/Lucas%E2%80%93Lehmer_primality_test>Lucas-Lehmer Primality Test</a>|

# Build
 Lastest Build by MacOS 12.7.2<br>
 ```gcc-13 -fopenmp mersenne.c -o mersenne```

# Run & Kill
```./mersenne option.config```
```./mersenne kill.config```
  
# mersenne database
  The mersenne.dat file contains unresolved Mersenne numbers.

# Question
 Kyojun.Kim ( kkjms2@gmail.com )

# License
  MIT License

# Reference
- https://t5k.org/
- https://en.wikipedia.org/wiki/Lucas%E2%80%93Lehmer_primality_test
- https://en.wikipedia.org/wiki/Mersenne_prime

#pragma once

#define _USE_MATH_DEFINES

#include <cmath>
#include <cstdint>
#include <random>
#include <iostream>

// This random number generator uses the xoshiro256++ algorithm. xoshiro256++ should be faster than the algorithms used by <random>: https://cplusplus.com/reference/random/
struct random_number_generator {
  uint64_t s[4] = {};
  uint64_t seed = 0;

  random_number_generator() {
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<uint64_t> dist;
    seed = dist(rng);
    std::cout << "Random Seed: " << std::hex << seed << std::dec << std::endl;
    init(seed);
  }

  random_number_generator(uint64_t seed) {
    this->seed = seed;
    init(seed);
  }

  void init(uint64_t seed) {
    s[0] = split_mix_64(seed);
    s[1] = split_mix_64(s[0]);
    s[2] = split_mix_64(s[1]);
    s[3] = split_mix_64(s[2]);
  }

  double uniform() { return double(next()) / double(UINT64_MAX); }

  int uniform_int(int min, int max) { return min + (next() % (max - min));}

  double normal(double mu = 0.0, double sigma = 1.0) {
      double u1 = uniform();
      double u2 = uniform();
      double x1 = sqrt(-2.0 * log(u1)) * cos(2 * M_PI * u2);
      double transformed = x1 * sigma + mu;
      return transformed;
  }

  /*  Written in 2018 by David Blackman and Sebastiano Vigna (vigna@acm.org)
  To the extent possible under law, the author has dedicated all copyright
  and related and neighboring rights to this software to the public domain
  worldwide. This software is distributed without any warranty.
  See <http://creativecommons.org/publicdomain/zero/1.0/>. */
  uint64_t next() {
    const uint64_t result = s[0] + s[3];

    const uint64_t t = s[1] << 17;

    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];

    s[2] ^= t;

    s[3] = rotl(s[3], 45);

    return result;
  }

  constexpr uint64_t rotl(const uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
  }

  constexpr uint64_t split_mix_64(uint64_t x) {
    uint64_t z = (x += 0x9e3779b97f4a7c15);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
    z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
    return z ^ (z >> 31);
  }

};

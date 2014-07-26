// -*-coding: mule-utf-8-unix; fill-column: 58; -*-
/**
 * @file
 * Calculates the "N-factor", it is yacoin/vertcoin concept
 * originally, see 
 * http://www.followthecoin.com/interview-creator-vertcoin/
 *
 * @author Sergei Lodyagin
 */

#include <chrono>
#include <cmath>
#include <stdexcept>
#include "types/time.h"
#include "n_factor.h"
#include "types/safe.h"
#include "pars.h"

//! Returns tuple with scrypt N,r,p parameters.
//! It differs from vertcoin. The assumption is doubling
//! a number of transistors (CPU cores) every 18 monthes
//! (Moore's law, N-factor) but increasing the speed of
//! CPU only by 10-20% (we use 15%).
//! Why increase N (amount of memory) but not p (a scrypt
//! parrallelism parameter)?  Because increasing of p not
//! protect us against ASIC/GPUs and increasing of CPU
//! cores means decreasing of block time checking (when a
//! node goes online and checks generated block in the
//! whole period from it was online last time). The block
//! checking period is a real limitation of scrypt memory
//! usage parameter (it is 128*N*r*p, the greater
//! parameter means better GPU protection but at the same
//! time greater block checking time on a single core
//! which is near proportional to N*r).
std::tuple<uint32_t, unsigned, unsigned>
GetNfactor(coin::times::block::time_point block_time) 
{
  using namespace coin::times::block;
  using namespace std::chrono;
  using namespace std;
  using namespace types;

  using days = chrono::duration<clock::rep, ratio<3600 * 24>>;
  using avg_months = chrono::duration
    <clock::rep, ratio<3600 * 24 * 30>>;

  using memory_size_t = safe<long long>;

  const memory_size_t MiB = 
    memory_size_t(1024LL) * memory_size_t(1024LL);

  //! The function birth time
  static const auto birth_time =
    coin::times::block::time_point(
      days(times::howard_hinnant::days_from_civil
        (2014, 05, 05))
    );

  if (block_time < birth_time)
    // invalid block time, use initial values
    block_time = birth_time;

  // he said "every 18 months"
  constexpr auto moore_period = avg_months(18);

  const int moore_steps = 
    duration_cast<avg_months>(block_time - birth_time)
    / moore_period; 

  const auto last_step_time =
    birth_time + moore_steps * moore_period;

  // the scrypt parameters
  memory_size_t p = 1; // increasing is not protect us
  const memory_size_t initial_mem = memory_size_t(1LL) * MiB;
  constexpr unsigned initial_r = 8;
  const auto moor_mem_multiplier = 
    memory_size_t(1ULL << moore_steps);
  memory_size_t mem = initial_mem * moor_mem_multiplier;
  assert(mem >= initial_mem);

  // today min amount of scrypt memory per device
  const memory_size_t birth_total_memory = 
    memory_size_t(128LL) * MiB;
  // today min amount of CPU cores
  constexpr size_t birth_cores = 2;

  // assume memory access speed will be not increased
  // much, so limit block load time to 1/6 of block period
  constexpr auto cell_access_time = nanoseconds(20);
  const memory_size_t max_mem = memory_size_t(
    std::min(
      pars::block_period_by_design.first, 
      pars::block_period_by_design.second
    ) / 6 / cell_access_time
  );

  while (mem / memory_size_t(p) > max_mem) {
    p *= 2;
  }

  unsigned long long moore_cores = (birth_cores << moore_steps);
  if (moore_cores < birth_cores)
    moore_cores = birth_cores;
  // memory used by scrypt when use moore_cores
  memory_size_t moore_total_memory =
    (birth_total_memory * moor_mem_multiplier);
  assert(moore_total_memory >= birth_total_memory);

  if (mem > moore_total_memory)
    mem = moore_total_memory; 
  // assert 2^x1

  if (mem / initial_mem > 
      memory_size_t(moore_cores)
     )
    mem = initial_mem * memory_size_t(moore_cores);
  // assert 2^x2

  // CPU speed grow 15% per moore_step, 
  // it is twiced in 5 steps
  const memory_size_t r0 = (initial_r << moore_steps / 5);
  const auto last_r_step_time = last_step_time - 
    moore_period * (moore_steps % 5);
  const auto r = r0 
    + r0 * (block_time - last_r_step_time) 
    / (5 * moore_period);

  // N must be power of 2
  const memory_size_t N = mem / (memory_size_t(128LL) * r0 * p);
  if (N > std::numeric_limits<uint32_t>::max())
    throw std::length_error("N factor is too large");
//  uint64_t N = 1 << (sizeof(N1) * 8 - __builtin_clz(N1) - 1);
  //if (N == 0) N = 1024;
  return make_tuple((uint32_t) N, (unsigned) r, (unsigned) p);
}

